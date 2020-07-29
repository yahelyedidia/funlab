import numpy as np
import pandas as pd
from itertools import islice
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.cm

GENES_B38 = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Homo_sapiens.GRCh38.98.gtf.gz"

GENES_B37 = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/genes/hg19.knownGene.gtf"

CSC = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/CSC/p_values_all_information_by_orig_vals.tsv"

DIR_CSC = DIR = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/CSC"

TAB = "\t"

BEDGRAPH_LINE_FORMAT = "s{i}\tchr{chr_name}\t{start}\t{number}\n"

COLUMNS = 1
REP_LIST = ['6', '10', '15']

"""number of cromosomes"""
NUM_OF_CHR = 24


def filter_data(filter, d, col, name):
    """
    A function that filter the data according do given arg and write it to new file
    :param filter: the arg to filtered by it
    :param d: the data
    :param col: the col to filter by
    :param name: the name of the filtered file
    :return: the data filters
    """
    data = pd.read_csv(d, sep=TAB)
    data = data.drop(data.columns[0], axis=COLUMNS)
    print("stop reading")
    header = data.columns.values.tolist()
    if header[0] != "chr":
        data = data.drop(data.columns[0], axis=1)
    print("stop drop")
    if filter >= 0:
        data = data[data[col] >= filter]
    else:
        data = data[data[col] <= filter]
    print("done filter")
    data = remove_duplicate(data, "chr", "start", "end", "change")
    print("no duplicate")
    pd.DataFrame(data=data).to_csv(name, sep=TAB)
    print("writing to file")
    return data


def remove_duplicate(data, chr_col, start, end, change=""):
    """
    get data and remove the duplicate site
    :param data: the data to fix as pd
    :param chr_col: a string that represent the name of the chrom column in the data
    :param start: the string that represent the name of the start column in the data
    :param end: the string that represent the name of the end column in the data
    :param change: the string that represent the name of the change column in the data
    :return: the fixed data
    """
    for ind, i in data.iterrows():
        for k, j in islice(data.iterrows(), 1, None):
            if ind >= k:
                continue
            if i[chr_col] == j[chr_col]:
                if abs(i[start] - j[start]) <= 10 and abs(i[end] - j[end]) <= 10:
                    data = data.drop(k)
                    continue
                if abs(i[start] - j[start]) <= 100 or abs(i[end] - j[end]) <= 100:
                    if change != "":
                        if i[change] == j[change]:
                            data = data.drop(k)
            continue
    return data


def read_genes_data(file, num_open_line=5, flag_38=False):
    """
    A function that reads the genes DB and filter it to the genes data
    :param file: the file to read
    :param num_open_line: the number of line to ignore in the begining of the file
    :return: the filtered data
    """
    if(flag_38):
        filter_by = 'gene'
        data = pd.read_csv(file, sep="\t", skiprows=[i for i in range(num_open_line)],
                           compression='gzip', header=None)
    else:
        filter_by = 'transcript'
        data = pd.read_csv(file, sep="\t", skiprows=[i for i in range(num_open_line)],
                            header=None)
    header = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    data.columns = header[:len(data.columns)]
    data = data[data['feature'] == filter_by]
    names, id = [], []
    data = data.drop_duplicates(['start', 'end'], keep='last')
    for row in data.iterrows():
        line = row[1]['attribute']
        sindex = line.find("gene_name")
        idinx = line.find("gene_id")

        # no gene name
        if sindex == -1:
            name_line = "no_name_found"
        else:
            name_line = line[sindex + 9:]
            name_line = name_line.split(";")[0]
        names.append(name_line)
        if(idinx == -1):
            id_line = "no_id_found"
        else:
            id_line = line[idinx + 7:]
            id_line = id_line.split(";")[0]
        id.append(id_line)
    data['attribute'] = names
    data['ids'] = id
    data['close_sites'] = [[] for i in range(data.shape[0])]
    if flag_38:
        data.to_csv("genes" + os.path.sep + "genes_19.csv", sep="\t", compression='gzip', index=False)
    else:
        data.to_csv("genes" + os.path.sep + "genes_19.csv", sep="\t", index=False)
    return data


def create_gene_data(flag_38):
    """
    creating an array of chromosomes with the genes data
    :param flag_38: flag to choose the genome build file
    :return: an array with genes data sorted by chromosomes
    """
    if flag_38:
        file = GENES_B38
    else:
        file = GENES_B37
    data = read_genes_data(file, flag_38)
    chroms = []
    for chr in range(NUM_OF_CHR - 2):
        if flag_38:
            char_name = chr + 1
        else:
            char_name = "chr{0}".format(chr + 1)
        chr_data = data[data['chr'] == char_name]
        chroms.append(chr_data)
    if flag_38:
        chr_x = 'X'
        chr_y = 'y'
    else:
        chr_x = 'chrX'
        chr_y = 'chrY'
    chr_data = data[data['chr'] == chr_x]
    chroms.append(chr_data)
    chr_data = data[data['chr'] == chr_y]
    chroms.append(chr_data)
    return chroms


def find_close_genes(filter, gene_data, site_file, name, i=False, csc=False,):
    """
    A function that finds which genes are close to the CTCF biding sites and count to how
    many site the gene is close.
    :param filter: the radius to look at
    :param gene_data: thr gene data
    :param site_file: the data of sites
    :param name: the final file's name
    :return: thr genes dict
    """
    gene_dict = {}
    data_sites = pd.read_csv(site_file, sep="\t")
    if csc:
        data_sites = data_sites[['ID_REF', 'chr', 'start', 'end', 'controls_{0}_month'.format(i), 'afters{0}_month'.format(i), 'p_values_{0}_month'.format(i)]]
        data_sites = data_sites[data_sites['p_values_{0}_month'.format(i)] <= 0.05]
    else:
        data_sites = data_sites[data_sites['p value'] <= 0.05]
    # print("pass")
    add_gene = []
    for site in data_sites.iterrows():
        genes = []
        fs = site[1]['start'] - filter
        fe = site[1]['end'] + filter
        chr = int(site[1]['chr'])
        for gene in gene_data[chr - 1].iterrows():
            if fs <= gene[1]['start'] and gene[1]['end'] <= fe:
                genes.append(gene[1]['attribute'])
                if gene[1]['attribute'] in gene_dict:
                    gene_dict[gene[1]['attribute']] += 1
                else:
                    gene_dict[gene[1]['attribute']] = 1
                gene[1]['close_sites'].append((chr, fs, fe))
        add_gene.append(genes)
    data_sites['close_genes'] = add_gene
    print()
    if csc:
        data_sites.to_csv("genes" + os.path.sep + "genes_close_to_sites_{1}_filter_{0}.csv".format(filter, name.format(i)), sep="\t", index=False)
        merge_genes_data = pd.concat(gene_data)
        merge_genes_data.to_csv("genes" + os.path.sep + "sites_close_to_genes_{1}_filter_{0}.csv".format(filter, name.format(i)), sep="\t", index=False)
    else:
        data_sites.to_csv("genes" + os.path.sep + "genes_close_to_sites_{1}_filter_{0}.csv".format(filter, name), sep="\t", index=False)
        merge_genes_data = pd.concat(gene_data)
        merge_genes_data.to_csv("genes" + os.path.sep + "sites_close_to_genes_{1}_filter_{0}.csv".format(filter, name), sep="\t", index=False)

    return gene_dict


def check_with_change_filter(list_of_filters, num_to_print, file_to_check, name, csc=False, flag_38=False):
    """
    A function that gets list of filter to look at, and print the most repetitive genes
    :param list_of_filters: the filters in list
    :param num_to_print: the num of repetitive elements to print
    """
    chroms = create_gene_data(flag_38)
    print("yay data")
    for f in list_of_filters:
        if csc:
            for label in REP_LIST:
                finds_and_print_genes(chroms, f, file_to_check, name, num_to_print, label, csc=True)
        else:
            finds_and_print_genes(chroms, f, file_to_check, name, num_to_print)


def finds_and_print_genes(chroms, f, file_to_check, name, num_to_print, p_label='p value', csc=False):
    d = find_close_genes(f, chroms, file_to_check, name, p_label, csc)
    print("dictionary after filter {0}".format(f))
    print("number of genes: {0}".format(len(d)))
    print_top_values(num_to_print, d)


def print_top_values(num_to_print, d):
    """
    A function that get dictionary and number of items to print and
    :param num_to_print:
    :param d: the dictionary to print
    """
    if num_to_print > len(d):
        num_to_print = len(d)
    counter = num_to_print
    while counter > 0:
        max_value = max(d.values())  # maximum value
        max_keys = [k for k, v in d.items() if v == max_value]
        if len(max_keys) > counter:
            break
        print("value is: {0}, genes with that value: {1}". format(max_value, max_keys))
        counter -= len(max_keys)
        for key in max_keys:
            del d[key]


def convert_csv_to_cn(file, s):
    """
    old function of converting csv files to IGV format
    """
    path = os.path.dirname(file)
    name = os.path.relpath(file)
    if name.endswith(".csv"):
        name = name.replace(".csv", ".cn")
    else:
        name = path + os.sep + name + ".cn"
    cn_file = open(name, 'w')
    csv_file = pd.read_csv(file, sep=s)
    csv_file = csv_file.drop(csv_file.columns[0], axis=1)
    csv_file = csv_file.drop(columns=['strand', 'no drugs avg', 'with drugs avg'])
    csv_file = csv_file.sort_values(by=['chr', 'start'])
    index = [f's{i}' for i in range(csv_file.shape[0])]
    csv_file = pd.merge(pd.DataFrame(index), csv_file)
    csv_file = csv_file.replace(23.0, 'X')
    csv_file = csv_file.replace(24.0, 'Y')
    csv_file.to_csv("temp.csv", index=False)
    csv_file = open("temp.csv", 'r')
    csv_file.readline()
    cn_file.write("sSNP\tchrChromosome\tPhysicalPosition\tctcfEnd\tcov\tchangeRate\n")
    for line in csv_file:
        x = line
        x = x.replace(",", '\t')
        cn_file.write(x)
    cn_file.close()


def convert_to_cn_2(file):
    """
    converting csv files to IGV files format
    :param file: the file to convert
    """
    path = os.path.dirname(file)
    name = os.path.relpath(file)
    if name.endswith(".csv"):
        name = name.replace(".csv", ".cn")
    else:
        name = path + os.sep + name + ".cn"
    csv_file = pd.read_csv(file, sep='\t')
    csv_file = csv_file.drop(csv_file.columns[0], axis=1)  # removing the index column
    csv_file = csv_file.sort_values(by=['chr', 'start', 'end'])
    csv_file = csv_file.iloc[1:]  # removing the first row with small values
    csv_file = csv_file.replace(23.0, 'X')  # replacing to X chromosome
    csv_file = csv_file.replace(24.0, 'Y')  # replacing to Y chromosome
    i = 0
    with open(name, "w") as output_file:
        output_file.write("sSNP\tchrChromosome\tPhysicalPosition\tchangeRate\n")
        for view in csv_file.iterrows():
            view = view[1]
            if view['chr'] != 'X' and view['chr'] != 'Y':
                chromosome = str(int(view['chr']))
            else:
                chromosome = view['chr']
            start = view['start']
            end = view['end']
            number = view['change']
            line = BEDGRAPH_LINE_FORMAT.format(i=i, chr_name=chromosome, start=start,
                                               number=number)
            output_file.write(line)
            i += 1


def creat_cns(dir):
    for file in os.listdir(dir):
        if file.endswith(".csv"):
            convert_to_cn_2(dir+os.path.sep+file)


def create_genes_files(up, down):
    for file in os.listdir("immortalization_result/by_window"):
        if file.endswith(".csv"):
            # f = os.path.abspath(file)
            # file = "imm_b1_filtered.csv"
            filter_data(up, "immortalization_result/by_window" + os.sep + file, "change", "immortalization_result/by_window/increase_" + file + "_{0}.csv".format(up))
            print("done1")
            filter_data(down, "immortalization_result/by_window" + os.sep + file, "change", "immortalization_result/by_window/decrease_" + file + "_{0}.csv".format(down))
            print("done2")
            check_with_change_filter([10000, 50000, 100000], 30, "immortalization_result/by_window/increase_" + file + "_{0}.csv".format(up), os.path.splitext(os.path.basename(file))[0])
            print("done increase")
            check_with_change_filter([10000, 50000, 100000], 30,  "immortalization_result/by_window/decrease_" + file + "_{0}.csv".format(down), os.path.splitext(os.path.basename(file))[0])
            print("done decrease")


def get_genes(file, window=500, flag_38=False, csc=False):
    if flag_38:
        check_with_change_filter([10000, 50000, 100000], 30, file, "imm_sagnificant_w_{0}".format(window), flag_38=True)
    else:
        if csc:
            check_with_change_filter([10000, 50000, 100000], 30, file, "csc_sgnificant_{0}", csc=True)
        else:
            check_with_change_filter([10000, 50000, 100000], 30, file, "imm_sgnificant_w_{0}".format(window))


def covnert_list_to_avg(data, col1, col2):
    control, treat = [], []
    for row in data.iterrows():
        temp = row[1][col1].split(",")
        s = [x for x in temp if x]
        c = [float(s[0][1:]), float(s[1]), float(s[2]), float(s[3][:len(s[3])-1])]
        temp = row[1][col2].split(",")
        s = [x for x in temp if x]
        t = [float(s[0][1:]), float(s[1]), float(s[2]), float(s[3][:len(s[3]) - 1])]
        control.append(sum(c) / len(c))
        treat.append(sum(t) / len(t))
    return control, treat


def get_output_gene_list(file, outputname, csc=False):
    data = pd.read_csv(file, sep="\t")
    no_name = "['no_name_found']"
    data = data[data['close_genes'] != '[]']
    data = data[data['close_genes'] != no_name]

    print(data.head())
    if csc:
        control_label = 'controls_{0}_month'.format(csc)
        after_label = 'afters{0}_month'.format(csc)
        data['metylation change'] = data[after_label] - data[control_label]
        # ref_data = pd.read_csv(DIR + os.path.sep + "after_fdr_correction_csc.tsv", sep="\t")
        # data['rejected'] = ref_data['rejected_{0}_month'.format(csc)]
    # data = data[data['rejected'] == True]
    genes_set = set()
    genes_lst = []
    for row in data.iterrows():
        at_site_lst = []
        gene = row[1]['close_genes'].split("'")
        for g in gene:
            if g != '[' and g != ']' and g != 'no_name_found' and g != ', ':
                g = g.strip('"')
                g = g.strip('" ')
                at_site_lst.append(g)
                genes_set.add(g)
        genes_lst.append(at_site_lst)
    data['close_genes'] = genes_lst
    with open(outputname, 'w') as f:
        for item in genes_set:
            f.write("%s\n" % item)

    # todo to csv ?
    return data

                # x.append(g)
                # ys.append(row[1]['metylation change'])
                # c.append(row[1]['chr'])
    # d = pd.DataFrame({'chr': c, 'genes': x, 'met': ys})
    # print(d.head())
    # d['genes'] = d['genes'].apply(lambda x: x.strip('" '))
    # d['genes'] = d['genes'].apply(lambda x: x.strip('"'))
    # x = d['genes']
    # return
    # x.to_csv("genes_at_{0}_month.txt".format(6), index=False)
    # d.to_csv("data_to_plot_gene_6_csc.tsv", sep="\t", index=False)
    # for i in range(1, 24):
    #     temp = d[d['chr'] == i]
    #     plt.scatter(x='genes', y='met', data=temp,
    #                 s=20, c=colors[i], label="chr " + str(i))
    # # i = 23
    # # temp = data[data['chr'] == "X"]
    # # plt.scatter(x='genes', y='met', data=temp,
    # #             s=20, c=colors[i], label="chr " + str(i))
    # # i = 24
    # # temp = data[data['chr'] == "Y"]
    # # plt.scatter(x='genes', y='met', data=temp,
    # #             s=20, c=colors[i], label="chr " + str(i))
    # plt.title(file)
    # plt.show()


if __name__ == '__main__':
    # file = sys.argv[1]
    # window = sys.argv[2]
    print("hi")
    # plot_sgnificant_genes("genes/genes_close_to_sites_imm_sgnificant_w_1000_filter_100000.csv", "genes/imm_1000_ref.txt")
    # plot_sgnificant_genes("genes/genes_close_to_sites_csc_sgnificant_15_filter_100000.csv", 'genes/genes_file_csc_15_ref.txt', 15)
    # plot_sgnificant_genes("genes/genes_close_to_sites_csc_sgnificant_10_filter_100000.csv", 'genes/genes_file_csc_10_ref.txt', 10)
    # plot_sgnificant_genes("genes/genes_close_to_sites_csc_sgnificant_6_filter_100000.csv", 'genes/genes_file_csc_6_ref.txt', 6)
    # get_genes(file, csc=True)
    # read_genes_data(GENES_B37)
    # read_genes_data(GENES_B38)

    # get_genes("corrected_t_test/t_test_by_site_with_population_all_w_500.csv")

# create_genes_files(0.2, -0.4)
# check_with_change_filter([50000], 30, "plass_result/filtered/increase_no_treatment_vs_with_dac.csv_0.6.csv", "increase_plass_no_treatment_vs_with_dac.csv_0.6.csv")
# check_with_change_filter([10000, 50000, 100000], 30, "plass_result/filtered/decrease_no_treatment_vs_with_dac_0.6.csv", "test")
# create_genes_files()
# creat_cns("plass_result")

# filter_data(0.001, "Compares files/after_dac_vs_after_dac_and_hdac.csv", "change", "Compares files/filtered/increase_mthylation_after_dac_vs_hdac_and_dac.csv")
# print("done1")
# filter_data(-0.1, "Compares files/after_dac_vs_after_dac_and_hdac.csv", "change", "Compares files/filtered/decrease_mthylation_after_dac_vs_hdac_and_dac.csv")
# print("done")
# check_with_change_filter([10000, 50000, 100000], 30, "plass_result/filtered/increase_no_treatment_vs_with_dac_0.6.csv", "test")
# remove_duplicate(pd.read_csv("probs_sort_and_uniq", sep='\t', header=None), 0, 1, 2)
