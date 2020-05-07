import pandas as pd
import lab
import numpy as np
import matplotlib.pyplot as plt
import os
import math
from scipy import stats
import matplotlib.cm as cm
import sys


COV = "cov"
T1 = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/adrenal_cells_try.txt"
DINFO = "Smoothed_Methylation_Level_H2_DMSO"

NODINFO = "Smoothed_Methylation_Level_H2_DAC"

PREVENT_INFO = "Smoothed_Methylation_Level_H2_SB939"

DAC_INFO = "Smoothed_Methylation_Level_H2_DAC_plus_SB939"

CONTROL = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz"

AFTER_TREATMENT = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz"

DAC_AND_HDAC = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/GSM2150387_H2_DAC_plus_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz"

HDAC_PREVENT = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/GSM2150389_H2_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz"

P3_control = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz"

P1_after_treatment = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz"

PLASS3 = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/ENCFF032DEW.bed.gz"

PLASS2 = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/ENCFF543VGD.bed.gz"

PLASS1 = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/files/Plass/ENCFF401ONY.bed.gz"

B1_ACTIVE = "/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/GSM1202804_B1_activation.csv.gz"

B2_ACTIVE = "/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/GSM1202805_B2_activation.csv.gz"

B3_ACTIVE = "/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/GSM1202806_B3_activation.csv.gz"

B1_TRANSFOR = "/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/GSM1202810_B1_transformed.csv.gz"

B2_TRANSFOR = "/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/GSM1202811_B2_transformed.csv.gz"

B3_TRANSFOR = "/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/GSM1202812_B3_transformed.csv.gz"

CHIP_B_CELLS_ACTIVE = "/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/ENCFF449NOT.bed"

TF_TRANS = "/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/ENCFF833FTF.bed"


def closest_to_peak(lst, peak, start):
    """
    find in list of probes i
    :param lst: the chromosome biding sits as we filtered from the data
    :param peak: the peak location at the site
    :param start: the biding site's start
    :return: the closest to peak value of methylation
    """
    first = 0
    last = len(lst) - 1
    val = start + peak
    mid = 0
    while first <= last:
        mid = (first + last) // 2
        if val == lst[mid][0]:
            break
        else:
            if val < lst[mid][0]:
                last = mid - 1
            else:
                first = mid + 1

    if val == lst[mid][0]:
        cov = calc_cov(mid, lst)
        return lst[mid][1], cov

    elif mid-1 < 0:
        cov = calc_cov(mid+1, lst)
        return lst[mid+1][1], cov

    elif mid+1 > last:
        cov = calc_cov(mid-1, lst)
        return lst[mid-1][1], cov

    cov = calc_cov(mid-1, lst)
    return max(lst[mid-1][1], lst[mid+1][1]), cov


def calc_cov(mid, chr_lst):
    # new_lst = [item[2] for item in chr_lst if chr_lst[mid][0] - 50 <= item[0] <= chr_lst[mid][0] + 50]
    l = np.array(chr_lst)
    l = l[(l[:, 0] >= chr_lst[mid][0] - 250) & (l[:, 0] <= chr_lst[mid][0] + 250)]
    if len(l) != 0:
        cov = np.mean(l[:, 2])
        # cov2 = sum(new_lst) / len(new_lst)
        # print("current cov: " + str(cov))
        return cov
    # print("current cov: 0")
    return 0


def read_gz_file(file1, file2, sep):
    """
    reading the file's and filter it if needed
    :param file1: the first file
    :param file2: the second fil
    e
    :return: the files
    """
    before = pd.read_csv(file1, sep=sep, low_memory=False, compression="gzip")
    after = pd.read_csv(file2, sep=sep, low_memory=False, compression="gzip")
    return before, after


def smooth_parse(data, level, chr_name, start):
    """
    parsing the smoothing file
    :param data: the data's file
    :param level: the file's name
    :return: an array with the data parsed from the file
    """
    chrom = [[] for i in range(24)]
    for index, loci, level, cov in zip(data[chr_name], data[start], data[level], data["Cov"]):
        if index == 'X' or index == "chrX":
            chrom[22].append([loci, level, cov])
        elif index == 'Y' or index == "chrY":
            chrom[23].append([loci, level, cov])
        else:
            index = int(''.join([s for s in list(index) if s.isdigit()]))
            chrom[index - 1].append([loci, level, cov])
    #     print(index)
    # print("before sorting")
    for i in chrom:
        i.sort()
    return chrom


# def find_position(lst, location):
#     """
#     binary search to find the CTCF biding site's start and end
#     :param lst: the chromosome biding sits as we filtered from the data
#     :param location: the chromosome's start position
#     :return: the start or end position
#     """
#     first = 0
#     last = len(lst) - 1
#     while first <= last:
#         mid = (first + last) // 2
#         if location == lst[mid][0]:
#             return mid
#         else:
#             if location < lst[mid][0]:
#                 last = mid - 1
#             else:
#                 first = mid + 1
#     return mid


def plot_cov(dir, graph_name, data, file=None):
    if file:
        data = pd.read_csv(file, sep='\t', low_memory=False, skiprows=[1])
    cov_data = data[COV]
    cov_dict = {}
    g = lambda x: int(x) + 0.5 if math.ceil(x) - x > 0.5 else int(x) - 0.5
    for i in cov_data:
        i = g(i)
        if i in cov_dict:
            cov_dict[i] += 1
        else:
            cov_dict[i] = 1
    # create coverage vs change graph
    xy = sorted(list(cov_dict.items()))
    x = list(map(lambda tup: tup[0], xy))
    y = list(map(lambda tup: tup[1], xy))
    fig, ax = plt.subplots()
    ax.plot(x, y)
    g = graph_name + "_coverage_vs_num_of_site"
    plt.title(g)
    plt.xlabel('coverage')
    plt.ylabel('number of sites')
    if dir != "":
        plt.savefig(dir + os.path.sep + g)
    else:
        plt.savefig(g)
    plt.show()
    fig, ax = plt.subplots()

    # create small coverage vs change graph
    sxy = list(filter(lambda tup: tup[0] < 20, xy))
    sx = list(map(lambda tup: tup[0], sxy))
    sy = list(map(lambda tup: tup[1], sxy))
    ax.scatter(sx, sy)
    plt.xlim(0, 20)
    g = graph_name + "small_coverage_vs_num_of_site"
    plt.title(g)
    plt.xlabel('coverage')
    plt.ylabel('number of sites')
    if dir != "":
        plt.savefig(dir + os.path.sep + graph_name + g)
    else:
        plt.savefig(graph_name + g)
    plt.show()

    # # create scatter graph of change vs coverage
    # change_data = data["change"]
    # colors = np.random.rand(change_data.size)
    # plt.scatter(cov_data, change_data, c=colors)
    # plt.xlim(0, 40)
    # g = graph_name + "_coverage_vs_change"
    # plt.title(g)
    # plt.xlabel('coverage')
    # plt.ylabel('change')
    # if dir != "":
    #     plt.savefig(dir + os.path.sep + graph_name + g)
    # else:
    #     plt.savefig(graph_name + g)
    # plt.show()

    # # create scatter graph of change vs coverage
    # change_data = data["change"]
    # print(1)
    # colors = np.random.rand(change_data.size)
    # xy = np.vstack([cov_data, change_data])
    # print(2)
    #
    # z = gaussian_kde(xy)(xy)
    # print(3)
    # plt.scatter(cov_data, change_data, c=z)
    # print(4)
    # plt.show()
    # plt.xlim(0, 40)
    # print(5)
    # g = graph_name + "_coverage_vs_change_try_gauss"
    # plt.title(g)
    # plt.xlabel('coverage')
    # plt.ylabel('change')
    # if dir != "":
    #     plt.savefig(dir + os.path.sep + graph_name + g)
    # else:
    #     plt.savefig(graph_name + g)
    # plt.show()

    y = data["change"]
    x = cov_data
    #histogram definition
    bins = [1000, 1000]  # number of bins

    # histogram the data
    hh, locx, locy = np.histogram2d(x, y, bins=bins)

    # Sort the points by density, so that the densest points are plotted last
    z = np.array([hh[np.argmax(a <= locx[1:]), np.argmax(b <= locy[1:])] for a, b in zip(x, y)])
    idx = z.argsort()
    x2, y2, z2 = x[idx], y[idx], z[idx]

    plt.figure(1, figsize=(8, 8)).clf()
    plt.xlim(0, 20)

    plt.scatter(x2, y2, c=z2, cmap='jet', marker='.')
    g = graph_name + "_coverage_vs_change"
    plt.title(g)
    plt.xlabel('coverage')
    plt.ylabel('change')
    if dir != "":
        plt.savefig(dir + os.path.sep + graph_name + g)
    else:
        plt.savefig(graph_name + g)
    plt.show()

    # create histogram of change distribution by different coverage value
    filter_range = [0, 2, 4, 6, 8, 10, 12, 14]
    for f in range(len(filter_range) - 1):
        d = data[data[COV] > filter_range[f]]
        d = d[data[COV] < filter_range[f+1]]
        cd = d["change"]
        plt.hist(cd, histtype='step', label="{0} - {1}".format(filter_range[f], filter_range[f+1]), density=True, stacked=True)
    plt.legend(title="coverage area")
    g = graph_name + "change_distribution"
    plt.title(g)
    plt.xlabel("change rate")
    plt.ylabel("num of sites")
    if dir != "":
        plt.savefig(dir + os.path.sep + graph_name + g)
    else:
        plt.savefig(graph_name + g)
    plt.show()
    #todo: check how the get the probability of this
    return


def make_box_plot(file, title, graph_name):
    """
    creating a box plot graph that shows the changes by chromosomes
    :param file: the file with the changing rate to show as graph
    """
    lst = pd.read_csv(file, sep="\t")  # not always is sep by tubs
    lst = lst.drop(lst.index[[0]])
    # lst = lst.drop(lst.index[[0]])

    # lst = pd.DataFrame(data=lst, columns=head)
    chroms = []
    for i in range(1, 25):
        this_chrom = lst[lst["chr"] == i]
        this_chrom = this_chrom.drop(columns=["chr", "start", "end", "no drugs avg", "with drugs avg", "cov"])
        this_chrom = this_chrom.drop(this_chrom.columns[0], axis=1)
        if not this_chrom.empty:
            chroms.append(np.array(this_chrom))
    plt.boxplot(chroms)
    plt.title(title)
    plt.xlabel("chromosomes, 23 = X, 24 = Y")
    plt.ylabel("change level")
    plt.grid()
    # plt.setp(plt, xticks=[c+1 for c in range(23)])
    plt.legend()
    plt.savefig(graph_name)
    plt.show()


def search(before, after, chip_data, name="in_progress.csv"):
    """
    searching if prob is at the area of an peak +- the buffer.
    :param before: the data before treatment
    :param after: the data after treatment
    :param chip_data: the data from the chip array experience
    :param name: the output file name
    :return: the ratio between sum of probes are is the buffer to the total amount of probes
    """
    head = ["chr", "start", "end", "no drugs avg", "with drugs avg", "cov"]
    strand_col = chip_data.drop(columns=["chromStart", "chromEnd", "chrom", "peak"])
    strand_col = np.vstack([strand_col, "."])
    # cov_col = chip_data.drop(columns=["chromStart", "chromEnd", "chrom", "peak"])
    print(strand_col.shape[0], strand_col.shape[1])
    lst = np.empty((1, len(head)))
    for start, end, chr, peak in zip(chip_data["chromStart"], chip_data["chromEnd"],
                                     chip_data["chrom"],  chip_data["peak"]):
        if chr == 'chrX':
            chrom = 23
        elif chr == 'chrY':
            chrom = 24
        else:
            chrom = int(chr[3:])
        if before[chrom-1] == []:
            line = np.array([chrom, start, end, 0, 0, -1])
            lst = np.vstack([lst, line])
            continue
        before_avg, cov1 = closest_to_peak(before[chrom - 1], peak, start)
        after_avg, cov2 = closest_to_peak(after[chrom - 1], peak, start)
        cov = min(cov1, cov2)
        line = np.array([chrom, start, end, before_avg, after_avg, cov])
        lst = np.vstack([lst, line])
        print("still alive at chr " + str(chr) + ", start site at " + str(start))

    #  adding column with the difference between the average with out / with the drag
    change = np.subtract(lst[:, 4], lst[:, 3])[np.newaxis]
    change = change.T
    lst = np.hstack([lst, strand_col])
    head.append("strand")
    lst = np.hstack([lst, change])
    head.append("change")
    lst = lst[lst[:, 0].argsort()]  # sorting by chromosomes
    #  writing the result
    pd.DataFrame(data=lst, columns=head).to_csv(name, sep="\t")


def main_plass(i, j):
    """
    main function to process plass files
    """
    plass = [PLASS1, PLASS2, PLASS3]
    data1 = []
    for p in plass:
        data1.append(lab.read_chip_file(p, 100))
    data = pd.concat(data1)
    data = data.drop_duplicates()
    print("done append data")
    plass_files = [(CONTROL, DINFO, "no_treatment"), (AFTER_TREATMENT, NODINFO, "with_dac"),
                   (DAC_AND_HDAC, DAC_INFO, "with_dac_and_hdac"), (HDAC_PREVENT, PREVENT_INFO, "with_hdac")]
    # for i in range(len(plass_files) - 1):
    #     for j in range(i + 1, len(plass_files)):
            # filters = [0.1, 0.3, 0.5]
            # for filter in filters:
            # print("results for filter: " + str(filter))
    ndrg, drg = read_gz_file(plass_files[i][0], plass_files[j][0], '\t')
    print("done reading the files")
    nodrag = smooth_parse(ndrg, plass_files[i][1], "Chromosome", "Start")
    drag = smooth_parse(drg, plass_files[j][1], "Chromosome", "Start")
    print("done parsing")
    name = "{0}_vs_{1}".format(plass_files[i][2], plass_files[j][2])
    search(nodrag, drag, data, name + ".csv")
    # plot_cov("", "{0}_covarge_histogram".format(name), data)
    # make_box_plot(name + ".csv", "mthylation level's changes at " + name, name)
    print("F I N I S H !")


def main_imm(i):
    """
    main function to process immortalization files
    """
    # immortalization = [CHIP_B_CELLS_ACTIVE, TF_TRANS]
    b_active = lab.read_chip_file(CHIP_B_CELLS_ACTIVE, 100)
    tf_trans = lab.read_chip_file(TF_TRANS, 100)
    data = pd.concat([b_active, tf_trans])
    data = data.drop_duplicates()
    print("done append data")
    imm_files = [(B1_ACTIVE, B1_TRANSFOR), (B2_ACTIVE, B2_TRANSFOR), (B3_ACTIVE, B3_TRANSFOR)]
    # imm_files = [(B1_ACTIVE, B1_TRANSFOR)]
    imm_active, imm_trans = [], []
    # for imm in imm_files[i-1]:
    imm = imm_files[i-1]
    a, b = read_gz_file(imm[0], imm[1], ',')
    imm_active.append(a)
    imm_trans.append(b)
    imm_active = pd.concat(imm_active)
    imm_trans = pd.concat(imm_trans)
    imm_active = imm_active.drop_duplicates()
    imm_trans = imm_trans.drop_duplicates()
    print("done reading the files")
    active = smooth_parse(imm_active, "smoothSmall", "chr", "pos")  # only small ! there is also large
    transformed = smooth_parse(imm_trans, "smoothSmall", "chr", "pos")
    print("done parsing")
    # output = "mthylation level's changes at b immortalization cells"
    # search(active, transformed, b_active, name.format(i+1, "active") + ".csv")  # active file
    search(active, transformed, data, "imm_result_b{0}_w_500.csv".format(i))  # transformed file
    print("D O N E")
    # make_box_plot("imm_result_b1.csv", output, "imm_result_graph_b1")
    # for i in range(len(imm_files)):
    #     act, trn = read_gz_file(imm_files[i][0], imm_files[i][1], ',')
    #     print("done reading the files")
    #     active = smooth_parse(act, "smoothSmall", "chr", "pos")  # only small ! there is also large
    #     transformed = smooth_parse(trn, "smoothSmall", "chr", "pos")
    #     print("done parsing")
    #     output = "mthylation level's changes at b{0} immortalization cells"
    #     name = "b{0}"
    #     # search(active, transformed, b_active, name.format(i+1, "active") + ".csv")  # active file
    #     search(active, transformed, data, name.format(i+1) + ".csv")  # transformed file
    #     make_box_plot(name.format(i+1) + ".csv", output.format(i+1), name.format(i+1) + "_graph")
        # make_box_plot(name.format(i+1, "trans") + ".csv", output.format(i+1, "transformed"), name.format(i+1, "transformed") + "_graph")
        # print("F I N I S H !")


def no_use():
    return

    # procecess the files
    # finding the start index of the methylation
    # startin = find_position(nodrags[chr], start)
    # endin = find_position(nodrags[chr], end)
    # nodrag_met = endin - startin
    # nodragcount += nodrag_met
    #  finding the end index of the methylation
    # startin = find_position(drags[chr], start)
    # endin = find_position(drags[chr], end)
    # drag_met = endin - startin
    # dragcount += drag_met
    #  calculation of the average
    # creating graphs
    # num_chr = 1
    # for chr in chroms:
        # x = chr['start']
        # start = int(min(x))
        # end = int(max(x))
        # no_d = chr['no drugs avg']
        # with_d = chr['with drugs avg']
        # # plt.xticks(range(start, end))
        # x_asix = np.linspace(start, end, chr.shape[0])
        # plt.plot(x_asix, no_d, label="No drug")
        # plt.plot(x_asix, with_d, label="With drug")
        # plt.title("change in matylation level at chromosome " + str(num_chr) + " after treatment")
        # plt.grid()
        # plt.legend()
        # plt.savefig("chr " + str(num_chr))
        # num_chr += 1
        # plt.clf()

    # iter in sarceh
# for chromStart, chromEnd, chrom, peak in zip(chip_data[i]["chromStart"],
#                             chip_data[i]["chromEnd"],
#                             chip_data[i]["chrom"],
#                                  chip_data[i]["peak"]):


def cut_by_filter(file):
    data = pd.read_csv(file, sep="\t")
    data = data[data["cov"] >= 5]
    data.to_csv("{0}_filtered.csv".format(file), sep="\t", index=False)


def remove_empty_sites(dir):
    for file in os.listdir(dir):
        if file.startswith("genes_"):
            df = pd.read_csv(dir + os.sep + file, sep="\t")
            df = df[df.astype(str)['close_genes'] != '[]']
            df.to_csv(dir + os.sep + file, index=False, sep="\t")


def sort_by_len(dir):
    for file in os.listdir(dir):
        if file.startswith("sites_"):
            df = pd.read_csv(dir + os.sep + file, sep="\t")
            count = [i.count(")") for i in df["close_sites"]]
            df = df.assign(count=count)
            df = df.sort_values(["count", "chr"])
            df = df.drop(columns=['count'])
            df.to_csv(dir + os.sep + file, index=False, sep="\t")


def get_index(index):
    if index == "chrX":
        return 22
    elif index == "chrY":
        return 23
    else:
        return int(''.join([s for s in list(index) if s.isdigit()])) - 1


def biding_vs_methylation(score=0):
    f = open("biding_vs_methylation.txt", 'w')
    f.write("starting")
    print("starting")
    range_dict = {(0, 0.2): 0, (0.2, 0.4): 0, (0.4, 0.6): 0, (0.6, 0.8): 0, (0.8, 1): 0}
    bed = pd.read_csv("/vol/sci/bio/data/yotam.drier/CTCF_and_DNAme/immortalization/ENCFF449NOT.bed", sep='\t', comment='t', header=None)
    header = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
              "signalVal", "pVal", "qVal", "peak"]
    bed.columns = header[:len(bed.columns)]
    #  optional : check by score
    # bed = bed[bed['score'] >= score]
    f.write("process data")
    print("process data")
    bed = bed.drop(columns=["name", "score", "strand", "signalVal", "pVal", "qVal", "peak"])
    chip = pd.read_csv(B1_ACTIVE, sep=",", low_memory=False, compression="gzip")
    chip = chip.drop(columns=["M", "Cov", "smoothLarge"])
    chip = chip.sort_values(["chr", "pos"])
    bed = bed.sort_values(["chrom", "chromStart"])

    chrom = []
    for i in range(1, 23):
        temp = bed[bed["chrom"] == "chr{0}".format(i)]
        chrom.append(temp)
    temp = bed[bed["chrom"] == "chrX"]
    chrom.append(temp)
    temp = bed[bed["chrom"] == "chrY"]
    chrom.append(temp)
    f.write("done process data")
    print("done process data")
    f.write("start checking")
    print("start checking")
    for row_data in chip.iterrows():
        inx = get_index(row_data[1]["chr"])
        for site in chrom[inx].iterrows():
            if site[1]["chromStart"] <= row_data[1]["pos"] <= site[1]["chromEnd"]:
                for key in range_dict:
                    if key[0] <= row_data[1]["smoothSmall"] < key[1]:
                        range_dict[key] += 1
                        break
    f.write("done checking")
    print("done checking")
    total = chip.shape()[0]
    total_range = {k[0]: k[1] / total for k in range_dict.items()}
    f.write("result:")
    f.write("range dict: " + str(range_dict))
    f.write("total range: " + str(total_range))
    f.close()
    print("result:")
    print(range_dict)
    print(total_range)
    print("yay")


def create_avg_file(b1, b2, b3):
    b1['source'] = 1
    b2['source'] = 2
    b3['source'] = 3
    all = pd.concat([b1, b2], ignore_index=True)
    all = pd.concat([all, b3], ignore_index=True)
    all = all.sort_values(by=['start', 'end', 'source'])
    all = all.reset_index(drop=True)
    all = all.drop(columns=['Unnamed: 0', 'strand'])
    inx = all.columns
    new_data = pd.DataFrame(columns=inx).T
    i = 0
    print("#starting")
    while i < len(all.index) - 2:
        if all['source'][i] != all['source'][i+1]:
            sum, count = -1, -1
            if all['start'][i] - 50 <= all['start'][i+1] <= all['start'][i] + 50 and all['end'][i] - 50 <= all['end'][i+1] <= all['end'][i] + 50:
                sum = all['change'][i+1] + all['change'][i]
                count = 2
                # print("we have the same !")
            if all['source'][i] != all['source'][i+2] and all['source'][i+1] != all['source'][i+2]:
                if all['start'][i] - 50 <= all['start'][i+2] <= all['start'][i] + 50 and all['end'][i] - 50 <= all['end'][i+2] <= all['end'][i] + 50:
                    sum += all['change'][i+2]
                    count += 1
            if sum != -1 and count != -1:
                avg = sum / count
                temp = pd.DataFrame(all.iloc[i]).T
                temp['change'] = avg
                # print("new line, before: i = {0}, data = \n{1}".format(i, new_data))
                new_data = pd.concat([temp, new_data], ignore_index=True, sort=False)
                # print("new line, after: i = {0}, data = \n{1}".format(i, new_data))
                i += count
                continue

        #add the old line
        # print("old line, before: i = {0}, data = \n{1}".format(i, new_data))
        new_data = pd.concat([pd.DataFrame(all.iloc[i]).T, new_data], ignore_index=True, sort=False)
        # print("old line, after: i = {0}, data = \n{1}".format(i, new_data))
        i += 1
    new_data = new_data.dropna()
    new_data = new_data.drop(columns=['source'])
    new_data.to_csv("imm_combine_with_avg_w500.csv", sep='\t')
    print("done !")


def read_methylation_file(file):
    print("a")
    f = pd.read_fwf(file)
    print("b")
    f.to_csv("try.csv")
    print("c")
    x = 0

def plot_change(b1, b2, b3):
    x1, y1 = get_uniq_rate(b1)
    x2, y2 = get_uniq_rate(b2)
    x3, y3 = get_uniq_rate(b3)
    # x1 = b1['start']
    # x2 = b2['start']
    # x3 = b3['start']
    # y1 = b1['change']
    # y2 = b2['change']
    # y3 = b3['change']
    # round = np.around(np.array(y1), 1)
    # unique_elements, counts_elements = np.unique(round, return_counts=True)
    plt.scatter(x1, y1, label='b1')
    plt.scatter(x2, y2, label='b2')
    plt.scatter(x3, y3, label='b3')
    # plt.hist(y1, alpha=0.5, label='b1')
    # plt.hist(y2, alpha=0.5, label='b2')
    # plt.hist(y3, alpha=0.5, label='b3')
    plt.xlabel("methylation rate")
    plt.ylabel("Amount of performances")
    plt.title("the diffrence of methylation rates by repeats, window :1000")
    plt.savefig("diffrence_in_repeats_scatter_2_points_w_1000")
    plt.legend()
    plt.show()


def get_uniq_rate(b):
    round = np.around(np.array(b['change']), 2)
    return np.unique(round, return_counts=True)


def t_test(b1, b2, b3):
    b1['source'] = 1
    b2['source'] = 2
    b3['source'] = 3
    all = pd.concat([b1, b2], ignore_index=True)
    all = pd.concat([all, b3], ignore_index=True)
    all = all.sort_values(by=['chr', 'start', 'end', 'source'])
    all = all.reset_index(drop=True)
    all = all.drop(columns=['Unnamed: 0', 'strand'])
    i = 0
    sites = np.empty((1, 7))
    print("#starting")
    while i < len(all.index) - 2:
        if all['source'][i] != all['source'][i+1]:
            if all['start'][i] - 50 <= all['start'][i+1] <= all['start'][i] + 50 and all['end'][i] - 50 <= all['end'][i+1] <= all['end'][i] + 50:
                if all['source'][i] != all['source'][i+2] and all['source'][i+1] != all['source'][i+2]:
                    if all['start'][i] - 50 <= all['start'][i+2] <= all['start'][i] + 50 and all['end'][i] - 50 <= all['end'][i+2] <= all['end'][i] + 50:
                        before = all['no drugs avg'][i:i+3]
                        after = all['with drugs avg'][i:i+3]
                        t_test = stats.ttest_ind(before, after, equal_var=False)
                        # if t_test.pvalue <= 0.05:
                        sites = np.vstack([sites, np.array([all['chr'][i], all['start'][i], all['end'][i], t_test.pvalue,
                                                                all['change'][i:i+3].mean(), np.array(before), np.array(after)])])
                        # sites = np.vstack([sites, np.array([int(all['chr'][i]), int(all['start'][i]), t_test.pvalue])])
                        # print("we have the same !")
                        i += 3
                        continue
        i += 1

    pd.DataFrame(sites, columns=['chr', 'start', 'end', 'p value','metylation change', 'control', 'after treatment']).to_csv("t_test_by_site_with_population_all_w_1000.csv", sep="\t")
    # pd.DataFrame(sites, columns=['chr', 'start', 'p value']).to_csv("t_test_by_site_with_population_all_w_500.csv", sep="\t")
    print(sites.shape)
    print("yay we finished")


def create_plot(t_file):
    x,ys,c = [],[],[]
    t = pd.read_csv(t_file, sep='\t')
    for row in t.iterrows():
        if row[1]['close_genes'] != '[]':
            gene = row[1]['close_genes'].split("'")
            for g in gene:
                if g != '[' and g != ']':
                    x.append(g)
                    ys.append(row[1]['metylation change'])
                    c.append(row[1]['chr'])
    d = pd.DataFrame({'chr': c, 'genes': x, 'met': ys})
    # sns.scatterplot(data=d, x='genes', y='met', hue='chr')
    print(d.head())
    # plt.scatter(x, y,edgecolors=c)
    # plt.xticks(rotation=90)
    # plt.show()


def filter_final_data(dir):
    for file in os.listdir(dir):
        if file.endswith(".csv"):
            data = pd.read_csv(dir + os.sep + file, sep="\t")
            # print("file : " + file)
            # print("max corrected p_value: " + str(max(data['correct_p'])))
            data['abs_met'] = data['metylation change'].abs()
            data = data[data['abs_met'] >= 0.2]
            data = data.drop(columns=['abs_met'])
            # data = data[data['metylation change'] <= 0.2]
            # data = data[data['metylation change'] >= -0.2]
            # data = data[data['correct_p'] <= 0.1]
            data.to_csv("genes/corrected/filtered" + os.sep + file, sep="\t", index=False)




if __name__ == '__main__':
    i = int(sys.argv[1])
    j = int(sys.argv[2])
    print(i,j)
    main_plass(i, j)
    print("done")
