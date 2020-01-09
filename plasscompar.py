import pandas as pd
import lab
import numpy as np
import matplotlib.pyplot as plt

DINFO = "Smoothed_Methylation_Level_H2_DMSO"

NODINFO = "Smoothed_Methylation_Level_H2_DAC"

PREVENT_INFO = "Smoothed_Methylation_Level_H2_SB939"

DAC_INFO = "Smoothed_Methylation_Level_H2_DAC_plus_SB939"

CONTROL = "PLASS/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

AFTER_TREATMENT = "PLASS/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

DAC_AND_HDAC = "PLASS/GSM2150387_H2_DAC_plus_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

HDAC_PREVENT = "PLASS/GSM2150389_H2_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

P3_control = "PLASS/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

P1_after_treatment = "PLASS/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

PLASS3 = "PLASS/ENCFF032DEW.bed"

PLASS2 = "PLASS/ENCFF543VGD.bed"

PLASS1 = "PLASS/ENCFF401ONY.bed"

B1_ACTIVE = "immortalization/B1_activation.csv.gz"

B2_ACTIVE = "immortalization/B2_activation.csv.gz"

B3_ACTIVE = "immortalization/B3_activation.csv.gz"

B1_TRANSFOR = "immortalization/B1_transformed.csv.gz"

B2_TRANSFOR = "immortalization/B2_transformed.csv.gz"

B3_TRANSFOR = "immortalization/B3_transformed.csv.gz"

CHIP_B_CELLS_ACTIVE = "immortalization/ENCFF449NOT.bed"

TF_TRANS = "imm/ENCFF833FTF.bed.gz"


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
            return lst[mid][1], lst[mid][2]
        else:
            if val < lst[mid][0]:
                last = mid - 1
            else:
                first = mid + 1
    if mid-1 < 0:
        return lst[mid+1][1], lst[mid+1][2]
    elif mid+1 > last:
        return lst[mid-1][1], lst[mid-1][2]

    return max(lst[mid-1][1], lst[mid+1][1]), lst[mid-1][2]


def read_gz_file(file1, file2, sep):
    """
    reading the file's and filter it if needed
    :param file1: the first file
    :param file2: the second file
    :return: the files
    """
    before = pd.read_csv(file1, sep=sep, low_memory=False)
    after = pd.read_csv(file2, sep=sep, low_memory=False)  # todo : change the "Chromosome" col to chr
    # before = before.rename(columns={'Chromosome': 'chr'}, axis='columns')
    # after = after.rename(columns={'Chromosome': 'chr'}, axis='columns')
    # nodrag = nodrag[nodrag[DINFO] >= filter]
    # drag = drag[drag[NODINFO] >= filter]
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
    for i in chrom:
        i.sort()
    return chrom


def find_position(lst, location):
    """
    binary search to find the CTCF biding site's start and end
    :param lst: the chromosome biding sits as we filtered from the data
    :param location: the chromosome's start position
    :return: the start or end position
    """
    first = 0
    last = len(lst) - 1
    while first <= last:
        mid = (first + last) // 2
        if location == lst[mid][0]:
            return mid
        else:
            if location < lst[mid][0]:
                last = mid - 1
            else:
                first = mid + 1
    return mid


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
    for i in range(1, 25):  # todo 25 ?
        this_chrom = lst[lst["chr"] == i]
        this_chrom = this_chrom.drop(columns=["chr", "start", "end", "strand", "no drugs avg", "with drugs avg", "cov"])
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
        befor_avg, cov1 = closest_to_peak(before[chrom - 1], peak, start)
        after_avg, cov2 = closest_to_peak(after[chrom - 1], peak, start)
        cov = min(cov1, cov2)
        line = np.array([chrom, start, end, befor_avg, after_avg, cov])
        lst = np.vstack([lst, line])
        print("still alive " + str(chr))

    x = 2
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


def main_plass():
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
    for i in range(len(plass_files) - 1):
        for j in range(i + 1, len(plass_files)):
            # filters = [0.1, 0.3, 0.5]
            # for filter in filters:
            # print("results for filter: " + str(filter))
            ndrg, drg = read_gz_file(plass_files[i][0],plass_files[j][0], '\t')
            print("done reading the files")
            nodrag = smooth_parse(ndrg, plass_files[i][1], "Chromosome", "Start")
            drag = smooth_parse(drg, plass_files[j][1], "Chromosome", "Start")
            print("done parsing")
            name = "{0}_vs_{1}".format(plass_files[i][2], plass_files[j][2])
            search(nodrag, drag, data, name + ".csv")
            make_box_plot(name + ".csv", "mthylation level's changes at " + name, name)
            print("F I N I S H !")


def main_imm():
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
    imm_active, imm_trans = [], []
    for imm in imm_files:
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
    output = "mthylation level's changes at b immortalization cells"
    name = "b{0}"
    # search(active, transformed, b_active, name.format(i+1, "active") + ".csv")  # active file
    search(active, transformed, data, "b_cells.csv")  # transformed file
    make_box_plot("b_cells.csv", output, "b_cells_graph")
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


# main_plass()
main_imm()
# make_box_plot("immortalization/b3_active.csv", "3", "active_test_original")
# make_box_plot("b3_active_test2.csv", "3", "active_test2_no_cov")
# make_box_plot("b3_trans_test2.csv", "3", "trans_test2_no_cov")

# make_box_plot("immortalization/b1_trans.csv", "3", "transformed")

# make_box_plot("Compares files/no_treatment_vs_hdac_only.csv", "1", "a")
# make_box_plot("test.csv", "1", "test1")
# t = "mthylation level's changes at b{0} {1} immortalization cells".format(number, state)
# g = "immortalization/b{0}_{1}_graph".format(number, state)



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
