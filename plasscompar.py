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

B1_ACTIVE = "immortalization/B1_activation.csv"

B2_ACTIVE = "immortalization/B2_activation.csv"

B3_ACTIVE = "immortalization/B3_activation.csv"

B1_TRANSFOR = "immortalization/B1_transformed.csv"

B2_TRANSFOR = "immortalization/B2_transformed.csv"

B3_TRANSFOR = "immortalization/B3_transformed.csv"

CHIP_B_CELLS_ACTIVE = "immortalization/ENCFF449NOT.bed"

TF_TRANS = "immortalization/ENCFF833FTF.bed"


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
            return lst[mid][1]
        else:
            if val < lst[mid][0]:
                last = mid - 1
            else:
                first = mid + 1
    if mid-1 < 0:
        return lst[mid+1][1]
    elif mid+1 > last:
        return lst[mid-1][1]

    return max(lst[mid-1][1], lst[mid+1][1])


def read_gz_file(file1, file2, sep):
    """
    reading the file's and filter it if needed
    :param file1: the first file
    :param file2: the second file
    :return: the files
    """
    before = pd.read_csv(file1, sep=sep, low_memory=False)
    after = pd.read_csv(file2, sep=sep, low_memory=False) # todo : change the "Chromosome" col to chr
    # before = before.rename(columns={'Chromosome': 'chr'}, axis='columns')
    # after = after.rename(columns={'Chromosome': 'chr'}, axis='columns')
    # nodrag = nodrag[nodrag[DINFO] >= filter]
    # drag = drag[drag[NODINFO] >= filter]
    return before, after


def smooth_parse(data, name):
    """
    parsing the smoothing file
    :param data: the data's file
    :param name: the file's name
    :return: an array with the data parsed from the file
    """
    chrom = [[] for i in range(24)]
    for index, loci, level in zip(data["Chromosome"], data["Start"], data[name]):
        if index == 'X' or index == "chrX":
            chrom[22].append([loci, level])
        elif index == 'Y' or index == "chrY":
            chrom[23].append([loci, level])
        else:
            index = int(''.join([s for s in list(index) if s.isdigit()]))
            chrom[index - 1].append([loci, level])
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


def make_box_plot(file):
    """
    creating a box plot graph that shows the changes by chromosomes
    :param file: the file with the changing rate to show as graph
    """
    lst = pd.read_csv(file)
    lst = lst.drop(lst.index[[0]])
    # lst = pd.DataFrame(data=lst, columns=head)
    chroms = []
    for i in range(1, 24):
        this_chrom = lst[lst['chr'] == i]
        this_chrom = this_chrom.drop(columns=["chr", "start", "end", "no drugs avg", "with drugs avg"])
        this_chrom = this_chrom.drop(this_chrom.columns[0], axis=1)
        chroms.append(np.array(this_chrom) * -1)
    plt.boxplot(chroms)
    plt.title("mthylation level's changes after treatment by DAC vs by DAC and HDAC")
    plt.xlabel("chromosomes, 22 = X, 23 = Y")
    plt.ylabel("change level")
    plt.grid()
    # plt.setp(plt, xticks=[c+1 for c in range(23)])
    plt.legend()
    plt.savefig("check")
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
    head = ["chr", "start", "end", "no drugs avg", "with drugs avg"]
    strand_col = chip_data.drop(columns=["chromStart", "chromEnd", "chrom", "peak"])
    strand_col = np.vstack([strand_col, "."])
    print(strand_col.shape[0], strand_col.shape[1])
    lst = np.empty((1, len(head)))
    for start, end, chr, peak in zip(chip_data["chromStart"], chip_data["chromEnd"],
                                     chip_data["chrom"],  chip_data["peak"]):
        if chr == 'chrX':
            chrom = 22
        elif chr == 'chrY':
            chrom = 23
        else:
            chrom = int(chr[3:])
        befor_avg = closest_to_peak(before[chrom - 1], peak, start)
        after_avg = closest_to_peak(after[chrom - 1], peak, start)
        line = np.array([chrom, start, end, befor_avg, after_avg])
        lst = np.vstack([lst, line])
        print("still alive " + str(chr))

    #  adding column with the difference between the average with out / with the drag
    change = np.subtract(lst[:, 3], lst[:, 4])[np.newaxis]
    change = change.T
    lst = np.hstack([lst, strand_col])
    head.append("strand")
    lst = np.hstack([lst, change])
    head.append("change")
    lst = lst[lst[:, 0].argsort()]  # todo change to 5
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
    # filters = [0.1, 0.3, 0.5]
    # for filter in filters:
    # print("results for filter: " + str(filter))
    ndrg, drg = read_gz_file(AFTER_TREATMENT, DAC_AND_HDAC, '\t')
    print("done reading the files")
    nodrag = smooth_parse(ndrg, NODINFO)
    drag = smooth_parse(drg, DAC_INFO)
    print("done parsing")
    search(nodrag, drag, data, "check_plass_with_strand.csv")
    print("F I N I S H !")


def main_imm():
    """
    main function to process immortalization files
    """
    # imm = [CHIP_B_CELLS_ACTIVE, TF_TRANS]
    b_active = lab.read_chip_file(CHIP_B_CELLS_ACTIVE, 100)
    tf_trans = lab.read_chip_file(TF_TRANS, 100)
    print("done append data")
    act, trn = read_gz_file(B1_ACTIVE, B1_TRANSFOR, ',')  # only B1
    print("done reading the files")
    active = smooth_parse(act, "smoothSmall")  # only small ! there is also large
    transformed = smooth_parse(trn, "smoothSmall")
    print("done parsing")
    search(active, transformed, b_active, "active.csv")  # active file
    search(active, transformed, tf_trans, "transformed.csv")  # transformed file
    print("F I N I S H !")
    return


main_imm()
# main_plass()
# make_box_plot("no_treatment_vs_dac_and_hdac.csv")
# make_box_plot("no_treatment_vs_hdac.csv")
# make_box_plot("Compares files/no_treatment_vs_hdac_only.csv")


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
