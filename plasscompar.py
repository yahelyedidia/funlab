import pandas as pd
import lab
import numpy as np
import matplotlib.pyplot as plt

DINFO = "Smoothed_Methylation_Level_H2_DMSO"

NODINFO = "Smoothed_Methylation_Level_H2_DAC"

CONTROL = "GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

AFTER_TREATMENT = "GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

PLASS3 = "ENCFF032DEW.bed"

PLASS2 = "ENCFF543VGD.bed"

PLASS1 = "ENCFF401ONY.bed"


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


def read_gz_file(file1, file2):
    """
    reading the file's and filter it if needed
    :param file1: the first file
    :param file2: the second file
    :return: the files
    """
    nodrag = pd.read_csv(file1, sep='\t', low_memory=False)
    drag = pd.read_csv(file2, sep='\t', low_memory=False)
    # nodrag = nodrag[nodrag[DINFO] >= filter]
    # drag = drag[drag[NODINFO] >= filter]
    return nodrag, drag


def smooth_parse(data, name):
    """
    parsing the smoothing file
    :param data: the data's file
    :param name: the file's name
    :return: an array with the data parsed from the file
    """
    chrom = [[] for i in range(24)]
    for index, loci, level in zip(data["Chromosome"], data["Start"], data[name]):
        if index == 'X':
            chrom[22].append([loci, level])
        elif index == 'Y':
            chrom[23].append([loci, level])
        else:
            chrom[int(index) - 1].append([loci, level])
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
    lst = pd.read_csv(file)
    lst = lst.drop(lst.index[[0]])
    # lst = pd.DataFrame(data=lst, columns=head)
    chroms = []
    for i in range(1, 24):
        this_chrom = lst[lst['chr'] == i]
        this_chrom = this_chrom.drop(columns=["chr", "start", "end", "no drugs avg", "with drugs avg"])
        this_chrom = this_chrom.drop(this_chrom.columns[0], axis=1)
        chroms.append(np.array(this_chrom))
    plt.boxplot(chroms)
    plt.title("changes at the mthylation level after treatment by chromosomes")
    plt.xlabel("chromosomes, 22 = X, 23 = Y")
    plt.ylabel("change level")
    plt.grid()
    # plt.setp(plt, xticks=[c+1 for c in range(23)])
    plt.legend()
    plt.savefig("boxplot_res")
    plt.show()

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


def search(nodrags, drags, chip_data):
    """
    searching if prob is at the area of an peak +- the buffer.
    :param chip_data: the data from the chip array experience
    :return: the ratio between sum of probes are is the buffer to the total amount of probes
    """
    # nodragcount = 0
    # dragcount = 0
    head = ["chr", "start", "end", "no drugs avg", "with drugs avg"]
    lst = np.empty((1, len(head)))
    for i in range(len(chip_data)):
        for start, end, chr, peak in zip(chip_data[i]["chromStart"],
                                    chip_data[i]["chromEnd"],
                                    chip_data[i]["chrom"],
                                         chip_data[i]["peak"]):
            if chr == 'chrX':
                chr = 22
            elif chr == 'chrY':
                chr = 23
            else:
                chr = int(chr[3:]) - 1
            no_drg_avg = closest_to_peak(nodrags[chr], peak, start)
            drg_avg = closest_to_peak(drags[chr], peak, start)
            line = np.array([chr + 1, start, end, no_drg_avg, drg_avg])
            lst = np.vstack([lst, line])

    #  adding column with the difference between the average with out / with the drag
    change = np.subtract(lst[:, 3], lst[:, 4])[np.newaxis]
    change = change.T
    lst = np.hstack([lst, change])
    head.append("change")
    lst = lst[lst[:, 0].argsort()] #todo change to 5
    #  writing the result
    pd.DataFrame(data=lst, columns=head).to_csv("in_progress.csv")


def main():
    plass = [PLASS1, PLASS2, PLASS3]
    data = []
    for p in plass:
        data.append(lab.read_chip_file(p, 100))
    print("done append data")
    # filters = [0.1, 0.3, 0.5]
    # for filter in filters:
    # print("results for filter: " + str(filter))
    ndrg, drg = read_gz_file(CONTROL, AFTER_TREATMENT)
    print("done reading the files")
    nodrag = smooth_parse(ndrg, DINFO)
    drag = smooth_parse(drg, NODINFO)
    print("done parsing")
    search(nodrag, drag, data)
    print("F I N I S H !")


make_box_plot("in_progress.csv")


def no_use():
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
    return