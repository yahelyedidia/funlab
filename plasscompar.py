import pandas as pd
import lab
import numpy as np
# import methylation
# from statistics import mean
import csv

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
    :param lst:
    :param peak:
    :param start:
    :param buffer:
    :return:
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
    nodrag = pd.read_csv(file1, sep='\t', low_memory=False)
    drag = pd.read_csv(file2, sep='\t', low_memory=False)
    # nodrag = nodrag[nodrag[DINFO] >= filter]
    # drag = drag[drag[NODINFO] >= filter]
    return nodrag, drag


def smooth_parse(data, name):
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


def find_start(lst, start):
    first = 0
    last = len(lst) - 1
    while first <= last:
        mid = (first + last) // 2
        if start == lst[mid][0]:
            return mid
        else:
            if start < lst[mid][0]:
                last = mid - 1
            else:
                first = mid + 1
    return mid


def search(nodrags, drags, chip_data, name):
    """
    searching if prob is at the area of an peak +- the buffer.
    :param parse_data: the probes
    :param chip_data: the data from the chip array experience
    :param buffer: the buffer of bases we look at around the peak
    :return: the ratio between sum of probes are is the buffer to the total amount of probes
    """
    nodragcount = 0
    dragcount = 0
    head = ["chr", "start", "end", "no drugs avg", "with drugs avg"]
    lst = np.empty((1, len(head)))
    # df = pd.DataFrame(columns=head)
    # f = open(name, "w+")
    # f.write("** chr ** \t ** start **\t ** end ** \t **no drugs avg** \t **with drugs avg** \n")
    # with open(name, mode='w') as res_file:
        # writer = csv.DictWriter(res_file, fieldnames=head)
        # writer.writeheader()
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
            startin = find_start(nodrags[chr], start)
            endin = find_start(nodrags[chr], end)
            #no_drg_lst = np.array(nodrags[chr][startin:endin])
            nodrag_met = endin - startin
            nodragcount += nodrag_met
            startin = find_start(drags[chr], start)
            endin = find_start(drags[chr], end)
            #drg_lst = np.array(drags[chr][startin:endin])
            drag_met = endin - startin
            dragcount += drag_met
            # if nodrag_met - drag_met > 1 or nodrag_met - drag_met < -1:
            no_drg_avg = closest_to_peak(nodrags[chr], peak, start)
            drg_avg = closest_to_peak(drags[chr], peak, start)
            # line = [str(chr+1), str(start), str(end), str(no_drg_avg), str(drg_avg)]
            # df.append(pd.DataFrame(data=line))
            line = np.array([chr + 1, start, end, no_drg_avg, drg_avg])
            lst = np.vstack([lst, line])

                # writer.writerow({"chr": str(chr+1), "start": str(start), "end": str(end),
                #                  "no drugs avg": str(no_drg_avg), "with drugs avg":  str(drg_avg), "change": str(no_drg_avg- drg_avg)})
                    # f.write(str(chr+1) + "\t" + str(start) + "\t" + str(end) + "\t" + str(no_drg_avg) + "\t" + str(drg_avg) + "\n")
    # f.close()
    # df.to_csv(name)
    change = np.subtract(lst[:, 3], lst[:, 4])[np.newaxis]
    change = change.T
    lst[:, :-1] = change
    head.append("change")
    print(lst)
    lst = lst[::, lst[5, ].argsort()]

    pd.DataFrame(lst).to_csv("changes.csv", header=head)
    # np.savetxt("changes.csv", int(lst.flatten()))
    # fmt = ",".join(["%s"] + ["%10.6e"] * (lst.shape[1] - 1))
    # np.savetxt("changes.csv", int(lst), fmt=fmt, header=str(head), comments='')
    print("no drags count :" + str(nodragcount) + "\nwith drags count : " + str(dragcount))
    print("the ratio :" + str(nodragcount - dragcount))


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
    nodrag = smooth_parse(ndrg, DINFO)
    drag = smooth_parse(drg, NODINFO)
    res = search(nodrag, drag, data, "change with no filter.csv")
    #     print(readin.head())
    #     chrom = lab.parse(readin,"chrom", "", chrom)


main()