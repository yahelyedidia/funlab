import pandas as pd
import lab
import methylation

DINFO = "Smoothed_Methylation_Level_H2_DMSO"

NODINFO = "Smoothed_Methylation_Level_H2_DAC"

DMC = "files/Plass/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

DMSO = "files/Plass/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv"

PLASS3 = "files/Plass/ENCFF032DEW.bed"

PLASS2 = "files/Plass/ENCFF543VGD.bed"

PLASS1 = "files/Plass/ENCFF401ONY.bed"


def read_gz_file(file1, file2, filter):
    nodrag = pd.read_csv(file1, sep='\t', low_memory=False)
    drag = pd.read_csv(file2, sep='\t', low_memory=False)
    # data = data[data['score'] >= score]
    nodrag = nodrag[nodrag[NODINFO] >= filter]
    drag = drag[drag[DINFO] >= filter]
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

def search(nodrags, drags, chip_data):
    """
    searching if prob is at the area of an peak +- the buffer.
    :param parse_data: the probes
    :param chip_data: the data from the chip array experience
    :param buffer: the buffer of bases we look at around the peak
    :return: the ratio between sum of probes are is the buffer to the total amount of probes
    """
    nodragcount = 0
    dragcount = 0
    for i in range(len(chip_data)):
        for start, end, chr in zip(chip_data[i]["chromStart"],
                                    chip_data[i]["chromEnd"],
                                    chip_data[i]["chrom"]):
            if chr == 'chrX':
                chr = 22
            elif chr == 'chrY':
                chr = 23
            else:
                chr = int(chr[3:]) - 1
            startin = find_start(nodrags[chr], start)
            endin = find_start(nodrags[chr], end)
            nodragcount += endin - startin
            startin = find_start(drags[chr], start)
            endin = find_start(drags[chr], end)
            dragcount += endin - startin
    print("no drags count :" + str(nodragcount) + "\nwith drags count : " + str(dragcount))
    print("the ratio :" + str(nodragcount - dragcount))









def main():
    plass = [PLASS1, PLASS2, PLASS3]
    data = []
    for p in plass:
        data.append(lab.read_chip_file(p, 100))
    filters = [0.1, 0.3, 0.5]
    for filter in filters:
        ndrg, drg = read_gz_file(DMSO, DMC, filter)
        nodrag = smooth_parse(ndrg, NODINFO)
        drag = smooth_parse(drg, DINFO)
        print("results for filter: " + str(filter))
        res = search(nodrag, drag, data)
    #     print(readin.head())
    #     chrom = lab.parse(readin,"chrom", "", chrom)


main()