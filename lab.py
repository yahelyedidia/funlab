import pandas as pd
#from sklearn import preprocessing as p
#import numpy as np
import matplotlib.pyplot as plt

Y_CHR = 23

X_CHR = 22

IMM2 = "ENCFF833FTF.bed.gz"

IMM1 = "ENCFF449NOT.bed.gz"

PLASS3 = "ENCFF032DEW.bed"

PLASS2 = "ENCFF543VGD.bed"

PLASS1 = "ENCFF401ONY.bed"

GES = "ges1_ctcf_12_fp-stylefactor.bed"

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 10000)


def read_chip_file(file, score, is_commpress=False):
    """
    A function that gets ChipSeq and convert the file to csv file while
    remove the row that ander the score treshold
    :param file: the file to convert
    :param score: the treshold
    :return: the file in csv format
    """
    if is_commpress:
        data = pd.read_csv(file, sep='\t', comment='t', header=None, compression="gzip")
    else:
        data = pd.read_csv(file, sep='\t', comment='t', header=None)
    if file == GES:
        data = data.drop([i for i in range(3, 6)], axis=1)
        header = ["chrom", "chromStart", "chromEnd"]
        data.columns = header[:len(data.columns)]
        peak = data['chromEnd'].sub(data['chromStart']).to_frame('peak')
        data.insert(3, "peak", peak//2)
    else:
        header = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
                  "signalVal", "pVal", "qVal", "peak"]
        data.columns = header[:len(data.columns)]
        #  check by score
        data = data[data['score'] >= score]
        # print(data["chromStart"])
        data = data.drop(
            columns=["name", "score", "signalVal", "pVal", "qVal"])  # todo : added the strand col to the output
        # data = data.sort_values(by=["chrom", "chromStart"])

    data = data.dropna()
    data = data.drop_duplicates()
    return data


def read_micro_info(file):
    """
    reading the information from the micro chip lab. filters by Genome_Build 36
    :param file: the micro chip csv file
    :return: pandas matrix with the filtered data
    """
    data = pd.read_csv(file)
    data = data.drop_duplicates()
    data = data.drop(
        columns=["IlmnID", "AddressA_ID", "AlleleA_ProbeSeq", "AddressB_ID",
                 "AlleleB_ProbeSeq",
                 "Next_Base", "Color_Channel", "Forward_Sequence", "SourceSeq",
                 "Chromosome_36",
                 "Coordinate_36", "Probe_SNPs", "Probe_SNPs_10",
                 "UCSC_RefGene_Group",
                 "UCSC_RefGene_Name", "UCSC_RefGene_Accession",
                 "UCSC_CpG_Islands_Name",
                 "DMR", "Enhancer", "Regulatory_Feature_Group",
                 "Methyl27_Loci", "Random_Loci",
                 "Regulatory_Feature_Name", "DHS", "HMM_Island",
                 "UCSC_CpG_Islands_Name"])
    data = data[data['Genome_Build'] >= 36]  # drop the empty genome build
    # data = data[data['CHR'] > 1] # drop the empty chrom
    return data
# def search(data, length):


def parse(data, charname, value, chrom=None, id=False):
    """
    parsing the data by the micro array test.
    :param data: the probes information
    :param charname: the number of name of chromosome column to parse by
    :param value:
    :param chrom: list of data arrange
    :return: array with the locations sorted by chromosomes
    """
    if (chrom == None):
        chrom = [[] for i in range(24)]
    if id:
        for index, info, id in zip(data[charname], data[value], data["IlmnID"]):
            if index == 'X':
                chrom[X_CHR].append((info, id))
            elif index == 'Y':
                chrom[Y_CHR].append((info, id))
            else:
                if pd.notna(index) and pd.notna(info):
                    chrom[int(index) - 1].append((info, id))
    else:
        for index, info, id in zip(data[charname], data[value]):
            if index == 'X':
                chrom[X_CHR].append(info)
            elif index == 'Y':
                chrom[Y_CHR].append(info)
            else:
                if pd.notna(index) and pd.notna(info):
                    chrom[int(index) - 1].append(info)
    for i in chrom:
        i.sort()
    return chrom


def search(parse_data, chip_data, buffer, to_write=False):
    """
    searching if prob is at the area of an peak +- the buffer.
    :param parse_data: the probes
    :param chip_data: the data from the chip array experience
    :param buffer: the buffer of bases we look at around the peak
    :return: the ratio between sum of probes are is the buffer to the total amount of probes
    """
    row_num = 0
    counter = 0
    id_to_select = []
    ind = []
    for i in range(len(chip_data)):
        for start, end, peak, chr in zip(chip_data[i]["chromStart"], chip_data[i]["chromEnd"],
                                    chip_data[i]["peak"],
                                    chip_data[i]["chrom"]):
            if chr == 'chrX':
                c = 22
                chr = 'X'
            elif chr == 'chrY':
                c = 23
                chr = 'Y'
            else:
                chr = chr[3:]
                c = int(chr) - 1
            id =  find_peak(parse_data[c], peak, start, buffer, to_write)
            if id:
                id_to_select.append(id)
                ind.append([id, chr, start, end])
                counter += 1
            row_num += 1
    if to_write:
        return id_to_select, ind
    return counter / row_num


def find_peak(lst, peak, start, buffer, id=False):
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
    while first <= last:
        mid = (first + last) // 2
        if id:
            if val - buffer <= lst[mid][0] <= val + buffer:
                return lst[mid][1]
            else:
                if val - buffer < lst[mid][0]:
                    last = mid - 1
                else:
                    first = mid + 1
        else:
            if val - buffer <= lst[mid] <= val + buffer:
                return mid
            else:
                if val - buffer < lst[mid]:
                    last = mid - 1
                else:
                    first = mid + 1
    return False


def create_graph(score, data):
    """
    not relevant, drawing bed graph
    """
    for (k, v) in dict.items():
        chip = []
        # print("***** results for " + str(l) + ": ******")
        for i in range(len(v)):
            chip.append(read_chip_file(v[i], score))
        # print(chip)
        buff = [250, 500, 750, 1000]
        y = []
        for b in buff:
            y.append(search(data, chip, b))
        plt.plot(buff, y, label=k)
    plt.title("results vs range of search, data filtered by score " + str(score))
    plt.grid()
    plt.legend()
    plt.savefig("results by different range for score " + str(score))
    plt.show()


if __name__ == '__main__':
    dict = {
        "all files": [IMM1, PLASS1, PLASS2,
             PLASS3, IMM2, GES],
        "Plass": [PLASS1, PLASS2, PLASS3],
        "immortalization": [IMM1, IMM2],
        "GES-1": [GES]
    }
    data = parse(read_micro_info("normal.csv"), "CHR", "MAPINFO")
    scores = [0, 600, 700, 800, 900, 1000]
    for score in scores:
        print("still alive")
        create_graph(score, data)


