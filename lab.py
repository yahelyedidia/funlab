import pandas as pd
from sklearn import preprocessing as p
import numpy as np
import matplotlib.pyplot as plt

IMM2 = "ENCFF833FTF.bed"

IMM1 = "ENCFF449NOT.bed"

PLASS3 = "ENCFF032DEW.bed"

PLASS2 = "ENCFF543VGD.bed"

PLASS1 = "ENCFF401ONY.bed"

GES = "ges1_ctcf_12_fp-stylefactor.bed"

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 10000)


def read_chip_file(file, score):
    data = pd.read_csv(file, sep='\t', comment='t', header=None)
    if file == GES:
        data = data.drop([i for i in range(3, 6)], axis=1)
        header = ["chrom", "chromStart", "chromEnd"]
        data.columns = header[:len(data.columns)]
        peak = data['chromEnd'].sub(data['chromStart']).to_frame('peak')
        data.insert(3, "peak", peak//2)
        data.drop_duplicates()
    else:
        header = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
                  "signalVal", "pVal", "qVal", "peak"]
        data.columns = header[:len(data.columns)]
        # לבדוק סינון לפי score
        data = data[data['score'] >= score]
        # print(data["chromStart"])
        data = data.drop_duplicates()
        data = data.drop(
            columns=["name", "score", "strand", "signalVal", "pVal", "qVal"])
        # data = data.sort_values(by=["chrom", "chromStart"])
    return data


def read_micro_info(file):
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


def parse(data):
    chrom = [[] for i in range(24)]
    for index, info in zip(data["CHR"], data["MAPINFO"]):
        if index == 'X':
            chrom[22].append(info)
        elif index == 'Y':
            chrom[23].append(info)
        else:
            chrom[int(index) - 1].append(info)
    for i in chrom:
        i.sort()
    return chrom


def search(parse_data, chip_data, buffer):
    row_num = 0
    counter = 0
    for i in range(len(chip_data)):
        for start, peak, chr in zip(chip_data[i]["chromStart"],
                                    chip_data[i]["peak"],
                                    chip_data[i]["chrom"]):
            if chr == 'chrX':
                chr = 22
            elif chr == 'chrY':
                chr = 23
            else:
                chr = int(chr[3:]) - 1
            if find_peak(parse_data[chr], peak, start, buffer):
                counter += 1
            row_num += 1
    return counter / row_num


def find_peak(lst, peak, start, buffer):
    first = 0
    last = len(lst) - 1
    val = start + peak
    while first <= last:
        mid = (first + last) // 2
        if val - buffer <= lst[mid] <= val + buffer:
            return True
        else:
            if val - buffer < lst[mid]:
                last = mid - 1
            else:
                first = mid + 1
    return False


def create_graph(score, data):
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
        "imm": [IMM1, IMM2],
        "GES-1": [GES]
    }
    data = parse(read_micro_info("normal.csv"))
    scores = [0, 600, 700, 800, 900, 1000]
    for score in scores:
        print("still alive")
        create_graph(score, data)


