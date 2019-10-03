import pandas as pd
from sklearn import preprocessing as p
import numpy as np

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 10000)

def readChipFile(file):
    # לבדוק סינון לפי score
    data = pd.read_csv(file, sep='\t', comment='t', header=None)
    header = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalVal", "pVal",
              "qVal", "peak "]
    data.columns = header[:len(data.columns)]
    print(data.describe())
    data = data.drop_duplicates()
    data = data.drop(columns=["name", "score", "strand",  "signalVal", "pVal", "qVal"])
    data = data.sort_values(["chrom", "chromStart"])
    return data



def readMicroInfo(file):
    data = pd.read_csv(file)
    data = data.drop(columns=["IlmnID", "AddressA_ID", "AlleleA_ProbeSeq", "AddressB_ID", "AlleleB_ProbeSeq",
                              "Next_Base","Color_Channel", "Forward_Sequence", "SourceSeq", "Chromosome_36" ,
                              "Coordinate_36", "Probe_SNPs", "Probe_SNPs_10","UCSC_RefGene_Group",
                              "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_CpG_Islands_Name",
                              "DMR", "Enhancer", "Regulatory_Feature_Group", "Methyl27_Loci", "Random_Loci",
                              "Regulatory_Feature_Name", "DHS", "HMM_Island", "UCSC_CpG_Islands_Name"])
    # print(data.size)
    data = data[data['Genome_Build'] >= 36] #drop the empty genome build
    # data = data[data[''] >= 36] # drop the empty genome build

    # print(data.size)
    # x_centered = p.scale(data, with_mean='True', with_std='False')
    # data = pd.DataFrame(x_centered, columns=data.columns)
    print(data.head())
# def search(data, length):


readMicroInfo("normal.csv")

# readChipFile("ENCFF401ONY.bed.gz")
# print(readChipFile().head())
