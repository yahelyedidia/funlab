import pandas as pd

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
    data = data.drop(columns=["IlmnID", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq", "Next_Base","Color_Channel","Forward_Sequence",  "SourceSeq", "UCSC_RefGene_Group"])
    
    print(data.head())

# def search(data, length):

readMicroInfo("normal.csv")

# readChipFile("ENCFF401ONY.bed.gz")
# print(readChipFile().head())
