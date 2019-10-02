import pandas as pd

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

def readfile():
    data = pd.read_csv("ENCFF401ONY.bed", sep='\t', comment='t', header=None)
    header = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalVal", "pVal", "qVal", "peak "]
    data.columns = header[:len(data.columns)]
    print(data.describe())
readfile()
