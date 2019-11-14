import pandas as pd
import lab
import methylation

PLASS3 = "files/Plass/ENCFF032DEW.bed"

PLASS2 = "files/Plass/ENCFF543VGD.bed"

PLASS1 = "files/Plass/ENCFF401ONY.bed"

def read_gz_file(file1, file2):
    nodrag = pd.read_csv(file1, sep='\t', low_memory=False)
    drag = pd.read_csv(file2, sep='\t', low_memory=False)
    return nodrag, drag


read_gz_file("files/Plass/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv", "files/Plass/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv")

def main():
    plass = [PLASS1, PLASS2, PLASS3]
