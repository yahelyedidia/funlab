import lab
import pandas as pd

#genome build 38
CHR_I = 3
MATRIX = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/site_matrix"
COLUMNS = ["chr", "start", "end"]
THRESHOLD = 50
PANCREAS = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/pancreas/ENCFF994QQT.bed"

CTCF_SITES = {}

def build_matrix(lst_of_site_files):
    counter = 1
    matrix = pd.DataFrame(columns=COLUMNS)
    for f in lst_of_site_files:
        data = pd.read_csv(f, sep='\t', skiprows=131630, header=None) #todo: deal with problemist rows
        data = data.drop(columns=[i for i in range(CHR_I, 10)])
        exist_site = False
        for chr, start, end in zip(data[0], data[1], data[2]):
            if chr[CHR_I:] in CTCF_SITES.keys():
                for site in CTCF_SITES[chr[CHR_I:]]:
                    if (site[0] - THRESHOLD <= start <= site[1]) and (site[0] <= end <= site[1] + THRESHOLD):
                        exist_site = True
                        break
            if not exist_site:
                if not chr[CHR_I:] in CTCF_SITES.keys():
                    CTCF_SITES[chr[CHR_I:]] = []
                CTCF_SITES[chr[CHR_I:]].append([start, end])
            exist_site = False
            matrix = matrix.append(pd.DataFrame([[chr[CHR_I:], start, end]], columns=COLUMNS))
        matrix = matrix.sort_values(COLUMNS)
        matrix.to_csv(MATRIX)
        return matrix


def add_cell(methylation_file, biding_file,matrix=MATRIX):
    pass

def new_site():
    pass

build_matrix([PANCREAS])