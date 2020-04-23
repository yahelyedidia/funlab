import lab
import pandas as pd

#genome build 38
PANCREAS = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/pancreas/ENCFF994QQT.bed"

CTCF_SITES = {}

def build_matrix(lst_of_site_files):
    counter = 1
    matrix = pd.DataFrame(columns=["chr", "start", "end"])
    for f in lst_of_site_files:
        data = pd.read_csv(f, sep='\t', skiprows=131630, header=None) #todo: deal with problemist rows
        data = data.drop(columns=[i for i in range(3,10)])
        exist_site = False
        for chr, start, end in zip(data[0], data[1], data[2]):
            if chr[3:] in CTCF_SITES.keys():
                for site in CTCF_SITES[chr[3:]]:
                    if (site[0] - 200 <= start <= site[1]) and (site[0] <= end <= site[1] + 200):
                        exist_site = True
                        break
            if not exist_site:
                if not chr[3:] in CTCF_SITES.keys():
                    CTCF_SITES[chr[3:]] = []
                CTCF_SITES[chr[3:]].append([start, end])
            exist_site = False
            matrix.append([[chr[3:], start, end]])
        x=1


def add_cell(methylation_file, biding_file):
    pass

def new_site():
    pass

build_matrix([PANCREAS])