import pandas as pd
#from sklearn import preprocessing as p
#import numpy as np
#import matplotlib.pyplot as plt
import lab

TRASHOLD = 0.5

def create_avg(file):
    data = pd.read_csv(file, sep='\t', comment='t', header=None, error_bad_lines=False)
    data = data.drop_duplicates()
    data = data.set_index("ID_REF")
    avg = []
    binar = []
    print("hi")
    for line in data.iterrows():
        lst = line[1]
        s = sum([float(pair[1]) for pair in lst])
        a = s/data.shape[1]
        avg.append(a)
        if a < TRASHOLD:
            binar.append(0)
        else:
            binar.append(1)
        print("end of line, going to the next line")
    data["avg"] = avg
    data["bin"] = binar
    return data


# print(create_avg("GES/GSE89269_series_matrix.txt.gz").head())
def try_on_file(file):
    c = pd.read_csv(file, sep="\t", skiprows=[0], headr=None)
    x=1

try_on_file("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/GSM1010981_UCSD.Adrenal_Gland.Bisulfite-Seq.STL003.wig")