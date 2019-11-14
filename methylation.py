import pandas as pd
from sklearn import preprocessing as p
import numpy as np
import matplotlib.pyplot as plt
import lab

TRASHOLD = 0.5


def create_avg(file):
    data = pd.read_csv(file, sep='\t', comment='t')
    data = data.drop_duplicates()
    data = data.set_index("ID_REF")
    avg = []
    binar = []
    for line in data.iterrows():
        lst = line[1]
        s = sum([float(pair[1]) for pair in lst])
        a = s/data.shape[1]
        avg.append(a)
        if a < TRASHOLD:
            binar.append(0)
        else:
            binar.append(1)
    data["avg"] = avg
    data["bin"] = binar
    return data

print(create_avg("GSE89269_series_matrix.txt").head())
