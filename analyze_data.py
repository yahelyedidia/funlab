import pandas as pd
# from sklearn import preprocessing as p
import numpy as np
import matplotlib.pyplot as plt
from itertools import islice

def filter_data(filter, d, col, name):
    """
    A function that filter the data according do given arg and write it to new file
    :param filter: the arg to filtered by it
    :param d: the data
    :param col: the col to filter by
    :param name: the name of the filtered file
    :return: the data filters
    """
    data = pd.read_csv(d)
    data = data.drop(data.columns[0], axis=1)
    if filter >= 0:
        data = data[data[col] >= filter]
    else:
        data = data[data[col] <= filter]
    data = remove_duplicate(data, "chr", "start", "end", "change")
    pd.DataFrame(data=data).to_csv(name)
    return data

def remove_duplicate(data, chr_col, start, end, change=""):
    """
    get data and remove the duplicate site
    :param data: the data to fix
    :param chr_col: a string that represent the name of the chrom column in the data
    :param start: the string that represent the name of the start column in the data
    :return: the fixed data
    """
    for ind, i in data.iterrows():
        for k, j in islice(data.iterrows(), 1, None):
            if ind >= k:
                continue
            if i[chr_col] == j[chr_col]:
                if abs(i[start] - j[start]) <= 10 and abs(i[end] - j[end]) <= 10:
                    data = data.drop(k)
                    continue
                if abs(i[start] - j[start]) <= 100 or abs(i[end] - j[end]) <= 100:
                    if change != "":
                        if i[change] == j[change]:
                            data = data.drop(k)
            continue
    return data

filter_data(0.6, "in_progress.csv", "change", "decrease_mthylation_plass")
filter_data(-0.6, "in_progress.csv", "change", "increase_mthylation_plass")