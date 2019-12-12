import pandas as pd
from sklearn import preprocessing as p
import numpy as np
import matplotlib.pyplot as plt

def filter_data(filter, d, col, name):
    """
    A function that filter the data according do given arg and write it to new file
    :param filter: the arg to filtered by it
    :param d: the data
    :param col: the col to filter by
    :param name: the name of the filtered file
    :return: the data filters
    """
    data = pd.read_csv(d, sep=',', comment='t')
    data = data[data[col] >= filter]
    pd.DataFrame(data=data).to_csv(name)
    return data

filter_data(0.6, "in_progress.csv", "change", "decrease_mthylation_plass")
filter_data(-0.6, "in_progress.csv", "change", "increase_mthylation_plass")