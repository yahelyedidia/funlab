import pandas as pd
import lab
import numpy as np
import matplotlib.pyplot as plt
import os

NUM_OF_FIRST_ROWS_AT_PROBS_FILE = 7

NUM_OF_FIRST_ROWS_AT_CSC = 72


def read_micro_data(f_data, f_probs, num_of_open_line_data=NUM_OF_FIRST_ROWS_AT_CSC,
                    num_of_open_line_probs=NUM_OF_FIRST_ROWS_AT_PROBS_FILE, score=0):
    """
    A function that gets file of micro array data and read it
    :param f_data: path to the data file
    :param f_probs: path to the probs file
    :param num_of_open_line_data: number of lines to ignore in f_data
    :param num_of_open_line_probs: number of lines to ignore in f_probs
    :return: the data as dataframe
    """
    array_data = pd.read_csv(f_data, sep="\t", skiprows=[i for i in range(num_of_open_line_data)], header=0)
    probs_data = pd.read_csv(f_probs, sep=",", compression="gzip", skiprows=[i for i in range(num_of_open_line_probs)],
                             header=0)
    chromosomes = lab.parse(probs_data, "CHR", "MAPINFO")
    # remove unnecessary probes
    cloumns = list(probs_data)
    plass = ["tehila/Plass/ENCFF032DEW.bed.gz", "tehila/Plass/ENCFF401ONY.bed.gz", "tehila/Plass/ENCFF543VGD.bed.gz"]
    imm = ["tehila/immortalization/ENCFF449NOT.bed", "tehila/immortalization/ENCFF833FTF.bed"]
    new_data = pd.DataFrame(columns=cloumns)
    for val in imm:
        data = pd.read_csv(val, sep="\t", header=None)
        for row in data.iterrows():
            row = row[1]
            #todo search in parse probes data and append to new data while true





    # parse the probes data:
    # chromosomes = lab.parse()


read_micro_data("tehila/CSC/GSE101673_series_matrix.txt", "tehila/CSC/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")
