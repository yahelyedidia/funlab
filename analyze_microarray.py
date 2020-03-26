import pandas as pd
import lab
import numpy as np
import matplotlib.pyplot as plt
import os

NUM_OF_FIRST_ROWS_AT_PROBS_FILE = 7

NUM_OF_FIRST_ROWS_AT_CSC = 72

M1_REP1 = ["GSM2711824", "GSM2711825", "control_vs_csc_after_1_month_rep1"]
M1_REP2 = ["GSM2711832", "GSM2711833", "control_vs_csc_after_1_month_rep2"]
M6_REP1 = ["GSM2711826", "GSM2711827", "control_vs_csc_after_6_month_rep1"]
M6_REP2 = ["GSM2711834", "GSM2711835", "control_vs_csc_after_6_month_rep2"]
M10_REP1 = ["GSM2711828", "GSM2711829", "control_vs_csc_after_10_month_rep1"]
M10_REP2 = ["GSM2711836", "GSM2711837", "control_vs_csc_after_10_month_rep2"]
M15_REP1 = ["GSM2711830", "GSM2711831", "control_vs_csc_after_15_month_rep1"]
M15_REP2 = ["GSM2711838", "GSM2711839", "control_vs_csc_after_15_month_rep2"]

LST_CSC = [M1_REP1, M1_REP2, M6_REP1, M6_REP2, M10_REP1, M10_REP2, M15_REP1, M15_REP2]


def read_micro_data(f_data, f_probs, num_of_open_line_data=NUM_OF_FIRST_ROWS_AT_CSC,
                    num_of_open_line_probs=NUM_OF_FIRST_ROWS_AT_PROBS_FILE, score=0, buffer=250):
                    #todo: decide the buffer and score sizes
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
    chromosomes = lab.parse(probs_data, "CHR", "MAPINFO", None, True)
    # remove unnecessary probes
    cloumns = list(probs_data)
    plass = ["tehila/Plass/ENCFF032DEW.bed.gz", "tehila/Plass/ENCFF401ONY.bed.gz", "tehila/Plass/ENCFF543VGD.bed.gz"]
    imm = ["tehila/immortalization/ENCFF449NOT.bed", "tehila/immortalization/ENCFF833FTF.bed"]
    chip = []
    for val in imm:
        chip.append(lab.read_chip_file(val, score))
    for val in plass:
        chip.append(lab.read_chip_file(val, score, True))
    ids, inds = lab.search(chromosomes, chip, buffer, True)
    new_data = probs_data.loc[probs_data["IlmnID"].isin(ids)]
    array_data = array_data.loc[array_data["ID_REF"].isin(ids)]
    array_data = array_data.sort_values(by="ID_REF")
    inds = pd.DataFrame(inds).sort_values(by=0)
    sites = []
    for id in ids:
        df = inds.loc[inds[0] == id]
        c = df[1].max()
        start = df[2].mean().round()
        end = df[3].mean().round()
        sites.append([id, c, start, end])
    sites = pd.DataFrame(sites).sort_values(by=0)
    sites = sites.drop_duplicates()
    for lst in LST_CSC:
        control = array_data[lst[0]]
        treatment = array_data[lst[1]]
        change = control - treatment
        temp = sites.loc[sites[0].isin(array_data["ID_REF"])].sort_values(by=0)
        chr = temp[1]
        start = temp[2]
        end = temp[3]
        result = pd.concat([chr, start, end], axis=1)
        result = pd.concat([result.set_index(temp[0]), pd.DataFrame(control).set_index(array_data["ID_REF"])], axis=1)
        result = pd.concat([result.set_index(temp[0]), pd.DataFrame(treatment).set_index(array_data["ID_REF"])], axis=1)
        result['cov'] = '.'
        result["strand"] = '.'
        result = pd.concat([result, pd.DataFrame(change).set_index(array_data["ID_REF"])], axis=1, ignore_index=True)
        result.columns = ["chr", "start", "end", "control", "treatment", "cov", "strand" ,"change"]
        result.to_csv(lst[2], sep="\t")
    # chromosomes = lab.parse()


read_micro_data("tehila/CSC/GSE101673_series_matrix.txt", "tehila/CSC/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")
