import lab
import pandas as pd
import numpy as np
import time
import os

#genome build 38
START = "start"
CHR = "chr"
END = "end"
MIN_COV = 5
COV = 3
METHYLATION = 4
MATRIX_SOURCE = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/CTCF.fimocentered200bpwherefound.min50.hg38.bed"
CHR_I = 3
MATRIX = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/site_matrix"
COLUMNS = ["chr", "start", "end"]
THRESHOLD = 50
PANCREAS_SITE = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/pancreas/ENCFF994QQT.bed"
PANCREAS_MET = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/pancreas/chrs_t"
STOMACH = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/stomach/edit_stomach"
STOMACH_SITE = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/stomach/ENCFF933XOI.bed"
t_file = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/chr1_test"
SPLEEN_SITE = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/spleen/ENCFF250RMB.bed"
SPLEEN_MET = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/spleen/edit_spleen_met"
TEST_MET = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/spleen/test_data"
TEST_BIND = SPLEEN_SITE




def build_matrix():
    """
    A function that initialize the site matrix and save it to csv file
    :param lst_of_site_files: list of BED file with sites
    :return: the matrix as DataFrame
    """
    # matrix = pd.DataFrame(columns=COLUMNS)
    # for f in lst_of_site_files:
    #     data = create_site_df(f)
    #     exist_site = False
    #     for chr, start, end in zip(data[0], data[1], data[2]):
    #         if chr[CHR_I:] in CTCF_SITES.keys():
    #             for site in CTCF_SITES[chr[CHR_I:]]:
    #                 if (site[0] - THRESHOLD <= start <= site[1]) and (site[0] <= end <= site[1] + THRESHOLD):
    #                     exist_site = True
    #                     break
    #             if not exist_site:
    #                 CTCF_SITES[chr[CHR_I:]].append([start, end])
    #                 exist_site = False
    #                 matrix = matrix.append(pd.DataFrame([[chr[CHR_I:], start, end]], columns=COLUMNS))
    #         else:
    #             CTCF_SITES[chr[CHR_I:]] = []
    #             CTCF_SITES[chr[CHR_I:]].append([start, end])
    #             exist_site = False
    #             matrix = matrix.append(pd.DataFrame([[chr[CHR_I:], start, end]], columns=COLUMNS))
    #     matrix = matrix.sort_values(COLUMNS)
    #     matrix.to_csv(MATRIX)
    #     return matrix
    matrix = pd.read_csv(MATRIX_SOURCE, sep='\t', header=None)
    matrix.rename(columns={0:"chr", 1:"start", 2:"end"}, inplace=True)
    matrix.to_csv(MATRIX)


def create_site_df(f, to_sort=False):
    data = pd.read_csv(f, sep='\t', skiprows=[131629], header=None)  # todo: deal with problemist rows
    data = data.drop(columns=[i for i in range(CHR_I, 10)])
    if to_sort:
        data = data.sort_values([0, 1, 2])
    return data


def add_cell(methylation_files_dir, binding_file, name, as_lst=True, matrix_as_df = False, matrix=MATRIX):
    """
    A function that add new cell column to the big matrix and save the matrix as file
    :param methylation_files_dir: the path to the methylation file
    :param binding_file: the path to the binding file
    :param name:the name of the new column
    :param as_lst: a boolean parameter that decide if the data came as list
    :param matrix_as_df: a boolean parameter that decide if the matrix came as dataframe
    :param matrix: the site matrix
    :return:the matrix as dataframe
    """
    if not matrix_as_df:
        matrix = pd.read_csv(matrix, index_col=0)
    biding = create_site_df(binding_file, True)
    matrix[name] = "."
    print("I am here")
    if as_lst:
        for c in os.listdir(methylation_files_dir):
            f = pd.read_csv(methylation_files_dir + os.sep + c, sep='\t')
            chrom_name = (list(f)[0].split("=")[1])[COV:]
            for chr, start, end in zip(matrix[CHR], matrix[START], matrix[END]):
                if chrom_name != chr:
                    continue
                else:
                    met = np.mean(f.loc[(start - THRESHOLD <= f.index) & (f.index <= end + THRESHOLD)])[0]
                    bind = biding[(biding[0] == CHR + chr) & (start - THRESHOLD <= biding[1]) & (biding[1] <= end)
                                  & (start <= biding[2]) & (biding[2] <= end + THRESHOLD)]
                    if bind.empty:
                        matrix.loc[(matrix[CHR] == chr) & (start == matrix[START]) &
                                   (matrix[END] == end), name] = str((met, 0))
                    else:
                        matrix.loc[(matrix[CHR] == chr) & (start == matrix[START]) &
                                   (matrix[END] == end), name] = str((met, 1))
    else:
        f = pd.read_csv(methylation_files_dir, sep='\t', header=None)
        print("after reading f")
        f = f[f[COV] >= MIN_COV]
        print("filtered by coverage")
        f[METHYLATION] = f[METHYLATION] / 100
        for c in range(1, 25):
            if c == 23:
                chrom_name = "X"
            elif c == 24:
                chrom_name = "Y"
            else:
                chrom_name = str(c)
            f = f[f[0] == CHR + chrom_name]
            for chr, start, end in zip(matrix[CHR], matrix[START], matrix[END]):
                if chrom_name != chr:
                    continue
                else:
                    met = np.mean(f[(CHR + chr == f[0]) & (start - THRESHOLD <= f[1]) & (f[1] <= end + THRESHOLD)])[METHYLATION]
                    bind = ((biding[0] == CHR + chr) & (start - THRESHOLD <= biding[1]) & (biding[1] <= end)
                            & (start <= biding[2]) & (biding[2] <= end + THRESHOLD)).any()
                    if bind:
                        matrix.loc[(matrix[CHR] == chr) & (start == matrix[START]) &
                                   (matrix[END] == end), name] = str((met, 0))
                    else:
                        matrix.loc[(matrix[CHR] == chr) & (start == matrix[START]) &
                                   (matrix[END] == end), name] = str((met, 1))
    matrix.to_csv(MATRIX, sep="\t") #todo: pay attention
    return matrix

if __name__ == '__main__':
    print("start runing")
    start = time.process_time()
    add_cell(SPLEEN_MET, TEST_BIND, "time_test", False, True)
    print("end running")
    print(time.process_time() - start)