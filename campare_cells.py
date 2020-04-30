import lab
import pandas as pd
import numpy as np
import os

#genome build 38
MATRIX_SOURCE = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/CTCF.fimocentered200bpwherefound.min50.hg38.bed"
CHR_I = 3
MATRIX = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/site_matrix"
COLUMNS = ["chr", "start", "end"]
THRESHOLD = 100
PANCREAS_SITE = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/pancreas/ENCFF994QQT.bed"
PANCREAS_MET = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/pancreas/chrs_t"
t_file = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/chr1_test"




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


def add_cell(methylation_files_dir, biding_file, name, as_lst=True, matrix=MATRIX):
    matrix = pd.read_csv(matrix, index_col=0)
    biding = create_site_df(biding_file, True)
    matrix[name] = "."
    if as_lst:
        for c in os.listdir(methylation_files_dir):
            f = pd.read_csv(methylation_files_dir + os.sep + c, sep='\t')
            chrom_name = (list(f)[0].split("=")[1])[3:]
            for chr, start, end in zip(matrix["chr"], matrix["start"], matrix["end"]):
                if chrom_name < chr:
                    continue
                elif chrom_name > chr:
                    continue
                else:
                    met = np.mean(f.loc[(start - THRESHOLD <= f.index) & (f.index <= end + THRESHOLD)])[0]
                    bind = biding.loc[(biding[0] == "chr" + chr) & (start - THRESHOLD <= biding[1]) &(biding[1] <= end)
                                      & (start <= biding[2]) & (biding[2] <= end + THRESHOLD)]
                    if bind.empty:
                        matrix.loc[(matrix["chr"] == chr) & (matrix["start"] == start) & (matrix["end"] == end), name] = str((met, 0))
                    else:
                        matrix.loc[(matrix["chr"] == chr) & (matrix["start"] == start) & (matrix["end"] == end), name] = str((met, 1))
    matrix.to_csv(MATRIX) #todo: pay attention
    return

if __name__ == '__main__':
    print("start runing")
    add_cell(PANCREAS_MET, PANCREAS_SITE, "pancreas")