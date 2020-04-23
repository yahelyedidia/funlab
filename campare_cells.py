import lab
import pandas as pd
import numpy as np

#genome build 38
CHR_I = 3
MATRIX = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/site_matrix"
COLUMNS = ["chr", "start", "end"]
THRESHOLD = 50
PANCREAS = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/pancreas/ENCFF994QQT.bed"
t_file = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/chr1_test"


CTCF_SITES = {}

def build_matrix(lst_of_site_files):
    """
    A function that initialize the site matrix and save it to csv file
    :param lst_of_site_files: list of BED file with sites
    :return: the matrix as DataFrame
    """
    matrix = pd.DataFrame(columns=COLUMNS)
    for f in lst_of_site_files:
        data = create_site_df(f)
        exist_site = False
        for chr, start, end in zip(data[0], data[1], data[2]):
            if chr[CHR_I:] in CTCF_SITES.keys():
                for site in CTCF_SITES[chr[CHR_I:]]:
                    if (site[0] - THRESHOLD <= start <= site[1]) and (site[0] <= end <= site[1] + THRESHOLD):
                        exist_site = True
                        break
                if not exist_site:
                    CTCF_SITES[chr[CHR_I:]].append([start, end])
                    exist_site = False
                    matrix = matrix.append(pd.DataFrame([[chr[CHR_I:], start, end]], columns=COLUMNS))
            else:
                CTCF_SITES[chr[CHR_I:]] = []
                CTCF_SITES[chr[CHR_I:]].append([start, end])
                exist_site = False
                matrix = matrix.append(pd.DataFrame([[chr[CHR_I:], start, end]], columns=COLUMNS))
        matrix = matrix.sort_values(COLUMNS)
        matrix.to_csv(MATRIX)
        return matrix


def create_site_df(f, to_sort=False):
    data = pd.read_csv(f, sep='\t', skiprows=131629, header=None)  # todo: deal with problemist rows
    data = data.drop(columns=[i for i in range(CHR_I, 10)])
    if to_sort:
        data = data.sort_values([0, 1, 2])
    return data


def add_cell(methylation_files, biding_file, name, as_lst=True, matrix=MATRIX):
    matrix = pd.read_csv(matrix)
    biding = create_site_df(biding_file, True)
    matrix[name] = "."
    if as_lst:
        for c in methylation_files:
            f = pd.read_csv(c, sep='\t')
            chrom_name = (list(f)[0].split("=")[1])[3:]
            for chr, start, end in zip(matrix["chr"], matrix["start"], matrix["end"]):
                if chrom_name < chr:
                    continue
                elif chrom_name > chr:
                    break
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


def new_site():
    pass

if __name__ == '__main__':
    # build_matrix([PANCREAS])
    add_cell([t_file], PANCREAS, "chr1")