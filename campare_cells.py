import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy.stats as st
from cells_dict import *

#genome build 38
DYNAMIC_STATE_TSV = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/sites_groups/dynamic_state.tsv"
STABLE_MET_TSV = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/sites_groups/stable_met.tsv"
BOUND_TSV = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/sites_groups/bound.tsv"
BOUND_STABLE_MET_TSV = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/sites_groups/bound_stable_met.tsv"
START = "start"
CHR = "chr"
END = "end"
MIN_COV = 5
COV = 3
METHYLATION = 4
MATRIX_SOURCE = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/CTCF.fimocentered200bpwherefound.min50.hg38.bed"
CHR_I = 3
MATRIX = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/the_big_matrix.tsv" #todo:change
MATRIX_FOR_PLAY = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/cell/site_&_bind_matrix.tsv"
COLUMNS = ["chr", "start", "end"]
THRESHOLD = 50
DIFFERENT = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/different_groups"


def build_matrix():
    """
    A function that initialize the site matrix and save it to csv file
    :param lst_of_site_files: list of BED file with sites
    :return: the matrix as DataFrame
    """
    matrix = pd.read_csv(MATRIX_SOURCE, sep='\t', header=None)
    matrix.rename(columns={0:"chr", 1:"start", 2:"end"}, inplace=True)
    matrix.to_csv(MATRIX, sep="\t")


def create_site_df(f, to_sort=False):
    """
    A function that get sites file and read it to dataframe
    :param f: the file to read
    :param to_sort: a boolean flag to sign if we need to sort the data
    :return: the file as dataframe
    """
    data = pd.read_csv(f, sep='\t', skiprows=[131629], header=None)  # todo: deal with problemist rows
    data = data.drop(columns=[i for i in range(CHR_I, 10)])
    if to_sort:
        data = data.sort_values([0, 1, 2])
    return data


def add_cell(methylation_files_dir, binding_file, name, as_lst=True, matrix_as_df = False, matrix=MATRIX):
    """
    A function that add new cell column to the big matrix and save the matrix as file
    PAY ATTENTION: at the end of the function you will save the update file instead of the
    original file. If you don't want to do it, please silence the relevant line.
    :param methylation_files_dir: the path to the methylation file
    :param binding_file: the path to the binding file
    :param name:the name of the new column
    :param as_lst: a boolean parameter that decide if the data came as list
    :param matrix_as_df: a boolean parameter that decide if the matrix came as dataframe
    :param matrix: the site matrix
    :return:the matrix as dataframe
    """
    if not matrix_as_df:
        matrix = pd.read_csv(matrix, index_col=0, sep="\t")
    print(matrix.describe())
    biding = create_site_df(binding_file, True)
    matrix[name + "_met"] = "."
    matrix[name + "_bind"] = "."
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
        f = f[f[COV] >= MIN_COV]
        f[METHYLATION] = f[METHYLATION] / 100
        for c in range(1, 25):
            if c == 23:
                chrom_name = CHR + "X"
            elif c == 24:
                chrom_name = CHR + "Y"
            else:
                chrom_name = CHR + str(c)
            file = f[f[0] == chrom_name]
            for chr, start, end in zip(matrix[CHR], matrix[START], matrix[END]):
                if chrom_name != chr:
                    continue
                else:
                    met = np.mean(file[(chr == file[0]) & (start - THRESHOLD <= file[1]) & (file[1] <= end + THRESHOLD)])[METHYLATION]
                    matrix.loc[(matrix[CHR] == chr) & (start == matrix[START]) & (matrix[END] == end), name + "_met"] = met
                    bind = ((biding[0] == chrom_name) & (((biding[1] <= start) & (end <= biding[2])) |
                                         ((start <= biding[1]) & (biding[2] <= end)) |
                                         ((biding[1] - THRESHOLD <= start) & (end <= biding[2]  + THRESHOLD)) )).any()
                    if bind:
                        matrix.loc[(matrix[CHR] == chr) & (start == matrix[START]) &
                                   (matrix[END] == end), name + "_bind"] = 1
                    else:
                        matrix.loc[(matrix[CHR] == chr) & (start == matrix[START]) &
                                   (matrix[END] == end), name + "_bind"] = 0
    matrix.to_csv(MATRIX, sep="\t") # WARNING: pay attention if you really want to replace the original file
    return matrix

def mann_whitney_and_fun(matrix):
    """
    A function that apply mann whitney test on each CTCF binding site
    :param matrix: the healthy cells matrix
    :return: save a tsv file with the results
    """
    matrix = pd.read_csv(matrix, sep="\t")
    col_name = list(matrix.columns)
    bind_col = [col_name[i] for i in range(5, len(col_name), 2)] # get the binding columns
    met_col = [col_name[i] for i in range(4, len(col_name), 2)] # get the methylation columns
    # calculate the average of methylation rate and binding rate for each binding site
    matrix["binding_rate"] = matrix[bind_col].mean(axis=1, skipna = True)
    matrix["met_rate"] = matrix[met_col].mean(axis=1, skipna = True)
    # create new dataframe of all the CTCF binding sites that always bound
    always_bound = matrix[matrix["binding_rate"] == 1]
    l = matrix.shape[0]
    # filter the matrix data to sites that bound at least in 5 different cell types
    matrix = matrix[(matrix["binding_rate"] > 5 / len(bind_col)) & (matrix["binding_rate"] < 1- (3 / len(bind_col)))]
    np_data = np.array(matrix[met_col])
    # calculate the variance of the methylation rate
    vars = np.var(np_data, axis=1)
    matrix['met_var'] = vars
    m = matrix[matrix["met_var"] != 0]
    # create new data frame of all the CTCF binding sites that there methylation variance is less then 0.01
    never_met = matrix[matrix["met_var"] <= 0.01]
    bound_avg = []
    unbound_avg = []
    p_val = []
    counter = 0
    # apply mann whitney test for each site
    for site in m.iterrows():
        bound = []
        unbound = []
        for cell in range(len(bind_col)):
            if site[1][bind_col[cell]] == 1:
                bound.append(float(site[1][met_col[cell]]))
            else:
                unbound.append(float(site[1][met_col[cell]]))
        if site[1]["met_rate"] == 0:
            p_val.append(None)
        else:
            p_val.append(st.mannwhitneyu(bound, unbound).pvalue)
        bound_avg.append((sum(bound))/len(bound))
        if len(unbound) != 0:
            unbound_avg.append((sum(unbound))/len(unbound))
        else:
            unbound_avg.append(0)
        counter += 1
    m["binded_avg"] = bound_avg
    # m["binded_avg"] = binded_avg
    m["unbinded_avg"] = unbound_avg
    m["p_val"] = p_val
    # create new dataframe of the significant values and not significant values
    sg_matrix = m[m["p_val"] <= 0.05]
    nsg_matrix = m[m["p_val"] > 0.05]
    sg_matrix.to_csv("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/significant_sites_all_chr_p=0.05.tsv", sep="\t")
    nsg_matrix.to_csv("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/not_significant_sites_all_chr_p=0.05.tsv", sep="\t")
    m.to_csv("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/significant_sites_all.tsv", sep="\t")

    # draw pie charts
    langs = ['Always bound', 'Stable methylation levels', 'Dynamic relationship']
    p = lambda x: (x/l) * 100
    a_r = p(always_bound.shape[0])
    n_r = p(never_met.shape[0])
    ratio = [a_r,n_r , 100 - a_r - n_r]
    print(ratio)
    pie(ratio, langs)
    s_r = p(sg_matrix.shape[0])

    langs.append("p value < 0.05")
    ratio2 = [a_r, n_r, 100 - a_r - n_r - s_r, s_r]
    print(ratio2)
    pie(ratio2, langs)

def pie(values, labels, color='RdPu'):
    """
    A function that draw pie chart
    :param values: the values for the chart
    :param labels: the labels for the chart
    :param color: the color of the chart
    :return: show the chart
    """
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    ax.axis('equal')
    c = plt.get_cmap(color)
    ax.set_prop_cycle("color", [c(1. * i/len(values)) for i in range(len(values))])
    ax.pie(values, labels = labels,autopct='%1.2f%%', labeldistance=1000)
    plt.legend()
    plt.show()


def play_with_data(matrix):
    """
    A function that split the big matrix to 4 groups:
    1. always bound (binding rate == 1) and dynamic methylation (met_var >= 0.01)
    2. always bound (binding rate == 1) and stable methylation (met_var < 0.01)
    3. stable methylation (met_var < 0.01) but not always bound (binding_rate != 1)
    4. dynamic state (met_var >=0.01, binding_rate != 1)
    :param matrix: the healthy cells matrix
    :return: save the groups as tsv
    """
    matrix = pd.read_csv(matrix, sep="\t")
    col_name = list(matrix.columns)
    bind_col = [col_name[i] for i in range(5, len(col_name), 2)]
    met_col = [col_name[i] for i in range(4, len(col_name), 2)]
    matrix["binding_rate"] = matrix[bind_col].mean(axis=1, skipna = True)
    matrix["met_rate"] = matrix[met_col].mean(axis=1, skipna = True)
    matrix["met_var"] = matrix[met_col].var(axis=1, skipna = True)
    matrix = matrix[matrix["binding_rate"] > 5/len(bind_col)]
    always_bound_and_stable_met = matrix[(matrix["binding_rate"] == 1) & (matrix["met_var"] <= 0.01)]
    always_bound_and_stable_met.to_csv(BOUND_STABLE_MET_TSV, sep="\t")
    always_bound = matrix[(matrix["binding_rate"] == 1) & (matrix["met_var"] > 0.01)]
    always_bound.to_csv(BOUND_TSV, sep="\t")
    stable_met = matrix[(matrix["binding_rate"] != 1) & (matrix["met_var"] <= 0.01)]
    stable_met.to_csv(STABLE_MET_TSV, sep="\t")
    dynamic_state = matrix[(matrix["binding_rate"] != 1) & (matrix["met_var"] > 0.01)]
    dynamic_state.to_csv(DYNAMIC_STATE_TSV, sep="\t")
    # is_it_the_same_distribution(matrix, met_col, ")
    all_cells = []
    for i in range(len(met_col)):
        all_cells.append(matrix[met_col[i]].tolist())
    s, p_val = st.kruskal(*zip(*all_cells))
    fig, axes = plt.subplots(2, 3, dpi=100, sharex=True, sharey=True)
    # fig, axes = plt.subplots(4, 5, dpi=100, sharex=True, sharey=True)
# colors = ['tab:red', 'tab:blue', 'tab:green', 'tab:pink', 'tab:olive']
    a = axes.flatten()
    cell_counter = 4
    axes_counter = 0
    bind_data = []
    unbind_data = []
    while cell_counter < len(col_name): #todo
        cell = (col_name[cell_counter].split("_"))[0]
        met = col_name[cell_counter]
        bind = col_name[cell_counter + 1]
        binded = matrix.loc[matrix[bind] == 1, met]
        unbinded = matrix.loc[matrix[bind] == 0, met]
        a[axes_counter].hist(binded, alpha=0.5, bins=50, density=True, stacked=True, label="methylation at the binded sites", color="darkturquoise")
        a[axes_counter].hist(unbinded, alpha=0.5, bins=50, density=True, stacked=True, label="methylation at the unbinded sutes", color="lightpink")
        # a[axes_counter].set_title(cell)
        plt.yscale("log")
        cell_counter = cell_counter + 2
        axes_counter = axes_counter + 1
        bind_data.append(binded.tolist())
        unbind_data.append(unbinded.tolist())
    # plt.title("methylation distribution at CTCF binding site in diffrent cells",  y=4.95, x=-2 ,size=16)
    plt.tight_layout()
    # plt.legend(loc="lower center)
    # plt.savefig("tests on healty data")
    plt.show()
    s, p_val = st.kruskal(*zip(*bind_data))
    print("all cell binding data")
    print("p value is {0}".format(p_val))
    s, p_val = st.kruskal(*zip(*unbind_data))
    print("all cell unbinding data")
    print("p value is {0}".format(p_val))

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')

def sort_by(l, m, h, val):
    if val > h:
        return 1
    if val > m:
        return 2
    if val > l:
        return 3
    return 4

def different_cuts(matrix):
    matrix = pd.read_csv(matrix, sep="\t")
    col_name = list(matrix.columns)
    bind_col = [col_name[i] for i in range(5, len(col_name), 2)]
    met_col = [col_name[i] for i in range(4, len(col_name), 2)]
    matrix["binding_rate"] = matrix[bind_col].mean(axis=1, skipna = True)
    matrix["met_rate"] = matrix[met_col].mean(axis=1, skipna = True)
    matrix["met_var"] = matrix[met_col].var(axis=1, skipna = True)
    matrix = matrix[matrix["binding_rate"] > 5/len(bind_col)]

    print(matrix["met_var"].quantile([.25, .5, .75]))
    ml = matrix["met_var"].quantile(.25)
    mm = matrix["met_var"].quantile(.5)
    mh = matrix["met_var"].quantile(.75)
    l_var = matrix[matrix["met_var"] <= ml]
    l_var.to_csv(DIFFERENT + os.sep + "methylation_variance_ander_{0:.2f}.tsv".format(ml), sep='\t')
    lm_var = matrix[(matrix["met_var"] > ml) & (matrix["met_var"] <= mm)]
    lm_var.to_csv(DIFFERENT + os.sep + "methylation_variance_between_{0:.2f}_to_{1:.2f}.tsv".format(ml, mm), sep='\t')
    mh_var = matrix[(matrix["met_var"] > mm) & (matrix["met_var"] <= mh)]
    mh_var.to_csv(DIFFERENT + os.sep + "methylation_variance_between_{0:.2f}_to_{1:.2f}.tsv".format(mm, mh), sep='\t')
    h_var = matrix[matrix["met_var"] > mh]
    h_var.to_csv(DIFFERENT + os.sep + "methylation_variance_above_{0:.2f}.tsv".format(mh), sep='\t')
    # methyaltion_matrix = pd.DataFrame(matrix["binding_rate"])
    # methyaltion_matrix["met_group"] = matrix["met_var"].apply(lambda x: sort_by(ml, mm, mh, x))
    # methyaltion_matrix.boxplot(by="met_group", grid=False)
    # plt.show()
    print(matrix["binding_rate"].quantile([.25, .5, .75]))
    bl = matrix["binding_rate"].quantile(.25)
    bm = matrix["binding_rate"].quantile(.5)
    bh = matrix["binding_rate"].quantile(.75)
    l_binding = matrix[matrix["binding_rate"] <= bl]
    l_binding.to_csv(DIFFERENT + os.sep + "binding_rate_ander_{0:.2f}.tsv".format(bl), sep='\t')
    lm_binding = matrix[(matrix["binding_rate"] > bl) & (matrix["binding_rate"] <= bm)]
    lm_binding.to_csv(DIFFERENT + os.sep + "binding_rate_between_{0:.2f}_to_{1:.2f}.tsv".format(bl, bm), sep='\t')
    mh_binding = matrix[(matrix["binding_rate"] > bm) & (matrix["binding_rate"] <= bh)]
    mh_binding.to_csv(DIFFERENT + os.sep + "binding_rate_between_{0:.2f}_to_{1:.2f}.tsv".format(bm, bh), sep='\t')
    h_binding = matrix[matrix["binding_rate"] > bh]
    h_binding.to_csv(DIFFERENT + os.sep + "binding_rate_above_{0:.2f}.tsv".format(bh), sep='\t')
    # l_binding.plot.scatter(x="met_rate", y="met_var")
    # lm_binding.plot.scatter(x="met_rate", y="met_var")
    # mh_binding.plot.scatter(x="met_rate", y="met_var")
    # h_binding.plot.scatter(x="met_rate", y="met_var")

    # binding_matrix = pd.DataFrame(matrix["met_rate"])
    # binding_matrix["binding_group"] = matrix["binding_rate"].apply(lambda x: sort_by(bl, bm, bh, x))
    # binding_matrix.boxplot(by="binding_group", grid=False)
    # plt.show()
    merge_matrix = pd.DataFrame(matrix["binding_rate"])
    merge_matrix["met_rate"] = matrix["met_rate"]
    merge_matrix["binding_group"] = matrix["binding_rate"].apply(lambda x: sort_by(bl, bm, bh, x))
    merge_matrix["met_group"] = matrix["met_var"].apply(lambda x: sort_by(ml, mm, mh, x))
    merge_matrix.boxplot(column="met_rate", by=["binding_group", "met_group"])
    plt.savefig("try save")
    plt.show()
    met_filter = [mh, 100]
    binding_filter = [0, bl, bm, bh, 1]
    for i in range(len(binding_filter) - 1):
        data =  matrix[(matrix["binding_rate"] > binding_filter[i]) & (matrix["binding_rate"] <= binding_filter[i+1])]
        for j in range(len(met_filter) - 1):
            d =  data[(data["met_var"] > met_filter[j]) & (data["met_var"] <= met_filter[j+1])]
            d.to_csv(DIFFERENT + os.sep + "binding_rate_between_{0:.2f}_to_{1:.2f}_methylation_variance_between_{2:.2f}_to_{3:.2f}.tsv".format(binding_filter[i], binding_filter[i + 1], met_filter[j], met_filter[j + 1]), sep='\t')
            print("i={0}, j={1}".format(i,j))




def compare_significant_sites(compare_to, num, significant_site):
    """
    A function that get a file that you want to compare and the file of significant sites
    and check if there overlap sites
    :param compare_to:
    :param num:
    :param significant_site:
    :return:
    """
    comp = pd.read_csv(compare_to, sep="\t")
    comp = comp[(comp["chr"] == 1) & (comp["p_values_{0}_month".format(num)] <= 0.1)]
    sg = pd.read_csv(significant_site, sep="\t")
    for chr, start, end in zip(sg[CHR], sg[START], sg[END]):
        s = ((((comp['start'] <= start) & (end <= comp['end'])) |
                                             ((start <= comp['start']) & (comp['end'] <= end)) |
                                             ((comp['start'] - THRESHOLD <= start) & (end <= comp['end']  + THRESHOLD)) )).any()
        if s:
            print("chr: {0}, start: {1}, end: {2}".format(chr, start, end))


def compare_at_significant(sg_file, title):
    sg = pd.read_csv(sg_file, sep="\t", index_col=0)
    sg = sg[sg["binding_rate"] != 1]
    col_name = list(sg.columns)[4:-6]
    bind = []
    unbind = []
    for sg_site in sg.iterrows():
        for cell in range(0, len(col_name), 2):
            if sg_site[1][col_name[cell + 1]] == 1:
                bind.append(sg_site[1][col_name[cell]])
            else:
                unbind.append(sg_site[1][col_name[cell]])
    plt.hist(bind, alpha=0.5, bins=50, density=True, stacked=True, label="methylation at the binded site")
    plt.hist(unbind, alpha=0.5, bins=50, density=True, stacked=True, label="methylation at the unbinded site")
    plt.yscale("log")
    plt.title(title)
    plt.legend()
    plt.show()





if __name__ == '__main__':
    # print("start runing")
    # # build_matrix()
    # first = True
    # matrix = None
    # for name, cell in cells_dict.items():
    #     if name == "A549":
    #         continue
    #     print("start ", name)
    #     if first:
    #         matrix = add_cell(cell[0], cell[1], name, False)
    #         first = False
    #     else:
    #         matrix = add_cell(cell[0], cell[1], name, False, True, matrix)
    #     print("end ", name)
    # print("end running")
    # add_cell(cells_dict["pancreas"][0], cells_dict["pancreas"][1], "pancreas", False)
    # play_with_data(MATRIX)
    # mann_witney_and_fun(MATRIX)
    # numbers = [6, 15, 10]
    # for num in numbers:
    #     print("for {0} month".format(num))
    #     compare_significant_sites("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/CSC/p                          
    # compare_at_significant("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/significant_sites_all_chr_p=0.05.tsv", "Methylation distribution at binding site with p value < 0.05")
    # compare_at_significant("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/not_significant_sites_all_chr_p=0.05.tsv", "Methylation distribution at binding site with p value >= 0.05")
    print("hi gal")
    # different_cuts(MATRIX)
    calculate_correlation(MATRIX)