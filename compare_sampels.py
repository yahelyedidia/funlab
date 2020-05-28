import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.spatial.distance import squareform
from scipy.spatial.distance import pdist




P_VALS = "p_values_all_information_by_orig_vals.tsv"

REP_LIST = ['1_month', '6_month', '10_month', '15_month']

DIR = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/CSC"
FIRST_FILE = "control_vs_csc_after_1_month_rep{0}"
ALL_NAME = "csc_all_rep{0}.tsv"
CHANGES_REP = "csc_all_changes_rep{0}.tsv"
ORIG_CHANGE = "all_changes_by_orig_data_rep_{0}.tsv"


def create_one_file(i, target_file_name, to_control=False):
    combined = pd.read_csv(DIR + os.path.sep + FIRST_FILE.format(i), sep="\t")
    combined = combined.drop(columns=['cov', 'strand', 'change'])
    if to_control:
        data = "change"
    else:
        data = "treatment"
    combined = combined.rename(columns={data: "1_month"})
    # a = os.listdir(DIR)
    for file in os.listdir(DIR):
        if file != FIRST_FILE and file.endswith(str(i)):
            time_interval = file.split("_")[4]
            print(time_interval)
            temp = pd.read_csv(DIR + os.path.sep + file, sep="\t")
            combined["{0}_month".format(time_interval)] = temp[data]
    combined = combined[['ID_REF',
                         'chr',
                         'start',
                         'end',
                         'control',
                         '1_month',
                         '6_month',
                         '10_month',
                         '15_month']]
    print(combined.head())
    combined.to_csv(DIR + os.path.sep + target_file_name, sep="\t", index=False)


def create_changes(k):
    combined = pd.read_csv(DIR + os.path.sep + FIRST_FILE.format(k), sep="\t")
    combined = combined.drop(columns=['cov', 'strand', 'treatment'])
    combined = combined.rename(columns={"change": "1_month"})
    data = pd.read_csv(DIR + os.path.sep + ALL_NAME.format(k), sep="\t")
    for i in range(1, len(REP_LIST)):
        combined[REP_LIST[i]] = data[REP_LIST[i-1]] - data[REP_LIST[i]]
        print(i)
    print(combined.head())
    combined.to_csv(DIR + os.path.sep + CHANGES_REP.format(k), sep="\t", index=False)


def get_uniq_rate(b, label):
    round = np.around(np.array(b[label]), 2)
    return np.unique(round, return_counts=True)


def plot_change(name):
    data = pd.read_csv(DIR + os.path.sep + name, sep="\t")
    xs_list, ys_list = [], []
    for col in REP_LIST:
        xs, ys = get_uniq_rate(data, col)
        xs_list += list(xs)
        ys_list += list(ys)
        y = data[col]
        plt.hist(y, alpha=0.5, label=col)
    # plt.scatter(xs_list, ys_list, color="black")
    # plt.plot(xs_list, ys_list, color="black")

    plt.xlabel("methylation change rate")
    plt.ylabel("Amount of performances")
    # plt.title("methylation change distribution according to time, replication : {0}".format(i))
    plt.title("compare between 1, 6 time point in both replictaions")
    plt.legend()
    # plt.savefig(DIR + os.path.sep + "change_distribution_graph_rep{0}".format(i))
    plt.show()


def statistic_test():
    # rep_1 = pd.read_csv(DIR + os.path.sep + ORIG_CHANGE.format(1), sep="\t")
    # rep_2 = pd.read_csv(DIR + os.path.sep + ORIG_CHANGE.format(2), sep="\t")
    all_vals_1 = pd.read_csv(DIR + os.path.sep + ALL_NAME.format(1), sep="\t")
    all_vals_2 = pd.read_csv(DIR + os.path.sep + ALL_NAME.format(2), sep="\t")

    result = all_vals_1
    result = result.drop(columns=REP_LIST)
    result = result.drop(columns=['control'])
    print("start to fill the data")
    for title in REP_LIST[1:]:
        c_list, a_list, p_vals = [], [], []
        for i in range(all_vals_1.shape[0]):
            first_month_1 = all_vals_1[REP_LIST[0]][i]
            first_month_2 = all_vals_2[REP_LIST[0]][i]
            control = [all_vals_1['control'][i], all_vals_2['control'][i], first_month_1, first_month_2]
            after = [all_vals_1[title][i], all_vals_2[title][i], all_vals_1[title][i] - first_month_1, all_vals_2[title][i] - first_month_2]
            t_test = stats.ttest_ind(control, after, equal_var=False)
            c_list.append(control)
            a_list.append(after)
            p_vals.append(t_test.pvalue)
        result["controls_{0}".format(title)] = c_list
        result["afters{0}".format(title)] = a_list
        result["p_values_{0}".format(title)] = p_vals
        print("done for " + title)
    print("done all")
    result.to_csv(DIR + os.path.sep + P_VALS, sep="\t", index=False)


def kde_scipy(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    # Note that scipy weights its bandwidth by the covariance of the
    # input data.  To make the results comparable to the other methods,
    # we divide the bandwidth by the sample standard deviation here.
    kde = stats.gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)


def smoothshow(i):
    data = pd.read_csv(DIR + os.path.sep + CHANGES_REP.format(i), sep="\t")
    data = data.drop(columns=['ID_REF', 'chr', 'start', 'end', 'control'])
    x = np.array(data)
    x_grid = np.linspace(-1, 1, 1000)
    pdf = kde_scipy(x, x_grid, bandwidth=0.2)
    plt.plot(x_grid, pdf, color='blue', alpha=0.5, lw=3)
    # plt.fill(x_grid, pdf_true, ec='gray', fc='gray', alpha=0.4)
    # plt.set_title(kde_funcnames[i])
    # plt.set_xlim(-4.5, 3.5)
    plt.show()


def data_to_plot():
    data = pd.read_csv(DIR + os.path.sep + CHANGES_REP.format(2), sep="\t")
    data = data.drop(columns=['ID_REF', 'chr', 'start', 'end', 'control'])
    data.to_csv("csc/changes_to_plot_1", index=False)


def compare_control_to_first_col(i, name):
    first = pd.read_csv(DIR + os.path.sep + FIRST_FILE.format(i), sep="\t")
    first = first.drop(columns=['cov', 'strand', 'change', 'treatment'])
    first = first.rename(columns={"control": "1_month"})
    # a= os.listdir(DIR)
    for file in os.listdir(DIR):
        if file != FIRST_FILE and file.endswith(str(i)):
            time_interval = file.split("_")[4]
            print(time_interval)
            temp = pd.read_csv(DIR + os.path.sep + file, sep="\t")
            first["{0}_month".format(time_interval)] = temp["control"]
    first = first[['ID_REF',
                         'chr',
                         'start',
                         'end',
                         '1_month',
                         '6_month',
                         '10_month',
                         '15_month']]
    print(first.head())
    first.to_csv(DIR + os.path.sep + name, sep="\t", index=False)


def compare_at_time():
    first = pd.read_csv(DIR + os.path.sep + ALL_NAME.format(1), sep="\t")
    first = first.drop(columns=['control', '10_month', '15_month'])
    second = pd.read_csv(DIR + os.path.sep + ALL_NAME.format(2), sep="\t")
    first['10_month'] = second['1_month']
    first['15_month'] = second['6_month']
    first.to_csv(DIR + os.path.sep + "compare_6_to_1.tsv", sep="\t", index=False)


def cov_1(rep):
    data = pd.read_csv(DIR + os.path.sep + ALL_NAME.format(rep), sep="\t")
    data = data.drop(columns=['ID_REF', 'chr', 'start', 'end'])
    pairwise = pd.DataFrame(squareform(pdist(data.loc[data.index])), columns=data.index, index=data.index)
    print("after create file")
    pairwise.to_csv("correlation csc.tsv", sep="\t")
    plt.matshow(pairwise)
    plt.savefig("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/correlation heat map by matshow")
    plt.close()
    plt.pcolor(pairwise)
    plt.savefig("/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/correlation heat map by pcolor")
    plt.show()


def mean_score():
    data = pd.read_csv("correlation csc.tsv", sep="\t")
    plt.pcolor(data)
    plt.show()


def divide_score():
    pass




if __name__ == '__main__':
    print("hi")
    mean_score()
    # compare_at_time()
    # plot_change("compare_6_to_1.tsv")
    # a = pd.read_csv(DIR + os.path.sep + P_VALS, sep="\t")
    # x = 2

