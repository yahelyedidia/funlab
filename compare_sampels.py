import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

REP_LIST = ['1_month', '6_month', '10_month', '15_month']

DIR = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/CSC"
FIRST_FILE = "control_vs_csc_after_1_month_rep{0}"
ALL_NAME = "csc_all_rep{0}.tsv"
CHANGES_REP = "csc_all_changes_rep{0}.tsv"


def create_one_file(target_file_name):
    combined = pd.read_csv(DIR + os.path.sep + FIRST_FILE, sep="\t")
    combined = combined.drop(columns=['cov', 'strand', 'change'])
    combined = combined.rename(columns={"treatment": "1_month"})
    # a = os.listdir(DIR)
    for file in os.listdir(DIR):
        if file != FIRST_FILE and file.endswith("2"):
            time_interval = file.split("_")[4]
            print(time_interval)
            temp = pd.read_csv(DIR + os.path.sep + file, sep="\t")
            combined["{0}_month".format(time_interval)] = temp["treatment"]
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


def plot_change(i):
    data = pd.read_csv(DIR + os.path.sep + CHANGES_REP.format(i), sep="\t")
    for col in REP_LIST:
        y = data[col]
        plt.hist(y, alpha=0.5, label=col)

    plt.xlabel("methylation change rate")
    plt.ylabel("Amount of performances")
    plt.title("methylation change distribution according to time, replication : {0}".format(i))
    plt.legend()
    plt.savefig(DIR + os.path.sep + "change_distribution_graph_rep{0}".format(i))
    plt.show()


def statistic_test(i):
    data = pd.read_csv(DIR + os.path.sep + CHANGES_REP.format(i), sep="\t")

    stat, p = stats.kruskal(data[REP_LIST[0]], data[REP_LIST[1]], data[REP_LIST[2]], data[REP_LIST[3]])
    print('stat=%.3f, p=%.3f' % (stat, p))
    if p > 0.05:
        print('Probably the same distribution')
    else:
        print('Probably different distributions')


if __name__ == '__main__':
    statistic_test(2)
