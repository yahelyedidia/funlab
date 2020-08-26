import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata


LOCATION = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/different_groups/genes/geneslist/biomart/10000/{0}/msigdb.txt"
output = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/different_groups/genes/geneslist/biomart/10000/compares.tsv"


def remove_first_row(files_dir):
    for folder in os.listdir(files_dir):
        for file in os.listdir(files_dir + os.path.sep + folder):
            if file.endswith(".txt"):
                with open(files_dir + os.path.sep + folder + os.path.sep + file, "r") as f:
                    rows = f.readlines()[1:]
                output = open(files_dir + os.path.sep + folder + os.path.sep + file, "w+")
                for line in rows:
                    output.write(line)
                output.close()
                print("done " + file)


def get_significant(file):
    data = pd.read_csv(file, sep='\t')
    # print(data.shape)
    data = data.sort_values(by=['logP'])
    reject, corrected_p = fdrcorrection(np.exp(data['logP']), 0.05)
    data['fdr'] = corrected_p
    # print(data.shape)
    # data = data[data['logP'] <= np.log(0.05)]
    data = data[data['fdr'] <= 0.05]
    # print(data.shape)
    # data.to_csv(LOCATION + os.path.sep + 'test.txt', sep="\t", index=False)
    # print(data.shape)
    # print("done")
    return data


def _ecdf(x):
    '''no frills empirical cdf used in fdrcorrection
    '''
    nobs = len(x)
    return np.arange(1, nobs + 1) / float(nobs)


def fdrcorrection(pvals, alpha=0.05):
    pvals = np.asarray(pvals)
    ecdffactor = _ecdf(pvals)
    reject = pvals <= ecdffactor * alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected > 1] = 1
    return reject, pvals_corrected


def compare_significant(names):
    data_frames = []
    path_set = set()
    for item in names:
        data = get_significant(LOCATION.format(item))
        data_frames.append(data)
        for row in data.iterrows():
            id = row[1]['TermID']
            path_set.add(id)
    z = np.zeros((len(data_frames), len(path_set)))
    pathes_data = pd.DataFrame(z, columns=path_set)
    for i in range(len(data_frames)):
        ids = data_frames[i]['TermID']
        for id in ids:
            pathes_data[id][i] = 1
    # pathes_data.to_csv(output, sep="\t")
    filtered = pathes_data.drop([col for col, val in pathes_data.sum().iteritems() if val > 1], axis=1)
    print(filtered.head())


if __name__ == '__main__':
    # get_significant(LOCATION + os.path.sep + "genes_close_to_sites_update_midhigh_binding_rateNhigh_met_rate_filter_50000/msigdb.txt")
    names = ["genes_close_to_sites_high_binding_rate_filter_10000",
             "genes_close_to_sites_low_binding_rate_filter_10000",
             "genes_close_to_sites_lowmid_binding_rate_filter_10000",
             "genes_close_to_sites_midhigh_binding_rate_filter_10000"]
    compare_significant(names)