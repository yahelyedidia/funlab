import os
import numpy as np
import pandas as pd
from scipy.stats import rankdata


LOCATION = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/different_groups/genes/geneslist/biomart/10000/{0}/msigdb.txt"
output = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/different_groups/genes/geneslist/biomart/10000/compares.tsv"
PATH = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/different_groups/genes/geneslist/biomart/10000/"
COL_NAMES = ['rank', 'motif_DB', 'motif_ID', 'motif_alt_ID', 'consensus', 'p-value', 'adj_p-value', 'E-value', 'tests', 'FASTA_max', 'pos', 'neg', 'PWM_min', 'TP', '%TP', 'FP', '%FP']
BAD_PATH_FILE = "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/different_groups/genes/geneslist/biomart/erelevant.txt"


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


def get_significant(file_name, fdr_thr=0.05, to_save=False):
    txt_file = open(BAD_PATH_FILE, "r")
    bad_lines = list(map(lambda x: x.strip("\n"), txt_file.readlines()))
    data = pd.read_csv(LOCATION.format(file_name), sep='\t')
    data = data.sort_values(by=['logP'])
    reject, corrected_p = fdrcorrection(np.exp(data['logP']), fdr_thr)
    data['fdr'] = corrected_p
    data = data[data['fdr'] <= fdr_thr]
    if to_save:
        data['bad_path'] = 0
        for i in range(data.shape[0]):
            id = data['TermID'][i]
            if id in bad_lines:
                pass
            else:
                data['bad_path'][i] = 1
        print(data.shape)
        data = data[data['bad_path'] == 1]
        print(data.shape)
        data.to_csv(PATH + os.path.sep + file_name + os.path.sep + 'after_fdr.txt', sep="\t", index=False)
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


def compare_significant(names, name, fltr, fdr):
    data_frames = []
    path_set = set()
    txt_file = open(BAD_PATH_FILE, "r")
    bad_lines = list(map(lambda x: x.strip("\n"), txt_file.readlines()))
    print(bad_lines)
    for item in names:
        data = get_significant(item, fdr)
        data_frames.append(data)
        for row in data.iterrows():
            id = row[1]['TermID']
            if id in bad_lines:
                pass
            else:
                path_set.add(id)
    z = np.zeros((len(data_frames), len(path_set)))
    pathes_data = pd.DataFrame(z, columns=path_set)
    for i in range(len(data_frames)):
        ids = data_frames[i]['TermID']
        for id in ids:
            if id in bad_lines:
                pass
            else:
                pathes_data[id][i] = 1
    print(pathes_data.shape)
    filtered = pathes_data.drop([col for col, val in pathes_data.sum().iteritems() if val > 1], axis=1)
    filtered.to_csv(PATH + os.sep + "{0}_pathways_filter_{1}_fdr_{2}.csv".format(name, fltr, fdr))
    print(filtered.shape)


def score_ctcf_sites(files_list):
    final = []
    for file in files_list:
        data = pd.read_csv(file, sep="\t")
        data = data[data['motif_ID'] == 'CTCF_HUMAN.H11MO.0.A']
        data['sorce'] = file
        final.append(data)
    for i in range(1, len(final)):
        final[0] = pd.concat([final[0], final[i]])
    data = final[0]
    data['p_ranked'] = data['adj_p-value'].shape[0] - rankdata(data['adj_p-value'].astype(float))
    data['TP_ranked'] = rankdata(data['%TP']).astype(float)
    data['total'] = data['p_ranked'] + data['TP_ranked']
    data = data.sort_values(by=['total'])
    print(data)
    return data
    # todo : save as file ?


if __name__ == '__main__':
    # get_significant(LOCATION + os.path.sep + "genes_close_to_sites_update_midhigh_binding_rateNhigh_met_rate_filter_50000/msigdb.txt")
    bind = ["genes_close_to_sites_high_binding_rate_filter_10000", "genes_close_to_sites_low_binding_rate_filter_10000",
            "genes_close_to_sites_lowmid_binding_rate_filter_10000", "genes_close_to_sites_midhigh_binding_rate_filter_10000"]
    vars = ["genes_close_to_sites_high_met_var_filter_10000", "genes_close_to_sites_low_met_var_filter_10000",
     "genes_close_to_sites_lowmid_met_var_filter_10000", "genes_close_to_sites_midhigh_met_var_filter_10000"]
    low_bind = ["genes_close_to_sites_low_binding_rateNlow_met_rate_filter_10000",
              "genes_close_to_sites_low_binding_rateNlowmid_met_rate_filter_10000",
              "genes_close_to_sites_low_binding_rateNmidhigh_met_rate_filter_10000",
              "genes_close_to_sites_update_low_binding_rateNhigh_met_rate_filter_10000"]
    low_mid_bind = ["genes_close_to_sites_lowmid_binding_rateNlow_met_rate_filter_10000",
                 "genes_close_to_sites_lowmid_binding_rateNlowmid_met_rate_filter_10000",
                 "genes_close_to_sites_lowmid_binding_rateNmidhigh_met_rate_filter_10000",
                 "genes_close_to_sites_update_lowmid_binding_rateNhigh_met_rate_filter_10000"]
    mid_high_bind = ["genes_close_to_sites_midhigh_binding_rateNlow_met_rate_filter_10000",
                     "genes_close_to_sites_midhigh_binding_rateNlowmid_met_rate_filter_10000",
                     "genes_close_to_sites_midhigh_binding_rateNmidhigh_met_rate_filter_10000",
                     "genes_close_to_sites_update_midhigh_binding_rateNhigh_met_rate_filter_10000"]
    high_bind = ["genes_close_to_sites_high_binding_rateNlow_met_rate_filter_10000",
                 "genes_close_to_sites_high_binding_rateNlowmid_met_rate_filter_10000",
                 "genes_close_to_sites_high_binding_rateNmidhigh_met_rate_filter_10000",
                 "genes_close_to_sites_update_high_binding_rateNhigh_met_rate_filter_10000"]
    correlation = ["neg_corr_filter_100000_genes_list", "pos_corr_filter_100000_genes_list", "low_corr_filter_100000_genes_list"]
    # compare_significant(bind, "bind", 10000, 0.05)
    for file in vars:
        get_significant(file, fdr_thr=0.01, to_save=True)
    # score_ctcf_sites(["/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/motives/bind_0.81-1.tsv",
    #                   "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/motives/bind_0.57-0.81.tsv",
    #                   "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/motives/bind_0.38-0.57.tsv",
    #                   "/vol/sci/bio/data/yotam.drier/Gal_and_Yahel/motives/bind_0-0.38.tsv"])
    print("done")
