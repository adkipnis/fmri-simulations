#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 20:15:10 2021

@author: alex
"""

import glob
import numpy as np
import pandas as pd

results_folder = "/home/alex/Downloads/cos_results/"

def generic_lables(results_folder = results_folder):
    csv_list = glob.glob(results_folder + "*.csv")
    csv_list.sort()
    full_csv = pd.read_csv(csv_list[0], index_col = 0)
    full_csv["nsf"] = 1/full_csv["snr_rel"]
    full_csv["n_stim"] = full_csv["pattern_subset"]
    full_csv["NN"] = full_csv["prec_type"].replace(
                    {'res-univariate': 'univariate', 'res-total': 'time-based'})
    perm_idx = full_csv['perm_num']==1
    labels = full_csv[perm_idx][["GT","comparison_method", "n_stim",
                                 "n_runs", "NN", "nsf"]]
    labels["GT"] = labels["GT"].replace(
                    {'_left': '-l', '_right': '-r'}, regex=True)
    labels.to_csv(results_folder+"labels.csv")
    print("Saved labels to", results_folder)
    return perm_idx

def extract_stats(results_folder = results_folder, perm_idx = None):
    npy_list = glob.glob(results_folder + "*.npy")
    npy_list.sort()
    means_superset = []
    stds_superset = []

    for npy_file in npy_list:
        means_list = []
        stds_list = []
        print("Loading", npy_file)
        full_npy_tmp = np.load(npy_file, allow_pickle=True)

        for i in range(len(full_npy_tmp)):
            means_list.append(np.nanmean(full_npy_tmp[i].evaluations, axis = 0))    
            stds_list.append(full_npy_tmp[i].variances)
        
        del full_npy_tmp
        means_stacked = np.stack(means_list, axis=0)
        stds_stacked = np.stack(stds_list, axis=0)
        
        if isinstance(perm_idx, pd.core.series.Series):
            means_superset.append(means_stacked[perm_idx])
            means_superset.append(means_stacked[~perm_idx])
            stds_superset.append(stds_stacked[perm_idx])
            stds_superset.append(stds_stacked[~perm_idx])
        else:
            means_superset.append(means_stacked[0:6000])
            means_superset.append(means_stacked[6000::])
            stds_superset.append(stds_stacked[0:6000])
            stds_superset.append(stds_stacked[6000::])
        
    means = np.stack(means_superset, axis=1)
    stds = np.stack(stds_superset, axis=1)
    np.save(results_folder + "means.npy", means)
    np.save(results_folder + "stds.npy", stds)


perm_idx = generic_lables()
extract_stats(perm_idx = perm_idx)

