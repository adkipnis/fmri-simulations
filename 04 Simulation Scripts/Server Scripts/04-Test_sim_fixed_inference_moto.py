#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for testing inference on average ground truth RDM

@author: alex
moto"""

def collect_RDMs(ds_dir, n_subs = 1, source = "PyRSA", beta_type = None):
    roi_rdms = []
    for sub in range(1, n_subs+1):
        fname = os.path.join(
            ds_dir, "derivatives", source, "rdms", "sub-" +
            str(sub).zfill(2), "RDM_" + beta_type)
        rdm = pyrsa.rdm.rdms.load_rdm(fname, file_type='hdf5')
        
        # Collect single RDMs
        if isinstance(roi_rdms, list):
            roi_rdms = rdm
        else:
            roi_rdms.append(rdm)
    return roi_rdms


def collect_model_rdms(rdms_full, rois, prec_type = 'none', method = None):
    model_rdms = []
    for roi in rois:
        rdms = rdms_full.subset('roi', roi).subset(
            'prec_type', prec_type)
        model_rdm = pyrsa.util.inference_util.pool_rdm(rdms, method=method)
        model_rdm.rdm_descriptors.update({'roi' : np.array([roi])})
        if isinstance(model_rdms, list):
            model_rdms = model_rdm
        else:
            model_rdms.append(model_rdm)         
    return model_rdms


def pattern_subset_rdms_sparse(model_rdms_full, data_rdms, n_stim):
    model_rdms_list = []
    data_rdms_list = []
    factor_list = []
    n_cond = model_rdms_full.n_cond
    
    for stim in n_stim:
        parts = partition_sets(n_cond, stim)
        model_rdms_list.append(model_rdms_full.subset_pattern(
             'index', parts[0]))
        data_rdms_list.append(data_rdms.subset_pattern(
             'index', parts[0]))
        factor_list.append(np.array([stim]))
    factors = np.stack(factor_list, axis=1).T
    return model_rdms_list, data_rdms_list, factors



def partition_sets(n_cond, stim, sampling = 'random'):
    parts = []
    lst = [i for i in range(n_cond)] 
    if sampling == 'random':
        random.shuffle(lst)
        
    for i in range(0, n_cond, stim):
        parts.append(lst[i:i + stim])
        
    return parts


def check_next_model_idx(winner_idx, model_names, gt_model, k):
    tie_winner = None
    if k == len(winner_idx):
        tie_winner = winner_idx[k]
    else:
        if model_names[winner_idx[k]][0] != gt_model:
            tie_winner = winner_idx[k]
        else:
            tie_winner = check_next_model_idx(winner_idx, model_names, gt_model, k+1)
    return tie_winner
            

def results_summary(fixed_results, roi_h):
    results = fixed_results.to_dict()
    evaluations = results['evaluations']
    if len(evaluations.shape) > 2:
        evaluations = np.squeeze(np.transpose(evaluations, (2, 1, 0)), axis=None)
    variances = results['variances']    
    dof = results['dof']
    noise_ceiling = results['noise_ceiling']
    noise_ceil_var = results['noise_ceil_var']
    
    # Names
    model_names_tmp = list(results['models'].keys())
    model_names = np.array([results['models'][i_model]['name']
                            for i_model in model_names_tmp])
    gt_model = roi_h
    gt_model_idx = np.where(model_names == roi_h)[0][0]
    
    # Determine winner model
    point_estimators = np.nanmean(evaluations, axis = 0)
    standard_deviations = np.nanstd(evaluations, axis = 0)
    best = np.max(point_estimators)
    winner_idx = np.where(point_estimators == best)
    
    if len(winner_idx) > 1: #handle ties
        winner_idx_tmp = check_next_model_idx(winner_idx, model_names, gt_model, 0)
        winner_idx = winner_idx_tmp
        
    winner_model = model_names[winner_idx]
    recovery = winner_model == gt_model
    
    # Significance testing
    p_values = pyrsa.util.inference_util.t_tests(evaluations, variances, dof)
    significance = fdr_control(p_values, alpha = 0.05)
    better = significance[gt_model_idx,:]
    
    # Noise ceiling tests
    if len(noise_ceiling.shape) > 1:
        noise_ceilings = np.nanmean(noise_ceiling, axis = 1)
    else:
        noise_ceilings = noise_ceiling
    above_nc = best > noise_ceilings[0]
    p = pyrsa.util.inference_util.t_test_nc(
        evaluations, variances, noise_ceilings[0], noise_ceil_var[:, 0], dof)
    nc_significance = p[gt_model_idx]
    
    
    # Putting everything together
    summary = {"GT": str(gt_model),
               "winner": str(winner_model[0]),
               "point_est": point_estimators[winner_idx],
               "std": standard_deviations[winner_idx],
               "recovered": int(recovery[0]),
               "n_sig_better": sum(better),
               "nc_low": noise_ceilings[0],
               "nc_high": noise_ceilings[1],
               "above_nc": int(above_nc),
               "dif_from_nc_sig": nc_significance}
    # for multi in [pe, stds, pw_sigs]:
    #     summary.update(multi)        
    return summary


def fdr_control(p_values, alpha = 0.05):
    ps = pyrsa.util.rdm_utils.batch_to_vectors(np.array([p_values]))[0][0]
    ps = np.sort(ps)
    criterion = alpha * (np.arange(ps.shape[0]) + 1) / ps.shape[0]
    k_ok = ps < criterion
    if np.any(k_ok):
        k_max = np.max(np.where(ps < criterion)[0])
        crit = criterion[k_max]
    else:
        crit = 0
    significant = p_values < crit
    
    return significant

###############################################################################

import glob
import os
import random
from time import gmtime, strftime
import mask_utils
import numpy as np
import pandas as pd
import pyrsa

# Set directories and specify ROIs
ds_dir          = '/moto/nklab/projects/ds001246/'
n_subs          = 5
#directory in which subject-specific volumetric ROI masks are saved by FS
freesurfer_mri  = "mri_glasser" 
mask_dir        = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" +
                               str(1).zfill(2), freesurfer_mri)
mask_dict       = mask_utils.load_dict(os.path.join(mask_dir, "sub-" +
                               str(1).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
roi_h_list      = list(mask_dict.keys())
mask_dict       = None
n_stim          = [5, 10, 20, 30, 50]
comp_methods    = ['cosine_cov', 'cosine']
###############################################################################
results_list = []
df = pd.DataFrame()
df_idx = -1

rdms = collect_RDMs(
    ds_dir, n_subs = n_subs, source = "PyRSA", beta_type = 'data_perm_mixed')
prec_types  = np.unique(rdms.rdm_descriptors['prec_type'])
run_subsets = np.unique(rdms.rdm_descriptors['n_runs'])
snr_range   = np.unique(rdms.rdm_descriptors['snr_rel'])
perms_range = np.unique(rdms.rdm_descriptors['perm'])
signal_rdms = collect_RDMs(
    ds_dir, n_subs = n_subs, source = "PyRSA_GT", beta_type = 'signal')


# method = 'cosine'
for method in comp_methods:
    # prec_type = 'instance-based'
    for prec_type in prec_types:
        # Collect full model RDMs (pooled according to precision type and comparison method)
        model_rdms_full = collect_model_rdms(signal_rdms, roi_h_list,
                                             prec_type = prec_type,
                                             method = method)
        # roi_h = 'V1_left'
        for roi_h in roi_h_list:        
            # snr = 1
            for snr in snr_range:
                # perm = 1
                for perm in perms_range:
                    # n_runs = 32
                    for n_runs in run_subsets: 
                        summary = {}            
                        model_rdms_list, data_rdms_list = [], []
            
                        # Collect data RDMs 
                        data_rdms = rdms.subset('prec_type', prec_type).subset(
                            'snr_rel', snr).subset('perm', perm).subset(
                                'n_runs', n_runs).subset('roi', roi_h)
    
                        # Pattern subset model RDMs
                        model_rdms_list, data_rdms_list, factors = pattern_subset_rdms_sparse(
                            model_rdms_full, data_rdms, n_stim)
                        n_subsets = len(model_rdms_list)
                        
                        # Do fixed inference for each subset
                        # comb = 0
                        for comb in range(n_subsets):
                            fixed_models = []
                            fixed_results = []
                            model_rdms = model_rdms_list[comb]
                            data_rdms_sub = data_rdms_list[comb]
                            
                            # Model selection
                            for i_model in roi_h_list:
                                fixed_models.append(pyrsa.model.ModelFixed(i_model,
                                    model_rdms.subset('roi', i_model)))
                            
                            # ModelFixed class rearranges pattern descriptor indices, and
                            # we need to redo this for data_rdms_sub, otherwise
                            # eval_bootstrap_pattern will throw up an exception
                            data_rdms_sub.pattern_descriptors['index'] = np.arange(
                                data_rdms_sub.n_cond)
                            
                            # Perform fixed inference 
                            if method == "cosine_cov":
                                fixed_results = pyrsa.inference.eval_fixed(
                                    fixed_models, data_rdms_sub, method=method)
                            else:   # with bootstrapping
                                fixed_results = pyrsa.inference.eval_bootstrap_pattern(
                                    fixed_models, data_rdms_sub, method=method)

                            
                            # pyrsa.vis.plot_model_comparison(fixed_results)
                            summary = results_summary(fixed_results, roi_h)
                            oe = dict(zip(["sub-" + str(i+1).zfill(2) + "_data_rdm_oe_rel"
                                           for i in range(data_rdms_sub.n_rdm)],
                                          data_rdms_sub.rdm_descriptors['oe_rel']))
                            roi_size = dict(zip(["sub-" + str(i+1).zfill(2) + "_GT_roi_size"
                                           for i in range(data_rdms_sub.n_rdm)],
                                          data_rdms_sub.rdm_descriptors['roi_size']))
                            for multi in [oe, roi_size]:
                                summary.update(multi) 
                            summary.update({"comparison_method" : method,
                                            "pattern_subset" : int(factors[comb]),
                                            "prec_type" : prec_type,
                                            "snr_rel" : snr,
                                            "perm_num" : perm,
                                            "n_runs" : n_runs})
                            df_idx += 1
                            summary_df = pd.DataFrame(summary, index=[df_idx])
                            df = df.append(summary_df)
                            results_list.append(fixed_results)
                        
    csv_fname = os.getcwd() + os.sep + "results_" + \
        strftime("%Y-%m-%d_%H-%M", gmtime()) + ".csv"           
    df.to_csv(csv_fname)
    npy_fname =  os.getcwd() + os.sep + "results_" + \
        strftime("%Y-%m-%d_%H-%M", gmtime()) + ".npy"   
    np.save(npy_fname, results_list)                
    # df_2 = pd.read_csv(csv_fname, index_col=0)
