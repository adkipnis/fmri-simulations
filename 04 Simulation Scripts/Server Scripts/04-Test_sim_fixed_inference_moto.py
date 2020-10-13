#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for testing inference on average ground truth RDM

@author: alex
"""

def collect_RDMs(ds_dir, n_subs = 1,  beta_type = None):
    roi_rdms = []
    for sub in range(1, n_subs+1):
        fname = os.path.join(
            ds_dir, "derivatives", "PyRSA", "rdms", "sub-" +
            str(sub).zfill(2), "RDM_" + beta_type)
        rdm = pyrsa.rdm.rdms.load_rdm(fname, file_type='hdf5')
        
        # Collect single RDMs
        if isinstance(roi_rdms, list):
            roi_rdms = rdm
        else:
            roi_rdms.append(rdm)
    return roi_rdms


def average_RDMs(rdms):
    # Create average RDM        
    rdms_matrices = rdms.get_matrices()
    average_matrix = np.mean(rdms_matrices, axis=0, keepdims=True)
    average_rdm = pyrsa.rdm.rdms.RDMs(average_matrix,
                          dissimilarity_measure=rdms.dissimilarity_measure,
                          descriptors=rdms.descriptors,
                          rdm_descriptors = {'roi' : np.array([rdms.rdm_descriptors['roi'][0] + '_mean']),
                                             'index' : np.array([rdms.rdm_descriptors['index'][0]])},
                          pattern_descriptors=rdms.pattern_descriptors)
    return average_rdm        


def collect_model_rdms(rdms, other_rois, snr = None, perm = None, n_runs = None):
    model_rdms = []
    for other_roi in other_rois:
        other_rdms = rdms.subset('snr_rel', snr).subset(
            'perm', perm).subset('n_runs', n_runs).subset(
            'roi', other_roi)
        model_rdm = average_RDMs(other_rdms)
        if isinstance(model_rdms, list):
            model_rdms = model_rdm
        else:
            model_rdms.append(model_rdm)         
    return model_rdms


def pattern_subset_rdms_sparse(model_rdms_full, data_rdms, n_stim, stim_sampling):
    model_rdms_list = []
    data_rdms_list = []
    factor_list = []
    n_cond = model_rdms_full.n_cond
    
    for stim in n_stim:
        for sampling in stim_sampling:
            if stim == n_cond and sampling == 'random':
                continue 
            parts = partition_sets(n_cond, stim, sampling = sampling)
            model_rdms_list.append(model_rdms_full.subset_pattern(
                 'index', parts[0]))
            data_rdms_list.append(data_rdms.subset_pattern(
                 'index', parts[0]))
            factor_list.append(np.array([stim, sampling]))
    factors = np.stack(factor_list, axis=1).T
    return model_rdms_list, data_rdms_list, factors



def partition_sets(n_cond, stim, sampling = 'serial'):
    parts = []
    lst = [i for i in range(n_cond)] 
    if sampling == 'random':
        random.shuffle(lst)
        
    for i in range(0, n_cond, stim):
        parts.append(lst[i:i + stim])
        
    return parts


def results_summary(fixed_results):
    results = fixed_results.to_dict()
    evaluations = results['evaluations']
    variances = results['variances']    
    dof = results['dof']
    noise_ceiling = results['noise_ceiling']
    noise_ceil_var = results['noise_ceil_var']
    
    # Names
    model_names_tmp = list(results['models'].keys())
    model_names = np.array([results['models'][i_model]['name']
                            for i_model in model_names_tmp])
    model_names[-1] = model_names[-1] + '_mean'
    gt_model_idx = len(model_names)-1
    gt_model = model_names[gt_model_idx]
    
    # Determine winner model
    point_estimators = np.nanmean(evaluations, axis = 0)
    standard_deviations = np.nanstd(evaluations, axis = 0)
    best = np.max(point_estimators)
    winner_idx = np.where(point_estimators == best)
    winner_model = model_names[winner_idx]
    recovery = winner_model == gt_model
    
    # Significance testing
    p_values = pyrsa.util.inference_util.t_tests(evaluations, variances, dof)
    significance = fdr_control(p_values, alpha = 0.05)
    better = significance[gt_model_idx,:]
    
    # Noise ceiling tests
    noise_ceilings = np.nanmean(noise_ceiling, axis = 1)
    above_nc = best > noise_ceilings[0]
    p = pyrsa.util.inference_util.t_test_nc(
        evaluations, variances, noise_ceilings[0], noise_ceil_var[:, 0], dof)
    nc_significance = p[gt_model_idx]
    
    # Column names for output arrays
    # pe_names = np.core.defchararray.add(
        # model_names, np.repeat("_point_est", len(point_estimators)))
    # pe = dict(zip(pe_names, point_estimators))
    # std_names = np.core.defchararray.add(
        # model_names, np.repeat("_std", len(standard_deviations)))
    # stds = dict(zip(std_names, standard_deviations))
    # pw_sig_names = np.core.defchararray.add(
        # model_names, np.repeat("_vs_GT", len(better)))
    # pw_sigs = dict(zip(pw_sig_names[:-1], better[:-1]))
    
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
beta_type       = "data_perm_mixed"
#directory in which subject-specific volumetric ROI masks are saved by FS
freesurfer_mri  = "mri_glasser" 
mask_dict        = mask_utils.load_dict(
    os.path.join(ds_dir, "derivatives","Masks", "sub-" + str(1).zfill(2) 
                 + "_mask_dict_EPI_disjoint.npy"))
roi_h_list      = list(mask_dict.keys())
mask_dict       = None
n_stim          = [5, 10, 25, 50]
stim_sampling   = ['serial', 'random']

###############################################################################

rdms = collect_RDMs(
    ds_dir, n_subs = n_subs, beta_type = beta_type)
prec_types  = np.unique(rdms.rdm_descriptors['prec_type'])
run_subsets = np.unique(rdms.rdm_descriptors['n_runs'])
snr_range   = np.unique(rdms.rdm_descriptors['snr_rel'])
perms_range = np.unique(rdms.rdm_descriptors['perm'])


# perm = 1
for perm in perms_range:
    results_list = []
    df = pd.DataFrame()
    df_idx = -1
    
    # prec_type = 'instance_based'
    for prec_type in prec_types:
        # snr = 1
        for snr in snr_range:
            # roi_h = 'V1_left'
            for roi_h in roi_h_list:    
                # n_runs = 8
                for n_runs in run_subsets: 
                    summary = {}            
                    model_rdms_list, data_rdms_list = [], []
        
                    # Collect data RDMs and get GT by averaging data for resp. ROI
                    other_rois = np.setdiff1d(roi_h_list, roi_h)
                    data_rdms = rdms.subset('prec_type', prec_type).subset(
                        'snr_rel', snr).subset('perm', perm).subset(
                            'n_runs', n_runs).subset('roi', roi_h)
                    ground_truth = average_RDMs(data_rdms)
                    
                    # Collect model RDMs (cross-subject averages for each ROI)
                    model_rdms_full = collect_model_rdms(
                        rdms, other_rois, snr = snr, perm = perm, n_runs = n_runs)  
                    model_rdms_full.append(ground_truth)
                    
                    # Pattern subset model RDMs
                    model_rdms_list, data_rdms_list, factors = pattern_subset_rdms_sparse(
                        model_rdms_full, data_rdms, n_stim, stim_sampling)
                    n_subsets = len(model_rdms_list)
                    
                    # Do flexible inference for each subset
                    # comb = 0
                    for comb in range(n_subsets):
                        fixed_models = []
                        fixed_results = []
                        model_rdms = model_rdms_list[comb]
                        data_rdms_sub = data_rdms_list[comb]
                        
                        # Model selection
                        for i_model in np.unique(other_rois):
                            i_model += '_mean'
                            fixed_models.append(pyrsa.model.ModelFixed(i_model,
                                model_rdms.subset('roi', i_model)))
                        fixed_models.append(pyrsa.model.ModelFixed(
                            roi_h, model_rdms.subset('roi', roi_h + '_mean'))) 
                        
                        # ModelFixed class rearranges pattern descriptor indices, and
                        # we need to redo this for data_rdms_sub, otherwise
                        # eval_bootstrap_pattern will throw up an exception
                        data_rdms_sub.pattern_descriptors['index'] = np.arange(
                            data_rdms_sub.n_cond)
                        
                        # Perform flexible inference with bootstrapping
                        fixed_results = pyrsa.inference.eval_bootstrap_pattern(
                            fixed_models, data_rdms_sub, method='corr')
                        
                        # pyrsa.vis.plot_model_comparison(fixed_results)
                        summary = results_summary(fixed_results)
                        oe = dict(zip(["sub-" + str(i+1).zfill(2) + "_data_rdm_oe_rel"
                                       for i in range(data_rdms_sub.n_rdm)],
                                      data_rdms_sub.rdm_descriptors['oe_rel']))
                        roi_size = dict(zip(["sub-" + str(i+1).zfill(2) + "_GT_roi_size"
                                       for i in range(data_rdms_sub.n_rdm)],
                                      data_rdms_sub.rdm_descriptors['roi_size']))
                        for multi in [oe, roi_size]:
                            summary.update(multi) 
                        summary.update({"pattern_subset" : int(factors[comb, 0]),
                                        "pattern_subset_type" : factors[comb, 1],
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
