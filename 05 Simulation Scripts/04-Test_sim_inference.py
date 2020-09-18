#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for testing inference on average ground truth RDM

@author: alex
"""

def collect_ROI_RDMs(n_subs, roi_h, ds_dir, beta_type, precision_type, perm, snr):
    roi_rdms = []
    for sub in range(1, n_subs+1): 
        rdm_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "rdms", "sub-"+
                               str(sub).zfill(2))
        rdm_filename = os.path.join(rdm_dir, "RDM_" + roi_h + "_" + precision_type
                                    + "_" + beta_type + "_" + str(perm).zfill(4)
                                    + "_snr_" + str(snr))
        rdm = pyrsa.rdm.rdms.load_rdm(rdm_filename, file_type='hdf5')
        
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
                          rdm_descriptors = {'ROI' : np.array([rdms.rdm_descriptors['ROI'][0] + '_mean']),
                                             'index' : np.array([rdms.rdm_descriptors['index'][0]])},
                          pattern_descriptors=rdms.pattern_descriptors)
    return average_rdm        


def collect_model_rdms(n_subs, ds_dir, beta_type, precision_type, perm,
                            snr, lures):
    model_rdms = []
    for i_model in lures:
        roi_rdms = collect_ROI_RDMs(n_subs, i_model, ds_dir, beta_type, 
                                    precision_type, perm, snr)
        lure_model = average_RDMs(roi_rdms)
        if isinstance(model_rdms, list):
            model_rdms = lure_model
        else:
            model_rdms.append(lure_model) 

    return model_rdms

def subset_rdms(model_rdms_full, data_rdms, n_stim, stim_sampling):
    model_rdms_list = []
    data_rdms_list = []
    factor_list = []
    n_cond = model_rdms_full.n_cond
    
    for stim in n_stim:
        for sampling in stim_sampling:
            if stim == n_cond and sampling == 'random':
                continue 
            parts = partition_sets(n_cond, stim, sampling = sampling)
            for part in parts:
                model_rdms_list.append(model_rdms_full.subset_pattern(
                     'index', part))
                data_rdms_list.append(data_rdms.subset_pattern(
                     'index', part))
                factor_list.append(str(stim) + '_' + str(sampling))
        
    return model_rdms_list, data_rdms_list, factor_list

def partition_sets(n_cond, stim, sampling = 'serial'):
    parts = []
    lst = [i for i in range(n_cond)] 
    if sampling == 'random':
        random.shuffle(lst)
        
    for i in range(0, n_cond, stim):
        parts.append(lst[i:i + stim])
        
    return parts


def results_summary(results_flex):
    results = results_flex.to_dict()
    evaluations = results['evaluations']
    
    # Names
    model_names_tmp = list(results['models'].keys())
    model_names = np.array([results['models'][i_model]['name'] for i_model in model_names_tmp])
    gt_model_idx = len(model_names)-1
    gt_model = model_names[gt_model_idx]
    
    # Determine winner model
    
    point_estimators = np.nanmean(evaluations, axis = (0,2))
    standard_deviations = np.nanstd(evaluations, axis = (0,2))
    # n = np.shape(evaluations)[0] * np.shape(evaluations)[2] 
    # standard_errors = standard_deviations / np.sqrt(n)
    best = np.max(point_estimators)
    
    winner_idx = np.where(point_estimators == best)
    winner_model = model_names[winner_idx]
    recovery = winner_model == gt_model
    
    # Significance testing
    p_values = pyrsa.util.inference_util.pair_tests(results['evaluations'])
    significance = fdr_control(p_values, alpha = 0.05)
    better = significance[gt_model_idx,:]
    
    # Determine Noise ceilings
    noise_ceilings = np.nanmean(results['noise_ceiling'], axis = 1)
    above_nc = best > noise_ceilings[0]
    

    # Putting everything together
    summary = {"Ground Truth": str(gt_model),
               "Point estimators": point_estimators,
               "Standard deviations": standard_deviations,
               "Winner estimator": best,
               "Winner standard deviation": standard_deviations[gt_model_idx],
               "Winner":str( winner_model[0]),
               "Recovered": recovery[0],
               "GT model comparisons": better,
               "Model comparison win count": sum(better),
               "Noise Ceilings": noise_ceilings,
               "Above lower NC": above_nc}
               
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
import copy
import mask_utils
import numpy as np
import pyrsa

# Set directories and specify ROIs
ds_dir          = "/home/alex/Datasets/ds001246/"
n_subs          = len(glob.glob(ds_dir + os.sep + "sub*"))
beta_type       = "data_perm_mixed"
#directory in which subject-specific volumetric ROI masks are saved by FS
freesurfer_mri  = "mri_glasser" 
mask_dir        = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" +
                               str(1).zfill(2), freesurfer_mri)
mask_dict       = mask_utils.load_dict(os.path.join(mask_dir, "sub-" +
                               str(1).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
roi_h_list      = list(mask_dict.keys())
precision_type  = 'instance-based' #opts: None, 'res-total', 'res-run-wise', 'instance-based'
snr_range       = [0.5, 1, 2]
n_perms         = 1
n_stim          = [5, 10, 25, 50]
stim_sampling   = ['serial', 'random']

###############################################################################
snr_summaries = []

# snr = 1
for snr in snr_range:
    perm_summaries = []
    
    # perm = 1
    for perm in range(1, n_perms+1):    
        summaries = []
        
        # roi_h = 'V1_left'
        for roi_h in roi_h_list:
            summary = {}            
            
            # Collect data RDMs and get GT by averaging data for resp. ROI
            other_rois = np.setdiff1d(roi_h_list, roi_h)
            data_rdms = collect_ROI_RDMs(
                n_subs, roi_h, ds_dir, beta_type, precision_type, perm, snr)
            ground_truth = average_RDMs(data_rdms)
            
            # Collect model RDMs (cross-subject averages for each ROI)
            model_rdms_full = collect_model_rdms(
                n_subs, ds_dir, beta_type, precision_type, perm, snr, other_rois)
            model_rdms_full.append(ground_truth)     
            
            # Subset model RDMs
            model_rdms_list, data_rdms_list, factor_list = subset_rdms(
                model_rdms_full, data_rdms, n_stim, stim_sampling)
            n_subsets = len(model_rdms_list)
            
            # Do flexible inference for each subset
            for comb in range(n_subsets):
                flexible_models = []
                results_flex = []
                model_rdms = model_rdms_list[comb]
                data_rdms = data_rdms_list[comb]
                
                # Model selection
                for i_model in np.unique(other_rois):
                    i_model += '_mean'
                    flexible_models.append(pyrsa.model.ModelSelect(i_model,
                        model_rdms.subset('ROI', i_model)))
                flexible_models.append(pyrsa.model.ModelSelect(
                    roi_h, model_rdms.subset('ROI', roi_h + '_mean'))) 
    
                # Perform flexible inference with bootstrapping
                if data_rdms.n_cond < 10:
                    k_pattern = 1
                else:
                    k_pattern = 2
                    
                results_flex = pyrsa.inference.bootstrap_crossval(
                    flexible_models, data_rdms, k_rdm=2, k_pattern=k_pattern,
                    method='corr', N=100)
                # pyrsa.vis.plot_model_comparison(results_flex)
                summary = results_summary(results_flex)
                summary.update({"RDM Subset" : factor_list[comb],
                                "Relative SNR" : snr,
                                "Permutation Num" : perm,
                                "CV pattern folds" : k_pattern,
                                "Precision Matrix" : precision_type})
                summaries.append(summary)
        perm_summaries.append(summaries) 
    snr_summaries.append(perm_summaries) 
            

                
                

