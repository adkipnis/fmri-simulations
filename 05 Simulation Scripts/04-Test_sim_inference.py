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


def collect_fixed_models(n_subs, ds_dir, beta_type, precision_type, perm, snr,
                         lures):
    models = []
    for i_model in lures:
        roi_rdms = collect_ROI_RDMs(n_subs, i_model, ds_dir, beta_type, 
                                    precision_type, perm, snr)
        lure_model = average_RDMs(roi_rdms)
        m = pyrsa.model.ModelFixed(i_model, lure_model)
        models.append(m)

    return models


def collect_flexible_models(n_subs, ds_dir, beta_type, precision_type, perm,
                            snr, lures):
    selected_models = []
    for i_model in lures:
        roi_rdms = collect_ROI_RDMs(n_subs, i_model, ds_dir, beta_type, 
                                    precision_type, perm, snr)
        lure_model = average_RDMs(roi_rdms)
        m = pyrsa.model.ModelSelect(i_model, lure_model)
        selected_models.append(m) 

    return selected_models


def mean_noise_celings(n_subs, ds_dir, beta_type, roi_h_list, precision_type,
                       perm, snr):
    noise_min_sup, noise_max_sup = [], []
    for roi_h in roi_h_list:
        # Collect ROI-specific RDM from each subject and average them
        roi_rdms = collect_ROI_RDMs(n_subs, roi_h, ds_dir, beta_type,
                                    precision_type, perm, snr)
        noise_min, noise_max = pyrsa.inference.boot_noise_ceiling(roi_rdms,
                               method='cosine', rdm_descriptor='index')
        noise_min_sup.append(noise_min)
        noise_max_sup.append(noise_max)
    median_lower = np.median(noise_min_sup)
    median_upper = np.median(noise_max_sup)
    return median_lower, median_upper, noise_min_sup, noise_max_sup

###############################################################################

import glob
import os
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
inference_type  = 'cv' #opts: 'fixed', 'cv'
snr_range       = [0.5, 1, 2]
n_perms         = 1

###############################################################################
for snr in snr_range:
    for perm in range(1, n_perms+1):
        # Get median noise ceilings for all ROIs
        # median_lower, median_upper, noise_min_sup, noise_max_sup = mean_noise_celings(
        #     n_subs, ds_dir, beta_type, roi_h_list, precision_type, perm, snr)
        
        # Fixed inference
        if inference_type == 'fixed':
            for roi_h in roi_h_list:
                # roi_h will be used as our ground truth
                other_rois = np.setdiff1d(roi_h_list, roi_h)
                data_rdm = collect_ROI_RDMs(n_subs, roi_h, ds_dir, beta_type,
                                            precision_type, perm, snr)
                ground_truth = average_RDMs(data_rdm)
                
                # average the data rdms for the remaining ROIs and use them as
                # model RDMs
                fixed_models = collect_fixed_models(n_subs, ds_dir, beta_type,
                                                    precision_type, perm, snr,
                                                    other_rois)
                fixed_models.append(pyrsa.model.ModelFixed(roi_h, ground_truth)) 
                
                # Fixed inference with correlational distance + bootstrapping
                # over patterns
                results_fixed = pyrsa.inference.eval_bootstrap_pattern(
                    fixed_models, data_rdm, method='corr')
                pyrsa.vis.plot_model_comparison(results_fixed)
            
        # Flexible inference    
        if inference_type == 'cv':
            for roi_h in roi_h_list:
                results_flex = []
                other_rois = np.setdiff1d(roi_h_list, roi_h)
                data_rdms = collect_ROI_RDMs(
                    n_subs, roi_h, ds_dir, beta_type, precision_type, perm, snr)
                ground_truth = average_RDMs(data_rdms)
                
                # Collect flexible models (cross-subject averages for each ROI)
                flexible_models = collect_flexible_models(
                    n_subs, ds_dir, beta_type, precision_type, perm, snr, other_rois)
                
                
                # add the "ground truth" model (the data RDM for roi_h)
                flexible_models.append(pyrsa.model.ModelSelect(roi_h, ground_truth)) 
                results_flex = pyrsa.inference.bootstrap_crossval(
                    flexible_models, data_rdms, k_pattern=2, k_rdm=2, method='corr',
                    N=100)
                pyrsa.vis.plot_model_comparison(results_flex)
                
                results = results_flex.to_dict()
                point_estimators = np.nanmean(results['evaluations'],
                                              axis = (0,2)) 

