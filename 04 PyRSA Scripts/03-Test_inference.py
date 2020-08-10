#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for testing inference on average ground truth RDM

@author: alex
"""

def collect_ROI_RDMs(n_subs, roi_h, ds_dir, beta_type):
    roi_rdms = []
    for sub in range(1, n_subs+1): 
        rdm_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "rdms", "sub-"+str(sub).zfill(2))
        rdm_filename = os.path.join(rdm_dir, beta_type+"_RDM_"+roi_h)
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
                          pattern_descriptors=rdms.pattern_descriptors)
    return average_rdm        


def collect_fixed_models(n_subs, ds_dir, beta_type, lures):
    models = []
    for i_model in lures:
        roi_rdms = collect_ROI_RDMs(n_subs, i_model, ds_dir, beta_type)
        lure_model = average_RDMs(roi_rdms)
        m = pyrsa.model.ModelFixed(i_model, lure_model)
        models.append(m)

    return models

def mean_noise_celings(n_subs, ds_dir, beta_type, roi_h_list):
    noise_min_sup, noise_max_sup = [], []
    for roi_h in roi_h_list:
        # Collect ROI-specific RDM from each subject and average them
        roi_rdms = collect_ROI_RDMs(n_subs, roi_h, ds_dir, beta_type)
        noise_min, noise_max = pyrsa.inference.boot_noise_ceiling(roi_rdms, method='cosine', rdm_descriptor='index')
        noise_min_sup.append(noise_min)
        noise_max_sup.append(noise_max)
    mean_lower = np.mean(noise_min_sup)
    mean_upper = np.mean(noise_max_sup)
    return mean_lower, mean_upper

###############################################################################

import glob
import os
import mask_utils
import numpy as np
import pyrsa

# Set directories and specify ROIs
ds_dir          = "/home/alex/Datasets/ds001246/"
n_subs          = len(glob.glob(ds_dir + os.sep + "sub*"))
beta_type       = 'SPM_3' 
freesurfer_mri  = "mri_glasser" #Name of the directory in which subject specific volumetric ROI masks are saved by FreeSurfer
mask_dir        = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" + str(1).zfill(2), freesurfer_mri)
mask_dict       = mask_utils.load_dict(os.path.join(mask_dir, "sub-" + str(1).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
roi_h_list      = list(mask_dict.keys())

# Get mean noise ceilings for all ROIs
mean_lower, mean_upper = mean_noise_celings(n_subs, ds_dir, beta_type, roi_h_list)
   
# Setup correlational inference
# roi_h = np.random.choice(roi_h_list, 1)[0]
# print('\n', roi_h, "was selected as the source for our ground truth model.")

for roi_h in roi_h_list:
    other_rois = np.setdiff1d(roi_h_list, roi_h)
    # lures = np.random.choice(other_rois, 5)
    # print(lures, "were selected as alternative candidates.")
    data_rdm = collect_ROI_RDMs(n_subs, roi_h, ds_dir, beta_type)
    ground_truth = average_RDMs(data_rdm)
    
    # Perform correlational inference with bootstrapping over stimuli
    fixed_models = collect_fixed_models(n_subs, ds_dir, beta_type, other_rois)
    fixed_models.append(pyrsa.model.ModelFixed(roi_h, ground_truth)) # add the "ground truth"
    results_fixed = pyrsa.inference.eval_bootstrap_pattern(fixed_models, data_rdm, method='corr')
    pyrsa.vis.plot_model_comparison(results_fixed)
    
