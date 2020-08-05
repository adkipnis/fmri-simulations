#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for testing inference on average ground truth RDM

@author: alex
"""

def collect_and_average_RDMs(n_subs, roi_h, ds_dir, beta_type):
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
    
    average_rdm = average_RDM(roi_rdms, rdm)
    
    return average_rdm


def average_RDM(rdms, rdm):
    # Create average RDM        
    rdms_matrices = rdms.get_matrices()
    average_matrix = np.mean(rdms_matrices, axis=0, keepdims=True)
    average_rdm = pyrsa.rdm.rdms.RDMs(average_matrix,
                          dissimilarity_measure=rdm.dissimilarity_measure,
                          descriptors=rdm.descriptors,
                          rdm_descriptors=rdm.rdm_descriptors,
                          pattern_descriptors=rdm.pattern_descriptors)
    return average_rdm        


def collect_fixed_models(sub, ds_dir, beta_type):
    rdm_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "rdms", "sub-"+str(sub).zfill(2))
    models = []
    for i_model in candidate_rois:
        rdm_filename = os.path.join(rdm_dir, beta_type+"_RDM_"+i_model)
        m = pyrsa.model.ModelFixed(i_model, pyrsa.rdm.rdms.load_rdm(rdm_filename, file_type='hdf5'))
        models.append(m)
    return models
###############################################################################

import glob
import os
import mask_utils
import numpy as np
import pyrsa

# Set directories and specify ROIs
ds_dir          = "/home/alex/Datasets/ds001246/"
n_subs          = len(glob.glob(ds_dir + os.sep + "sub*"))
beta_type       = 'SPM_0' 
freesurfer_mri  = "mri_glasser" #Name of the directory in which subject specific volumetric ROI masks are saved by FreeSurfer
mask_dir        = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" + str(1).zfill(2), freesurfer_mri)
mask_dict       = mask_utils.load_dict(os.path.join(mask_dir, "sub-" + str(1).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
roi_h_list      = list(mask_dict.keys())
roi_h           = np.random.choice(roi_h_list, 1)[0]
other_rois      = np.setdiff1d(roi_h_list, roi_h)
print('\n', roi_h, "was selected as the source for our ground truth model.")

# Collect ROI-specific RDM from each subject and average them
average_rdm = collect_and_average_RDMs(n_subs, roi_h, ds_dir, beta_type)
   
# Setup correlational inference
lures = np.random.choice(other_rois, 5)
print(lures, "were selected as alternative candidates.")
candidate_rois = np.union1d(roi_h, lures)

# Perform correlational inference with bootstrapping over stimuli
for sub in range(1, n_subs+1): 
    fixed_models = collect_fixed_models(sub, ds_dir, beta_type)
    # Perform correlational inference
    # results_1 = pyrsa.inference.eval_fixed(fixed_models, average_rdm, method='corr')
    # pyrsa.vis.plot_model_comparison(results_1)
    results_2 = pyrsa.inference.eval_bootstrap_pattern(fixed_models, average_rdm, method='corr')
    pyrsa.vis.plot_model_comparison(results_2)
    
