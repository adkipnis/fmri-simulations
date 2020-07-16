#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline and utility functions for generating pyrsa datasets and corresponding residuals for volumetric ROIs

@author: alex
"""


def beta_id_to_label(beta_ids, n_stim, label_dict, crop=False):
    '''
    Transform beta_id to object labels
    
    Inputs:
    - beta_ids (list / df): sorted list of all GLM predictors
    - n_stim: The first n_stim GLM predictors which represent stimulus identities of interest
    - label_dict (dict): {synset ID : object label}
    - crop (bool): Use only first term in object label before the comma

    Output:
    - labels (list): sorted list of the stimulus object labels corresponding to GLM predictors
    '''
    labels = []
    for i in range(n_stim):
        synset = 'n0'+beta_ids['RegressorNames'][:n_stim][i].split('.')[0]
        if crop:
            labels.append(label_dict[synset].split(',')[0])
        else:
            labels.append(label_dict[synset])
    
    return labels


def mask_dict_pipeline(t1w_mask_path, glm_dir_example, merge_list = None, merged_names = None, overwrite = False, threshold = 0):   
    # Load original atlas and betas
    atlas_o = nifti1.load(t1w_mask_path)
    betas_example = nifti1.load(glm_dir_example + os.sep + "beta_0001.nii")
    
    # Get ROI IDs and merge tuples of masks in T1w space
    roi_ids = mask_utils.get_roi_ids_glasser(txt_dir, target_ROIs)
    if merge_list is not None:
        atlas, roi_ids = mask_utils.merge_masks(atlas_o, roi_ids, merge_list, merged_names, rh_constant = 200)  
    else:
        atlas = atlas_o
    
    # Load or create mask dictionary, then make all masks disjunct
    mask_dict_d_path = os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_mask_dict_T2*w_disjunct.npy")
    if os.path.isfile(mask_dict_d_path) and not overwrite:
        # mask_dict_o = np.load(mask_dict_o_path,allow_pickle='TRUE').item()
        mask_dict_d = np.load(mask_dict_d_path,allow_pickle='TRUE').item()
        intersect_bad, intersect_good = mask_utils.test_roi_id_overlap(mask_dict_d, merge_list, merged_names)
        
        if len(intersect_bad)>0 or len(intersect_good)<len(merged_names):
            overwrite = True
            print("Loaded mask dictionary contains masks that still require merging. Let's create it anew...")
           
    if not os.path.isfile(mask_dict_d_path) or overwrite:
        mask_dict_o, mask_dict_s, mask_dict_r, mask_dict_t = mask_utils.create_mask_dict(atlas, betas_example, roi_ids, fwhm = fwhm, interpolation='nearest', threshold = threshold, rh_constant = 200)
        mask_dict_d = mask_utils.remove_mask_overlap(mask_dict_r, mask_dict_t)
        # np.save(mask_dict_o_path, mask_dict_o) 
        np.save(mask_dict_d_path, mask_dict_d) 
    return mask_dict_d

def superset_pipeline(spm_dir, voxel_ids_dict, mask_dict_d, roi_h, ses_type, n_ses=1, label_dict = None, processing_mode = 'both', mask_method = 'nilearn'):   
    measurements_superset = []
    residuals_superset = [] 
    session_desc_superset = []
    run_desc_superset = []
    stim_desc_superset = []
          
    for ses in range(1, n_ses+1):
        # get number of runs in this session    
        run_dir = os.path.join(spm_dir, "sub-"+str(sub).zfill(2), ses_type + str(ses).zfill(2))           
        n_runs = len(glob.glob(run_dir + os.sep + 'run*'))

        for run in range(1, n_runs+1):
            glm_dir = os.path.join(run_dir, "run-"+str(run).zfill(2))
            if processing_mode in ['both', 'datasets']:
                # 1. Apply mask to all beta coefficient images in this run     
                measurements = mask_utils.mask_and_get_SPM_measurements(mask_dict_d, roi_h, glm_dir, n_stim, method = mask_method)
                measurements_superset.append(measurements)
                
                # 2. Save stimulus ID's to stim_desc_superset
                beta_ids_dir = os.path.join(glm_dir,"sub-"+str(sub).zfill(2) + "-" + ses_type \
                                + str(ses).zfill(2) + "-run-" + str(run) + "_stim_ids.txt")
                beta_ids = pd.read_csv(beta_ids_dir, sep=" ")
                if type(label_dict) == dict:
                    labels = beta_id_to_label(beta_ids, n_stim, label_dict, crop=True)
                    stim_desc_superset.append(labels)
                else:
                    stim_desc_superset.append(beta_ids["RegressorNames"][:n_stim].to_list())
                        
    
            if processing_mode in ['both', 'residuals']:
                # 3. Apply mask to all residuals images in this run   
                residuals = mask_utils.mask_and_get_SPM_residuals(mask_dict_d, roi_h, glm_dir, n_res=n_res, method = mask_method)
                residuals_superset.append(residuals)
            
             
            # 4. If not already done in a previous run: save mask-specific voxel IDs
            if len(mask_dict_d.keys()) > len(voxel_ids_dict.keys()):
                voxel_ids_dict.update({roi_h: mask_utils.get_voxel_ids(mask_dict_d, roi_h)})
            run_desc = np.repeat(run, n_stim)
            run_desc_superset.append(run_desc)
                            
        session_desc = np.repeat(ses, n_stim*n_runs)
        session_desc_superset.append(session_desc)
    return measurements_superset, residuals_superset, session_desc_superset, run_desc_superset, stim_desc_superset, voxel_ids_dict


def remove_empty_voxels(pooled_observations, voxel_ids = None, nan_threshold = 0):
    '''
    After mask smoothing, resampling and thresholding, some mask voxels end up outside of the brain and need to be removed
    
    Args:
        pooled_observations (ndarray):
            n_obs resp. n_res x n_voxels (either measurements or residuals)
        voxel_ids (ndarray):
            n_voxel x 3 - dimensional array of the X,Y,Z coordinates of voxels that are contained in mask
        nan_threshold (float OR int):
            if between 0 and 1 - maximal percentage of allowed NaN entries
            if int > 1 - maximal number of allowed NaN entries 
   
    Returns:
        pooled_observations_cleaned (ndarray):
            n_obs resp. n_res x n_brain_voxels
        voxel_ids_cleaned (ndarray):
            n_brain_voxels x 3 - dimensional array of the X,Y,Z coordinates of voxels that are contained in mask
        bad_cols (list):
            indices of critical voxels
            
    '''
    assert nan_threshold >= 0, "nan_threshold must be a nonnegative number"
    
    obs_shape = np.shape(pooled_observations)
    
    if voxel_ids is not None:
        ROI_shape = np.shape(voxel_ids)
        assert obs_shape[1] == ROI_shape[0], "Dimension mismatch between data array and number of voxels"
    
    if nan_threshold<1:
        max_nans = np.around(obs_shape[0]*nan_threshold).astype(int)
    else:
        max_nans = nan_threshold.astype(int)
        
    nan_rows = np.sum(np.isnan(pooled_observations).astype(int), axis = 0)
    bad_cols = np.where(nan_rows>max_nans)
    
    if len(bad_cols[0]) > 0:
        print("Out of", obs_shape[0], "observations, more than", max_nans, "were NaN-values for the voxels:", bad_cols[0]) 
        pooled_observations_cleaned = np.nan_to_num(np.delete(pooled_observations, bad_cols, axis = 1))
    else:
        pooled_observations_cleaned = pooled_observations
        
    if voxel_ids is not None:
        voxel_ids_cleaned = np.delete(voxel_ids, bad_cols, axis = 0)
    else:
        voxel_ids_cleaned = None
    
    return pooled_observations_cleaned, voxel_ids_cleaned, bad_cols[0]

def remove_empty_voxels_pipeline(measurements_superset, residuals_superset, voxel_ids = None, nan_threshold = 0):
    '''
    Pipeline for remove_empty_voxels for every processing_mode
    
    Args:
        measurements_superset (list):
            list (n_ses * n_run) of measurement ndarrays (n_obs x n_voxels)
        residuals_superset (list):
            list (n_ses * n_run) of residual ndarrays (n_res x n_voxels)
        voxel_ids (ndarray):
            n_voxel x 3 - dimensional array of the X,Y,Z coordinates of voxels that are contained in mask
        nan_threshold (float OR int):
            if between 0 and 1 - maximal percentage of allowed NaN entries
            if int > 1 - maximal number of allowed NaN entries 
   
    Returns:
        measurements_cleaned (ndarray):
            n_obs x n_brain_voxels
        residuals_cleaned (ndarray):
            n_res x n_brain_voxels
        voxel_ids_cleaned (ndarray):
            n_brain_voxels x 3 - dimensional array of the X,Y,Z coordinates of voxels that are contained in mask
            
    '''
    
    assert len(measurements_superset) > 0 or len(residuals_superset) > 0, "Provide at least one nonempty superset"
    measurements_cleaned = []
    residuals_cleaned = []
    voxel_ids_cleaned = []
    
    if len(measurements_superset) > 0:
        measurements = np.vstack((measurements_superset[:]))
        measurements_cleaned, voxel_ids_cleaned, bad_cols_m = remove_empty_voxels(measurements, voxel_ids = voxel_ids, nan_threshold = nan_threshold)
    
    if len(residuals_superset) > 0:
        residuals = np.vstack((residuals_superset[:]))
        if len(voxel_ids_cleaned) == 0:
             residuals_cleaned, voxel_ids_cleaned, bad_cols_r = remove_empty_voxels(residuals, voxel_ids = voxel_ids, nan_threshold = nan_threshold)
        else:
             residuals_cleaned, _, bad_cols_r = remove_empty_voxels(residuals, voxel_ids = None, nan_threshold = nan_threshold)
    
    if len(measurements_superset) > 0 and len(residuals_superset) > 0 and not np.array_equal(bad_cols_m, bad_cols_r):
        print("There was a mismatch between bad voxels in measurements and residuals, bad voxel indices will be pooled")
        bad_cols = np.union1d(bad_cols_m, bad_cols_r).astype(int)
        measurements_cleaned = np.nan_to_num(np.delete(measurements, bad_cols, axis = 1))
        residuals_cleaned = np.nan_to_num(np.delete(residuals, bad_cols, axis = 1))
        voxel_ids_cleaned = np.delete(voxel_ids, bad_cols, axis = 0)
        
    return measurements_cleaned, residuals_cleaned, voxel_ids_cleaned

###############################################################################

import os
import glob
from nibabel import nifti1
import numpy as np
import pandas as pd
import mask_utils
import pyrsa

# Data analysis parameters
processing_mode = 'both' # alternatively: 'datasets' or 'both'
beta_type = 'SPM_s' # 'SPM' or 'SPM_s'
ses_type = 'perceptionTest'
n_stim = 50 # Use first n_stim beta coefficients

# Set directories, specify ROIs and load dictionary for labels
ds_dir = "/home/alex/Datasets/ds001246/"
txt_dir = "/home/alex/Datasets/templateflow/tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt" #directory of mask descriptors
spm_dir = os.path.join(ds_dir, "derivatives", beta_type)
mask_dir = os.path.join(ds_dir, "derivatives", "ROI_masks")
target_ROIs = ['V%d' % i for i in range(1,5)] + ['VMV%d' % i for i in range(1,4)] + ['PHA%d' % i for i in range(1,4)] + ['VVC', 'FFC', 'TF', 'PeEc', 'MT', 'MST']
label_dict = np.load(os.path.join(ds_dir, "stimulus_label_dictionary.npy"),allow_pickle='TRUE').item()
n_subs = len(glob.glob(ds_dir + os.sep + "sub*"))

# Mask parameters
fwhm = np.array([3,3,3]) # For mask smoohing (the functional EPI images had a voxel size of 3 × 3 × 3 mm)
mask_threshold = 0.4 # For mask thresholding
mask_merging = True
if mask_merging:
    merge_list = [tuple(['PeEc', 'TF']), tuple('PHA%d' % i for i in range(1,4)), tuple('VMV%d' % i for i in range(1,4)), tuple(['MT', 'MST'])]
    merged_names = ['IT','PHA', 'VMV', 'MT+MST']

##############################################################################

for sub in range(1, n_subs+1):     
    # Set respective paths to atlas, betas and residuals
    t1w_mask_path = os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_Mask_T1w_Glasser.nii.gz")
    glm_dir_example = os.path.join(spm_dir, "sub-"+str(sub).zfill(2), ses_type + '01', "run-01")
    ds_output_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "datasets", "sub-"+str(sub).zfill(2))
    if not os.path.isdir(ds_output_dir):
        os.makedirs(ds_output_dir)
    res_output_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "noise", "sub-"+str(sub).zfill(2))
    if not os.path.isdir(res_output_dir):
        os.makedirs(res_output_dir)
     
    # Infer number of sessions and run-wise residuals 
    n_ses = len(glob.glob(os.path.join(spm_dir, "sub-"+str(sub).zfill(2), ses_type+"*")))
    n_res = len(glob.glob(glm_dir_example + os.sep + "Res_*"))
    
    # Initialize and fetch dictionaries
    voxel_ids_dict = {} 
    mask_dict_d = mask_dict_pipeline(t1w_mask_path, glm_dir_example,
                                     merge_list = merge_list, merged_names = merged_names, overwrite = False, threshold = mask_threshold)
    
    # Generate and save pyrsa dataset and pooled residuals for each ROI
    for roi_h in mask_dict_d.keys():
        # Collect measurements, residuals or both, as well as respective descriptors for each run in each session
        measurements_superset, residuals_superset, session_desc_superset, \
            run_desc_superset, stim_desc_superset, voxel_ids_dict = \
                superset_pipeline(spm_dir, voxel_ids_dict, mask_dict_d, roi_h, ses_type, n_ses=n_ses,
                                  label_dict = None, processing_mode = processing_mode, mask_method = 'custom')
        
        # Remove voxels that are outside of the brain (artifact of mask smoothing, resampling and thresholding)
        measurements_cleaned, residuals_cleaned, voxel_ids_cleaned = \
            remove_empty_voxels_pipeline(measurements_superset, residuals_superset, voxel_ids = voxel_ids_dict[roi_h], nan_threshold = 0.01)
        voxel_ids_dict[roi_h] = voxel_ids_cleaned
        
        # Create pyrsa dataset    
        if processing_mode in ['both', 'datasets']:
            dataset = pyrsa.data.Dataset(measurements_cleaned,
                                 descriptors = {'ROI':roi_h}, 
                                 obs_descriptors = {'stim': np.hstack((stim_desc_superset[:])),
                                                    #'session': np.hstack((session_desc_superset[:])),
                                                    'run': np.hstack((run_desc_superset[:]))+ 100*np.hstack((session_desc_superset[:]))},                             
                                 channel_descriptors = {'positions':voxel_ids_dict[roi_h]})
            
            # Save dataset and residuals array
            dataset_filename = os.path.join(ds_output_dir,"RSA_dataset_"+roi_h+"_"+beta_type)
            dataset.save(dataset_filename, file_type='hdf5', overwrite=True)
            print("Created pyrsa dataset:", dataset_filename)
        
        # Save pooled residuals
        if processing_mode in ['both', 'residuals']:
            residuals_filename = os.path.join(res_output_dir,"Residuals_"+roi_h+"_"+beta_type)
            np.save(residuals_filename, residuals_cleaned)
            print("Created residuals dataset:", residuals_filename)
    
    # Save dictionary of voxel IDs    
    np.save(os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_mask-wise_voxel_ids.npy"), voxel_ids_dict)