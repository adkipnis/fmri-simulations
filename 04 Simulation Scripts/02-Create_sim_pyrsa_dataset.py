#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline and utility functions for generating simulated pyrsa datasets
for volumetric ROIs (from freesurfer)  
"""

def get_perm_range(glm_dir):
    files = glob.glob(glm_dir + os.sep + "*perm*")
    p_num = [int(file.split('_')[-4]) for file in files]
    return np.unique(p_num)

def stim_id_to_label(stim_ids, label_dict, n_stim=None, crop=False):
    '''
    Transform beta_id to object labels
    
    Inputs:
    - stim_ids (list / df): sorted list of all GLM predictors
    - n_stim: The first n_stim GLM predictors which represent
              stimulus identities of interest
    - label_dict (dict): {synset ID : object label}
    - crop (bool): Use only first term in object label before the comma

    Output:
    - labels (list): sorted list of the stimulus object labels
                     corresponding to GLM predictors
    '''
    labels = []
    if n_stim == None:
       n_stim = len(stim_ids) 
        
    if stim_ids['RegressorNames'][0][0:2] != 'n0':
        for i in range(n_stim):
            synset = 'n0'+stim_ids['RegressorNames'].\
                astype(str)[:n_stim][i].split('.')[0]
            if crop:
                labels.append(label_dict[synset].split(',')[0])
            else:
                labels.append(label_dict[synset])
    else:
        for i in range(n_stim):
            synset = stim_ids['RegressorNames'].astype(str)[:n_stim][i]
            if crop:
                labels.append(label_dict[synset].split(',')[0])
            else:
                labels.append(label_dict[synset])
    
    return labels 
    

def get_stim_desc(descriptors_path, label_dict = None):
     
    stim_ids = pd.read_csv(descriptors_path, sep='[, _]', header=None,
                           engine='python')
    stim_ids.columns = ['row', 'RegressorNames', 'run']
    runs = stim_ids['run'].to_list()

    if isinstance(label_dict, dict):
        labels = stim_id_to_label(stim_ids, label_dict, n_stim=None, crop=True)
    else:
        stim_ids['synset_id'] = 'n0'+stim_ids['synset_id'].astype(str)
        labels = stim_ids['RegressorNames'] .to_list()
    return labels, runs


def cronbachs_alpha(stim_measurements):
    cmatrix = np.corrcoef(stim_measurements, rowvar=True)
    N = len(cmatrix)
    rs = list(cmatrix[np.triu_indices(N, k=1)])
    mean_r = np.mean(rs)
    alpha = (N * mean_r) / (1 + (N - 1) * mean_r)
    return alpha

def batch_cronbachs_alpha(measurements, stimulus_descriptors):
    cronbach_alphas = []
    stimulus_descriptors = np.array(stimulus_descriptors)
    unique_stimuli = np.unique(stimulus_descriptors)
    for stimulus in unique_stimuli:
        stim_rows = np.where(stimulus_descriptors == stimulus)
        stim_measurements = measurements_cleaned[stim_rows]
        cronbach_alphas.append(cronbachs_alpha(stim_measurements))
        # plt.figure(dpi=600)
        # plt.matshow(stim_measurements)  
    return cronbach_alphas

def remove_empty_voxels(pooled_observations, voxel_ids = None,
                        nan_threshold = 0):
    '''
    After mask smoothing, resampling and thresholding,
    some mask voxels end up outside of the brain and need to be removed
    
    Args:
        pooled_observations (ndarray):
            n_obs resp. n_res x n_voxels (either measurements or residuals)
        voxel_ids (ndarray):
            n_voxel x 3 - dimensional array of the X,Y,Z coordinates
            of voxels that are contained in mask
        nan_threshold (float OR int):
            if between 0 and 1 - maximal percentage of allowed NaN entries
            if int > 1 - maximal number of allowed NaN entries 
   
    Returns:
        pooled_observations_cleaned (ndarray):
            n_obs resp. n_res x n_brain_voxels
        voxel_ids_cleaned (ndarray):
            n_brain_voxels x 3 - dimensional array of the X,Y,Z coordinates
            of voxels that are contained in mask
        bad_cols (list):
            indices of critical voxels
            
    '''
    assert nan_threshold >= 0, "nan_threshold must be a nonnegative number"
    
    obs_shape = np.shape(pooled_observations)
    
    if voxel_ids is not None:
        ROI_shape = np.shape(voxel_ids)
        assert obs_shape[1] == ROI_shape[0], "Dimension mismatch between data \
            array and number of voxels"
    
    if nan_threshold<1:
        max_nans = np.around(obs_shape[0]*nan_threshold).astype(int)
    else:
        max_nans = nan_threshold.astype(int)
        
    nan_rows = np.sum(np.isnan(pooled_observations).astype(int), axis = 0)
    bad_cols = np.where(nan_rows>max_nans)
    
    if len(bad_cols[0]) > 0:
        print("\nOut of", obs_shape[0], "observations",
              np.round(len(nan_rows)/obs_shape[0] * 100, 1),
              "% contained at least one NaN-value.")
        print("The following voxels contained at least",
              max_nans, "NaN-values:\n", bad_cols[0], "\n")
        print("These", len(bad_cols[0]), "voxels made up",
              np.round(len(bad_cols[0])/obs_shape[1] * 100, 1),
              "% of all voxels in this ROI and will be removed"+
              " from further analyses.\n")
        pooled_observations_cleaned = np.nan_to_num(
            np.delete(pooled_observations, bad_cols, axis = 1))
    else:
        pooled_observations_cleaned = pooled_observations
        
    if voxel_ids is not None:
        voxel_ids_cleaned = np.delete(voxel_ids, bad_cols, axis = 0)
    else:
        voxel_ids_cleaned = None
    
    return pooled_observations_cleaned, voxel_ids_cleaned, bad_cols[0]


def remove_empty_voxels_pipeline(measurements = None, residuals = None,
                                 voxel_ids = None, nan_threshold = 0):
    '''
    Pipeline for remove_empty_voxels for every processing_mode
    
    Args:
        measurements (ndarray):
            numpy array (pattern vector x ROI size) 
        residuals (ndarray):
            numpy array (residual vector x ROI size) 
        voxel_ids (ndarray):
            n_voxel x 3 - dimensional array of the X,Y,Z coordinates 
            of voxels that are contained in mask
        nan_threshold (float OR int):
            if between 0 and 1 - maximal percentage of allowed NaN entries
            if int > 1 - maximal number of allowed NaN entries 
   
    Returns:
        measurements_cleaned (ndarray):
            n_obs x n_brain_voxels
        residuals_cleaned (ndarray):
            n_res x n_brain_voxels
        voxel_ids_cleaned (ndarray):
            n_brain_voxels x 3 - dimensional array of the X,Y,Z coordinates
            of voxels that are contained in mask
            
    '''
    
    assert isinstance(measurements, np.ndarray) or isinstance(residuals, np.ndarray), \
        "Provide at least one nonempty array."
    
    measurements_cleaned = []
    residuals_cleaned = []
    voxel_ids_cleaned = []
    
    if isinstance(measurements, np.ndarray):
        measurements_cleaned, voxel_ids_cleaned, bad_cols_m = \
            remove_empty_voxels(measurements, voxel_ids = voxel_ids,
                                nan_threshold = nan_threshold)
    
    if isinstance(residuals, np.ndarray):
        if len(voxel_ids_cleaned) == 0:
             residuals_cleaned, voxel_ids_cleaned, bad_cols_r = \
                 remove_empty_voxels(residuals, voxel_ids = voxel_ids,
                                     nan_threshold = nan_threshold)
        else:
             residuals_cleaned, _, bad_cols_r = remove_empty_voxels(
                 residuals, voxel_ids = None, nan_threshold = nan_threshold)
    
    if isinstance(measurements, np.ndarray) and isinstance(residuals, np.ndarray) \
        and not np.array_equal(bad_cols_m, bad_cols_r):
        print("There was a mismatch between bad voxels in measurements" +
              "and residuals, bad voxel indices will be pooled")
        bad_cols = np.union1d(bad_cols_m, bad_cols_r).astype(int)
        measurements_cleaned = np.nan_to_num(np.delete(measurements,
                                                       bad_cols, axis = 1))
        residuals_cleaned = np.nan_to_num(np.delete(residuals,
                                                    bad_cols, axis = 1))
        voxel_ids_cleaned = np.delete(voxel_ids, bad_cols, axis = 0)
        
    return measurements_cleaned, residuals_cleaned, voxel_ids_cleaned


###############################################################################

import os
import glob
import shutil
import numpy as np
import pandas as pd
from nibabel import nifti1
import mask_utils
import pyrsa

# Data analysis parameters
processing_mode   = 'both' # Options: 'datasets', 'residuals' or 'both'
spm_type          = 'Data_perm'
task              = 'perception'
stimulus_set      = 'Test'
ses_type          = 'ses-' + task + stimulus_set
estimate_cronbach = False
cronbachs_list    = []
save_dataset      = True
snr_range         = [0.1, 1, 10]
delete_inputs     = False


# Set directories, specify ROIs and load dictionary for labels
ds_dir           = "/Users/heiko/fmri-simulations/ds001246/"
#directory of mask descriptors
txt_dir          = "/Users/heiko/fmri-simulations/tpl-Glasser/"+ \
                    "HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt"                    
spm_dir          = os.path.join(ds_dir, "derivatives", spm_type)
#directory in which subject-specific volumetric ROI masks are saved by FS
freesurfer_mri   = "mri_glasser"                                                
label_dict       = np.load(os.path.join(ds_dir, "custom_synset_dictionary.npy"),
                           allow_pickle='TRUE').item()
n_subs           = len(glob.glob(spm_dir + os.sep + "sub*"))


##############################################################################

for sub in range(1, n_subs+1):     
    
    # Set output directories    
    ds_output_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "datasets",
                                 "sub-"+str(sub).zfill(2))
    if not os.path.isdir(ds_output_dir):
        os.makedirs(ds_output_dir)
    res_output_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "noise",
                                  "sub-"+str(sub).zfill(2))
    if not os.path.isdir(res_output_dir):
        os.makedirs(res_output_dir)

    # Load mask dictionaries and 4d image
    mask_dir = os.path.join(
        ds_dir, "derivatives", "freesurfer","sub-" + 
        str(sub).zfill(2), freesurfer_mri)
    mask_dict = mask_utils.load_dict(
        os.path.join(mask_dir, "sub-" +str(sub).zfill(2)
                     + "_mask_dict_EPI_disjoint.npy"))
    voxel_ids_dict = mask_utils.load_dict(
        os.path.join(mask_dir, "sub-" + str(sub).zfill(2)
                     + "_voxel_IDs_dict_EPI_disjoint.npy"))
    voxel_ids_dict_tmp = voxel_ids_dict.copy()
    
    # Collect measurements as well as respective descriptors
    glm_dir = os.path.join(spm_dir, "sub-"+str(sub).zfill(2)) 
    
    perm_range = get_perm_range(glm_dir)
    
    for snr in snr_range:
        for perm in perm_range:
            if processing_mode in ['datasets', 'both']:
                signal_4d = None
                signal_path = os.path.join(glm_dir, "sub-"+str(sub).zfill(2)+ "_" +
                                           task +"_" + stimulus_set +
                                           "_data_perm_mixed_" + str(perm).zfill(4)+
                                           "_snr_" + str(snr) + "_signal.nii.gz")
                descriptors_path = os.path.join(glm_dir, "sub-"+str(sub).zfill(2)+
                                                "_" + task + "_" + stimulus_set +
                                                "_signal.csv")
                signal_4d = nifti1.load(signal_path)
            
            if processing_mode in ['residuals', 'both']:
                noise_4d = None
                noise_path = os.path.join(glm_dir, "sub-"+str(sub).zfill(2)+ "_" +
                                          task +"_" + stimulus_set +
                                           "_data_perm_mixed_" + str(perm).zfill(4)+
                                           "_snr_" + str(snr) + "_noise.nii.gz")
                noise_4d = nifti1.load(noise_path)
                   
            # Generate and save pyrsa dataset and pooled residuals for each ROI
            for roi_h in mask_dict.keys():
                measurements = None
                residuals = None
                
                if processing_mode in ['datasets', 'both']:
                    measurements = mask_utils.apply_roi_mask_4d(
                        signal_4d, mask_dict[roi_h], method = 'custom')
                    stimulus_descriptors, run_descriptors = get_stim_desc(
                        descriptors_path, label_dict = label_dict)
                
                if processing_mode in ['residuals', 'both']:
                    residuals = mask_utils.apply_roi_mask_4d(
                        noise_4d, mask_dict[roi_h], method = 'custom')
                
                
                # Remove voxels that are outside of the brain
                # (artifact of mask smoothing, resampling and thresholding)
                measurements_cleaned, residuals_cleaned, voxel_ids_cleaned = [],[],[]
                measurements_cleaned, residuals_cleaned, voxel_ids_cleaned = \
                    remove_empty_voxels_pipeline(measurements = measurements,
                                                 residuals = residuals,
                                                 voxel_ids = voxel_ids_dict[roi_h],
                                                 nan_threshold = 0.01)
                voxel_ids_dict_tmp[roi_h] = voxel_ids_cleaned
                measurements = None
                residuals = None
                
                # View measurements for one stimulus and get their internal consistency
                if estimate_cronbach and processing_mode in ['datasets', 'both']:
                    cronbach_alphas = batch_cronbachs_alpha(measurements_cleaned,
                                                            stimulus_descriptors)
                    cronbachs_list.append(cronbach_alphas)
                    print("For", roi_h, "the mean internal consistency of all"+
                          " stimuli is:", np.mean(cronbach_alphas))
                    
                if save_dataset:    
                    if processing_mode in ['datasets', 'both']:
                        # Create pyrsa dataset    
                        dataset = pyrsa.data.Dataset(measurements_cleaned,
                                  descriptors = {'ROI':roi_h}, 
                                  obs_descriptors = {'stim': stimulus_descriptors,
                                                     'run': run_descriptors},                             
                                  channel_descriptors = {'positions':
                                                         voxel_ids_dict_tmp[roi_h]})
                        
                        # Save dataset and residuals array
                        dataset_filename = os.path.join(ds_output_dir, 
                                            "RSA_dataset_" + roi_h +
                                            "_data_perm_mixed_" + str(perm).zfill(4)
                                            + "_snr_" + str(snr))
                        dataset.save(dataset_filename, file_type='hdf5',
                                     overwrite = True)
                        print("Created pyrsa dataset:", dataset_filename)
                        
                    if processing_mode in ['both', 'residuals']:
                        residuals_filename = os.path.join(res_output_dir,
                                              "Residuals_" + roi_h 
                                              + "_data_perm_mixed_"
                                              + str(perm).zfill(4) + "_snr_" +
                                              str(snr))
                        np.save(residuals_filename, residuals_cleaned)
                        print("Created residuals dataset:", residuals_filename)
                measurements_cleaned = None
                residuals_cleaned = None
    
    if delete_inputs:
        shutil.rmtree(glm_dir)
            