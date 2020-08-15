#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for pooling the residuals of the noise-regressors-only GLM
(block averaging approach)

@author: alex
"""

def get_path_to_tsv(run_dir, spm_type):
    tsv_dir = run_dir
    substring = "derivatives" + os.sep + spm_type + os.sep
    tsv_dir = tsv_dir.replace(substring, "", 1)
    return tsv_dir


def load_events(event_file_path, TR):
    events = pd.read_csv(event_file_path, sep="\t")
    events['stim_id'] = np.round(events['stim_id'], decimals = 0)
    events['stim_id'] = "n0"+events['stim_id'][1:-1].astype(int).astype(str)
    events['scan'] = (events['onset']/TR + 1).astype(int)
    return events


def get_residual_scans_dict(events, TR):
    residual_scans_dict = {}
    stim_ids = np.unique(events['stim_id'].to_list())
    stim_ids = np.delete(stim_ids, stim_ids == 'nan')
    for stim in stim_ids:
        stim_rows = events['stim_id']==stim
        durations = events[stim_rows].duration.to_list()
        scans_per_stim = sum(np.array(durations)/TR).astype(int)
        first_scan = events[stim_rows].scan.to_list()
        scans = [i for i in range(first_scan[0], first_scan[0]+scans_per_stim)]
        residual_scans_dict.update({stim: scans})
    return residual_scans_dict


def get_scan_numbers(event_file_path, TR):
    events = load_events(event_file_path, TR)
    residual_scans_dict = get_residual_scans_dict(events, TR)
    return residual_scans_dict


def load_and_average_block_residuals(glm_dir, scan_numbers,
                                     drop_first_res = False):
    if drop_first_res:
        scan_numbers = scan_numbers[1:]
    
    residual_superset = []
    for scan_num in scan_numbers:
        res_path = os.path.join(glm_dir, "Res_" + str(scan_num).zfill(4) +
                                ".nii")
        residual = nifti1.load(res_path)
        residual_array = residual.get_fdata()
        # residual_array = np.nan_to_num(residual_array)
        residual_superset.append(residual_array)
        
    mean_residual_array = np.mean(np.stack(residual_superset), axis=0)
    # mean_residual = nifti1.Nifti1Image(mean_residual_array,
    #                                    residual.affine.copy())
    return mean_residual_array, residual.affine.copy()


def pooled_residuals_pipeline(sub, spm_dir, spm_type, ses_type, TR,
                              drop_first_res):
    run_counter = 0
    fourth_dimension_descriptors = []
    mean_residual_array_superset = []
    n_ses = len(glob.glob(os.path.join(spm_dir, "sub-"+str(sub).zfill(2),
                                       ses_type+"*")))      
    
    for ses in range(1, n_ses+1):      
        run_dir = os.path.join(spm_dir, "sub-"+str(sub).zfill(2), ses_type
                               + str(ses).zfill(2))           
        n_runs = len(glob.glob(run_dir + os.sep + 'run*'))                      # relative number of runs   
        event_file_dir = get_path_to_tsv(run_dir, spm_type)
        
        for run in range(1, n_runs+1):
            run_counter += 1
            glm_dir = os.path.join(run_dir, "run-"+str(run).zfill(2))
            event_file_path = os.path.join(event_file_dir, "func", "sub-" + 
                             str(sub).zfill(2) + "_" + ses_type +
                             str(ses).zfill(2) + "_task-" + task + "_run-" +
                             str(run).zfill(2) +"_events.tsv")
            
            # Get dictionary with stim_id as key and a list of corresponding scan numbers as value
            residual_scans_dict = get_scan_numbers(event_file_path, TR)
            
            for condition in residual_scans_dict.keys():
                # load and average these scans, add them to superset
                mean_residual_array, _ = load_and_average_block_residuals(
                                         glm_dir, residual_scans_dict[condition],
                                         drop_first_res = drop_first_res)
                mean_residual_array_superset.append(mean_residual_array)
                fourth_dimension_descriptors.append(condition + "_" +
                                                    str(run_counter))    
    
    # Get affine matrix
    _, generic_affine = load_and_average_block_residuals(glm_dir, 
                        residual_scans_dict[condition],
                        drop_first_res = drop_first_res)
    # Pool the arrays
    pooled_residuals_array = np.stack(mean_residual_array_superset, axis=3)
    
    return pooled_residuals_array, generic_affine, fourth_dimension_descriptors

import os
import glob
import numpy as np
import pandas as pd
from nibabel import nifti1


# Data analysis parameters
processing_mode   = 'both' # Options: 'datasets', 'residuals' or 'both'
spm_type          = 'Res_0'
task              = 'perception'
stimulus_set      = 'Test'
ses_type          = 'ses-' + task + stimulus_set
drop_first_res    = True
save_dataset      = False

# Set directories, specify ROIs and load dictionary for labels
TR               = 3
ds_dir           = "/home/alex/Datasets/ds001246/"
spm_dir          = os.path.join(ds_dir, "derivatives", spm_type)
n_subs           = len(glob.glob(ds_dir + os.sep + "sub*"))

##############################################################################

for sub in range(1, n_subs+1):     
     
    img_output_dir = os.path.join(spm_dir, "sub-"+str(sub).zfill(2))  
    
    # For each run, get averaged condition-specific residuals 
    pooled_residuals_array, generic_affine, fourth_dimension_descriptors = \
        pooled_residuals_pipeline(sub, spm_dir, spm_type, ses_type, TR,
                                  drop_first_res)
        
    # Make subject-specific nifti image of block-averaged residuals
    pooled_residuals = nifti1.Nifti1Image(pooled_residuals_array,
                                          generic_affine)
    nifti_filename = os.path.join(img_output_dir, "sub-"+str(sub).zfill(2) +
                                  "_V_" + stimulus_set + ".nii.gz")
    nifti1.save(pooled_residuals, nifti_filename)
    csv_filename = os.path.join(img_output_dir, "sub-"+str(sub).zfill(2) +
                                  "_uq_conditions_" + stimulus_set + ".csv")
    df = pd.DataFrame({'descriptor':fourth_dimension_descriptors})
    df.to_csv(csv_filename, header=False)
       
            
            