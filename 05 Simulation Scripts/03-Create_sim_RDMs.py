#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for generating pyrsa RDMs with precision matrices
for crossnobis distance estimates

@author: alex
"""
def sort_invert_and_numerate_dict(dictionary):
    inv_dict = {}
    i = 0
    values = list(dictionary.values())
    values.sort()
    for value in values:  
        inv_dict.update({value: i})
        i += 1
    return inv_dict

def get_n_perm(dataset_dir):
    files = glob.glob(dataset_dir + os.sep + "*perm*")
    p_num = [int(file.split('_')[-3]) for file in files]
    return np.unique(p_num)

def runwise_split_residuals(residuals, n_runs=35):
    res_per_run = int(len(residuals)/n_runs)
    run_splits = []
    
    # Split residual vectors into separate runs
    for i in range(n_runs):
        run_splits.append(residuals[i*res_per_run:(i+1)*res_per_run])
        
    return run_splits


def oe_split_residuals(residuals, n_runs=35):
    run_splits = runwise_split_residuals(residuals, n_runs=n_runs)
    
    # Split odd and even runs
    odd_residuals_list = run_splits[0::2]
    even_residuals_list = run_splits[1::2]
    
    # List to matrix
    odd_residuals = np.concatenate(odd_residuals_list, axis=0)
    even_residuals = np.concatenate(even_residuals_list, axis=0)
    
    return odd_residuals, even_residuals, odd_residuals_list, \
        even_residuals_list


def oe_split_reliability(dataset, residuals=None, l1_obs_desc='stim',
                         l2_obs_desc='run', n_runs=35, get_precision='res-total'):
    # Split measurements
    odd_dataset, even_dataset = pyrsa.data.dataset.nested_odd_even_split(
        dataset, l1_obs_desc, l2_obs_desc)
    
    # Split residuals and get precision matrices
    if not isinstance(residuals, np.ndarray):
        odd_precision = None
        even_precision = None
        get_precision = None
    else:
        odd_residuals, even_residuals, odd_residuals_list, even_residuals_list = \
            oe_split_residuals(residuals, n_runs = n_runs)
        odd_precision = calc_precision(
            dataset = odd_dataset, residuals = odd_residuals,
            get_precision = get_precision, n_runs = n_runs)
        even_precision = calc_precision(
            dataset = even_dataset, residuals = even_residuals,
            get_precision = get_precision, n_runs = n_runs)
        
    # Calculate respective rdms
    odd_rdm = pyrsa.rdm.calc.calc_rdm(odd_dataset, method='crossnobis',
                                      descriptor='stim', cv_descriptor='run',
                                      noise=odd_precision)
    even_rdm = pyrsa.rdm.calc.calc_rdm(even_dataset, method='crossnobis',
                                       descriptor='stim', cv_descriptor='run',
                                       noise=even_precision)
    
    # Calculate Pearson's product moment correlation coefficient
    # between vectorized rdms
    odd_vector = odd_rdm.get_vectors()
    even_vector = even_rdm.get_vectors()
    pearson_r = np.corrcoef(odd_vector, even_vector, rowvar=True)[0,1] 
    return pearson_r

def subset_dataset(dataset_full, n_runs):
    dataset_list = []
    
    for runs_goal in n_runs:
        run_list = [i+1 for i in range(runs_goal)]
        dataset_list.append(dataset_full.subset_obs('run', run_list))

    return dataset_list

def subset_residuals(residuals, n_runs, per_run = 1):
    residuals_list = []
    n_res = residuals.shape[0]
    n_runs_full = n_res / per_run
    assert n_runs_full.is_integer(), "odd number of residuals per run"
    
    for runs_goal in n_runs:
        residuals_list.append(residuals[0:runs_goal*per_run, :])
    
    return residuals_list

def calc_precision(dataset = None, residuals = None, get_precision = None,
                   obs_desc = 'stim', n_runs = 35):
    precision = None
    if get_precision == 'res-total':                            
        precision = pyrsa.data.noise.prec_from_residuals(
            residuals, dof=None)
    if get_precision == 'res-univariate':
        sample_cov, _ = pyrsa.data.noise.sample_covariance(residuals)
        precision = np.linalg.inv(sample_cov)
        precision = np.multiply(
               precision, np.identity(len(precision)))
    elif get_precision == 'res-run-wise':
        runwise_residuals = runwise_split_residuals(
            residuals, n_runs = n_runs)
        precision = pyrsa.data.noise.prec_from_residuals(
            runwise_residuals, dof=None)
    elif get_precision == 'instance-based':
        precision = pyrsa.data.noise.prec_from_measurements(
            dataset, obs_desc = obs_desc)
    return precision
       

###############################################################################

import glob
import os
import shutil
import mask_utils
import numpy as np
import pyrsa

# Set directories and specify ROIs
ds_dir           = "/home/alex/Datasets/ds001246/"
n_subs           = 1#len(glob.glob(ds_dir + os.sep + "sub*"))
beta_type        = "data_perm_mixed"
estimate_rel     = True
oe_reliabilities = []
precision_types  = [None, 'instance-based', 'res-total', 'res-univariate'] #opts: None, 'res-total', 'res-run-wise', 'instance-based', 'res-univariate'
calculate_rdm    = True
remove_ds        = False
freesurfer_mri   = "mri_glasser" #Name of the directory in which subject specific volumetric ROI masks are saved by FreeSurfer
mask_dir         = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" +
                               str(1).zfill(2), freesurfer_mri)
mask_dict        = mask_utils.load_dict(os.path.join(mask_dir, "sub-" +
                               str(1).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
roi_h_list       = list(mask_dict.keys())
label_dict       = np.load(os.path.join(ds_dir, "custom_synset_dictionary.npy"), allow_pickle='TRUE').item()
label2num        = sort_invert_and_numerate_dict(label_dict)
order            = ["n01443537", "n01943899", "n01976957", "n02071294",                              # water animals
                    "n01621127", "n01846331", "n01858441", "n01677366", "n02190790", "n02274259",    # air and land animals (non-mammals)
                    "n02128385", "n02139199", "n02416519", "n02437136", "n02437971",                 # land-Mammals
                    "n02951358", "n03272010", "n03482252", "n03495258",                              # humans in the picture
                    "n04254777", "n03237416", "n03124170", "n03379051", "n04572121",                 # clothing 
                    "n02824058", "n02882301", "n03345837", "n04387400", "n03716966", "n03584254", "n04533802", "n03626115", "n03941684", "n03954393", "n04507155", # small, handy objects 
                    "n02797295", "n02690373", "n02916179", "n02950256", "n03122295", "n04252077",    # machines
                    "n03064758", "n04210120", "n04554684", "n03452741", "n03761084",                 # large indoor objects
                    "n03710193", "n03455488", "n03767745", "n04297750"]                              # landmarks


p                = [] # Permutation vector
for i in range(len(order)):
    p.append(label2num[label_dict[order[i]]])
p                = np.array(p)
snr_range        = [0.5, 1, 2]
run_subsets      = [2**(i+1) for i in range(5)]

###############################################################################
#sub = 5
for sub in range(1, n_subs+1): 
    rdms = []
    
    # Set subject-specific paths
    dataset_dir = os.path.join(
        ds_dir, "derivatives", "PyRSA", "datasets", "sub-"+str(sub).zfill(2))
    res_dir = os.path.join(
        ds_dir, "derivatives", "PyRSA", "noise", "sub-"+str(sub).zfill(2)) 
    rdm_output_dir = os.path.join(
        ds_dir, "derivatives", "PyRSA", "rdms", "sub-"+str(sub).zfill(2))
    if not os.path.isdir(rdm_output_dir):
        os.makedirs(rdm_output_dir)
    
    perm_range = get_n_perm(dataset_dir)
    
    # Load datasets
    # snr = 1
    for snr in snr_range:
        # perm = 1
        for perm in perm_range:
            # get_precision = 'instance-based'
            for get_precision in precision_types:
                # roi_h = 'V1_left'
                for roi_h in roi_h_list:        
                    # Load dataset
                    dataset_filename = os.path.join(
                        dataset_dir, "RSA_dataset_" + roi_h + "_" + beta_type +
                        "_" + str(perm).zfill(4) + "_snr_" + str(snr))
                    dataset_full = pyrsa.data.dataset.load_dataset(
                        dataset_filename, file_type='hdf5')    
                    
                    # Load residuals
                    residuals_filename = os.path.join(
                        res_dir, "Residuals_" + roi_h + "_" + beta_type +
                        "_" + str(perm).zfill(4) + "_snr_" + str(snr)) + ".npy"   
                    residuals_full = np.load(residuals_filename)
                    
                    # Subset datasets
                    dataset_list = subset_dataset(dataset_full, run_subsets)
                    n_subsets = len(dataset_list)
                    residuals_list = subset_residuals(
                        residuals_full, run_subsets, per_run = 178)
                    
                    for ds_num in range(n_subsets):
                        dataset = dataset_list[ds_num]
                        n_runs = run_subsets[ds_num]
                        residuals = residuals_list[ds_num]
                        
                        # Precision matrix
                        precision = calc_precision(
                            dataset = dataset, residuals = residuals,
                            get_precision = get_precision, n_runs = n_runs)
                        
                        # Estimate odd-even reliability
                        oe_reliability = np.nan
                        if estimate_rel and n_runs > 2:
                            oe_reliability = oe_split_reliability(
                                dataset, residuals=precision, l1_obs_desc='stim',
                                l2_obs_desc='run', n_runs=n_runs, get_precision=get_precision)
                        
                        # Calculate RDM with crossnobis distance estimates    
                        if calculate_rdm:
                            rdm = pyrsa.rdm.calc.calc_rdm(
                                dataset, method='crossnobis', descriptor='stim',
                                cv_descriptor='run', noise=None)
                            rdm_p = pyrsa.rdm.rdms.permute_rdms(rdm, p = p)
                            rdm_p.rdm_descriptors = {'index': np.array([0]),
                                                     'sub': np.array([sub]),
                                                     'prec_type': np.array([get_precision]),
                                                     'roi': np.array([roi_h]),
                                                     'roi_size': np.array(
                                                         [dataset.channel_descriptors['positions'].shape[0]]),
                                                     'oe_rel': np.array([oe_reliability]),
                                                     'perm': np.array([perm]),
                                                     'snr_rel': np.array([snr]),
                                                     'n_runs': np.array([n_runs])}
                                
                            # Collect single RDMs
                            if isinstance(rdms, list):
                                rdms = rdm_p
                            else:
                                rdms.append(rdm_p)
    if remove_ds:
        shutil.rmtree(dataset_dir)
        shutil.rmtree(res_dir)
            
    if calculate_rdm:
        # Save subject RDM
        rdm_filename = os.path.join(rdm_output_dir, "RDM_" + beta_type)
        rdms.save(rdm_filename, file_type='hdf5', overwrite = True)
        print("Created subject RDM:", rdm_filename)
        
     