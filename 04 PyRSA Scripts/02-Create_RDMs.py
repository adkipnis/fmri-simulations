#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for generating pyrsa RDMs with precision matrices for crossnobis distance estimates

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
    
    return odd_residuals, even_residuals, odd_residuals_list, even_residuals_list


def oe_split_reliability(dataset, residuals, l1_obs_desc='stim', l2_obs_desc='run', n_runs=35, get_precision='total'):
    # Split measurements
    odd_dataset, even_dataset = pyrsa.data.dataset.nested_odd_even_split(dataset, l1_obs_desc, l2_obs_desc)
    
    # Split residuals and get precision matrices
    odd_residuals, even_residuals, odd_residuals_list, even_residuals_list = oe_split_residuals(residuals, n_runs=35)
    
    if get_precision == 'total':
        odd_precision = pyrsa.data.noise.prec_from_residuals(odd_residuals, dof=None)
        even_precision = pyrsa.data.noise.prec_from_residuals(even_residuals, dof=None)
    elif get_precision == 'run-wise':
        odd_precision = pyrsa.data.noise.prec_from_residuals(odd_residuals_list, dof=None)
        even_precision = pyrsa.data.noise.prec_from_residuals(even_residuals_list, dof=None)
           
    # Calculate respective rdms
    odd_rdm = pyrsa.rdm.calc.calc_rdm(odd_dataset, method='crossnobis', descriptor='stim', cv_descriptor='run', noise=odd_precision)
    even_rdm = pyrsa.rdm.calc.calc_rdm(even_dataset, method='crossnobis', descriptor='stim', cv_descriptor='run', noise=even_precision)
    
    # Calculate Pearson's product moment correlation coefficient between vectorized rdms
    odd_vector = odd_rdm.get_vectors()
    even_vector = even_rdm.get_vectors()
    pearson_r = np.corrcoef(odd_vector, even_vector, rowvar=True)[0,1]
    
    return pearson_r


###############################################################################

import glob
import os
import mask_utils
import numpy as np
import pyrsa

# Set directories and specify ROIs
ds_dir           = "/home/alex/Datasets/ds001246/"
n_subs           = len(glob.glob(ds_dir + os.sep + "sub*"))
beta_type        = 'SPM_6' 
estimate_rel     = True
oe_reliabilities = []
get_precision    = 'run-wise' #opts: 'run-wise', 'total'
calculate_rdm    = False
freesurfer_mri   = "mri_glasser" #Name of the directory in which subject specific volumetric ROI masks are saved by FreeSurfer
label_dict       = np.load(os.path.join(ds_dir, "custom_synset_dictionary.npy"),allow_pickle='TRUE').item()
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


p = [] # Permutation vector
for i in range(len(order)):
    p.append(label2num[label_dict[order[i]]])
p = np.array(p)


#sub = 1
for sub in range(1, n_subs+1): 
    
    # Set subject-specific paths
    mask_dir = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" + str(sub).zfill(2), freesurfer_mri)
    mask_dict = mask_utils.load_dict(os.path.join(mask_dir, "sub-" + str(sub).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
    dataset_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "datasets", "sub-"+str(sub).zfill(2))
    res_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "noise", "sub-"+str(sub).zfill(2)) 
    rdm_output_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "rdms", "sub-"+str(sub).zfill(2))
    if not os.path.isdir(rdm_output_dir):
        os.makedirs(rdm_output_dir)

    rdms = []
    
    # roi_h = 'FFC_right'
    # roi_h = 'V1_left'
    # Load datasets
    for roi_h in mask_dict.keys():        
        
        # Load residuals and estimate precision matrix
        residuals_filename = os.path.join(res_dir,"Residuals_"+roi_h+"_"+beta_type+".npy")
        residuals = np.load(residuals_filename)
        
        if get_precision == 'total':
            precision = pyrsa.data.noise.prec_from_residuals(residuals, dof=None)
        elif get_precision == 'run-wise':
            runwise_residuals = runwise_split_residuals(residuals, n_runs=35)
            precision = pyrsa.data.noise.prec_from_residuals(runwise_residuals, dof=None)
        
        
        # Load dataset
        dataset_filename = os.path.join(dataset_dir,"RSA_dataset_"+roi_h+"_"+beta_type)
        dataset = pyrsa.data.dataset.load_dataset(dataset_filename, file_type='hdf5')    
        
        # Estimate odd-even reliability
        if estimate_rel:
            oe_reliability = oe_split_reliability(dataset, residuals, l1_obs_desc='stim', l2_obs_desc='run', n_runs=35, get_precision=get_precision)
            oe_reliabilities.append(oe_reliability)
        if calculate_rdm:
            # Calculate RDM with crossnobis distance estimates    
            rdm = pyrsa.rdm.calc.calc_rdm(dataset, method='crossnobis', descriptor='stim', cv_descriptor='run', noise=precision)
            rdm_p = pyrsa.rdm.rdms.permute_rdms(rdm, p = p)
            rdm_p.rdm_descriptors = {'index':np.array([0]), 'ROI':np.array([roi_h])}
             
            # Save ROI RDM
            
            rdm_filename = os.path.join(rdm_output_dir,beta_type+"_RDM_"+roi_h)
            rdm_p.save(rdm_filename, file_type='hdf5')
            print("Created ROI RDM:", rdm_filename)
                
            # Collect single RDMs
            if isinstance(rdms, list):
                rdms = rdm_p
            else:
                rdms.append(rdm_p)
    
    if calculate_rdm:
        # Save subject RDM
        rdm_filename = os.path.join(rdm_output_dir,beta_type+"_RDM")
        rdms.save(rdm_filename, file_type='hdf5')
        print("Created subject RDM:", rdm_filename)
        
     