#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline for generating pyrsa RDMs with precision matrices for crossnobis distance estimates

@author: alex
"""

def stim_to_label(stim_list, label_dict, crop=False):
    '''
    Transform beta_id to object labels
    
    Inputs:
    - stim_list (list): sorted list of all GLM predictors
    - label_dict (dict): {synset ID : object label}
    - crop (bool): Use only first term in object label before the comma

    Output:
    - label_list (list): sorted list of the stimulus object labels corresponding to GLM predictors
    '''
    label = []
    
    if crop:
        for i in range(len(stim_list)):
            synset = 'n0'+stim_list[i].split('.')[0]
            label.append(label_dict[synset].split(',')[0])
    else:
        for i in range(len(stim_list)):
            synset = 'n0'+stim_list[i].split('.')[0]    
            label.append(label_dict[synset])
    label_list = np.array(label)
    return label_list


# def permute_rdm(rdm, p = None):
#         """ Permute rows, columns and corresponding pattern descriptors of RDM matrices according to a permutation vector
        
#         Args:
#             p (numpy array):
#                permutation vector (values must be unique integers from 0 to n_cond of RDM matrix).
#                If p = None, a random permutation vector is created.
               
#         Returns:
#             rdm_p(pyrsa.rdm.RDMs): the rdm object with a permuted matrix and pattern descriptors

#         """
#         if p is None:
#             p = np.random.permutation(rdm.n_cond)
#             print("No permutation vector specified, performing random permutation.")
        
#         assert p.dtype == 'int', "permutation vector must have integer entries."
#         assert min(p) == 0 and max(p) == rdm.n_cond-1, "permutation vector must have entries ranging from 0 to n_cond"
#         assert len(np.unique(p)) == rdm.n_cond, "permutation vector must only have unique integer entries"
        
#         rdm_mats = rdm.get_matrices()
#         descriptors = rdm.descriptors.copy()
#         rdm_descriptors = rdm.rdm_descriptors.copy()
#         pattern_descriptors = rdm.pattern_descriptors.copy()

#         p_inv = np.arange(len(p))[np.argsort(p)] # To easily reverse permutation later
#         descriptors.update({'p_inv': p_inv})
#         rdm_mats = rdm_mats[:,p,:]
#         rdm_mats = rdm_mats[:,:,p]
#         stims = np.array(pattern_descriptors['stim'])
#         pattern_descriptors.update({'stim': list(stims[p].astype(np.str_))})
                
#         rdms_p = pyrsa.rdm.RDMs(dissimilarities=rdm_mats,
#                     descriptors=descriptors,
#                     rdm_descriptors=rdm_descriptors,
#                     pattern_descriptors=pattern_descriptors)
#         return rdms_p

###############################################################################

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pyrsa

# Set directories and specify ROIs
ds_dir = "/home/alex/Datasets/ds001246/"
txt_dir = "/home/alex/Datasets/templateflow/tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt" #directory of mask descriptors
mask_dir = os.path.join(ds_dir, "derivatives", "ROI_masks")
label_dict = np.load(os.path.join(ds_dir, "custom_synset_dictionary.npy"),allow_pickle='TRUE').item()
n_subs = len(glob.glob(ds_dir + os.sep + "sub*"))
beta_type = 'SPM_s' 
rename_stim = True
p = [0] # Permutation vector

sub = 1
for sub in range(1, n_subs+1): 
    # Set respective paths
    mask_dict_d_path = os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_mask_dict_T2*w_disjunct.npy")
    mask_dict_d = np.load(mask_dict_d_path,allow_pickle='TRUE').item()
    dataset_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "datasets", "sub-"+str(sub).zfill(2))
    res_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "noise", "sub-"+str(sub).zfill(2))
    
    # Load datasets and 
    roi_h = 'FFC_left'
    for roi_h in mask_dict_d.keys():        
       
        # Load residuals and estimate precision matrix
        residuals_filename = os.path.join(res_dir,"Residuals_"+roi_h+"_"+beta_type+".npy")
        residuals = np.load(residuals_filename)
        precision_matrix = pyrsa.data.noise.prec_from_residuals(residuals, dof=None)
        
        # Load dataset
        dataset_filename = os.path.join(dataset_dir,"RSA_dataset_"+roi_h+"_"+beta_type)
        dataset = pyrsa.data.dataset.load_dataset(dataset_filename, file_type='hdf5')    
        
        # Optionally rename stimulus IDs
        if rename_stim:
            stim_list = dataset.obs_descriptors['stim']
            label_list = stim_to_label(stim_list, label_dict, crop=True)
            dataset.obs_descriptors['stim'] = label_list
        
        # Calculate RDM with crossnobis distance estimates    
        rdm = pyrsa.rdm.calc.calc_rdm(dataset, method='crossnobis', descriptor='stim', cv_descriptor='run', noise=precision_matrix)
        rdm_p = pyrsa.rdm.rdms.permute_rdms(rdm, p = p)
        x = rdm_p.get_matrices()
        
        # TODO: Remember to pull request change in calc_rdm 
        pyrsa.vis.rdm_plot.show_rdm(rdm, do_rank_transform=False, pattern_descriptor='stim',
             cmap=None, rdm_descriptor=None)
        pyrsa.vis.rdm_plot.show_rdm(rdm_p, do_rank_transform=False, pattern_descriptor='stim',
             cmap=None, rdm_descriptor=None)
        
        
        
     