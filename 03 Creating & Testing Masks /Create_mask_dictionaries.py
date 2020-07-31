#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create subject-specific dictionaries containing binary volumetric masks for each ROI (based on FreeSurfer output)

@author: alex
"""

import os
import glob
from nibabel import nifti1
import numpy as np
import mask_utils #mask_utils.py must be in your working directory


# Set directories, specify ROIs and load dictionary for labels
ds_dir           = "/home/alex/Datasets/ds001246/"
txt_dir          = "/home/alex/Datasets/templateflow/tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt" #directory of mask descriptors
spm_dir          = os.path.join(ds_dir, "derivatives", 'SPM_0')
freesurfer_mri   = "mri_glasser" #Name of the directory in which subject specific volumetric ROI masks are saved by FreeSurfer
target_ROIs      = ['V%d' % i for i in range(1,5)] + ['VMV%d' % i for i in range(1,4)] + ['PHA%d' % i for i in range(1,4)] + ['VVC', 'FFC', 'TF', 'PeEc', 'MT', 'MST']
n_subs           = len(glob.glob(ds_dir + os.sep + "sub*"))

# Mask parameters
voxel_dimensions = {'T1w': [1,1,1], 'EPI': [3,3,3]}
fwhm             = np.array(voxel_dimensions['EPI']) # For mask smoohing prior to downsampling
mask_threshold   = 'adaptive' # Options: 'adaptive', a natural number (absolute threshold) or a real number between 0 and 1 
mask_merging     = True # Merge several small ROIs
rh_constant      = 200 #right hemisphere constant (roi_id in rh := roi_id in lh + rh_constant)
if mask_merging:
    merge_list   = [tuple(['PeEc', 'TF']), tuple('PHA%d' % i for i in range(1,4)), tuple('VMV%d' % i for i in range(1,4)), tuple(['MT', 'MST'])]
    merged_names = ['IT','PHA', 'VMV', 'MT+MST']

##############################################################################

for sub in range(1, n_subs+1):     
    
    # Set paths
    mask_dir = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" + str(sub).zfill(2), freesurfer_mri)
    glm_dir_example = os.path.join(spm_dir, "sub-"+str(sub).zfill(2), 'ses-perceptionTest' + '01', "run-01")
    betas_example = nifti1.load(glm_dir_example + os.sep + "beta_0001.nii")
    roi_ids = mask_utils.get_roi_ids_glasser(txt_dir, target_ROIs)
    
    # Create T1w atlas from freesurfer single ROI nifti files
    atlas_o = mask_utils.atlas_from_freesurfer_masks(sub, mask_dir, roi_ids, rh_constant = rh_constant, overwrite = True)
    atlas, roi_ids = mask_utils.merge_masks(atlas_o, roi_ids, merge_list, merged_names, rh_constant = rh_constant)  
    
    # Create mask dictionaries (T1w -> smooth -> resample -> threshold -> make disjoint) and get voxel IDs for the final one
    mask_dict_o, mask_dict_s, mask_dict_r, mask_dict_t = mask_utils.create_mask_dicts(atlas, betas_example, roi_ids, fwhm = fwhm, interpolation='nearest',
                                                 threshold = mask_threshold, voxel_dimensions = voxel_dimensions, rh_constant = rh_constant)
    mask_dict_d = mask_utils.remove_mask_overlap(mask_dict_r, mask_dict_t)
    voxel_ids_dict = mask_utils.get_voxel_ids_from_dict(mask_dict_d)
    
    # Save the final dictionary to mask directory
    mask_utils.save_dict(mask_dict_d, mask_dir, "sub-" + str(sub).zfill(2) + "_mask_dict_EPI_disjoint")
    mask_utils.save_dict(voxel_ids_dict, mask_dir, "sub-" + str(sub).zfill(2) + "_voxel_IDs_dict_EPI_disjoint")
    