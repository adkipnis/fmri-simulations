#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility funtions for creating ROI masks from an atlas in native space and applying them to fMRI images

@author: alex
"""

import os
import glob
from nibabel import nifti1
import numpy as np
import pandas as pd
from nilearn import image
from nilearn import masking


def get_roi_ids_glasser(txt_dir, target_ROIs):
    '''
    Load label txt file to find the index for each target ROI in target_ROIs

    Args:
        csv_dir (str):
            path to Glasser txt file
        target_ROIs (list):
            strings of target ROI names

    Returns:
        roi_ids (dict):
            {target_ROI: Left Hemisphere ROI index}
            Note: Right Hemisphere ROI index := Left Hemisphere ROI index + 200

    '''
    
    # Load table
    labels = pd.read_csv(txt_dir, sep=" ", header=None)
    labels.columns = ["ID", "ROI"]

    # Transform target ROIs
    targets = str()
    for roi in target_ROIs:
        targets += str(roi)+'_|'
    label_bool = labels['ROI'].str.contains(targets[:-1])
    roi_ids_tab = labels[label_bool].copy()
    roi_ids_tab['ROI'] = roi_ids_tab['ROI'].str.replace('L_', '').str.replace('_ROI', '')
    roi_ids = roi_ids_tab.set_index('ROI').T.to_dict('list')
    return roi_ids


def merge_masks(atlas_o, roi_ids, merge_list, merged_names, rh_constant = 200):
    '''
    Merge prespecified tuples of ROIs in original atlas for both hemispheres.
    This function assumes that IDs for right hemisphere ROIs are a function of IDs for left hemisphere ROIs

    Args:
        atlas_o (Nifti1Image):
            Original T1w Atlas with voxel values indicative of ROI identity
        roi_ids (dict):
            {target_ROI: Left Hemisphere ROI index}
        merge_list (list):
            List of tuples of ROI names. All ROIs within a Tuple will be merged.
        rh_constant (int):
            Right Hemisphere ROI index := Left Hemisphere ROI index + rh_constant
        
    Returns:
        atlas (Nifti1Image):
            T1w Atlas with voxel values indicative of ROI identity
        roi_ids_m (dict):
            {target_ROI: Left Hemisphere ROI index}, where specified target_ROI's are merged

    '''
    
    
    # Extract atlas and target nii data
    atlas_array = atlas_o.get_fdata()
    atlas_affine = atlas_o.affine.copy()
    s = 0
    
    for sublist in merge_list:
        # take first ROI in tuple as reference
        ref_id = roi_ids[sublist[0]][0]
        for other_roi in sublist[1:]:
            # Remap each remaining ROI to the reference ROI
            rest_id = roi_ids[other_roi][0]
            atlas_array[atlas_array==rest_id] = ref_id
            if rh_constant is not None:
                atlas_array[atlas_array==rest_id+rh_constant] = ref_id+rh_constant # repeat for right hemisphere
            roi_ids.pop(other_roi)
        roi_ids[merged_names[s]] = roi_ids.pop(sublist[0])
        s+=1
    atlas = nifti1.Nifti1Image(atlas_array, atlas_affine)
    roi_ids_m = roi_ids.copy()
    return atlas, roi_ids_m




def make_roi_mask(atlas, betas_example, roi_id, fwhm=None, interpolation='nearest'):
    '''
    Extract ROI-specific mask from Atlas in T1w space and sample it down to have
    the dimensions of the to-be-masked image of beta coefficients

    Args:
        atlas (Nifti1Image):
            Atlas with voxel values indicative of ROI identity
        betas_example (Nifti1Image):
            Example image of GLM coefficients which will be masked later (only needed for shape and affine)
        roi_id (int):
            ROI-specific integer used in "atlas"
        fwhm (float OR np.ndarray OR None):
            FWHM in mm (along each axis, axis-specific OR skip smoothing)
        interpolation (str):
            Interpolation method by nilearn.image.resample_img (consult nilearn documentation for other options)
    
    Returns:
        mask_nifti_r (Nifti1Image):
            resampled ROI-specific mask image (target dimensions)
        mask_nifti_s (Nifti1Image):
            smoothed ROI-specific mask image (original dimensions)
        mask_nifti (Nifti1Image):
            original binary ROI-specific mask image

    '''    

    # Extract atlas and target nii data
    atlas_array = atlas.get_fdata()
    atlas_affine = atlas.affine.copy()
    betas_array = betas_example.get_fdata()
    betas_affine = betas_example.affine.copy()

    # Extract ROI-specific mask and resample it
    mask = (atlas_array == roi_id).astype(float)    
    mask_nifti = nifti1.Nifti1Image(mask, atlas_affine)
    
    # Gaussian smoothing of mask
    if fwhm is not None:
        mask_nifti_s = image.smooth_img(mask_nifti, fwhm)
    else:
        mask_nifti_s = mask_nifti 
    
    # Resample mask
    mask_nifti_r = image.resample_img(mask_nifti_s, target_affine=betas_affine,
                                      target_shape=betas_array.shape, interpolation=interpolation, copy=True) 
    return mask_nifti_r, mask_nifti_s, mask_nifti

def threshold_mask(mask, threshold = 0):
    '''
    After smoothing, the downsampled masks are not binary anymore. Make a mask binary again.
    
    Args:
        mask (Nifti1Image):
            resampled ROI-specific mask image 
        threshold (float):
            between 0 and 1
    
    Returns:
        mask_t (Nifti1Image):
            binary resampled ROI-specific mask image after thresholding

    '''
    mask_array = mask.get_fdata()
    mask_thresholded = (mask_array >= threshold).astype(float)
    mask_t = nifti1.Nifti1Image(mask_thresholded, mask.affine.copy())
    return mask_t


def create_mask_dict(atlas, betas_example, roi_ids, fwhm = None, interpolation='nearest', threshold = 0, rh_constant = 200):
    '''
    Create a dictionary of mask nifti files for every target ROI
    
    Args:
        atlas (Nifti1Image):
            Atlas with voxel values indicative of ROI identity
        betas_example (Nifti1Image):
            Example image of GLM coefficients which will be masked later (only needed for shape and affine)
            Note: The created masks are not specific to the predictor whose weight is estimated by the beta coefficients
            - only the metadata (like shape) of 'betas' are needed
        roi_ids (dict): {target_ROIs_g:
                         Left Hemisphere ROI index}
        fwhm (float or np.ndarray):
            FWHM in mm (along each axis, axis specif or none)
        threshold (float):
            between 0 and 1 for thresholding after smoothing and resampling
        rh_constant (int OR None):
            Right Hemisphere ROI index := Left Hemisphere ROI index + rh_constant
            If 0 or None, this function assumes that all ROI indices are bilateral

    Returns:
        mask_dict_o (dict):
            {target_ROI+hemisphere : original mask}
        mask_dict_o (dict):
            {target_ROI+hemisphere : original mask after smoothing}
        mask_dict_r (dict):
            {target_ROI+hemisphere : smoothed and resampled mask}
        mask_dict_t (dict):
            {target_ROI+hemisphere : smoothed, resampled and thresholded mask}
    
    '''
    
    mask_dict_o = {}
    mask_dict_s = {}
    mask_dict_r = {}
    mask_dict_t = {}
    
    if rh_constant > 0:
        hemi_dict = {'left': 0, 'right':rh_constant}
    else:
        hemi_dict = {'bilateral': 0}
    
    for roi in roi_ids.keys(): 
        for side in hemi_dict.keys():
            roi_id = roi_ids[roi][0] + hemi_dict[side]
            mask_nifti_r, mask_nifti_s , mask_nifti_o = make_roi_mask(atlas, betas_example, roi_id,
                                                                      fwhm=fwhm, interpolation=interpolation)
            mask_nifti_t = threshold_mask(mask_nifti_r, threshold=threshold)
            mask_dict_o.update({roi+"_"+side : mask_nifti_o})
            mask_dict_s.update({roi+"_"+side : mask_nifti_s})
            mask_dict_r.update({roi+"_"+side : mask_nifti_r})
            mask_dict_t.update({roi+"_"+side : mask_nifti_t})

    return mask_dict_o, mask_dict_s, mask_dict_r, mask_dict_t


def test_roi_id_overlap(mask_dict, merge_list, merged_names):
    '''
    Test if to-be-merged ROIs are in the keys of the loaded mask dictionary
    
    Args:
        mask_dict (dict):
            {target_ROI+hemisphere : mask} (preferably mask_dict_d)
        merge_list (list):
            List of tuples of ROI names. All ROIs within a Tuple will be merged.
        merged_names (list):
            List of ROI names for each merged tuple
    
    Returns:
        intersect_bad (array):
            List of ROI names in merge_list that occur in mask_dict
        intersect_bad (array):
            List of ROI names in merged_names that occur in mask_dict
    '''
    intersect = []
    key_list = list(mask_dict.keys())
    substring = '_'+key_list[0].split('_')[1]
    
    # Flatten merge_list and append each item with an occurring substring in the keylist
    test_list = []
    for i in range(len(merge_list)): 
        test_list.append(list(merge_list[i]))
    test_list = [item for sublist in test_list for item in sublist]    
    test_list = [s + substring for s in test_list]
    should_be_insed_list = [s + substring for s in merged_names]
    
    # Intersection list
    intersect_bad = np.intersect1d(test_list, key_list)
    intersect_good = np.intersect1d(should_be_insed_list, key_list)
    return intersect_bad, intersect_good


def remove_mask_overlap(mask_dict_r, mask_dict_t):
    '''
    Find overlapping voxels in thresholded masks and only keep the voxel with highest pre-threshold voxel intensity
    
    Args:
        mask_dict_r (dict):
            {ROI_name : smoothed and resampled mask}
        mask_dict_t (dict):
            The above after thresholding

    Returns:
        mask_dict_d (dict):
            mask_dict_t after removal of overlapping voxels
    '''
    macro_mask = []
    mask_dict_d = []
    overlap_indices = []
    mask_dict_d = mask_dict_t.copy()
    
   # Add up all thresholded masks 
    for roi_h in mask_dict_t.keys():
        mask_array = mask_dict_t[roi_h].get_fdata()
   
        if len(macro_mask) == 0:
            macro_mask = mask_array.copy()
        else:
            macro_mask += mask_array.copy()
            
    # Find voxels with mask overlap
    overlap_indices = np.array(np.where(macro_mask > 1)).T

    # For each overlapping voxel: Create list of overlapping ROIs 
    for idx_number in range(len(overlap_indices)):
        idx = tuple(overlap_indices[idx_number])
        overlap_list = []
        voxel_intensities = []
              
        for roi_h in mask_dict_r.keys():
            voxel = mask_dict_d[roi_h].get_fdata()[idx]
            if voxel == 1:
                overlap_list.append(roi_h)
        print("The ROIs", overlap_list, " overlap on voxel: ", idx)
        
        # Check which voxel value is largest on smoothed versions of the overlapping masks 
        for overlap_mask in overlap_list:
            voxel_intensities.append(mask_dict_r[overlap_mask].get_fdata()[idx])
        most_intense = np.amax(voxel_intensities)
        most_intense_idx = np.where(voxel_intensities==most_intense)
        delete_idx = most_intense_idx[0][0]
        
        # Remove the ROI with higher voxel intensity from overlap list and set the voxel to 0 for all remaining entries
        del(overlap_list[delete_idx])
        print("Setting said voxel's intensity to 0 for the mask of: ", overlap_list) 
        for roi_d in overlap_list:
            mask = mask_dict_d[roi_d]
            mask_array = mask.get_fdata()
            mask_array[idx] = 0
            mask_nifti = nifti1.Nifti1Image(mask_array, mask.affine.copy())
            mask_dict_d[roi_d] = mask_nifti
    
    return mask_dict_d    


def apply_roi_mask_SPM(glm_dir, img_num, mask, target = 'betas', method = 'nilearn'):
    '''
    Load 3-D image for prespecified beta-coefficient or residual and apply mask to get a pattern vector
    This function assumes that beta and residual images are in the same directory and follow the SPM naming convention
    
    Args:
        glm_dir (str):
            Path to SPM outputs
        img_num (int):
            Stimulus number
        mask (Nifti1Image):
            To be applied mask (must have same dimensions as the GLM betas image)
        method (str):
            Method for applying the mask 

    Returns:
        beta_vector (1-D array OR Nifti Image):
            Vector of Beta-coefficients for each voxel in mask OR Nifti Image after masking
    '''
    # Check args
    if target not in ['betas', 'residuals']:
        print(target, "is not a valid SPM output descriptor. Please use 'betas' or 'residuals'")
    if method not in ['nilearn', 'custom', '3d']:
        print(method, "is not a valid method descriptor. Please use 'nilearn', 'custom', or '3d'")
    
    # Load image
    if target == 'betas':
        betas_path = os.path.join(glm_dir,"beta_" + str(img_num).zfill(4) + ".nii")
        image = nifti1.load(betas_path)
    elif target == 'residuals':
        residuals_path = os.path.join(glm_dir,"Res_" + str(img_num).zfill(4) + ".nii")
        image = nifti1.load(residuals_path)
    
    # Apply mask
    if method == 'nilearn':
        masked = masking.apply_mask(image, mask, dtype='f', smoothing_fwhm=None, ensure_finite=True)            
    elif method == 'custom':
    # Alternative to nilearn's implementation (allows for plotting of masked 3d image)
        mask_array = mask.get_fdata()
        image_array = image.get_fdata()
        masked = image_array[mask_array.astype(bool)]
    elif method == '3d':
        mask_array = mask.get_fdata()
        image_array = image.get_fdata()
        image_masked = np.multiply(mask_array, image_array)
        masked = nifti1.Nifti1Image(image_masked, mask.affine.copy())
    
    return masked


def mask_and_get_SPM_measurements(mask_dict_d, roi_h, glm_dir, n_stim=None, method = 'nilearn'):
    '''
    Apply mask to all beta images of interest and store them in a measurements array (later input for pyrsa.dataset())
    
    Args:
        mask_dict_d (dict):
            {roi_h : disjunct, smoothed and resampled mask}
        roi_h (str):
            ROI name
        glm_dir (str):
            Path to GLM directory
        n_stim (int):
            number of stimuli (the first n_stim beta images to loop over)

    Returns:
        measurements (2D-array):
            array of beta patterns (stimulus x channels)
    '''
    if n_stim == None:
        n_stim = len(glob.glob(glm_dir + os.sep + "beta_*"))
    
    mask_tmp = mask_dict_d[roi_h]
    mask_array = mask_tmp.get_fdata()
    mask_size = sum(sum(sum(mask_array)))
    measurements = np.empty((n_stim, int(mask_size))) 
   
    for stim_num in range(n_stim):
        beta_vector = apply_roi_mask_SPM(glm_dir, stim_num+1, mask_tmp, target = 'betas', method = method)
        measurements[stim_num,:] = beta_vector.copy() 
        
    return measurements


def mask_and_get_SPM_residuals(mask_dict_d, roi_h, glm_dir, n_res=None, method = 'nilearn'):
    '''
    Apply mask to all residual images and store them in an array (used for estimating the precision matrix for crossnobis distance)
    
    Args:
        mask_dict_d (dict):
            {roi_h : disjunct, smoothed and resampled mask}
        roi_h (str):
            ROI name
        glm_dir (str):
            Path to GLM directory
        n_res (int):
            Number of residual images per run
    
    Returns:
        run_noise (2D-array):
            array of residual patterns (n_residuals x channels)
    '''
    if n_res == None:
         n_res = len(glob.glob(glm_dir + os.sep + "Res_*"))
         
    mask_tmp = mask_dict_d[roi_h]
    mask_array = mask_tmp.get_fdata()
    mask_size = sum(sum(sum(mask_array)))
    residuals = np.empty((n_res, int(mask_size))) 
   
    for res_num in range(n_res):
        res_vector = apply_roi_mask_SPM(glm_dir, res_num+1, mask_tmp, target = 'residuals', method = method)
        residuals[res_num,:] = res_vector.copy() 
        
    return residuals


def get_voxel_ids(mask_dict, roi):
    '''
    Get Voxel positions of the applied mask (useful for channel descriptors)
    
    Args:
        mask_dict (dict):
            {ROI_name : mask} (Preferably mask_dict_d)
        roi (str):
            key in mask_dict
    
    Returns:
        voxel_ids (array):
            n_voxel x 3 - dimensional array of the X,Y,Z coordinates of voxels that are contained in mask
                
    '''
    mask_tmp = mask_dict[roi]
    mask_array = mask_tmp.get_fdata()
    voxel_ids = np.array(np.where(mask_array.astype(bool)==True)).T
    return voxel_ids