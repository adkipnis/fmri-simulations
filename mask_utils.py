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








#-------------------------- Mask alterations

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


def merge_freesurfer_masks(mask_dict_raw, merge_list, merged_names): 
    '''
    Merge prespecified tuples of ROIs in mask_dict_raw for both hemispheres.

    Args:
        mask_dict_raw (dict):
            {hemi_roi : nifti object of the corresponding binary mask}
        merge_list (list):
            List of tuples of ROI names. All ROIs within a Tuple will be merged.
        merged_names (list):
            Corresponding new names for merged masks.
        
    Returns:
        atlas (Nifti1Image):
            T1w Atlas with voxel values indicative of ROI identity
        roi_ids_m (dict):
            {target_ROI: Left Hemisphere ROI index}, where specified target_ROI's are merged

    '''
    assert len(merge_list) == len(merged_names), "Every tuple in 'merge_list' requires a unique new name, provided in merged_names."
    
    mask_dict = mask_dict_raw.copy()
    s = 0 # counter for merged_names
    for sublist in merge_list:
        for hemi in ['L', 'R']:
            mask_raw = mask_dict[sublist[0] + '_' + hemi]
            mask_array = mask_raw.get_fdata()
            mask_affine = mask_raw.affine.copy()
            
            for m in range(1, len(sublist)):
                mask_array += mask_dict[sublist[m] + '_' + hemi].get_fdata() # add remaining masks
            mask_array[mask_array > 1] = 1 # set overlapping voxels to 1
            mask_merged = nifti1.Nifti1Image(mask_array, mask_affine)
            mask_dict.update({merged_names[s] + '_' + hemi : mask_merged})
        s += 1
    
    # Remove original parts of now merged masks
    m = 0
    for sublist in merge_list:     
        for roi in sublist:
            for hemi in ['L', 'R']:
                mask_dict.pop(roi + '_' + hemi, None)
                m += 1
                            
    assert len(mask_dict) == len(mask_dict_raw) - m + 2*len(merged_names), "Unexpected number of masks in final mask_dict. Check, if merged_names are different from strings in merge_list."
                          
    return mask_dict


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


def threshold_mask(mask_EPI, threshold = 0.5, voxel_dimensions = None, mask_T1w = None):
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
    mask_array = mask_EPI.get_fdata()
    
    assert type(threshold) == float or threshold == 'adaptive', "threshold must be floating point number or 'adaptive'."
    if type(threshold) == float:
        assert (threshold >= 0 and threshold <= 1), "numeric threshold must be between 0 and 1."
        assert np.unique(mask_array)[1] < threshold, "numeric threshold set too low, all values will pass."
    else:
        assert type(voxel_dimensions) == dict, "for adaptive thresholding you must provide the voxel dimensions in mm."

    
    if type(threshold) == float:
        mask_thresholded = (mask_array >= threshold).astype(float)
    else:
        mask_array_T1w = mask_T1w.get_fdata()
        geom_vol_T1w = sum(sum(sum(mask_array_T1w))) * np.prod(voxel_dimensions['T1w'])
        mask_thresholded = optimize_threshold(mask_array, geom_vol_T1w, voxel_dimensions)
        
    mask_t = nifti1.Nifti1Image(mask_thresholded, mask_EPI.affine.copy())
    return mask_t


def optimize_threshold(mask_array, geom_vol_T1w, voxel_dimensions, threshold = 0.5, stepsize = 0.01, epsilon = 10, max_iter = 100):
    mask_thresholded = (mask_array >= threshold).astype(float)
    geom_vol_EPI = sum(sum(sum(mask_thresholded))) * np.prod(voxel_dimensions['EPI'])
    delta = geom_vol_EPI - geom_vol_T1w  
    threshold_tmp = threshold
    i = 0
    change = []
    while abs(delta) > epsilon and i <= max_iter:
        if delta > 0:
            threshold_tmp += stepsize
        else:
            threshold_tmp -= stepsize        
        mask_thresholded = (mask_array >= threshold_tmp).astype(float)
        geom_vol_EPI = sum(sum(sum(mask_thresholded))) * np.prod(voxel_dimensions['EPI'])
        delta_tmp = geom_vol_EPI - geom_vol_T1w
        change.append(delta - delta_tmp)
        if i > 2 and change[-1] > change[-2]:
            stepsize *= 0.9
        delta = delta_tmp
        if i > 20 and len(np.unique(np.absolute(change[-10:-1]))) == 1:
           break 
        i += 1
    print("Converged to threshold:", threshold_tmp)
    return mask_thresholded


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







#-------- Mask dictionaries

def create_mask_dicts(atlas, betas_example, roi_ids, fwhm = None, interpolation='nearest', threshold = 0.5, voxel_dimensions = None, rh_constant = 200):
    '''
    Create a series of dictionaries containing mask nifti files for every target ROI
    
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
            mask_nifti_t = threshold_mask(mask_nifti_r, threshold = threshold, voxel_dimensions = voxel_dimensions, mask_T1w = mask_nifti_o)
            mask_dict_o.update({roi+"_"+side : mask_nifti_o})
            mask_dict_s.update({roi+"_"+side : mask_nifti_s})
            mask_dict_r.update({roi+"_"+side : mask_nifti_r})
            mask_dict_t.update({roi+"_"+side : mask_nifti_t})

    return mask_dict_o, mask_dict_s, mask_dict_r, mask_dict_t


def save_dict(dict_object, dict_dir, dict_name):
    '''
    Self-explanatory...
    '''
    dict_path = os.path.join(dict_dir, dict_name + '.npy')
    np.save(dict_path, dict_object) 
    return


def load_dict(dict_path, merge_list = None, merged_names  = None):
    '''
    Load a previously saved mask dictionary. Optionally test if it still matches your target ROIs (in case oyu have changed your merge_list)
    
    Args:
        dict_path (str):
            Absolute path to where you have saved your target (mask) dictionary
        merge_list (list):
            List of tuples of ROI names. All ROIs within a Tuple will be merged.
        merged_names (list):
            List of ROI names for each merged tuple
            
    Returns:
        loaded_dict (dict):
            Dictionary object
    '''
    assert os.path.isfile(dict_path), "No file found matching the provided path."    
    
    loaded_dict = np.load(dict_path, allow_pickle = 'TRUE').item()
    
    if merge_list is not None and merged_names is not None:
        intersect_bad, intersect_good = test_roi_id_overlap(loaded_dict, merge_list, merged_names)
        if len(intersect_bad)>0 or len(intersect_good)<len(merged_names):
            print("Warning: Loaded mask dictionary contains masks that still require merging. You create a new one.")
    return loaded_dict




   




#-------- Mask application (masking)
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


def mask_and_get_SPM_measurements(mask, glm_dir, n_stim=None, method = 'nilearn'):
    '''
    Apply mask to all beta images of interest and store them in a measurements array (later input for pyrsa.dataset())
    
    Args:
        mask (Nifti1Image):
            To be applied mask (must have same dimensions as the GLM betas image)
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
    
    mask_array = mask.get_fdata()
    mask_size = sum(sum(sum(mask_array)))
    measurements = np.empty((n_stim, int(mask_size))) 
   
    for stim_num in range(n_stim):
        beta_vector = apply_roi_mask_SPM(glm_dir, stim_num+1, mask, target = 'betas', method = method)
        measurements[stim_num,:] = beta_vector.copy() 
        
    return measurements


def mask_and_get_SPM_residuals(mask, glm_dir, n_res=None, method = 'nilearn'):
    '''
    Apply mask to all residual images and store them in an array (used for estimating the precision matrix for crossnobis distance)
    
    Args:
        mask (Nifti1Image):
            To be applied mask (must have same dimensions as the GLM betas image)
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
         
    mask_array = mask.get_fdata()
    mask_size = sum(sum(sum(mask_array)))
    residuals = np.empty((n_res, int(mask_size))) 
   
    for res_num in range(n_res):
        res_vector = apply_roi_mask_SPM(glm_dir, res_num+1, mask, target = 'residuals', method = method)
        residuals[res_num,:] = res_vector.copy() 
        
    return residuals







#-------------------------- Miscellaneous

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


def get_voxel_ids(mask):
    '''
    Get Voxel positions of the mask (useful for channel descriptors)
    
    Args:
        mask (Nifti1Image):
            Mask object
    
    Returns:
        voxel_ids (array):
            n_voxel x 3 - dimensional array of the X,Y,Z coordinates of voxels that are contained in mask
                
    '''
    mask_array = mask.get_fdata()
    voxel_ids = np.array(np.where(mask_array.astype(bool)==True)).T
    return voxel_ids


def get_voxel_ids_from_dict(mask_dict):
    '''
    Apply get_voxel_ids() to each mask in a dictionary and create a dictionary for respective voxel IDs.
    '''
    voxel_ids_dict = {}
    for roi in mask_dict.keys():
        voxel_ids = get_voxel_ids(mask_dict[roi])
        voxel_ids_dict.update({roi: voxel_ids})

    return voxel_ids_dict

        
def atlas_from_freesurfer_masks(sub, mask_dir, roi_ids, rh_constant = 200, overwrite = False):
    '''
    Create atlas from individual binary masks produced by freesurfer's 'mri_label2vol' tool

    Args:
        sub (int):
            subject number
        mask_dir (str):
            path to subject's folder with freesurfer output
        roi_ids (dict):
            {target_ROI: Left Hemisphere ROI index}
        rh_constant (int):
            Right Hemisphere ROI index := Left Hemisphere ROI index + rh_constant
        overwrite (bool):
            If true: delete old atlas and create new one
            

    Returns:
        atlas (Nifti1Image):
            T1w Atlas with voxel values indicative of ROI identity (as defined in 'roi_ids')
    '''
    atlas_path = os.path.join(mask_dir, "sub-" + str(sub).zfill(2) + "_Mask_T1w_Glasser.nii.gz")
    
    if os.path.isfile(atlas_path) and not overwrite:
        atlas = nifti1.load(atlas_path)
    else:
        # create atlas
        hemi_dict = {'L':'lh.L_', 'R':'rh.R_'}
        hemi_add = {'L':0, 'R':rh_constant}
        mask_array = []
        mask_array_tmp = []
        for roi in roi_ids.keys():    
            for hemi in ['L', 'R']:
                mask_raw = nifti1.load(os.path.join(mask_dir, hemi_dict[hemi] + roi + '_ROI.fbr.nii.gz'))
                if len(mask_array) == 0:
                    mask_array_tmp = mask_raw.get_fdata()
                    mask_array_tmp[mask_array_tmp==1] = roi_ids[roi][0] + hemi_add[hemi] # convert binary mask by assigning roi_id to masked voxels
                    mask_array = mask_array_tmp.copy()
                    freesurfer_affine = mask_raw.affine.copy()
                else:
                    mask_array_tmp = mask_raw.get_fdata()
                    mask_array[mask_array_tmp==1] = roi_ids[roi][0] + hemi_add[hemi]
        atlas = nifti1.Nifti1Image(mask_array, freesurfer_affine)
        nifti1.save(atlas, atlas_path)
    return atlas






