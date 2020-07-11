def get_roi_ids_glasser(txt_dir, target_ROIs):
    '''
    load label txt to find the index for each target ROI

    Inputs:
    - csv_dir (str): path to Glasser txt file
    - target_ROIs (list): strings of target ROI names

    Output:
    - roi_ids (dict): {target_ROIs_g: Left Hemisphere ROI index}. RH ROI index := left counterpart + 100

    '''
    import pandas as pd

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

    # #add right hemisphere
    # roi_ids_tab['Hemisphere'] = 'L'
    # roi_ids_tab_tmp = roi_ids_tab.copy()
    # roi_ids_tab_tmp['ID'] = roi_ids_tab['ID']+100
    # roi_ids_tab_tmp['Hemisphere'] = 'R'

    roi_ids = roi_ids_tab.set_index('ROI').T.to_dict('list')
    return roi_ids

def merge_masks(atlas_o, roi_ids, merge_list, merged_names):
    '''
    Merge prespecified tuples of ROIs in original atlas

    Inputs:
    - atlas_o (Nifti1Image): Original T1w Atlas with voxel values indicative of ROI identity
    - roi_ids (dict): {target_ROIs_g: Left Hemisphere ROI index}
    - merge_list (list): List of tuples of ROI names. All ROIs within a Tuple will be merged.
        
    Output:
    - atlas (Nifti1Image):  T1w Atlas with voxel values indicative of ROI identity
    - roi_ids (dict): {target_ROIs_g: Left Hemisphere ROI index}, with remapped ROI indices

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
            atlas_array[atlas_array==rest_id+200] = ref_id+200 # repeat for right hemisphere
            roi_ids.pop(other_roi)
        roi_ids[merged_names[s]] = roi_ids.pop(sublist[0])
        s+=1
    atlas = nifti1.Nifti1Image(atlas_array, atlas_affine)
    return atlas, roi_ids
            
def make_roi_mask(atlas, betas_example, roi_id, fwhm):
    '''
    Extract ROI-specific mask from Atlas in T1w space and sample it down to have the dimensions of the to-be-masked image of betas

    Inputs:
    - atlas (Nifti1Image): Atlas with voxel values indicative of ROI identity
    - betas_example (Nifti1Image): Example image of GLM coefficients which will be masked later (only needed for shape and affine)
    - roi_id (int): ROI-specific integer used in "atlas"
    - fwhm (float or np.ndarray): FWHM in mm (along each axis, axis specif or none)
    
    Output:
    - mask_nifti_r (Nifti1Image): resampled ROI-specific mask image (target dimensions)
    - mask_nifti_s (Nifti1Image): smoothed ROI-specific mask image (original dimensions)
    - mask_nifti (Nifti1Image): original binary ROI-specific mask image

    '''
    from nilearn import image

    # Extract atlas and target nii data
    atlas_array = atlas.get_fdata()
    atlas_affine = atlas.affine.copy()
    betas_array = betas_example.get_fdata()
    betas_affine = betas_example.affine.copy()

    # Extract ROI-specific mask and resample it
    mask = (atlas_array == roi_id).astype(float)    
    mask_nifti = nifti1.Nifti1Image(mask, atlas_affine)
    mask_nifti_s = image.smooth_img(mask_nifti, fwhm)
    mask_nifti_r = image.resample_img(mask_nifti_s, target_affine=betas_affine, target_shape=betas_array.shape, interpolation='nearest', copy=True)
    
    return mask_nifti_r, mask_nifti_s, mask_nifti

def threshold_mask(mask, threshold):
    '''
    After smoothing, the downsampled masks are not binary any more. Make a mask binary again.
    
    Inputs:
    - mask (Nifti1Image): resampled ROI-specific mask image 
    - threshold (float): between 0 and 1
    
    Output:
    - mask_t (Nifti1Image): binary resampled ROI-specific mask image after thresholding

    '''
    mask_array = mask.get_fdata()
    mask_thresholded = (mask_array >= threshold).astype(float)
    mask_t = nifti1.Nifti1Image(mask_thresholded, mask.affine.copy())
    return mask_t

def create_mask_dict(atlas, betas_example, roi_ids, fwhm, threshold):
    '''
    Create a dictionary of mask nifti files for every target ROI
    
    Inputs:
    - atlas (Nifti1Image): Atlas with voxel values indicative of ROI identity
    - betas_example (Nifti1Image): Example image of GLM coefficients which will be masked later (only needed for shape and affine)
    - roi_ids (dict): {target_ROIs_g: Left Hemisphere ROI index}
    - fwhm (float or np.ndarray): FWHM in mm (along each axis, axis specif or none)
    - threshold (float): between 0 and 1 for thresholding after smoothing and resampling

    Output:
    - mask_dict_o (dict): {target_ROI+hemisphere : original mask}
    - mask_dict_o (dict): {target_ROI+hemisphere : original mask after smoothing}
    - mask_dict_r (dict): {target_ROI+hemisphere : smoothed and resampled mask}
    - mask_dict_t (dict): {target_ROI+hemisphere : smoothed, resampled and thresholded mask}
    
    '''
    hemi_dict = {'left': 0, 'right':200}
    mask_dict_o = {}
    mask_dict_s = {}
    mask_dict_r = {}
    mask_dict_t = {}
    
    
    for roi in roi_ids.keys(): 
        for side in hemi_dict.keys():
            roi_id = roi_ids[roi][0] + hemi_dict[side]
            mask_nifti_r, mask_nifti_s , mask_nifti_o = make_roi_mask(atlas, betas_example, roi_id, fwhm) # Mask is not specific to the stimulus of the beta coefficients - only the metadata (like shape) of 'betas' are needed
            mask_nifti_t = threshold_mask(mask_nifti_r, threshold)
            mask_dict_o.update({roi+"_"+side : mask_nifti_o})
            mask_dict_s.update({roi+"_"+side : mask_nifti_s})
            mask_dict_r.update({roi+"_"+side : mask_nifti_r})
            mask_dict_t.update({roi+"_"+side : mask_nifti_t})

    return mask_dict_o, mask_dict_s, mask_dict_r, mask_dict_t

def remove_mask_overlap(mask_dict_r, mask_dict_t):
    '''
    Find overlapping voxels in thresholded masks and only keep the voxel with highest pre-threshold voxel intensity
    
    Inputs:
    - mask_dict_r (dict): {ROI_name : smoothed and resampled mask}
    - mask_dict_t (dict): The above after thresholding

    Output:
    - mask_dict_d (dict): mask_dict_t after removal of overlapping voxels
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
            
def apply_roi_mask(glm_dir, img_num, mask, target = 'betas'):
    '''
    Load GLM pattern for prespecified stimulus and apply mask to get a vectorized version
    
    Inputs:
    - glm_dir (str): Path to GLM directory (for SPM outputs)
    - img_num (int): Stimulus number
    - mask (Nifti1Image): To be applied mask (must have same dimensions as the GLM betas image)

    Output:
    - beta_vector (1-D array): Vector of Beta-coefficients for each voxel in mask
    '''
    
    from nilearn import masking
    
    if target == 'betas':
        betas_path = os.path.join(glm_dir,"beta_" + str(img_num).zfill(4) + ".nii")
        image = nifti1.load(betas_path)
    elif target == 'residuals':
        residuals_path = os.path.join(glm_dir,"Res_" + str(img_num).zfill(4) + ".nii")
        image = nifti1.load(residuals_path)
        
    # Alternative to nilearn's implementation (allows for plotting of masked 3d image)
    # mask_array = mask.get_fdata()
    # betas_array = betas.get_fdata()
    # betas_masked = np.multiply(mask_array, betas_array)
    # betas_masked_nifti = nifti1.Nifti1Image(betas_masked, mask.affine.copy())
    # betas_masked_vector = betas_array[mask_array.astype(bool)]
    
    beta_vector = masking.apply_mask(image, mask, dtype='f', smoothing_fwhm=None, ensure_finite=True)
   
    return beta_vector

def get_voxel_ids(mask_dict_d, roi_h):
    '''
    Get Voxel positions of the applied mask (useful for channel descriptors)
    
    Inputs:
    - mask_dict_d (dict): {ROI_name : smoothed, resampled and disjunct mask}
    - roi_h (str): key in mask_dict_d
    '''
    mask_tmp = mask_dict_d[roi_h]
    mask_array = mask_tmp.get_fdata()
    voxel_ids = np.array(np.where(mask_array.astype(bool)==True)).T
    return voxel_ids

def generate_measurements(mask_dict_d, roi_h, glm_dir, n_stim):
    '''
    Apply mask to all beta images of interest and store them in a measurements array (later input for pyrsa.dataset())
    
    Inputs:
    - mask_dict_d (dict): {roi_h : disjunct, smoothed and resampled mask}
    - roi_h (str): ROI name
    - glm_dir (str): Path to GLM directory
    - n_stim (int): number of stimuli (the first n_stim beta images to loop over)

    Output:
    - measurements (2D-array): array of beta patterns (stimulus x channels)
    '''
    mask_tmp = mask_dict_d[roi_h]
    mask_array = mask_tmp.get_fdata()
    mask_size = sum(sum(sum(mask_array)))
    measurements = np.empty((n_stim, int(mask_size))) 
   
    for stim_num in range(n_stim):
        beta_vector = apply_roi_mask(glm_dir, stim_num+1, mask_tmp, target = 'betas')
        measurements[stim_num,:] = beta_vector.copy() 
        
    return measurements

def generate_residuals(mask_dict_d, roi_h, glm_dir, n_res=None):
    '''
    Apply mask to all residual images and store them in an array (used for estimating the precision matrix for crossnobis distance)
    
    Inputs:
    - mask_dict_d (dict): {roi_h : disjunct, smoothed and resampled mask}
    - roi_h (str): ROI name
    - glm_dir (str): Path to GLM directory
    - n_res (int): Number of residual images per run
    

    Output:
    - run_noise (2D-array): array of residual patterns (n_residuals x channels)
    '''
    if n_res == None:
         n_res = len(glob.glob(glm_dir + os.sep + "Res_*"))
         
    mask_tmp = mask_dict_d[roi_h]
    mask_array = mask_tmp.get_fdata()
    mask_size = sum(sum(sum(mask_array)))
    residuals = np.empty((n_res, int(mask_size))) 
   
    for res_num in range(n_res):
        res_vector = apply_roi_mask(glm_dir, res_num+1, mask_tmp, target = 'residuals')
        residuals[res_num,:] = res_vector.copy() 
        
    return residuals


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

#####################################################################################

import glob
import os
from nibabel import nifti1
import numpy as np
import pandas as pd
import pyrsa

# Set directories and specify ROIs
ds_dir = "/home/alex/ds001246/"
txt_dir = "/home/alex/templateflow/tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt" #directory of mask descriptors
beta_type = "SPM"
ses_type = "perceptionTest"
target_ROIs = ['V%d' % i for i in range(1,5)] + ['VMV%d' % i for i in range(1,4)] + ['PHA%d' % i for i in range(1,4)] + ['VVC', 'FFC', 'TF', 'PeEc', 'MT', 'MST']
n_stim = 50 # number of unique stimuli
label_dict = np.load(os.path.join(ds_dir, "stimulus_label_dictionary.npy"),allow_pickle='TRUE').item()
spm_dir = os.path.join(ds_dir, "derivatives", beta_type)
mask_dir = os.path.join(ds_dir, "derivatives", "ROI_masks")
n_subs = len(glob.glob(ds_dir + os.sep + "sub*"))
processing_mode = 'residuals' # alternatively: 'datasets' or 'both'

# Set mask parameters
mask_merging = True
fwhm = np.array([3,3,3]) # For mask smoohing (the functional EPI images had a voxel size of 3 × 3 × 3 mm)
threshold = 0.4 # For mask thresholding
if mask_merging:
    merge_list = [tuple('PHA%d' % i for i in range(1,4)), tuple('VMV%d' % i for i in range(1,4)), tuple(['MT', 'MST'])]
    merged_names = ['PHA', 'VMV', 'MT']


#####################################################################################

# sub = 1
for sub in range(2, n_subs+1): 
    
    ### 1. Make Subject-specific ROI dictionary
    # Set respective paths to atlas and betas
    ses, run = 1, 1
    t1w_mask = os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_Mask_T1w_Glasser.nii.gz")
    run_dir = os.path.join(spm_dir, "sub-"+str(sub).zfill(2), ses_type + str(ses).zfill(2))
    glm_dir = os.path.join(run_dir, "run-"+str(run).zfill(2))  
    betas_path_example = os.path.join(glm_dir, "beta_0001.nii")
    # mask_dict_o_path = os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_mask_dict_T1w.npy")
    mask_dict_d_path = os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_mask_dict_T2*w_disjunct.npy")
    roi_ids = get_roi_ids_glasser(txt_dir, target_ROIs)
    ds_output_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "datasets", "sub-"+str(sub).zfill(2))
    if not os.path.isdir(ds_output_dir):
        os.makedirs(ds_output_dir)
    res_output_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "noise", "sub-"+str(sub).zfill(2))
    if not os.path.isdir(res_output_dir):
        os.makedirs(res_output_dir)
    
    # Infer number of sessions and run-wise residuals
    n_ses = len(glob.glob(os.path.join(spm_dir, "sub-"+str(sub).zfill(2), ses_type+"*")))
    n_res = len(glob.glob(glm_dir + os.sep + "Res_*"))
    
    # Load original atlas and betas
    atlas_o = nifti1.load(t1w_mask)
    betas_example = nifti1.load(betas_path_example)
    
    # Merge tuples of masks in T1w space
    if mask_merging:
        atlas, roi_ids = merge_masks(atlas_o, roi_ids, merge_list, merged_names)     
    else:
        atlas = atlas_o
    
    # Load or create mask dictionary, then make all masks disjunct
    if os.path.isfile(mask_dict_d_path):
        # mask_dict_o = np.load(mask_dict_o_path,allow_pickle='TRUE').item()
        mask_dict_d = np.load(mask_dict_d_path,allow_pickle='TRUE').item()
    else:
        mask_dict_o, mask_dict_s, mask_dict_r, mask_dict_t = create_mask_dict(atlas, betas_example, roi_ids, fwhm, threshold)
        mask_dict_d = remove_mask_overlap(mask_dict_r, mask_dict_t)
        # np.save(mask_dict_o_path, mask_dict_o) 
        np.save(mask_dict_d_path, mask_dict_d) 
        
           
    ### 2. Save masked beta patterns for each session and run as pyrsa Dataset
    voxel_ids_dict = {}
    #roi_h = 'V1_left'
    for roi_h in mask_dict_d.keys():
        # We will collect measurements and respective descriptors for each run in each session in a superset
        measurements_superset = []
        session_desc_superset = []
        run_desc_superset = []
        stim_desc_superset = []
        residuals_superset = []      
        # ses = 1
        for ses in range(1, n_ses+1):
            run_dir = os.path.join(spm_dir, "sub-"+str(sub).zfill(2), ses_type + str(ses).zfill(2))           
            n_runs = len(glob.glob(run_dir + os.sep + 'run*'))

            for run in range(1, n_runs+1):
                if processing_mode in ['both', 'datasets']:
                    # 1. Apply mask to all beta coefficient images in this run  
                    glm_dir = os.path.join(run_dir, "run-"+str(run).zfill(2))
                    measurements = generate_measurements(mask_dict_d, roi_h, glm_dir, n_stim)
                    measurements_superset.append(measurements)
                    
                    # 2. Save stimulus ID's to stim_desc_superset
                    beta_ids_dir = os.path.join(glm_dir,"sub-"+str(sub).zfill(2) + "-" + ses_type \
                                    + str(ses).zfill(2) + "-run-" + str(run) + "_stim_ids.txt")
                    beta_ids = pd.read_csv(beta_ids_dir, sep=" ")
                    labels = beta_id_to_label(beta_ids, n_stim, label_dict, crop=True)
                    stim_desc_superset.append(labels)
                    
                
                if processing_mode in ['both', 'residuals']:
                    # 3. Apply mask to all residuals images in this run   
                    residuals = generate_residuals(mask_dict_d, roi_h, glm_dir, n_res)
                    residuals_superset.append(residuals)
                
                 
                # 4. If not already done in a previous run: save mask-specific voxel IDs
                if len(mask_dict_d.keys()) > len(voxel_ids_dict.keys()):
                    voxel_ids_dict.update({roi_h: get_voxel_ids(mask_dict_d, roi_h)})
                run_desc = np.repeat(run, n_stim)
                run_desc_superset.append(run_desc)
                
                
            session_desc = np.repeat(ses, n_stim*n_runs)
            session_desc_superset.append(session_desc)
        
        # Create pyrsa dataset    
        if processing_mode in ['both', 'datasets']:
            dataset = pyrsa.data.Dataset(np.vstack((measurements_superset[:])),
                                 descriptors = {'ROI':roi_h}, 
                                 obs_descriptors = {'stim': np.hstack((stim_desc_superset[:])),
                                                    #'session': np.hstack((session_desc_superset[:])),
                                                    'run': np.hstack((run_desc_superset[:]))+ 100*np.hstack((session_desc_superset[:]))},                             
                                 channel_descriptors = {'positions':voxel_ids_dict[roi_h]}, 
                                 )  
            
            # Save dataset and residuals array
            dataset_filename = os.path.join(ds_output_dir,"RSA_dataset_"+roi_h+"_"+beta_type)
            dataset.save(dataset_filename, file_type='hdf5', overwrite=True)
            print("Created pyrsa dataset:", dataset_filename)
        
        if processing_mode in ['both', 'residuals']:
            residuals_filename = os.path.join(res_output_dir,"Residuals_"+roi_h+"_"+beta_type)
            np.save(residuals_filename, np.vstack((residuals_superset[:])))
            print("Created residuals dataset:", residuals_filename)
        
    np.save(os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_mask-wise_voxel_ids.npy"), voxel_ids_dict)

#####################################################################################
# Plotting different masks 
from nilearn import plotting
bg_img = os.path.join(ds_dir, "sub-"+str(sub).zfill(2), "ses-anatomy", "anat", "sub-"+str(sub).zfill(2)+"_ses-anatomy_T1w.nii.gz")
plot_cords = (-7, -48, 0)
plot_cords = (-41, -37, 16) 
roi_plot = "MT_left"


mask_plot = plotting.view_img(mask_dict_o[roi_plot], bg_img = bg_img , cut_coords=plot_cords, title=roi_plot +" Mask (HCP-MMP1.0)")
mask_plot.open_in_browser()
mask_plot.save_as_html("Sub-0"+str(sub)+"_"+roi_plot+'_T1w.html')

atlas_plot = plotting.view_img(atlas, bg_img = bg_img , cut_coords=plot_cords, title="HCP-MMP1.0 full")
atlas_plot.open_in_browser()
atlas_plot.save_as_html("Sub-0"+str(sub)+"_HCP_T1w.html")

# v1_mask = plotting.view_img(mask_dict_d[roi_h], bg_img = bg_img , cut_coords=plot_cords, title="Left V1 Mask (HCP-MMP1.0)")
# v1_mask_smoothed = plotting.view_img(mask_nifti_s, bg_img = bg_img , cut_coords=plot_cords, title="Left V1 Mask smoothed")
# v1_mask_resampled = plotting.view_img(mask_nifti_r, bg_img = bg_img , cut_coords=plot_cords, title="Left V1 Mask resampled to T2W* resolution")
# v1_mask_resampled_t = plotting.view_img(mask_nifti_t, bg_img = bg_img , cut_coords=plot_cords, title="Left V1 Mask resampled and thresholded")

# # v2_betas = plotting.view_img(betas_masked_nifti, bg_img = bg_img , cut_coords=plot_cords, title="Beta coefficients after V2 masking")
# v1_mask.save_as_html('v1_mask.html')
# v1_mask_smoothed.save_as_html('v1_mask_smoothed.html')
# v1_mask_resampled.save_as_html('v1_mask_resampled.html')
# v1_mask_resampled_t.save_as_html('v1_mask_resampled_threshold='+ str(threshold)+'.html')
# # v2_betas.save_as_html('v2_betas.html')
