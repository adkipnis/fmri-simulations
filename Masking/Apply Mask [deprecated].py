def get_roi_ids_glasser(txt_dir, target_ROIs):
    '''
    load label txt to find the index for each target ROI

    Inputs:
    - csv_dir (str): path to Glasser txt file
    - target_ROIs (list): strings of target ROI names

    Output:
    - roi_ids (dict): target_ROIs_g as keys and [RH Index, LH Index] as values

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

def remove_mask_overlap(mask_dict_s, mask_dict_t):
    '''
    Find overlapping voxels in thresholded masks and only keep the voxel with highest pre-threshold voxel intensity
    
    Inputs:
    - mask_dict_s (dict): Each ROI's smoothed and resampled mask
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
              
        for roi_h in mask_dict_s.keys():
            voxel = mask_dict_d[roi_h].get_fdata()[idx]
            if voxel == 1:
                overlap_list.append(roi_h)
        print("The ROIs", overlap_list, " overlap on voxel: ", idx)
        
        # Check which voxel value is largest on smoothed versions of the overlapping masks 
        for overlap_mask in overlap_list:
            voxel_intensities.append(mask_dict_s[overlap_mask].get_fdata()[idx])
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
            

def apply_roi_mask(glm_dir, stim, mask):
    from nilearn import masking
    betas_path = glm_dir+"beta_"+str(stim).zfill(4)+".nii"
    betas = nifti1.load(betas_path)
    beta_vector = masking.apply_mask(betas, mask, dtype='f', smoothing_fwhm=None, ensure_finite=True)
    return beta_vector

#####################################################################################

from nilearn import plotting

from nibabel import nifti1
import numpy as np
import pandas as pd

# Set directories and specify ROIs
ds_path = "/home/alex/ds001246/"
mask_dir = ds_path + "derivatives/ROI_masks/"
txt_dir = "/home/alex/templateflow/tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt"
target_ROIs = ['V%d' % i for i in range(1,5)] + ['VMV%d' % i for i in range(1,4)] + ['PHA%d' % i for i in range(1,4)] + ['VVC', 'FFC', 'TF', 'PeEc', 'MT', 'MST']
hemi_dict = {'left': 0, 'right':100}
roi_ids = get_roi_ids_glasser(txt_dir, target_ROIs)
n_stim = 50 # number of unique stimuli

# Set mask parameters
fwhm = np.array([3,3,3]) # For mask smoohing (the functional EPI images had a voxel size of 3 × 3 × 3 mm)
threshold = 0.4 # For mask thresholding


#####################################################################################

# Set runner variables and respective paths to atlas and betas
sub = 1
ses = 1
run = 1


t1w_mask = mask_dir + "T1w/sub-"+str(sub).zfill(2)+"_Mask_T1w_Glasser.nii.gz"
glm_dir = ds_path + "derivatives/SPM/sub-"+str(sub).zfill(2)+"/perceptionTest" \
    +str(ses).zfill(2)+"/run-"+str(run).zfill(2)+"/"
beta_ids_path = glm_dir+"sub-"+str(sub).zfill(2)+"-perceptionTest"+str(ses).zfill(2)+"-run-"+str(run)+"_stim_ids.txt"
beta_ids = pd.read_csv(beta_ids_path, sep=" ")
betas_path = glm_dir+"beta_0001.nii"

# Load atlas and betas
atlas = nifti1.load(t1w_mask)
betas_example = nifti1.load(betas_path)

# Create mask dictionary for later use
mask_dict_s = {}
mask_dict_t = {}
for roi in target_ROIs:
    for side in hemi_dict.keys():
        roi_id = roi_ids[roi][0] + hemi_dict[side]
        mask_nifti_r, _, _ = make_roi_mask(atlas, betas_example, roi_id, fwhm) # Mask is not specific to the stimulus of the beta coefficients - only the metadata (like shape) of 'betas' are needed
        mask_nifti_t = threshold_mask(mask_nifti_r, threshold)
        mask_dict_s.update({roi+"_"+side : mask_nifti_r})
        mask_dict_t.update({roi+"_"+side : mask_nifti_t})

# Make all masks disjunct
mask_dict_d = remove_mask_overlap(mask_dict_s, mask_dict_t)


#####################################################################################
## Apply mask to all stimuli in this run
# for roi_h in mask_dict_d.keys():
    mask_tmp = mask_dict_d[roi_h]
    mask_size = sum(sum(sum(mask_tmp.get_fdata())))
    measurements = np.empty((n_stim, int(mask_size))) # observation (here: stimulus x channels)
   
    for stim_num in range(50):
        beta_vector = apply_roi_mask(glm_dir, stim_num+1, mask_tmp)
        measurements[stim_num,:] = beta_vector.copy() 
      

# mask_r = mask_nifti_r.get_fdata()
# betas_masked = np.multiply(mask_r, betas_array)
# betas_masked_nifti = nifti1.Nifti1Image(betas_masked, mask_nifti_r.affine)
# betas_masked_vector = betas_array[mask_r.astype(bool)]
betas_masked_vector = masking.apply_mask(betas, mask_nifti_r, dtype='f', smoothing_fwhm=None, ensure_finite=True)

bg_img = ds_path + "sub-01/ses-anatomy/anat/sub-01_ses-anatomy_T1w.nii.gz"
plot_cords = (-7, -48, 0)
v1_mask = plotting.view_img(mask_nifti, bg_img = bg_img , cut_coords=plot_cords, title="Left V1 Mask (HCP-MMP1.0)")
v1_mask_smoothed = plotting.view_img(mask_nifti_s, bg_img = bg_img , cut_coords=plot_cords, title="Left V1 Mask smoothed")
v1_mask_resampled = plotting.view_img(mask_nifti_r, bg_img = bg_img , cut_coords=plot_cords, title="Left V1 Mask resampled to T2W* resolution")
v1_mask_resampled_t = plotting.view_img(mask_nifti_t, bg_img = bg_img , cut_coords=plot_cords, title="Left V1 Mask resampled and thresholded")

# v2_betas = plotting.view_img(betas_masked_nifti, bg_img = bg_img , cut_coords=plot_cords, title="Beta coefficients after V2 masking")
v1_mask.save_as_html('v1_mask.html')
v1_mask_smoothed.save_as_html('v1_mask_smoothed.html')
v1_mask_resampled.save_as_html('v1_mask_resampled.html')
v1_mask_resampled_t.save_as_html('v1_mask_resampled_threshold='+ str(threshold)+'.html')
# v2_betas.save_as_html('v2_betas.html')
