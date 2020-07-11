def get_roi_ids(csv_dir, target_ROIs):
    '''
    load label csv to find the LH and RH index for each target ROI

    Inputs:
    - csv_dir (str): path to CerebrA csv file (altered for pandas compatibility)
    - target_ROIs (str): list of Mindboggle ROI names

    Output:
    - roi_ids (dict): target_ROIs as keys and [RH Index, LH Index] as values
    '''
    import pandas as pd

    labels = pd.read_csv(csv_dir)
    targets = str()
    for roi in target_ROIs:
        targets += str(roi)+'|'

    label_bool = labels['Label Name'].str.contains(targets[:-1])
    roi_ids_tab = labels[['Label Name', 'RH Label','LH Label']][label_bool]
    roi_ids = roi_ids_tab.set_index('Label Name').T.to_dict('list')
    return roi_ids


def get_shape_dif(mni_dir, cerebra_dir):
    '''

    Inputs:
    - mnidir (str): path to MNI nifti file
    - cerebradir (str): path to CerebrA nifti file

    Output:
    - shape_diff (tuple): dimension ratios between MNI space and Atlas space
    '''
    from nibabel import nifti1
    mni = nifti1.load(mni_dir)
    mni_shape = mni.dataobj.shape
    atlas = nifti1.load(cerebra_dir)
    atlas_shape = atlas.dataobj.shape
    shape_diff = tuple(np.array(mni_shape)/np.array(atlas_shape))
    return shape_diff


def adjust_affine(affine, shape_diff):
    '''
    Adjust diagonal elements of the affine transform of the CerebrA nifti file
    '''
    affine_tmp = affine.copy()
    for k in range(3):
        affine_tmp[k, k] = 1 / np.array(shape_diff)[k]
    affine_n = affine_tmp.copy()
    return affine_n


def gen_mask(cerebra_dir, roi_id, shape_diff, plot_u, plot_n):
    '''
    Generate unilateral mask file for pre-specified anatomical region

    Inputs:
    #- atlas_array (np array): CerebrA Atlas dataobject
    #- affine (np array): CerebrA Atlas affine transformation matrix
    - roi_id (int): ROI index
    - shape_diff (tuple):  dimension ratios between MNI space and Atlas space
    - plot_u (bool): view html-plot of mask in its initial dimensions
    - plot_n (bool): view html-plot of normalized mask

    Output:
    - mask_n_nifti (nifti image): ROI mask in target dimension
    '''
    from nibabel import nifti1
    from scipy.ndimage import zoom
    from nilearn import plotting

    # Load in atlas
    atlas = nifti1.load(cerebra_dir)
    atlas_array = atlas.get_fdata()
    affine = atlas.affine.copy()

    # Extract ROI-specific mask and normalize it
    mask = (atlas_array == roi_id).astype(float)
    mask_n_tmp = zoom(mask, shape_diff)  # normalize dimensions
    mask_n = np.rint(mask_n_tmp)  # round to nearest integer (to get binary array again)
    affine_n = adjust_affine(affine, shape_diff)
    mask_n_nifti = nifti1.Nifti1Image(mask_n, affine_n)

    if plot_u:
        mask_nifti = nifti1.Nifti1Image(mask, affine)
        plotting.view_img(mask_nifti, cut_coords=(6, -80, 10)).open_in_browser()

    if plot_n:
        plotting.view_img(mask_n_nifti, cut_coords=(6, -80, 10)).open_in_browser()

    return mask_n_nifti

#####################################################################################

import numpy as np
import nibabel as nib
# Specify Direcories
template_dir = "/home/alex/templateflow/tpl-MNI152NLin2009cAsym/"
cerebra = "mni_icbm152_CerebrA_tal_nlin_sym_09c.nii"
mni_template = "tpl-MNI152NLin2009cAsym_res-02_desc-brain_T1w.nii.gz"
labels_dir = "CerebrA_LabelDetails_altered.csv"
output_dir = "/home/alex/Documents/10. Semester - CU/Master Thesis/Python Code/ds001246/derivatives/ROI_masks/MNI152NLin2009cAsym/"

# Specify ROIs
target_ROIs = ['Pericalcarine', 'Lateral Occipital', 'Inferior temporal',
               'Fusiform','Parahippocampal','Superior Parietal', 'Rostral Middle Frontal']
hemi_dict = {'right': 0, 'left':1}
hemi_dict_inv = {0 : 'right', 1 : 'left'}
# r = 0; hemi = 0

# Get ROI IDs and difference in the dimensions of MNI template and CerebrA Atlas
roi_ids = get_roi_ids(template_dir+labels_dir, target_ROIs)
shape_diff = get_shape_dif(template_dir+mni_template, template_dir+cerebra)

# Create mask for each ROI for each hemisphere
for r in range(len(target_ROIs)):
    for hemi in range(2):
        target = roi_ids[target_ROIs[r]][hemi]
        mask_n_nifti = gen_mask(template_dir+cerebra, target, shape_diff, True, True)
        # Note: If you plot both the raw and the normalized mask, the latter will only have intermediate values because of the affine transform.
        # The mask array itself has binary values only (cf. gen_mask() ).
        fname = output_dir + "Mask_MNI152NLin2009cAsym_" + target_ROIs[r] + '_' + hemi_dict_inv[hemi] + '.nii.gz'
        nib.save(mask_n_nifti, fname)

######################################################################################################################


## interpolate atlas to target dimensions
# data = atlas.get_fdata()
# data_n = zoom(data, shape_diff) #normalize dimensions
# atlas_n = nifti1.Nifti1Image(data_n, affine_tmp)
# atlas_n = nifti1.Nifti1Image(data_n, atlas.affine)
# plotting.view_img(atlas, cut_coords=(6,-80,10)).open_in_browser()
# plotting.view_img(atlas_n, cut_coords=(6,-80,10)).open_in_browser()

# This would be more sparse but interpolation "mixes" the numeric ROI labels

# Plot some subject data
from nilearn import plotting
from nibabel import nifti1
import numpy as np

mask_dir = "/home/alex/Documents/10. Semester - CU/Master Thesis/Python Code/ds001246/derivatives/ROI_masks/"
mni_mask = "MNI152NLin2009cAsym/Mask_MNI152NLin2009cAsym_Pericalcarine_left.nii.gz"
t1w_mask = "T1w/Sub-01/Mask_T1w_Pericalcarine_left.nii.gz"
t1w_mask = "T1w/Mask_T1w_Glasser.nii.gz"
ref = '/home/alex/Documents/10. Semester - CU/Master Thesis/Python Code/ds001246/derivatives/fmriprep/sub-01/ses-perceptionTraining01/func/sub-01_ses-perceptionTraining01_task-perception_run-01_space-T1w_boldref.nii.gz'
bg_img = '/home/alex/Documents/10. Semester - CU/Master Thesis/Python Code/ds001246/sub-01/ses-anatomy/anat/sub-01_ses-anatomy_T1w.nii.gz'
mni_mask_html = plotting.view_img(mask_dir+mni_mask, bg_img = 'MNI152' , cut_coords=(6, -80, 10))
t1w_mask_html = plotting.view_img(mask_dir+t1w_mask, bg_img = bg_img , cut_coords=(6, -80, 10))

#Save html viewers
mni_mask_html.save_as_html('mni_mask.html')
t1w_mask_html.save_as_html('t1w_mask.html')

# Check out new mask
new_mask = nifti1.load(mask_dir+t1w_mask)
new_mask_array = new_mask.get_fdata()
new_mask_array.shape
np.unique(new_mask_array)
plotting.view_img(new_mask, bg_img = bg_img, cut_coords=(0, 0, 0)).open_in_browser()

# Check out bg_img
bg = nifti1.load(ref)
bg_array = bg.get_fdata()
bg_array.shape
np.unique(bg_array)
plotting.view_img(new_mask, bg_img = bg_img, cut_coords=(0, 0, 0)).open_in_browser()
