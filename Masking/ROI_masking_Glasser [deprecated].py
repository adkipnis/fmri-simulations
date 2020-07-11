#######################
# Unfinished / Buggy! #
#######################

def get_roi_ids_glasser(txt_dir, target_ROIs):
    '''
    load label txt to find the index for each target ROI

    Inputs:
    - csv_dir (str): path to Glasser txt file
    - target_ROIs (list): strings of target ROI names

    Output:
    - target_ROIs_g (list): strings of original ROI names as found in the txt file (corresponding to "target_ROIs")
    - roi_ids (dict): target_ROIs_g as keys and [RH Index, LH Index] as values

    '''
    import pandas as pd

    labels = pd.read_csv(txt_dir, sep=" ", header=None)
    labels.columns = ["ID", "ROI"]
    targets = str()
    for roi in target_ROIs:
        targets += str(roi)+'_|'

    label_bool = labels['ROI'].str.contains(targets[:-1])
    roi_ids_tab = labels[label_bool]
    roi_ids = roi_ids_tab.set_index('ROI').T.to_dict('list')
    target_ROIs_g = list(roi_ids.keys())
    return roi_ids, target_ROIs_g


def get_shape_dif(mni_dir, glasser_dir):
    '''

    Inputs:
    - mnidir (str): path to MNI nifti file
    - glasser_dir (str): path to Glasser nifti file

    Output:
    - shape_diff (tuple): dimension ratios between MNI space and Atlas space
    '''
    from nibabel import nifti1
    mni = nifti1.load(mni_dir)
    mni_shape = mni.dataobj.shape
    atlas = nifti1.load(glasser_dir)
    atlas_shape = atlas.dataobj.shape
    shape_diff = tuple(np.array(mni_shape)/np.array(atlas_shape))
    return shape_diff

def check_affine_order(affine):
    '''
    Check whether the first minor of the affine matrix has the shape of the identity matrix and
    '''
    affine_tmp = affine.copy()
    # check whether he dimensions are messed up (thx Chris Gorgolewsky)
    main_diag = np.diagonal(affine_tmp, offset=0)
    if 0 in main_diag:
        # find positions of nonzero elements
        col_1_nz = (np.where(affine_tmp[0, 0:3] != 0))[0][0]
        col_2_nz = (np.where(affine_tmp[1, 0:3] != 0))[0][0]
        col_3_nz = (np.where(affine_tmp[2, 0:3] != 0))[0][0]
        order = [col_1_nz, col_2_nz, col_3_nz]
    else:
        order = [0, 1, 2]
    return order

def adjust_affine(affine, shape_diff, order, rectify):
    '''
    Adjust diagonal elements of the affine transform of the atlas nifti file
    '''
    affine_tmp = affine.copy()
    for k in range(3):
        affine_tmp[k, order[k]] = affine_tmp[k, order[k]] * 1 / np.array(shape_diff)[k]

    if rectify and order != [0,1,2]:
        affine_rect = affine_tmp.copy()
        for k in range(3):
            affine_rect[:,k] = affine_tmp[:, order[k]]
            affine_rect[k, 3] = affine_tmp[order[k],3]
        affine_tmp = affine_rect.copy()

    affine_n = affine_tmp.copy()
    return affine_n


def gen_mask(glasser_dir, roi_id, shape_diff, plot_u, plot_n):
    '''
    Generate unilateral mask file for pre-specified anatomical region

    Inputs:
    - glasser_dir (str): Path to Glasser atlas nii.gz file
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
    atlas = nifti1.load(glasser_dir)
    atlas_array_tmp = atlas.get_fdata()
    atlas_array = np.rint(atlas_array_tmp)

    # Load in affine matrix
    affine = atlas.affine.copy()
    order = check_affine_order(affine)
    if order != [0, 1, 2]:
        atlas_array_tmp = np.transpose(atlas_array, tuple(order))
        atlas_array = atlas_array_tmp
        shape_diff_tmp = (shape_diff[order[0]], shape_diff[order[1]], shape_diff[order[2]])
        shape_diff = shape_diff_tmp
    affine_n = adjust_affine(affine, shape_diff, order, True)


    # full_atlas = nifti1.Nifti1Image(atlas_array, affine)
    # plotting.view_img(full_atlas, threshold=0, cut_coords=(6, -80, 10)).open_in_browser()

    # Extract ROI-specific mask and normalize it
    mask = (atlas_array == roi_id).astype(float)
    mask_n_tmp = zoom(mask, shape_diff)  # normalize dimensions
    mask_n = np.rint(mask_n_tmp)  # round to nearest integer (to get binary array again)
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
template_dir = "/home/alex/templateflow/"
# glasser = "tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii.gz" # wrong MNI template
glasser = "tpl-Glasser/MMP_in_MNI_corr.nii.gz"
mni_template = "tpl-MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_res-02_desc-brain_T1w.nii.gz"
labels_txt = "tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt"
glasser_dir = template_dir+glasser
mni_dir = template_dir+mni_template
output_dir = "/home/alex/Documents/10. Semester - CU/Master Thesis/Python Code/ds001246/derivatives/ROI_masks/MNI152NLin2009cAsym/"

# Specify ROIs
target_ROIs = ['V%d' % i for i in range(1,5)] + ['VMV%d' % i for i in range(1,4)] + ['PHA%d' % i for i in range(1,4)] + ['VVC', 'FFC', 'TF', 'PeEc', 'MT', 'MST']
hemi_dict = {'right': 0, 'left':1}
hemi_dict_inv = {0 : 'right', 1 : 'left'}


# Get ROI IDs and difference in the dimensions of MNI template and Glasser Atlas
roi_ids, target_ROIs_g = get_roi_ids_glasser(template_dir+labels_txt, target_ROIs)
shape_diff = get_shape_dif(mni_dir, glasser_dir)
# r = 0; hemi = 0

# Create mask for each ROI for each hemisphere
for r in range(len(target_ROIs)):
    for hemi in range(2):
        target = roi_ids[target_ROIs_g[r]][hemi]
        mask_n_nifti = gen_mask(template_dir+glasser, target, shape_diff, True, True)
        # Note: If you plot both the raw and the normalized mask, the latter will only have intermediate values because of the affine transform.
        # The mask array itself has binary values only (cf. gen_mask() ).
        fname = output_dir + "Mask_MNI152NLin2009cAsym_" + target_ROIs[r] + '_' + hemi_dict_inv[hemi] + '.nii.gz'
        nib.save(mask_n_nifti, fname)

######################################################################################################################
