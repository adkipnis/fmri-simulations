#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test if ROIs contain significant voxels

@author: alex
"""

import os
import glob
import csv
import mask_utils
import numpy as np
from nibabel import nifti1
from nilearn import plotting

ds_dir              = "/home/alex/Datasets/ds001246/"
session_type        = "perceptionTest"
method              = "pooled"
output_type         = "SPM_6" + "_" + method
# image_type          = "significant_clusters_Effects-of-interest_FWE"
image_type          = "significant_clusters_cosine_basis_functions_FWE"
freesurfer_mri      = "mri_glasser" #Name of the directory in which subject specific volumetric ROI masks are saved by FreeSurfer

# Optional plotting
plot_img            = False
resample_plot       = 'nearest' #alternative: 'continuous' or 'linear')
cmap_name           = 'hot' #color map, see https://matplotlib.org/3.3.0/gallery/color/colormap_reference.html
plot_cords          = (-7, -48, 0)
n_subs              = len(glob.glob(ds_dir + os.sep + "sub*"))
sig_dict_superset   = []

#sub=1       
for sub in range(1, n_subs+1):  
    
    # Load statistical image (plot each run if bg_img is specified in function)
    stat_img = []
    sig_dict = {}    
    
    sub_dir = os.path.join(ds_dir, "derivatives", output_type, "sub-"+str(sub).zfill(2))
    ses_dir = sub_dir + os.sep + session_type + "-results" 
    
    # Load table footer to get FWER conrolling F-threshold
    # footer = []
    # with open(ses_dir + os.sep + "TabDat_Effects-of-interest_footer.txt", newline='') as inputfile:
    #     for row in csv.reader(inputfile):
    #         footer.append(row)
    # stat_thresh = float(footer[4][0].split(' ')[1])
        
    # Load mask dictionaries
    mask_dir = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" + str(sub).zfill(2), freesurfer_mri)
    mask_dict = mask_utils.load_dict(os.path.join(mask_dir, "sub-" + str(sub).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
    
    # Load image
    img_path = os.path.join(ses_dir, image_type+".nii")
    stat_img = nifti1.load(img_path)
    stat_img_array = stat_img.get_fdata()
    
    # Apply mask and statistical threshold
    for roi_h in mask_dict.keys():    
        mask = mask_dict[roi_h]
        mask_array = mask.get_fdata()
        stat_img_array_masked = np.multiply(mask_array, stat_img_array)
        # stat_img_array_thresholded = (stat_img_array_masked>=stat_thresh).astype(int)
        
        # Save number of significant voxels in ROI
        sig_dict.update({roi_h: sum(sum(sum(stat_img_array_masked)))})
        
        # Plotting
        if plot_img:
            stat_img_thresholded = nifti1.Nifti1Image(stat_img_array_masked, stat_img.affine.copy())
            bg_img = os.path.join(ds_dir, "sub-"+str(sub).zfill(2), "ses-anatomy", "anat", "sub-"+str(sub).zfill(2)+"_ses-anatomy_T1w.nii.gz")
            plot_title = "Sub-0"+str(sub)+"_"+output_type+"_"+image_type[:-4]+"_"+method
            img_plot = plotting.view_img(stat_img_thresholded, bg_img = bg_img , cut_coords = plot_cords, title = plot_title,
                                         resampling_interpolation = resample_plot, vmax = 1, opacity = 1, symmetric_cmap = False, cmap = cmap_name)
            img_plot.open_in_browser()
            img_plot.save_as_html(plot_title+'.html')
        
    sig_dict_superset.append(sig_dict)
np.save(os.path.join(ds_dir, "derivatives", "ROI_wise_" + image_type + "_" + output_type + ".npy"), sig_dict_superset)