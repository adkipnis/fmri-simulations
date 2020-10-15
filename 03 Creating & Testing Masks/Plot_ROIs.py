#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot masks

@author: alex
"""

import os
import glob
import mask_utils
from nilearn import plotting

ds_dir              = "/home/alex/Datasets/ds001246/"
freesurfer_mri      = "mri_glasser" #Name of the directory in which subject specific volumetric ROI masks are saved by FreeSurfer

# Optional plotting
plot_img            = False
resample_plot       = 'nearest' #alternative: 'continuous' or 'linear')
cmap_name           = 'hot' #color map, see https://matplotlib.org/3.3.0/gallery/color/colormap_reference.html
plot_cords          = (-7, -48, 0)
n_subs              = len(glob.glob(ds_dir + os.sep + "sub*"))


#sub=1       
for sub in range(1, n_subs+1):  
     
    # Load mask dictionaries
    mask_dir = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" + str(sub).zfill(2), freesurfer_mri)
    mask_dict = mask_utils.load_dict(os.path.join(mask_dir, "sub-" + str(sub).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
      
    # Plot ROI Mask
    for roi_h in mask_dict.keys():    
        mask = mask_dict[roi_h]
        bg_img = os.path.join(ds_dir, "sub-"+str(sub).zfill(2), "ses-anatomy", "anat", "sub-"+str(sub).zfill(2)+"_ses-anatomy_T1w.nii.gz")
        plot_title = "Sub-0"+str(sub)+"_"+roi_h
        img_plot = plotting.view_img(mask, bg_img = bg_img , cut_coords = plot_cords, title = plot_title,
                                     resampling_interpolation = resample_plot, vmax = 1, opacity = 1, symmetric_cmap = False, cmap = cmap_name)
        img_plot.open_in_browser()
        # img_plot.save_as_html(plot_title+'.html')
        
