#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot F-Values

@author: alex
"""
import os
from nibabel import nifti1
from nilearn import plotting

sub = 1
ds_dir = "/home/alex/Datasets/ds001246/"
run_dir = os.path.join(ds_dir, "derivatives", "SPM_3", "sub-"+str(sub).zfill(2),"ses-perceptionTest01","run-01")



# Plot 

bg_img = os.path.join(ds_dir, "sub-"+str(sub).zfill(2), "ses-anatomy", "anat", "sub-"+str(sub).zfill(2)+"_ses-anatomy_T1w.nii.gz")
plot_cords = (-7, -48, 0)
plot_cords = (-41, -37, 16) 
plot_title = "spmF_0001.nii"
img_path = os.path.join(run_dir, "spmF_0001.nii")
img = nifti1.load(img_path)
img.header
img_plot = plotting.view_img(img, bg_img = bg_img , cut_coords=plot_cords, title=plot_title)
img_plot.open_in_browser()
img_plot.save_as_html("Sub-0"+str(sub)+"_"+roi_plot+'_T1w.html')