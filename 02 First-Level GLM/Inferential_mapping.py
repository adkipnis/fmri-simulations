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



# Plot F-Map
bg_img = os.path.join(ds_dir, "sub-"+str(sub).zfill(2), "ses-anatomy", "anat", "sub-"+str(sub).zfill(2)+"_ses-anatomy_T1w.nii.gz")
plot_cords = (-7, -48, 0)
plot_title = "F-Map"
img_path = os.path.join(run_dir, "spmF_0001.nii")
img = nifti1.load(img_path)
img_plot = plotting.view_img(img, bg_img = bg_img , cut_coords=plot_cords, title=plot_title)
img_plot.open_in_browser()
img_plot.save_as_html("Sub-0"+str(sub)+"_"+plot_title+'_T1w.html')


# Plot significant Voxels
plot_title = "significant voxels"
img_path = os.path.join(run_dir, "sig_Effects-of-interest.nii")
img = nifti1.load(img_path)
img_plot = plotting.view_img(img, bg_img = bg_img , cut_coords=plot_cords, title=plot_title)
img_plot.open_in_browser()
img_plot.save_as_html("Sub-0"+str(sub)+"_"+plot_title+'_T1w.html')

# Plot significant clusters
plot_title = "significant clusters"
img_path = os.path.join(run_dir, "sigclus_Effects-of-interest.nii")
img = nifti1.load(img_path)
img_plot = plotting.view_img(img, bg_img = bg_img , cut_coords=plot_cords, title=plot_title)
img_plot.open_in_browser()
img_plot.save_as_html("Sub-0"+str(sub)+"_"+plot_title+'_T1w.html')