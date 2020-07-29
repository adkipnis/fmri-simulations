#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot inferential maps created by SPM

@author: alex
"""

def get_statistical_img(sub, ds_dir, session_type, FWHM = None, image_type = "spmF_0001.nii", method = 'run-wise', bg_img = None):    
    """
    Load or create statistical parametric map of F-values, clusters or significance
    
    Args:
        sub (int):
            subject ID 
        ds_dir (str):
            path to BIDS directory
        session_type (str):
            only runs for sessions that contain this string in their name will be loaded
        FWHM (int):
            Used to generate folder name of SPM outouts in the derivatives directory of your dataset
        image_type (str):
            relative filename of the statistical image of interest
        method (str):
            run-wise (averages each run's statistical image) or pooled (takes image from '-pooled' folder)
        bg_img(Nifti):
            Nifti object used as background image for single run plots. If none is provided, no single runs will be plotted.
                      
    Returns:
        stat_img(Nifti):
            Nifti object for the specified type of statistical map
        array_superset (list):
            Set of the img_arrays for each run 
    """
    assert method in ['run-wise', 'pooled'], "Non-existent method specified."
    output_type = "SPM"
    
    if FWHM is not None:
        output_type += "_"+str(FWHM)
    if method == 'pooled':
        output_type += "_"+method
    stat_img = []
    array_superset = []    
    sub_dir = os.path.join(ds_dir, "derivatives", output_type, "sub-"+str(sub).zfill(2))
    
    if method == 'run-wise':     
        sessions = glob.glob(sub_dir + os.sep + "*" + session_type + "*" )
        sessions.sort()
        n_ses = len(sessions)
        array_superset = []
        
        for ses in range(1, n_ses+1):   
            ses_dir = sessions[ses-1]
            n_runs = len(glob.glob(ses_dir + os.sep + "run*"))
            
            for run in range(1, n_runs+1):
                run_dir = os.path.join(ses_dir, "run-"+str(run).zfill(2))
                img_path = os.path.join(run_dir, image_type)
                img = nifti1.load(img_path)
                img_array = img.get_fdata()
                array_superset.append(img_array)
                if bg_img is not None:    
                    plot_title = "Sub-0"+str(sub)+"_"+output_type+"_"+image_type[:-4]+"_"+"ses-"+str(ses).zfill(2)+"_"+"run-"+str(run).zfill(2)
                    plotting.view_img(img, bg_img = bg_img , cut_coords=plot_cords, title=plot_title,
                                 resampling_interpolation=resample_plot, vmax=vmax, opacity=1, symmetric_cmap=False, cmap=plt.get_cmap(cmap_name)).save_as_html(plot_title+'.html')
                    # plotting.plot_stat_map(img, bg_img = bg_img , display_mode='z', cut_coords=5, title=plot_title,
                    #              resampling_interpolation=resample_plot, vmax=vmax, colorbar=False, symmetric_cbar=False, output_file=plot_title+'.png')
                    
        img_affine = img.affine.copy()
        img_array_mean = np.average(array_superset, axis=0)
        stat_img = nifti1.Nifti1Image(img_array_mean, img_affine)
    
    elif method == 'pooled':
        ses_dir = sub_dir + os.sep + session_type + "-results" 
        img_path = os.path.join(ses_dir, image_type)
        stat_img = nifti1.load(img_path)
        
    
    
    return stat_img, array_superset, output_type

###################
import os
import glob
import numpy as np
from nibabel import nifti1
from nilearn import plotting
import matplotlib.pyplot as plt

ds_dir = "/home/alex/Datasets/ds001246/"
session_type = "perceptionTest"
FWHM = 3
# method = "pooled"
method = "run-wise"
vmax = 8 #max statistical value (for standard color representations)
image_type = "spmF_0001.nii"
# image_type = "sigclus_Effects-of-interest_FWE.nii"
# image_type = "sig_Effects-of-interest_FWE.nii"
resample_plot = 'nearest' #alternative: 'continuous' or 'linear')
cmap_name = 'hot' #color map, see https://matplotlib.org/3.3.0/gallery/color/colormap_reference.html
plot_cords = (-7, -48, 0)
n_subs = len(glob.glob(ds_dir + os.sep + "sub*"))

       
for sub in range(1, n_subs+1):  
    # Load statistical image (plot each run if bg_img is specified in function)
    bg_img = os.path.join(ds_dir, "sub-"+str(sub).zfill(2), "ses-anatomy", "anat", "sub-"+str(sub).zfill(2)+"_ses-anatomy_T1w.nii.gz")
    stat_img, array_superset, output_type = get_statistical_img(sub, ds_dir, session_type, FWHM = FWHM, image_type = image_type, method = method, bg_img = None)
    
    # Plotting
    plot_title = "Sub-0"+str(sub)+"_"+output_type+"_"+image_type[:-4]+"_"+method
    img_plot = plotting.view_img(stat_img, bg_img = bg_img , cut_coords=plot_cords, title=plot_title,
                                 resampling_interpolation=resample_plot, vmax=vmax, opacity=1, symmetric_cmap=False, cmap=plt.get_cmap(cmap_name))
    img_plot.open_in_browser()
    img_plot.save_as_html(plot_title+'.html')