#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 11 15:11:54 2020

@author: alex
"""


import mask_utils
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
