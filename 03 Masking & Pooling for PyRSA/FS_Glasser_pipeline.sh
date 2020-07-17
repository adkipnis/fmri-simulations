#!/bin/bash
# Pipeline for mapping the Glasser surface atlas in standard space to several ROI masks in T1w space
# scripts inspired by https://github.com/solleo

# Prerequisite: fmriprep's freesurfer output is in the derivatives directory and the Glasser surface atlas (rh/lh.HCP-MMP1.annot) in the .../fsaverage/label/ directory

./01-FS_register_Glasser_annots.sh
./02-FS_transform_annots_to_labels.sh
./03-FS_map_multi_labels_to_T1w_GM_fill.sh -p
