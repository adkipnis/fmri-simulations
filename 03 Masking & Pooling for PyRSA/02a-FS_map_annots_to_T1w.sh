#!/bin/bash
# Map Glasser surface atlas in native space to volumetric space

# Prerequisite: Successfully run FS_register_Glasser_annots.sh

# Options:  -f - fill ROI voxels between white matter and pial surface
#           -p - plot atlas on sub-01's T1 template

cd $SUBJECTS_DIR 

if  [[ $1 = "-f" ]]; then
    for (( sub=1; sub <= 5; sub=sub+1 )); do
        mkdir -p sub-0$sub/mri_glasser/
        for hemi in lh rh; do
            mri_label2vol \
                --annot sub-0$sub/label_glasser/$hemi.HCP-MMP1.annot \
                --subject sub-0$sub \
                --hemi $hemi \
                --identity \
                --temp sub-0$sub/mri/T1.mgz \
                --o sub-0$sub/mri_glasser/$hemi.HCP-MMP1.nii.gz \
                --proj frac 0 1 0.01
        done
    done
else
    for (( sub=1; sub <= 5; sub=sub+1 )); do
        mkdir -p sub-0$sub/mri_glasser/
        for hemi in lh rh; do
            mri_label2vol \
                --annot sub-0$sub/label_glasser/$hemi.HCP-MMP1.annot \
                --subject sub-0$sub \
                --hemi $hemi \
                --identity \
                --temp sub-0$sub/mri/T1.mgz \
                --o sub-0$sub/mri_glasser/$hemi.HCP-MMP1.nii.gz
        done
    done
fi

# Optionally plot atlas on T1 template
if  [[ $2 = "-p" ]]; then
    freeview sub-01/mri/T1.mgz sub-01/mri_glasser/lh.HCP-MMP1.nii.gz 
fi
