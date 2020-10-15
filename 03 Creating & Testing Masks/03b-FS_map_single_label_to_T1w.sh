#!/bin/bash
# Map a single ROI surface label to native volumetric space.
# Default ROI is lh V1, rename this variable if needed.

# Prerequisite: FS_transform_annots_to_labels.sh

# Options:  -f - fill ROI voxels between white matter and pial surface
#           -p - plot sub-01's left V1 on volume and surface

ROI=V1_ROI

cd $SUBJECTS_DIR 

if  [[ $1 = "-f" ]]; then
    for (( sub=1; sub <= 5; sub=sub+1 )); do
        mkdir -p sub-0$sub/mri_glasser
        for hemi in lh.L_ rh.R_; do
            mri_label2vol \
            --label sub-0$sub/label_glasser/$hemi$ROI.label \
            --subject sub-0$sub \
            --hemi lh \
            --identity \
            --temp sub-0$sub/mri/T1.mgz \
            --o sub-0$sub/mri_glasser/$hemi$ROI.nii.gz \
            --proj frac 0 1 0.01
        done
    done
else
    for (( sub=1; sub <= 5; sub=sub+1 )); do
        mkdir -p sub-0$sub/mri_glasser
        for hemi in lh.L_ rh.R_; do
            mri_label2vol \
            --label sub-0$sub/label_glasser/$hemi$ROI.label \
            --subject sub-0$sub \
            --hemi lh \
            --identity \
            --temp sub-0$sub/mri/T1.mgz \
            --o sub-0$sub/mri_glasser/$hemi$ROI.nii.gz
        done
    done
fi

# Optionally plot ROI on volume
if  [[ $2 = "-p" ]]; then
    freeview sub-01/mri/T1.mgz sub-01/mri_glasser/lh.L_$ROI.nii.gz \
    -f sub-01/surf/lh.white:label=sub-01/label_glasser/lh.L_$ROI.label \
    -f sub-01/surf/lh.pial:label=sub-01/label_glasser/lh.L_$ROI.label
fi
