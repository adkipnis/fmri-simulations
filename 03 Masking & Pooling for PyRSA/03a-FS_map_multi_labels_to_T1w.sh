#!/bin/bash
# Map a list of ROI surface labels to native volumetric space.
# Rename ROIs according to your goals.

# Prerequisite: FS_transform_annots_to_labels.sh

ROIs="V1 V2 V3 V4 VMV1 VMV2 VMV3 PHA1 PHA2 PHA3 VVC FFC TF PeEc MT MST"
S="_ROI"

cd $SUBJECTS_DIR

for (( sub=1; sub <= 5; sub=sub+1 )); do
    mkdir -p sub-0$sub/mri_glasser
    for hemi in lh.L_ rh.R_; do
        for ROI in $ROIs; do
            mri_label2vol \
            --label sub-0$sub/label_glasser/$hemi$ROI$S.label \
            --subject sub-0$sub \
            --hemi lh \
            --identity \
            --temp sub-0$sub/mri/T1.mgz \
            --o sub-0$sub/mri_glasser/$hemi$ROI$S.nii.gz
        done
    done
done

