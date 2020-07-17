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
        hemi_short=${hemi:0:2}
        for ROI in $ROIs; do   
            mri_label2vol \
                --label sub-0$sub/label_glasser/$hemi$ROI$S.label \
                --subject sub-0$sub \
                --hemi $hemi_short \
                --identity \
                --temp sub-0$sub/mri/T1.mgz \
                --o sub-0$sub/mri_glasser/$hemi$ROI$S.f.nii.gz \
                --proj frac 0 1 0.01
            
            mri_binarize \
                --dilate 1 \
                --erode 1 \
                --i sub-0$sub/mri_glasser/$hemi$ROI$S.f.nii.gz \
                --o sub-0$sub/mri_glasser/$hemi$ROI$S.fb.nii.gz \
                --min 1         
            
            mris_calc \
                -o sub-0$sub/mri_glasser/$hemi$ROI$S.fbr.nii.gz \
                sub-0$sub/mri_glasser/$hemi$ROI$S.fb.nii.gz \
                mul sub-0$sub/mri/$hemi_short.ribbon.mgz    
            
        done
    done
done

# Optionally plot ROI on T1 template
if  [[ $1 = "-p" ]]; then
    freeview sub-01/mri/T1.mgz sub-01/mri_glasser/lh.L_V1_ROI.fbr.nii.gz \
    -f sub-01/surf/lh.white:label=sub-01/label_glasser/lh.L_V1_ROI.label \
    -f sub-01/surf/lh.pial:label=sub-01/label_glasser/lh.L_V1_ROI.label
fi

