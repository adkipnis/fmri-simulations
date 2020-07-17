#!/bin/bash
# Export Glasser surface atlas in native space to ROI labels

# Prerequisite: Successfully run FS_register_Glasser_annots.sh

# Options: -p - plot V1 on sub-01's surface


cd $SUBJECTS_DIR 

for (( sub=1; sub <= 5; sub=sub+1 )); do
    for hemi in lh rh; do
        mri_annotation2label \
            --subject sub-0$sub \
            --hemi $hemi \
            --annotation sub-0$sub/label_glasser/$hemi.HCP-MMP1.annot \
            --outdir sub-0$sub/label_glasser
    done
done

# Optionally plot V1 on surface
if  [[ $1 = "-p" ]]; then
	freeview -f sub-01/surf/lh.white:label=sub-01/label_glasser/lh.L_V1_ROI.label
fi
