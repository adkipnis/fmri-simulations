#!/bin/bash
# Apply ANTs registration of Glasser Atlas to individual T1w space 
      
data_path="/home/alex/ds001246/derivatives/"
input_path="/home/alex/templateflow/tpl-Glasser/MMP_in_MNI_corr.nii.gz"

for (( sub=1; sub <= 5; sub=sub+1 )); do
	
	anat_path=$data_path"fmriprep/sub-0"$sub"/anat/"
	ref_path=$anat_path"sub-0"$sub"_desc-brain_mask.nii.gz"
	transform_path=$anat_path"sub-0"$sub"_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5"
	output_path=$data_path"ROI_masks/T1w/sub-0"$sub"_Mask_T1w_Glasser.nii.gz"



	antsApplyTransforms \
		-i $input_path \
		-r $ref_path \
		-t $transform_path \
		-n NearestNeighbor -o \
		-o $output_path \
		-v 1
done

