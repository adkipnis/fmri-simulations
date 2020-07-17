#!/bin/bash
# Register Glasser surface atlas from standard space to native space (fsaverage -> fsnative)

# Prerequisite: fmriprep's freesurfer output is in the derivatives directory and the Glasser surface atlas (rh/lh.HCP-MMP1.annot) in the .../fsaverage/label/ directory

# Options: -p - plot atlas on fsaverage surface
#          -s - plot atlas on native surface of sub-01
#          -ps - first -p then -s

cd $SUBJECTS_DIR     

# Register
for (( sub=1; sub <= 5; sub=sub+1 )); do
    mkdir -p sub-0$sub/label_glasser/
    for hemi in lh rh; do
        mri_surf2surf \
            --srcsubject fsaverage \
            --trgsubject sub-0$sub \
            --hemi lh \
            --sval-annot fsaverage/label/$hemi.HCP-MMP1.annot \
            --tval  sub-0$sub/label_glasser/$hemi.HCP-MMP1.annot
    done
done


# Optionally plot annots on surface
if  [[ $1 = "-p" ]]; then
	freeview -f fsaverage/surf/lh.white:annot=fsaverage/label/lh.HCP-MMP1.annot
fi

if  [[ $1 = "-s" ]]; then
	freeview -f sub-01/surf/lh.white:annot=sub-01/label_glasser/lh.HCP-MMP1.annot
fi

if  [[ $1 = "-ps" ]]; then
	freeview -f fsaverage/surf/lh.white:annot=fsaverage/label/lh.HCP-MMP1.annot
	freeview -f sub-01/surf/lh.white:annot=sub-01/label_glasser/lh.HCP-MMP1.annot
fi



