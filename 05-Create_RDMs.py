def stim_to_label(stim_list, label_dict, crop=False):
    '''
    Transform beta_id to object labels
    
    Inputs:
    - stim_list (list): sorted list of all GLM predictors
    - label_dict (dict): {synset ID : object label}
    - crop (bool): Use only first term in object label before the comma

    Output:
    - label_list (list): sorted list of the stimulus object labels corresponding to GLM predictors
    '''
    label = []
    
    if crop:
        for i in range(len(stim_list)):
            synset = 'n0'+stim_list[i].split('.')[0]
            label.append(label_dict[synset].split(',')[0])
    else:
        for i in range(len(stim_list)):
            synset = 'n0'+stim_list[i].split('.')[0]    
            label.append(label_dict[synset])
    label_list = np.array(label)
    return label_list

###############################################################################

import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pyrsa

# Set directories and specify ROIs
ds_dir = "/home/alex/ds001246/"
txt_dir = "/home/alex/templateflow/tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt" #directory of mask descriptors
mask_dir = os.path.join(ds_dir, "derivatives", "ROI_masks")
label_dict = np.load(os.path.join(ds_dir, "stimulus_label_dictionary.npy"),allow_pickle='TRUE').item()
n_subs = len(glob.glob(ds_dir + os.sep + "sub*"))
beta_type = 'SPM_s' 
rename_stim = True

sub = 1
for sub in range(1, n_subs+1): 
    # Set respective paths
    mask_dict_d_path = os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_mask_dict_T2*w_disjunct.npy")
    mask_dict_d = np.load(mask_dict_d_path,allow_pickle='TRUE').item()
    dataset_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "datasets", "sub-"+str(sub).zfill(2))
    res_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "noise", "sub-"+str(sub).zfill(2))
    
    # Load datasets and 
    roi_h = 'FFC_left'
    for roi_h in mask_dict_d.keys():        
       
        # Load residuals and estimate precision matrix
        residuals_filename = os.path.join(res_dir,"Residuals_"+roi_h+"_"+beta_type+".npy")
        residuals = np.load(residuals_filename)
        precision_matrix = pyrsa.data.noise.prec_from_residuals(residuals, dof=None)
        
        # Load dataset
        dataset_filename = os.path.join(dataset_dir,"RSA_dataset_"+roi_h+"_"+beta_type)
        dataset = pyrsa.data.dataset.load_dataset(dataset_filename, file_type='hdf5')    
        
        # Optionally rename stimulus IDs
        if rename_stim:
            stim_list = dataset.obs_descriptors['stim']
            label_list = stim_to_label(stim_list, label_dict, crop=True)
            dataset.obs_descriptors['stim'] = label_list
        
        # Calculate RDM with crossnobis distance estimates    
        rdm = pyrsa.rdm.calc.calc_rdm(dataset, method='crossnobis', descriptor='stim', cv_descriptor='run', noise=precision_matrix)
        # TODO: Remember to pull request change in calc_rdm 
        pyrsa.vis.rdm_plot.show_rdm(rdm, do_rank_transform=False, pattern_descriptor='stim',
             cmap=None, rdm_descriptor=None)
        
        
     