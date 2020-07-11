import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import pyrsa

# Set directories and specify ROIs
ds_dir = "/home/alex/ds001246/"
txt_dir = "/home/alex/templateflow/tpl-Glasser/HCP-MMP1_on_MNI152_ICBM2009a_nlin.txt" #directory of mask descriptors
mask_dir = os.path.join(ds_dir, "derivatives", "ROI_masks")
n_subs = len(glob.glob(ds_dir + os.sep + "sub*"))


residuals_filename = os.path.join(res_output_dir,"Residuals_"+roi_h+"_"+beta_type)
            np.save(residuals_filename, np.vstack((residuals_superset[:])))



sub = 1
for sub in range(1, n_subs+1): 
    ### 1. Make Subject-specific ROI dictionary
    # Set respective paths to atlas and betas
    mask_dict_d_path = os.path.join(mask_dir, "Native","sub-" + str(sub).zfill(2) + "_mask_dict_T2*w_disjunct.npy")
    dataset_dir = os.path.join(ds_dir, "derivatives", "PYRSA", "datasets", "sub-"+str(sub).zfill(2))
    
    
    # Load mask dictionary
    mask_dict_d = np.load(mask_dict_d_path,allow_pickle='TRUE').item()
    
    # Load RSA data set
    roi_h = 'FFC_left'
    for roi_h in mask_dict_d.keys():        
        dataset = pyrsa.data.dataset.load_dataset(os.path.join(dataset_dir,"RSA_dataset_"+roi_h), file_type='hdf5')
        rdm = pyrsa.rdm.calc.calc_rdm(dataset, method='crossnobis', descriptor='stim', cv_descriptor='run', noise=None)
        # TODO: Remember to pull request change in calc_rdm 
        pyrsa.vis.rdm_plot.show_rdm(rdm, do_rank_transform=False, pattern_descriptor='stim',
             cmap=None, rdm_descriptor=None)
        
        
        # rdm.getmatrix()
        # 