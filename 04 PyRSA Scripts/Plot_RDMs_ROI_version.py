#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ROI RDMs (single plots or collage per subject) 

@author: alex
"""

def plot_rdm(rdm, title=None, do_rank_transform=False, pattern_descriptor=None,
             cmap=None, rdm_descriptor=None, dpi=300, filename=None):
    """shows an rdm object

    Parameters
    ----------
    rdm : pyrsa.rdm.RDMs
        RDMs object to be plotted
    do_rank_transform : bool
        whether we should do a rank transform before plotting
    pattern_descriptor : String
        name of a pattern descriptor which will be used as an axis label
    rdm_descriptor : String
        name of a rdm descriptor which will be used as a title per RDM
    cmap : color map
        colormap or identifier for a colormap to be used
        conventions as for matplotlib colormaps
    dpi : int
        dots per inch (determines visual resolution of plots)
    filename : str
        relative path to which the plot will be saved
        (if None: do not save plot)

    """
    
    plt.figure(dpi=dpi)
    plt.suptitle(title, fontsize=8)
    if do_rank_transform:
        rdm = pyrsa.rdm.rank_transform(rdm)
    rdm_mat = rdm.get_matrices()
    if rdm.n_rdm  > 1:
        m = 3
        n = 2
        for idx in range(rdm.n_rdm):
            plt.subplot(n, m, idx + 1)
            plt.imshow(rdm_mat[idx], cmap=cmap)
            _add_descriptor_labels(rdm, pattern_descriptor)
            if rdm_descriptor in rdm.rdm_descriptors:
                plt.title(rdm.rdm_descriptors[rdm_descriptor][idx], fontsize=4)
            elif isinstance(rdm_descriptor, str):
                plt.title(rdm_descriptor)
        plt.subplot(n, m, n * m)
        plt.imshow(np.mean(rdm_mat, axis=0), cmap=cmap)
        _add_descriptor_labels(rdm, pattern_descriptor)
        plt.title('Average', fontsize=4)
    elif rdm.n_rdm == 1:
        plt.imshow(rdm_mat[0], cmap=cmap)
        _add_descriptor_labels(rdm, pattern_descriptor)
        if rdm_descriptor in rdm.rdm_descriptors:
            plt.title(rdm.rdm_descriptors[rdm_descriptor][0], fontsize=4)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None,
                            wspace=0.4, hspace=0.5)
    if isinstance(filename, str):    
        fig = plt.gcf()
        fig.savefig(filename, bbox_inches='tight')

def _add_descriptor_labels(rdm, descriptor, ax=None):
    """ adds a descriptor as ticklabels """
    if ax is None:
        ax = plt.gca()
    if descriptor is not None:
        
        desc = rdm.pattern_descriptors[descriptor]
        ax.set_xticks(np.arange(rdm.n_cond))
        ax.set_xticklabels(desc, {'fontsize': 2,
                     'fontweight': 'normal',
                     'verticalalignment': 'center',
                     'horizontalalignment':'center'})
        ax.set_yticks(np.arange(rdm.n_cond))
        ax.set_yticklabels(desc, {'fontsize': 2,
                     'fontweight': 'normal',
                     'verticalalignment': 'center',
                     'horizontalalignment':'right'})
        ax.tick_params(length=2, width=0.5)
        plt.ylim(rdm.n_cond - 0.5, -0.5)
        plt.xlim(-0.5, rdm.n_cond - 0.5)
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
                 rotation_mode="anchor")
    else:
        ax.set_xticks([])
        ax.set_yticks([])


###############################################################################
import glob
import os
import mask_utils
import numpy as np
import matplotlib.pyplot as plt
import pyrsa

# Set directories and specify ROIs
ds_dir          = "/home/alex/Datasets/ds001246/"
n_subs          = len(glob.glob(ds_dir + os.sep + "sub*"))
beta_types      = ['SPM_6', 'R1_GLM', 'Res_0']
freesurfer_mri  = "mri_glasser" 
mask_dir        = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" +
                               str(1).zfill(2), freesurfer_mri)
mask_dict       = mask_utils.load_dict(
                    os.path.join(mask_dir, "sub-" +str(1).zfill(2) +
                                 "_mask_dict_EPI_disjoint.npy"))
roi_h_list      = list(mask_dict.keys())

for beta_type in beta_types:
    for roi_h in roi_h_list:
        rdms = [] 
        
        for sub in range(1, n_subs+1): 
            
            # Set subject-specific paths
            rdm_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "rdms",
                                   "sub-"+str(sub).zfill(2))
            
            # Load RDM
            rdm_filename = os.path.join(rdm_dir, beta_type+"_RDM_"+roi_h)
            rdm = pyrsa.rdm.rdms.load_rdm(rdm_filename, file_type='hdf5')
            rdm.rdm_descriptors.update({"sub": np.array(["sub-" +
                                                         str(sub).zfill(2)])})    
            # Collect single RDMs
            if isinstance(rdms, list):
                rdms = rdm
            else:
                rdms.append(rdm)
             
        # Plot RDMs for subject
        plot_rdm(rdms, dpi=1200, do_rank_transform=True, title=beta_type+" "+
                 roi_h, pattern_descriptor='stim', cmap=None,
                 rdm_descriptor='sub', filename=beta_type+"_"+roi_h)
        
        
        
     