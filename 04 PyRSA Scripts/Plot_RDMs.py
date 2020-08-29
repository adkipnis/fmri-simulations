#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ROI RDMs (single plots or collage per subject) 

@author: alex
"""

def plot_rdm(rdm, do_rank_transform=False, pattern_descriptor=None,
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
        relative path to which the plot will be saved (if None: do not save plot)

    """
    
    plt.figure(dpi=dpi)
    if do_rank_transform:
        rdm = pyrsa.rdm.rank_transform(rdm)
    rdm_mat = rdm.get_matrices()
    if rdm.n_rdm  > 1:
        plot_dims = factorization(rdm.n_rdm)
        m = plot_dims[0]
        n = plot_dims[1]  
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
                            wspace=-0.75, hspace=1.1)
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
        ax.set_xticklabels(desc, {'fontsize': 0.8,
                     'fontweight': 'normal',
                     'verticalalignment': 'center',
                     'horizontalalignment':'center'})
        ax.set_yticks(np.arange(rdm.n_cond))
        ax.set_yticklabels(desc, {'fontsize': 0.8,
                     'fontweight': 'normal',
                     'verticalalignment': 'center',
                     'horizontalalignment':'right'})
        ax.tick_params(length=0, width=0)
        plt.ylim(rdm.n_cond - 0.5, -0.5)
        plt.xlim(-0.5, rdm.n_cond - 0.5)
        plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
                 rotation_mode="anchor")
    else:
        ax.set_xticks([])
        ax.set_yticks([])


def factorization(n):
    from math import gcd
    factors = []
    def get_factor(n):
        x_fixed = 2
        cycle_size = 2
        x = 2
        factor = 1
        while factor == 1:
            for count in range(cycle_size):
                if factor > 1: break
                x = (x * x + 1) % n
                factor = gcd(x - x_fixed, n)
            cycle_size *= 2
            x_fixed = x
        return factor

    while n > 1:
        next = get_factor(n)
        factors.append(next)
        n //= next

    return factors




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
beta_type       = 'R1_GLM' 
freesurfer_mri  = "mri_glasser" #Name of the directory in which subject specific volumetric ROI masks are saved by FreeSurfer
plot_roi        = False
plot_subject    = True


for sub in range(1, n_subs+1): 
    
    # Set subject-specific paths
    mask_dir = os.path.join(ds_dir, "derivatives", "freesurfer","sub-" + str(sub).zfill(2), freesurfer_mri)
    mask_dict = mask_utils.load_dict(os.path.join(mask_dir, "sub-" + str(sub).zfill(2) + "_mask_dict_EPI_disjoint.npy"))
    rdm_dir = os.path.join(ds_dir, "derivatives", "PyRSA", "rdms", "sub-"+str(sub).zfill(2))
    rdms = []   
    
    for roi_h in mask_dict.keys():        
        
        # Load RDM
        rdm_filename = os.path.join(rdm_dir, beta_type+"_RDM_"+roi_h)
        rdm = pyrsa.rdm.rdms.load_rdm(rdm_filename, file_type='hdf5')
        
        # Plot RDM for single ROI
        if plot_roi:
            plot_rdm(rdm, dpi=1200, do_rank_transform=True, pattern_descriptor='stim',
                                        cmap=None, rdm_descriptor='ROI', filename=beta_type+"_"+"sub-"+str(sub).zfill(2)+"_"+roi_h)
        
        # Collect single RDMs
        if isinstance(rdms, list):
            rdms = rdm
        else:
            rdms.append(rdm)
    
        
    # Plot RDMs for subject
    if plot_subject:
        plot_rdm(rdms, dpi=1200, do_rank_transform=True, pattern_descriptor='stim',
                                    cmap=None, rdm_descriptor='ROI', filename=beta_type+"_"+"sub-"+str(sub).zfill(2))
        
        
        
     