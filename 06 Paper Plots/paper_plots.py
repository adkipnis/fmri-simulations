#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 10:49:28 2021

@author: alex
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

results_folder= '/home/alex/Code/fmri-sim/06 Paper Plots/'

def plot_paper(simulation_folder=results_folder, savefig=True):
    labels = pd.read_csv(os.path.join(simulation_folder, 'labels.csv'),
                         index_col = 0)
    means = np.load(os.path.join(simulation_folder, 'means.npy'))
    stds = np.load(os.path.join(simulation_folder, 'stds.npy'))
    stds = stds[:,:,0:20, 0:20] #remove noise ceiling stds
    means = means[:len(labels)]
    # remove nan entries
    idx_nan = ~np.any(np.isnan(means[:, :, 0]), axis=1)
    labels = labels[list(idx_nan)]
    means = means[idx_nan]
    stds = stds[idx_nan]
    # compute statistics
    true_std = np.nanstd(means, axis=1)
    std_mean = np.nanmean(stds, axis=1)
    std_mean = np.array([np.diag(i) for i in std_mean])
    std_mean = np.sqrt(std_mean)  # those are actually variances!
    std_var = np.nanvar(stds, axis=1)
    std_var = np.array([np.diag(i) for i in std_var])
    std_relative = std_mean / true_std
    std_std = np.sqrt(std_var)
    snr = (np.var(np.mean(means, axis=1), axis=1)
           / np.mean(np.var(means, axis=1), axis=1))
    # seaborn based plotting
    # create full data table
    data_df = pd.DataFrame()
    
    for i_model in range(means.shape[2]):
        labels['true_std'] = true_std[:, i_model]
        labels['std_mean'] = std_mean[:, i_model]
        labels['std_var'] = std_var[:, i_model]
        labels['std_relative'] = std_relative[:, i_model]
        labels['std_std'] = std_std[:, i_model]
        labels['model_layer'] = i_model
        labels['snr'] = snr
        labels['log-snr'] = np.log10(snr)
    data_df = data_df.append(labels)
    data_df_sans_50 = data_df[data_df["n_stim"]!=50]


    with sns.axes_style('ticks'):
        sns.set_context('paper', font_scale=2)
        g1_m = sns.catplot(data=data_df, legend=False,
                           x='GT', y='log-snr', hue='nsf', kind='point',
                           ci='sd', palette='ch:s=.25,rot=-.25_r',
                           aspect = 2, dodge=.8)
        g1_m.add_legend(
            frameon=False, title='noise scaling factor',
            bbox_to_anchor=(1.0, 1.0), loc=2)
        
        g1_m.set_xlabels('ground truth model')
        g1_m.set_ylabels('signal to noise ratio')
        g1_m.set_xticklabels(rotation=90)
        plt.ylim([-2, 3.5])
        plt.yticks([-2, -1, 0, 1, 2, 3],
                   ['10^-2', '10^-1', '10^0', '10^1', '10^2', '10^3'])
        sns.despine(trim=True)


        g2_m = sns.catplot(data=data_df, legend=False, x='NN', y='log-snr',
                           order=["none", "univariate", "time-based",
                                  "instance-based",], hue='nsf', kind='point',
                           ci='sd', palette='ch:s=.25,rot=-.25_r', dodge=.2)
        g2_m.add_legend(
            frameon=False, title='noise scaling factor',
            bbox_to_anchor=(1.0, 1.0), loc=2)
        g2_m.set_xlabels('noise normalization method')
        g2_m.set_ylabels('signal to noise ratio')
        g2_m.set_xticklabels(rotation=90, size = 12)
        plt.xlabel('noise normalization method', labelpad=10)
        plt.ylim([-2, 3.5])
        plt.yticks([-2, -1, 0, 1, 2, 3],
                   ['10^-2', '10^-1', '10^0', '10^1', '10^2', '10^3'])
        sns.despine(trim=True)


        g3_m = sns.catplot(data=data_df, legend=False, x='n_stim', y='log-snr',
                           col = 'nsf', hue='n_runs', kind='point',
                           ci='sd', palette='Blues_d', dodge=.4)
        g3_m.add_legend(
            frameon=False, title='# of runs',
            bbox_to_anchor=(1.0, 1.0), loc=2)
        g3_m.set_xlabels('# of stimuli')
        g3_m.set_ylabels('signal to noise ratio')
        # plt.xlabel('noise normalization method', labelpad=10)
        plt.ylim([-2, 4.5])
        plt.yticks([-2, -1, 0, 1, 2, 3, 4],
                   ['10^-2', '10^-1', '10^0', '10^1', '10^2', '10^3', '10^4'])
        sns.despine(trim=True)
    
        n_stim_levels = data_df_sans_50["n_stim"].unique()
        expected_inflation = np.sqrt(1 / (1 - n_stim_levels/50))
        g4_m = sns.catplot(data=data_df_sans_50, legend=False, x='n_stim',
                           y='std_relative', col = 'nsf', hue='n_runs',
                           kind='point', ci='sd', palette='Blues_d', dodge=.4)
        for j in range(3):
            g4_m.axes[0, j].plot(np.linspace(0,3,4), expected_inflation, 'r--')
        g4_m.add_legend(
            frameon=False, title='# of runs',
            bbox_to_anchor=(1.0, 1.0), loc=2)
        g4_m.set_ylabels(r'relative uncertainty $[\sigma_{boot}/\sigma_{true}]$')
        g4_m.set_xlabels('# of stimuli')
        sns.despine(trim=True)

        if savefig:
            g1_m.fig.savefig(results_folder+'SNR_roi.svg',
                             bbox_inches='tight')
            g2_m.fig.savefig(results_folder+'SNR_nn.svg',
                             bbox_inches='tight')
            g3_m.fig.savefig(results_folder+'SNR_stim_runs.svg',
                             bbox_inches='tight')
            g4_m.fig.savefig(results_folder+'std_rel_stim_runs.svg',
                             bbox_inches='tight')

plot_paper()
