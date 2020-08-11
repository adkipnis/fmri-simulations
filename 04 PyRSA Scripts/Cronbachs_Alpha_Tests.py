#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cronbach's Alpha Tests

@author: alex
"""
import numpy as np

def cronbachs_alpha(stim_measurements):
    cmatrix = np.corrcoef(stim_measurements, rowvar=True)
    N = len(cmatrix)
    rs = list(cmatrix[np.triu_indices(N, k=1)])
    mean_r = np.mean(rs)
    alpha = (N * mean_r) / (1 + (N - 1) * mean_r)
    return alpha


def simulate_spatial_covariance(d_cors, data_sd):
    """ Simulate covariance matrix such that the cov of neighboring columns
    decreases with their distance
    
    Args:
        d_cors (np.array):
            vector with decreasing positive pearson correlation coefficients
            (starting with 1) - the index of this vector array corresponds to
            the distance between column indices in our data.
            
            e.g. if d_cors[2] = 0.3, data from column c will be be correlated
            to 0.3 with data from column c-2 and c+2.
            
        data_sd (float):
            standard deviation of signal or noise
    """
    n = len(d_cors)
    d_vars = d_cors * np.square(data_sd)
    
    # Initialize Covariance matrix
    cov_mat = np.zeros((n, n))
    
    # Fill its diagonals
    for i in range(n):
        diagonal = np.ones((1, n-i))[0]*d_vars[i]
        if i == 0:
            cov_mat += np.diag(diagonal, i)
        else:
            cov_mat += np.diag(diagonal, i) + np.diag(diagonal, -i)
    
    return cov_mat

"""
# SNR = Power(S) / Power(N) = E[S²] / E[N²] - Since noise is simulated with
# zero mean, we can use its variance = E[N²] - E[N]² (with E[N]² = 0) to
# calculate power. Similarly, E[S²] = var(S) + E(S)², hence:
# SNR = (var(S) + E(S)²) / var(N)
"""

# Simulation Parameters
SNR = 0.05                  # positive influence on calpha
signal_mean = 20            # negative influence on calpha (Note: 20 is an empirically observed value) 
signal_sd = 16              # does not influence calpha (Note: 16 is an empirically observed value) 
roi_size = 100              # decreases variance of calpha 
n_observations = 35         # positive influence on calpha
multivariate_noise = True   # does not influence calpha
d_cors = np.array([np.power(0.5,i) for i in range(roi_size)]) # does not influence calpha
n_sim = 100

# Initialize list variables
SNR_realized = []
empirical_alphas = []
noise_alphas = []

# Start simulation
for i in range(n_sim):
    # Simulate Signal under univariate normal distribution
    ground_truth = np.random.normal(signal_mean, signal_sd, roi_size)
    true_patterns = np.tile(ground_truth,(n_observations,1))
    
    # Simulate Noise 
    power_s = np.var(ground_truth) + np.square(np.mean(ground_truth))
    noise_sd = np.sqrt(power_s/SNR)
    if multivariate_noise:
        # ...under multivariate normal distribution
        noise_cov = simulate_spatial_covariance(d_cors, noise_sd)
        additive_noise = np.random.multivariate_normal(np.zeros(roi_size), noise_cov, n_observations)
        pure_noise_cov = simulate_spatial_covariance(d_cors, signal_sd)
        pure_noise = np.random.multivariate_normal(np.ones(roi_size)*signal_mean, pure_noise_cov, n_observations) + additive_noise
    else:
        # ...under univariate normal distribution
        additive_noise = np.random.normal(0, noise_sd, true_patterns.shape)
        pure_noise = np.random.normal(signal_mean, signal_sd, (n_observations, roi_size)) + additive_noise
    
    # Simulate data and pure noise with similar characteristics
    data = true_patterns + additive_noise
    
    # Check SNR
    power_n = np.var(additive_noise)
    SNR_realized.append(power_s/power_n)
    
    # Compare estimates of internal consistency 
    empirical_alphas.append(cronbachs_alpha(data))
    noise_alphas.append(cronbachs_alpha(pure_noise))
    

print("SNR of simulated data: \nμ =", np.mean(SNR_realized), "\nσ =", np.std(SNR_realized))
print("\nCronbach's alpha of simulated data: \nμ =", np.mean(empirical_alphas), "\nσ =", np.std(empirical_alphas))
print("\nCronbach's alpha of simulated noise: \nμ =", np.mean(noise_alphas), "\nσ =", np.std(noise_alphas))
