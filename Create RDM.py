import numpy as np
from scipy import io
import matplotlib.pyplot as plt
import pyrsa

dataset = pyrsa.data.dataset(measurements #(numpy.ndarray): n_obs x n_channel 2d-array,
                             descriptors #(dict):           descriptors (metadata)
                             obs_descriptors #(dict):       observation descriptors (all are array-like with shape = (n_obs,...))
                             channel_descriptors #(dict):   channel descriptors (all are array-like with shape = (n_channel,...))
                             )