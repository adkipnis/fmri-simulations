# fMRI Simulations for Testing RSA 
1. Analyses of fMRI data from Horikawa & Kamitani (2017), available on [openneuro](https://openneuro.org/datasets/ds001246/versions/1.2.1)
2. Based on (1), synthetic data with known noise structure are simulated
3. Simulations are used to test statistical inference of [RSA 3.0](https://github.com/rsagroup/pyrsa) under different conditions:
   - Signal-to-Noise ratio
   - Noise Normalization procedure
   - Number of fMRI runs
   - Number of unique representations

### Process Flowchart
<img src="https://github.com/adkipnis/fmri-simulations/blob/master/flowchart.png" width="700">
