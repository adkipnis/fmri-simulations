#!/bin/bash
#
#
#SBATCH --account=nklab          # The account name for the job.
#SBATCH --job-name=full_sim    # The job name.
#SBATCH -c 8                     # The number of cpu cores to use.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=120:00:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=8gb        # The memory the job will use per cpu core.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ak4572@columbia.edu

module load matlab
module load anaconda/3-5.3.1
source activate fmri-sim
    
for (( time=1; time <= 3; time=time+1 )); do 

    matlab -nodisplay -nosplash - nodesktop -r "cd /moto/home/ak4572/; try, run ('/moto/home/ak4572/Noise_shuffling_moto.m'); catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit"
    matlab -nodisplay -nosplash - nodesktop -r "cd /moto/home/ak4572/; try, run ('/moto/home/ak4572/GLM_on_sim_moto.m'); catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit"
    python 01-Pool_simulation_results_moto.py
    python 02-Create_sim_pyrsa_dataset_moto.py
    python 03-Create_sim_RDMs_moto.py
    python 04-Test_sim_fixed_inference_moto.py

done
