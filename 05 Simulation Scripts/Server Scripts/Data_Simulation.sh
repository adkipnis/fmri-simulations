#!/bin/bash
#
#
#SBATCH --account=nklab          # The account name for the job.
#SBATCH --job-name=Noise_shuffling    # The job name.
#SBATCH -c 8                     # The number of cpu cores to use.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=8gb        # The memory the job will use per cpu core.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ak4572@columbia.edu


module load matlab
matlab -nodisplay -nosplash - nodesktop -r "cd /moto/home/ak4572/; try, run ('/moto/home/ak4572/Noise_shuffling_moto.m'); catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit"
matlab -nodisplay -nosplash - nodesktop -r "cd /moto/home/ak4572/; try, run ('/moto/home/ak4572/GLM_on_sim_moto.m'); catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit"
