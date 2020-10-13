#!/bin/bash
#
#
#SBATCH --account=nklab          # The account name for the job.
#SBATCH --job-name=Sim_Test      # The job name.
#SBATCH -c 8                     # The number of cpu cores to use.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00              # The time the job will take to run.
#SBATCH --mem-per-cpu=8gb        # The memory the job will use per cpu core.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ak4572@columbia.edu

module load anaconda/3-5.3.1
source activate fmri-sim
python 03-Create_sim_RDMs_moto.py
python 04-Test_sim_fixed_inference_moto.py
