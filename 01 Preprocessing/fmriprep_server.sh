#!/bin/bash
#
#SBATCH --account=nklab          # The account name for the job.
#SBATCH --job-name=b-fmriprep    # The job name.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8gb        # The memory the job will use per cpu core.
#SBATCH --output=fmriprep_%A-%a.out    # Standard output and error log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ak4572@columbia.edu
#SBATCH --output=slurm-output/slurm_eco_%a.out


for (( sub=1; sub <= 5; sub=sub+1 )); do
	module load singularity
	singularity run --cleanenv /moto/nklab/projects/singularity_images/fmriprep-stable.simg \
	    /moto/nklab/projects/ds001246/ /moto/nklab/projects/ds001246/derivatives \
	    participant \
	    --participant_label $sub \
	    --fs-license-file /moto/home/ak4572/fs_license.txt \
	    --output-spaces T1w \
	    --task perception\
	    --stop-on-first-crash
done

 
