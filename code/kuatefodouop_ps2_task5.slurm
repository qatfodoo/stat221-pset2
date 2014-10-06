#!/bin/bash
#SBATCH -J mcmc-mispecified # name for job array
#SBATCH -o task5.out #Standard output
#SBATCH -e task5.err #Standard error
#SBATCH -p general #Partition
#SBATCH -t 02:00:00 #Running time of 20 mins.
#SBATCH --mem-per-cpu 3000 #Memory request
#SBATCH -n 1 #Number of cores
#SBATCH -N 1 #All cores on one machine
#SBATCH --mail-user=kuatefodouop@fas.harvard.edu
#SBATCH --mail-type=ALL    # When to send mail

# argument = job id

Rscript kuatefodouop_ps2_task5.R $SLURM_ARRAY_TASK_ID
