#!/bin/sh
#                     # lines starting with #SBATCH is an instruction to the job scheduler
#SBATCH --job-name=hailey   # Job name
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=haileyjp@stanford.edu # Where to send mail  
#SBATCH --nodes=1                   # Use one node
#SBATCH --ntasks=1                  # Run a single task
#SBATCH --mem-per-cpu=50gb           # Memory per processor
#SBATCH --qos long
#SBATCH --time=7-00:00:00            # Time limit hrs:min:sec
#SBATCH --output=array_%A-%a.out    # Standard output and error log
#SBATCH --array=1               # Array range
#SBATCH -n 1                # Array range
#SBATCH -C CPU_GEN:MLN          # Allocate a specific type of node

date
hostname

ml R/4.2.0

Rscript run-mcmc.R ${SLURM_ARRAY_TASK_ID}

wait