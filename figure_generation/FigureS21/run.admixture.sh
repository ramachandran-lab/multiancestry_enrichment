#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --time=120:00:00

#SBATCH --array=1,4,67,89,54,43,66,42,99,23
mkdir ADMIXTURE_K_3_SEED_$SLURM_ARRAY_TASK_ID
cd ADMIXTURE_K_3_SEED_$SLURM_ARRAY_TASK_ID/
.././admixture -s $SLURM_ARRAY_TASK_ID african_pruned.bed 3
