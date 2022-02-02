#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=samuel_smith1@brown.edu
module load R/3.4.3
module load plink/2.00
Rscript genee_script.R british 1 Height_norm
Rscript genee_script_nobuffer.R british 1 Height_norm