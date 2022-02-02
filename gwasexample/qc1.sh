#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=samuel_smith1@brown.edu
module load R/3.4.0
module load plink/2.00
plink2 --pgen european_chr1_C.pgen --psam european_chr1_C.psam --pvar european_chr1_C.pvar --make-bed --out ../../6.HESS/UKB/plink_files/european_chr1_C
