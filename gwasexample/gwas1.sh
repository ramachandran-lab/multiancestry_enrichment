#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --mail-user=samuel_smith1@brown.edu
module load R/3.4.0
module load plink/2.00
plink2 --pgen european_plink_files/european_chr1.pgen --pvar  european_plink_files/european_chr1.pvar --psam european_plink_files/european_chr1.psam --glm hide-covar --pheno european_continuous_phenos.txt --pheno-name Height_norm --chr 1 --covar european_covars.txt --covar-name Sex,Age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --covar-variance-standardize --out /users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/redefining_replication/2.imputed_ukb_data/1.gwas/outputs/european/Height_norm.1