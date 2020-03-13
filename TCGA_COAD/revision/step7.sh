#!/bin/bash
#SBATCH --job-name=tcga_m
#SBATCH --time=06:00:00
#SBATCH --nodes=20
#SBATCH -A newcomb_p

ml R/3.3.3-foss-2016b-fh2

mpirun -n 1 R CMD BATCH step7_microbiomeGWAS_22ct_kernel_L2_and_cor.R Step7_tcga_22ct_L2.Rout
