#!/bin/bash
#SBATCH --job-name=UTRSoyFACEopt
#SBATCH --time=120:00:00
#SBATCH -p lowmem
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=END
#SBATCH --mail-user=xpiao2@illinois.edu

Rscript multiyear_partitioning_senescence_coefs_optimization_lsrp.R
