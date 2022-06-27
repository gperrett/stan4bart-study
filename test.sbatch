#!/bin/bash
#
#SBATCH --job-name=myJobarrayTest
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=12GB
#SBATCH --output=block1_%A_%a.out
#SBATCH --error=block1_%A_%a.err
module purge
module load r/gcc/4.1.2

Rscript test.R
