#!/bin/bash
#SBATCH --job-name=RStanJob
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --partition=short

# Load R and any required modules (check module list)
module load R
# module load gcc # if needed for C++ compilation

# Run the R script
Rscript model.R


