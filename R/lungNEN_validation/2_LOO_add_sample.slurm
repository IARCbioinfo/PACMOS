#!/bin/bash
#SBATCH --job-name=L_all_add_sample
#SBATCH --partition=high_p
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=logs/L_all_add_sample.out
#SBATCH --error=logs/L_all_add_sample.err

echo "======================================"
echo "Running LOO MOFA pipeline (XXX)"
echo "======================================"

# -------------------------------
# Load environment
# -------------------------------
source ~/.bashrc

# Activate conda env with MOFA + R deps
conda activate mofa_env

echo "Working directory: $(pwd)"

# -------------------------------
# Run R pipeline script
# -------------------------------
Rscript 2_LOO_add_sample.R

echo "======================================"
echo "Finished"
echo "======================================"

