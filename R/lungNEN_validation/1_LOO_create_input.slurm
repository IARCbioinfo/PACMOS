#!/bin/bash
#SBATCH --job-name=L_create_input
#SBATCH --partition=high_p
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=logs/L_create_input.out
#SBATCH --error=logs/L_create_input.err

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
Rscript 1_LOO_create_input.R

echo "======================================"
echo "Finished"
echo "======================================"

