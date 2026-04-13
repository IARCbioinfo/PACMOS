#!/bin/bash
#SBATCH --job-name=M_M_s3
#SBATCH --partition=high_p
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=logs/M_M_s3.out
#SBATCH --error=logs/M_M_s3.err

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
Rscript 5_run_MOFA_projection_LOO.R

echo "======================================"
echo "Finished"
echo "======================================"

