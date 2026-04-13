#!/bin/bash
#SBATCH --job-name=L_GM_run_mofa
#SBATCH --partition=high_p
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=logs/L_GM_run_mofa.out
#SBATCH --error=logs/L_GM_run_mofa.err

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
Rscript 3_LOO_run_mofa.R

echo "======================================"
echo "Finished"
echo "======================================"

