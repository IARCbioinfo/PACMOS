#!/bin/bash
#SBATCH --job-name=L_M
#SBATCH --partition=high_p
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --output=logs/L_M_plot_cluster.out
#SBATCH --error=logs/L_M_plot_cluster.err

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
Rscript 4_LOO_plot_query_sample_cluster.R

echo "======================================"
echo "Finished"
echo "======================================"

