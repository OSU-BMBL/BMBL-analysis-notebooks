#!/bin/bash
#SBATCH --account=PAS1475
#SBATCH --job-name=expr_analysis
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --output=slurm-%j.out

echo "Starting analysis at $(date)"
echo "Working directory: $(pwd)"

# Load R module (requires gcc on Ascend)
module load gcc/12.3.0
module load R/4.4.0

echo "R version:"
R --version | head -1

# Run analysis
Rscript analysis.R

echo "Analysis completed at $(date)"
