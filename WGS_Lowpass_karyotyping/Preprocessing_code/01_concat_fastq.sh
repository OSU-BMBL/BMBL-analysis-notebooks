#!/bin/bash
#
# LP-WGS Karyotyping Pipeline - FASTQ Concatenation
# Handles multi-lane samples by concatenating FASTQs from different flowcells
#
# Usage: sbatch 01_concat_fastq.sh
#
# Prerequisites: Run 00_setup.sh first
#

###############################################
# Load configuration
###############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

#SBATCH --account ${SLURM_ACCOUNT}
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16GB
#SBATCH --job-name=lpwgs_concat
#SBATCH --output=logs/concat_%A_%a.out
#SBATCH --error=logs/concat_%A_%a.err
#SBATCH --array=1-${NUM_SAMPLES}

set -e

###############################################
# Get sample name for this array task
###############################################
cd ${WD}
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

if [ -z "$SAMPLE" ]; then
    echo "ERROR: Could not get sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "=== Processing sample: ${SAMPLE} ==="
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"

###############################################
# Check if output already exists
###############################################
OUT_R1="${WD}/merged_fastq/${SAMPLE}_R1.fq.gz"
OUT_R2="${WD}/merged_fastq/${SAMPLE}_R2.fq.gz"

if [ -f "${OUT_R1}" ] && [ -f "${OUT_R2}" ]; then
    echo "Output files already exist. Skipping."
    echo "  ${OUT_R1}"
    echo "  ${OUT_R2}"
    exit 0
fi

###############################################
# Find FASTQ files for this sample
###############################################
SAMPLE_DIR="${DATA_DIR}/${SAMPLE}"

if [ ! -d "${SAMPLE_DIR}" ]; then
    echo "ERROR: Sample directory not found: ${SAMPLE_DIR}"
    exit 1
fi

# Find R1 and R2 files (pattern: *_1.fq.gz and *_2.fq.gz)
R1_FILES=(${SAMPLE_DIR}/*_1.fq.gz)
R2_FILES=(${SAMPLE_DIR}/*_2.fq.gz)

echo "Found ${#R1_FILES[@]} R1 file(s):"
printf '  %s\n' "${R1_FILES[@]}"
echo "Found ${#R2_FILES[@]} R2 file(s):"
printf '  %s\n' "${R2_FILES[@]}"

# Validate
if [ ${#R1_FILES[@]} -eq 0 ] || [ ! -f "${R1_FILES[0]}" ]; then
    echo "ERROR: No R1 FASTQ files found in ${SAMPLE_DIR}"
    exit 1
fi

if [ ${#R2_FILES[@]} -eq 0 ] || [ ! -f "${R2_FILES[0]}" ]; then
    echo "ERROR: No R2 FASTQ files found in ${SAMPLE_DIR}"
    exit 1
fi

if [ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]; then
    echo "ERROR: Mismatched number of R1 (${#R1_FILES[@]}) and R2 (${#R2_FILES[@]}) files"
    exit 1
fi

###############################################
# Concatenate or symlink
###############################################
if [ ${#R1_FILES[@]} -gt 1 ]; then
    # Multiple lanes - concatenate
    echo "Multiple lanes detected. Concatenating..."

    echo "Concatenating R1 files..."
    cat "${R1_FILES[@]}" > "${OUT_R1}"

    echo "Concatenating R2 files..."
    cat "${R2_FILES[@]}" > "${OUT_R2}"

    echo "Concatenation complete."
else
    # Single lane - create symlink to save space
    echo "Single lane detected. Creating symlinks..."

    ln -sf "${R1_FILES[0]}" "${OUT_R1}"
    ln -sf "${R2_FILES[0]}" "${OUT_R2}"

    echo "Symlinks created."
fi

###############################################
# Verify output
###############################################
echo ""
echo "Output files:"
ls -lh "${OUT_R1}" "${OUT_R2}"

echo ""
echo "=== Completed: ${SAMPLE} ==="
echo "End time: $(date)"
