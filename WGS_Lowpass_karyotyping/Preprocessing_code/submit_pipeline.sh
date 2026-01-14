#!/bin/bash
#
# LP-WGS Karyotyping Pipeline - Master Submission Script
# Submits all pipeline steps with proper dependencies
#
# Usage: bash submit_pipeline.sh
#
# This script will:
#   1. Run setup (00_setup.sh)
#   2. Submit FASTQ concatenation job (01_concat_fastq.sh)
#   3. Submit processing job with dependency (02_process_sample.sh)
#

set -e

###############################################
# Load configuration
###############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

echo "=============================================="
echo "LP-WGS Karyotyping Pipeline"
echo "=============================================="
echo ""

# Validate configuration
if ! validate_config; then
    exit 1
fi

echo "Configuration validated successfully."
echo ""

###############################################
# Step 1: Run setup
###############################################
echo "=== Step 1: Running setup ==="
bash ${SCRIPT_DIR}/00_setup.sh

###############################################
# Step 2: Submit FASTQ concatenation
###############################################
echo ""
echo "=== Step 2: Submitting FASTQ concatenation jobs ==="

# Update array range dynamically
ACTUAL_SAMPLES=$(wc -l < ${SAMPLE_LIST})
echo "Submitting array job for ${ACTUAL_SAMPLES} samples..."

CONCAT_JOB=$(sbatch \
    --account=${SLURM_ACCOUNT} \
    --array=1-${ACTUAL_SAMPLES} \
    --output=${WD}/logs/concat_%A_%a.out \
    --error=${WD}/logs/concat_%A_%a.err \
    ${SCRIPT_DIR}/01_concat_fastq.sh \
    | awk '{print $NF}')

echo "FASTQ concatenation job submitted: ${CONCAT_JOB}"

###############################################
# Step 3: Submit processing with dependency
###############################################
echo ""
echo "=== Step 3: Submitting processing jobs (with dependency) ==="

PROCESS_JOB=$(sbatch \
    --account=${SLURM_ACCOUNT} \
    --array=1-${ACTUAL_SAMPLES} \
    --dependency=afterok:${CONCAT_JOB} \
    --output=${WD}/logs/process_%A_%a.out \
    --error=${WD}/logs/process_%A_%a.err \
    ${SCRIPT_DIR}/02_process_sample.sh \
    | awk '{print $NF}')

echo "Processing job submitted: ${PROCESS_JOB}"

###############################################
# Summary
###############################################
echo ""
echo "=============================================="
echo "Pipeline submitted successfully!"
echo "=============================================="
echo ""
echo "Job IDs:"
echo "  FASTQ concatenation: ${CONCAT_JOB}"
echo "  Sample processing:   ${PROCESS_JOB} (depends on ${CONCAT_JOB})"
echo ""
echo "Monitor jobs with:"
echo "  squeue -u \$USER"
echo ""
echo "Check logs in:"
echo "  ${WD}/logs/"
echo ""
echo "After completion, proceed with karyotype analysis:"
echo "  1. Copy counts files to analysis directory"
echo "  2. Run 1_digital_karyotype_analysis.rmd"
