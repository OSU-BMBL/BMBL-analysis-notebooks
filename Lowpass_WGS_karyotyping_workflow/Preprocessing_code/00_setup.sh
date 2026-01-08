#!/bin/bash
#
# LP-WGS Karyotyping Pipeline - Setup Script
# Creates directory structure, sample list, and genome bins
#
# Usage: bash 00_setup.sh
#
# This script should be run once before starting the pipeline.
#

set -e  # Exit on error

###############################################
# Load configuration
###############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

# Validate configuration
if ! validate_config; then
    exit 1
fi

# Load required modules
load_modules

###############################################
# Create directory structure
###############################################
echo "=== Creating directory structure ==="
mkdir -p ${WD}/{merged_fastq,fastp_out,alignment,counts,logs,reference,scripts}
echo "Created directories in ${WD}"

###############################################
# Generate sample list from data directory
###############################################
echo "=== Generating sample list ==="
ls -d ${DATA_DIR}/*/ 2>/dev/null | xargs -n1 basename > ${SAMPLE_LIST}

# Count samples
FOUND_SAMPLES=$(wc -l < ${SAMPLE_LIST})
echo "Found ${FOUND_SAMPLES} samples:"
cat ${SAMPLE_LIST}
echo ""

# Reminder to update config
if [ ${FOUND_SAMPLES} -ne ${NUM_SAMPLES} ]; then
    echo "NOTE: Found ${FOUND_SAMPLES} samples but config.sh has NUM_SAMPLES=${NUM_SAMPLES}"
    echo "      Please update NUM_SAMPLES in config.sh if needed."
    echo ""
fi

###############################################
# Create genome bins
###############################################
if [ -f "${BINS_FILE}" ]; then
    echo "=== Genome bins file already exists: ${BINS_FILE} ==="
else
    echo "=== Creating ${BIN_SIZE}bp genome bins ==="

    # Check if reference genome index exists
    if [ ! -f "${REF_GENOME}.fai" ]; then
        echo "Creating reference genome index..."
        samtools faidx ${REF_GENOME}
    fi

    # Create genome file for bedtools (chr, length)
    # Filter to standard chromosomes only (1-22, X, Y)
    echo "Filtering to standard chromosomes (1-22, X, Y)..."
    grep -E '^[0-9]+\s|^X\s|^Y\s' ${REF_GENOME}.fai | \
        awk -v OFS='\t' '{print $1, $2}' > ${WD}/reference/genome_sizes.txt

    # Create bins using bedtools
    echo "Creating ${BIN_SIZE}bp windows..."
    bedtools makewindows \
        -g ${WD}/reference/genome_sizes.txt \
        -w ${BIN_SIZE} > ${BINS_FILE}

    # Count bins
    NUM_BINS=$(wc -l < ${BINS_FILE})
    echo "Created ${NUM_BINS} genomic bins"
fi

###############################################
# Summary
###############################################
echo ""
echo "=== Setup Complete ==="
echo "Working directory: ${WD}"
echo "Sample list: ${SAMPLE_LIST}"
echo "Genome bins: ${BINS_FILE}"
echo ""
echo "Next steps:"
echo "  1. Review sample_list.txt"
echo "  2. Update NUM_SAMPLES in config.sh if needed"
echo "  3. Run: sbatch 01_concat_fastq.sh"
echo "  4. Run: sbatch --dependency=afterok:<JOB_ID> 02_process_sample.sh"
