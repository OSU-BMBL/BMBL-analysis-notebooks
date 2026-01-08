#!/bin/bash
#
# LP-WGS Karyotyping Pipeline - Configuration File
#
# INSTRUCTIONS:
#   1. Copy this file and edit the values below for your environment
#   2. All other scripts source this file for configuration
#
# For OSC users: Update SLURM_ACCOUNT and paths under /fs/
# For other clusters: Update paths and module names accordingly
#

###############################################
# SLURM Configuration
###############################################
# Your SLURM account/allocation
SLURM_ACCOUNT="YOUR_ACCOUNT"

###############################################
# Directory Paths
###############################################
# Working directory (where outputs will be written)
# Should be on scratch/fast storage for large files
WD="/path/to/working/directory"

# Raw FASTQ data directory
# Should contain subdirectories named by sample, each with *_1.fq.gz and *_2.fq.gz files
DATA_DIR="/path/to/raw/fastq/data"

###############################################
# Reference Genome
###############################################
# GRCh38/hg38 reference genome FASTA
REF_GENOME="/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# BWA index prefix (without file extension)
BWA_INDEX="/path/to/bwa_index"

###############################################
# Tool Paths
###############################################
# fastp executable path
# Set to "fastp" if installed in PATH via conda/module
FASTP="/path/to/fastp"

###############################################
# Pipeline Parameters
###############################################
# Genome bin size for read counting (default: 50kb = 50000)
BIN_SIZE=50000

# Number of threads for parallel processing
THREADS=16

# Number of samples (update based on your sample_list.txt)
# Used for SLURM array job range
NUM_SAMPLES=14

###############################################
# Module Configuration (OSC-specific)
###############################################
# Uncomment and modify for your cluster's module system
# If using conda, you can comment out module loads and activate your environment instead
#
# load_modules() {
#     module load samtools/1.21
#     module load bedtools2/2.31.0
#     module load bwa/0.7.17
#     module load R/4.3.0
# }

# OSC Ascend defaults:
load_modules() {
    module load samtools/1.21
    module load bedtools2/2.31.0
    module load bwa/0.7.17
}

###############################################
# Derived Paths (usually don't need to change)
###############################################
BINS_FILE="${WD}/reference/hg38_${BIN_SIZE}bp_bins.bed"
SAMPLE_LIST="${WD}/sample_list.txt"

###############################################
# Validation
###############################################
validate_config() {
    local errors=0

    if [ "${SLURM_ACCOUNT}" == "YOUR_ACCOUNT" ]; then
        echo "ERROR: SLURM_ACCOUNT not configured"
        errors=$((errors + 1))
    fi

    if [ ! -d "${WD}" ]; then
        echo "WARNING: Working directory does not exist: ${WD}"
        echo "  Will be created by 00_setup.sh"
    fi

    if [ ! -d "${DATA_DIR}" ]; then
        echo "ERROR: Data directory not found: ${DATA_DIR}"
        errors=$((errors + 1))
    fi

    if [ ! -f "${REF_GENOME}" ]; then
        echo "ERROR: Reference genome not found: ${REF_GENOME}"
        errors=$((errors + 1))
    fi

    if [ ${errors} -gt 0 ]; then
        echo ""
        echo "Please fix the above errors in config.sh before running the pipeline."
        return 1
    fi

    return 0
}
