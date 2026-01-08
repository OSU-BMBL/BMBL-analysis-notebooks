#!/bin/bash
#
# LP-WGS Karyotyping Pipeline - Main Processing Script
# fastp QC -> BWA-MEM alignment -> sorted BAM -> bedtools coverage -> .counts.bed
#
# Usage: sbatch --dependency=afterok:<CONCAT_JOB_ID> 02_process_sample.sh
#
# Prerequisites: Run 00_setup.sh and 01_concat_fastq.sh first
#

###############################################
# Load configuration
###############################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

#SBATCH --account ${SLURM_ACCOUNT}
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=${THREADS}
#SBATCH --mem=64GB
#SBATCH --job-name=lpwgs_process
#SBATCH --output=logs/process_%A_%a.out
#SBATCH --error=logs/process_%A_%a.err
#SBATCH --array=1-${NUM_SAMPLES}

set -e

# Load required modules
load_modules

###############################################
# Get sample name
###############################################
cd ${WD}
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_LIST})

if [ -z "$SAMPLE" ]; then
    echo "ERROR: Could not get sample name for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "=============================================="
echo "Processing sample: ${SAMPLE}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"
echo "=============================================="

###############################################
# Define file paths
###############################################
INPUT_R1="${WD}/merged_fastq/${SAMPLE}_R1.fq.gz"
INPUT_R2="${WD}/merged_fastq/${SAMPLE}_R2.fq.gz"

TRIM_R1="${WD}/fastp_out/${SAMPLE}_R1.trimmed.fq.gz"
TRIM_R2="${WD}/fastp_out/${SAMPLE}_R2.trimmed.fq.gz"
FASTP_HTML="${WD}/fastp_out/${SAMPLE}.fastp.html"
FASTP_JSON="${WD}/fastp_out/${SAMPLE}.fastp.json"

SORTED_BAM="${WD}/alignment/${SAMPLE}.sorted.bam"
BAM_INDEX="${WD}/alignment/${SAMPLE}.sorted.bam.bai"

COUNTS_FILE="${WD}/counts/${SAMPLE}.counts.bed"

###############################################
# Check if final output exists
###############################################
if [ -f "${COUNTS_FILE}" ]; then
    echo "Final output already exists: ${COUNTS_FILE}"
    echo "Skipping this sample."
    exit 0
fi

###############################################
# Validate inputs
###############################################
if [ ! -f "${INPUT_R1}" ] || [ ! -f "${INPUT_R2}" ]; then
    echo "ERROR: Input FASTQ files not found:"
    echo "  R1: ${INPUT_R1}"
    echo "  R2: ${INPUT_R2}"
    exit 1
fi

if [ ! -f "${BINS_FILE}" ]; then
    echo "ERROR: Genome bins file not found: ${BINS_FILE}"
    exit 1
fi

###############################################
# Step 1: fastp - Quality control and trimming
###############################################
echo ""
echo "=== Step 1: fastp QC and trimming ==="
echo "Time: $(date)"

if [ -f "${TRIM_R1}" ] && [ -f "${TRIM_R2}" ]; then
    echo "Trimmed files already exist. Skipping fastp."
else
    ${FASTP} \
        -w ${THREADS} \
        -i ${INPUT_R1} \
        -I ${INPUT_R2} \
        -o ${TRIM_R1} \
        -O ${TRIM_R2} \
        -h ${FASTP_HTML} \
        -j ${FASTP_JSON} \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36

    echo "fastp complete."
fi

###############################################
# Step 2: BWA-MEM alignment
###############################################
echo ""
echo "=== Step 2: BWA-MEM alignment ==="
echo "Time: $(date)"

if [ -f "${SORTED_BAM}" ]; then
    echo "Sorted BAM already exists. Skipping alignment."
else
    # BWA-MEM with read group, pipe to samtools for sorting
    bwa mem -t ${THREADS} \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}" \
        ${BWA_INDEX} \
        ${TRIM_R1} \
        ${TRIM_R2} \
        2> ${WD}/logs/${SAMPLE}.bwa.log \
        | samtools view -@ 4 -bS - \
        | samtools sort -@ 4 -m 4G -o ${SORTED_BAM} -

    echo "Alignment complete."
fi

###############################################
# Step 3: Index BAM
###############################################
echo ""
echo "=== Step 3: Index BAM ==="
echo "Time: $(date)"

if [ -f "${BAM_INDEX}" ]; then
    echo "BAM index already exists. Skipping."
else
    samtools index -@ ${THREADS} ${SORTED_BAM}
    echo "Indexing complete."
fi

###############################################
# Step 4: Calculate alignment stats
###############################################
echo ""
echo "=== Step 4: Alignment statistics ==="
echo "Time: $(date)"

TOTAL_READS=$(samtools view -c ${SORTED_BAM})
MAPPED_READS=$(samtools view -c -F 4 ${SORTED_BAM})
PROPER_PAIRS=$(samtools view -c -f 2 ${SORTED_BAM})

echo "Total reads: ${TOTAL_READS}"
echo "Mapped reads: ${MAPPED_READS}"
echo "Properly paired: ${PROPER_PAIRS}"

# Estimate coverage (assuming 150bp reads, 3Gb genome)
COVERAGE=$(echo "scale=2; ${MAPPED_READS} * 150 / 3000000000" | bc -l)
echo "Estimated coverage: ${COVERAGE}x"

###############################################
# Step 5: Count reads in bins
###############################################
echo ""
echo "=== Step 5: Count reads in genomic bins ==="
echo "Time: $(date)"

bedtools coverage \
    -a ${BINS_FILE} \
    -b ${SORTED_BAM} \
    -counts \
    > ${COUNTS_FILE}

NUM_BINS=$(wc -l < ${COUNTS_FILE})
echo "Read counting complete. Output: ${COUNTS_FILE}"
echo "Total bins: ${NUM_BINS}"

###############################################
# Summary
###############################################
echo ""
echo "=============================================="
echo "=== Completed: ${SAMPLE} ==="
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -lh ${FASTP_HTML} ${FASTP_JSON}
ls -lh ${SORTED_BAM} ${BAM_INDEX}
ls -lh ${COUNTS_FILE}
echo ""
echo "Alignment stats:"
echo "  Total reads: ${TOTAL_READS}"
echo "  Mapped reads: ${MAPPED_READS}"
echo "  Coverage: ${COVERAGE}x"
