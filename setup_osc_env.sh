#!/bin/bash
# BMBL Workflow Environment Setup for OSC (Ohio Supercomputer Center)
# Usage: source setup_osc_env.sh
#
# This script loads all required modules for BMBL bioinformatics workflows
# Works on both Pitzer and Ascend clusters

echo "=========================================="
echo "BMBL Workflow Environment Setup"
echo "Ohio Supercomputer Center"
echo "=========================================="
echo ""

# Determine cluster name
if hostname | grep -q "pitzer"; then
    CLUSTER="Pitzer"
elif hostname | grep -q "ascend"; then
    CLUSTER="Ascend"
else
    CLUSTER="Unknown"
    echo "Warning: Could not determine cluster. Attempting general setup..."
fi

echo "Detected cluster: $CLUSTER"
echo ""

# Purge default modules to avoid conflicts
echo "Loading modules..."
module purge

# Core languages
module load gcc-compatibility/8.5.0
module load R/4.3.0
module load python/3.12
module load perl/5.34.0

# Bioinformatics tools - Quality Control
module load fastqc/0.12.1
module load multiqc/1.14

# Bioinformatics tools - Alignment
module load samtools/1.21
module load bedtools2/2.31.0
module load bwa/0.7.17
module load star/2.7.10b
module load bowtie2/2.5.1

# Bioinformatics tools - Peak calling and analysis
module load macs2/2.2.9.1
module load homer/4.11
module load meme/5.4.1

# Bioinformatics tools - Utilities
module load sratoolkit/3.0.2
module load trim_galore/0.6.10
module load subread/2.0.6

# Pipeline runners
module load nextflow/24.10.4

# Visualization
module load graphviz/7.2.0

echo ""
echo "=========================================="
echo "Module Status"
echo "=========================================="
echo ""

# Display versions of key tools
echo "R version:"
R --version | head -1

echo ""
echo "Python version:"
python --version

echo ""
echo "Key bioinformatics tools:"
echo "  samtools: $(samtools --version 2>&1 | head -1)"
echo "  bedtools: $(bedtools --version 2>&1 | head -1)"
echo "  fastqc: $(fastqc --version 2>&1 || echo 'version check not available')"
echo "  multiqc: $(multiqc --version 2>&1 || echo 'version check not available')"
echo "  nextflow: $(nextflow --version 2>&1 | head -1)"

echo ""
echo "=========================================="
echo "Environment Variables"
echo "=========================================="

# Set up environment variables
export BMBL_ROOT=$(pwd)
export R_LIBS_USER="$HOME/R/bmbl-workflows"
export PYTHONPATH="$BMBL_ROOT:$PYTHONPATH"
export OMP_NUM_THREADS=8

# Create R library directory if needed
mkdir -p $R_LIBS_USER

echo "  BMBL_ROOT: $BMBL_ROOT"
echo "  R_LIBS_USER: $R_LIBS_USER"
echo "  PYTHONPATH: $PYTHONPATH"
echo "  OMP_NUM_THREADS: $OMP_NUM_THREADS"

echo ""
echo "=========================================="
echo "Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Install R packages: Rscript install_r_packages.R"
echo "  2. Browse workflows: ls *_workflow/"
echo "  3. Check workflow README: cat <workflow>/README.md"
echo ""
echo "For help, see AGENTS.md or contact shaopeng.gu@osumc.edu"
