# BMBL scATACseq General Workflow - Docker Image
# Optimized for single-cell ATAC-seq analysis

FROM rocker/rstudio:4.3.0

LABEL maintainer="BMBL Lab <shaopeng.gu@osumc.edu>"
LABEL description="Docker environment for scATACseq_general_workflow"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    pandoc \
    git \
    # Bioinformatics tools
    samtools \
    bedtools \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install Seurat and Signac (for ATAC)
RUN R -e "install.packages(c('Seurat', 'SeuratObject', 'Signac'), repos='https://cloud.r-project.org/')"

# Install ATAC-specific Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'GenomicRanges', \
    'IRanges', \
    'rtracklayer', \
    'Biostrings', \
    'Rsamtools', \
    'TFBSTools', \
    'JASPAR2020', \
    'chromVAR', \
    'motifmatchr', \
    'ArchR', \
    'SingleCellExperiment' \
), ask=FALSE)"

# Install annotation databases
RUN R -e "BiocManager::install(c( \
    'EnsDb.Hsapiens.v86', \
    'BSgenome.Hsapiens.UCSC.hg38', \
    'BSgenome.Hsapiens.UCSC.hg19' \
), ask=FALSE)"

# Install visualization packages
RUN R -e "install.packages(c( \
    'patchwork', \
    'ggplot2', \
    'RColorBrewer', \
    'pheatmap' \
), repos='https://cloud.r-project.org/')"

# Install utility packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'rmarkdown', \
    'knitr' \
), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /home/rstudio/workspace

# Expose RStudio port
EXPOSE 8787

# Default command
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize", "0"]
