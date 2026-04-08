# BMBL RNAseq nfcore Workflow - Docker Image
# Optimized for bulk RNA-seq analysis

FROM rocker/rstudio:4.3.0

LABEL maintainer="BMBL Lab <shaopeng.gu@osumc.edu>"
LABEL description="Docker environment for RNAseq_nfcore_workflow"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    pandoc \
    git \
    default-jdk \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install core RNA-seq packages
RUN R -e "BiocManager::install(c( \
    'DESeq2', \
    'edgeR', \
    'limma', \
    'tximport', \
    'tximeta', \
    'SummarizedExperiment', \
    'GenomicFeatures', \
    'rtracklayer' \
), ask=FALSE)"

# Install pathway analysis packages
RUN R -e "BiocManager::install(c( \
    'clusterProfiler', \
    'DOSE', \
    'enrichplot', \
    'pathview', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'ReactomePA', \
    'fgsea' \
), ask=FALSE)"

# Install visualization packages
RUN R -e "install.packages(c( \
    'ggplot2', \
    'pheatmap', \
    'RColorBrewer', \
    'EnhancedVolcano', \
    'PCAtools', \
    'ggrepel' \
), repos='https://cloud.r-project.org/')"

# Install EnhancedVolcano from Bioconductor
RUN R -e "BiocManager::install('EnhancedVolcano', ask=FALSE)"

# Install utility packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'rmarkdown', \
    'knitr', \
    'DT', \
    'here' \
), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /home/rstudio/workspace

# Expose RStudio port
EXPOSE 8787

# Default command
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize", "0"]
