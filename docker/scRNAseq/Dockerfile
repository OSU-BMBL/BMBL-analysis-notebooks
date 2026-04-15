# BMBL scRNAseq General Workflow - Docker Image
# Optimized for single-cell RNA-seq analysis with Seurat

FROM rocker/rstudio:4.3.0

LABEL maintainer="BMBL Lab <shaopeng.gu@osumc.edu>"
LABEL description="Docker environment for scRNAseq_general_workflow"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libmagick++-dev \
    pandoc \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install Seurat and single-cell packages
RUN R -e "install.packages(c('Seurat', 'SeuratObject', 'sctransform'), repos='https://cloud.r-project.org/')"

# Install Bioconductor single-cell packages
RUN R -e "BiocManager::install(c( \
    'SingleCellExperiment', \
    'scater', \
    'scran', \
    'scuttle', \
    'SingleR', \
    'celldex', \
    'limma', \
    'edgeR', \
    'DESeq2', \
    'batchelor' \
), ask=FALSE)"

# Install annotation databases
RUN R -e "BiocManager::install(c( \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'AnnotationDbi' \
), ask=FALSE)"

# Install visualization packages
RUN R -e "install.packages(c( \
    'patchwork', \
    'ggplot2', \
    'RColorBrewer', \
    'viridis', \
    'pheatmap', \
    'ggrepel', \
    'plotly', \
    'clustree' \
), repos='https://cloud.r-project.org/')"

# Install utility packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'rmarkdown', \
    'knitr', \
    'devtools', \
    'remotes', \
    'here', \
    'DT', \
    'shiny' \
), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /home/rstudio/workspace

# Expose RStudio port
EXPOSE 8787

# Default command
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize", "0"]
