# BMBL Analysis Notebooks - Base Docker Image
# Provides a reproducible environment for bioinformatics workflows

FROM rocker/r-ver:4.3.0

LABEL maintainer="BMBL Lab <shaopeng.gu@osumc.edu>"
LABEL description="Reproducible environment for BMBL bioinformatics workflows"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=/usr/local/lib/R/site-library
ENV WORKSPACE=/home/rstudio/workspace

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Build tools
    build-essential \
    gfortran \
    cmake \
    # Version control
    git \
    wget \
    curl \
    # Python (for mixed workflows)
    python3 \
    python3-pip \
    python3-dev \
    # Bioinformatics tools
    samtools \
    bedtools \
    bowtie2 \
    # R Markdown dependencies
    pandoc \
    pandoc-citeproc \
    # XML/HTML parsing
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    # Image processing
    libmagick++-dev \
    # Java (for some bioinformatics tools)
    default-jdk \
    # File compression
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    # Clean up
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for scanpy-based workflows
RUN pip3 install --no-cache-dir \
    scanpy \
    anndata \
    pandas \
    numpy \
    scipy \
    matplotlib \
    seaborn \
    leidenalg \
    scrublet

# Install commonly used R packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install CRAN packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'rmarkdown', \
    'knitr', \
    'devtools', \
    'remotes', \
    'patchwork', \
    'RColorBrewer', \
    'pheatmap', \
    'ggrepel', \
    'plotly', \
    'shiny', \
    'DT', \
    'here', \
    'fs', \
    'glue', \
    'readr', \
    'dplyr', \
    'ggplot2' \
), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'SingleCellExperiment', \
    'SummarizedExperiment', \
    'MultiAssayExperiment', \
    'DESeq2', \
    'edgeR', \
    'limma', \
    'GenomicRanges', \
    'IRanges', \
    'Biobase', \
    'S4Vectors', \
    'DelayedArray', \
    'HDF5Array', \
    'rhdf5', \
    'org.Hs.eg.db', \
    'org.Mm.eg.db', \
    'AnnotationDbi', \
    'GO.db', \
    'KEGG.db', \
    'reactome.db', \
    'clusterProfiler', \
    'DOSE', \
    'enrichplot', \
    'fgsea', \
    'singscore', \
    'GSEABase', \
    'GSVA', \
    'BioNet' \
), ask=FALSE)"

# Install single-cell packages (CRAN and Bioconductor)
RUN R -e "install.packages(c( \
    'Seurat', \
    'SeuratObject', \
    'Signac' \
), repos='https://cloud.r-project.org/')"

RUN R -e "BiocManager::install(c( \
    'scater', \
    'scran', \
    'scuttle', \
    'SingleR', \
    'celldex', \
    'monocle', \
    'slingshot', \
    'tradeSeq', \
    'TrajectoryUtils', \
    'uwot', \
    'Rtsne', \
    'batchelor' \
), ask=FALSE)"

# Install visualization packages
RUN R -e "install.packages(c( \
    'viridis', \
    'circlize', \
    'ComplexHeatmap', \
    'ggpubr', \
    'ggsci', \
    'dittoSeq' \
), repos='https://cloud.r-project.org/')"

# Create workspace directory
RUN mkdir -p ${WORKSPACE}
WORKDIR ${WORKSPACE}

# Copy repository contents
COPY . ${WORKSPACE}

# Set permissions
RUN chown -R rstudio:rstudio ${WORKSPACE}

# Expose port for RStudio Server
EXPOSE 8787

# Default command
CMD ["R"]
