# BMBL Spatial Transcriptomics Workflow - Docker Image
# Optimized for spatial transcriptomics analysis

FROM rocker/rstudio:4.3.0

LABEL maintainer="BMBL Lab <shaopeng.gu@osumc.edu>"
LABEL description="Docker environment for ST_general_workflow"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libtiff-dev \
    libfftw3-dev \
    pandoc \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install Seurat (includes spatial support in v4+)
RUN R -e "install.packages(c('Seurat', 'SeuratObject'), repos='https://cloud.r-project.org/')"

# Install Spatial-specific packages
RUN R -e "BiocManager::install(c( \
    'SpatialExperiment', \
    'SingleCellExperiment', \
    'scater' \
), ask=FALSE)"

# Install SeuratData for example datasets
RUN R -e "remotes::install_github('satijalab/seurat-data')" || \
    echo "SeuratData installation optional"

# Install Giotto (alternative spatial framework)
RUN R -e "install.packages(c('BiocManager', 'remotes'), repos='https://cloud.r-project.org/')" && \
    R -e "remotes::install_github('drieslab/Giotto')" || \
    echo "Giotto installation optional, continuing..."

# Install visualization packages
RUN R -e "install.packages(c( \
    'patchwork', \
    'ggplot2', \
    'RColorBrewer', \
    'viridis', \
    'pheatmap' \
), repos='https://cloud.r-project.org/')"

# Install utility packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'rmarkdown', \
    'knitr', \
    'shiny' \
), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /home/rstudio/workspace

# Expose RStudio port
EXPOSE 8787

# Default command
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize", "0"]
