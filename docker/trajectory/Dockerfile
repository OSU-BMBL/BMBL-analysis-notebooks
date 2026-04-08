# BMBL Trajectory Analysis (Slingshot) - Docker Image
# Optimized for pseudotime and trajectory analysis

FROM rocker/rstudio:4.3.0

LABEL maintainer="BMBL Lab <shaopeng.gu@osumc.edu>"
LABEL description="Docker environment for scRNAseq_trajectory_Slingshot"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    pandoc \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install Seurat (base)
RUN R -e "install.packages(c('Seurat', 'SeuratObject'), repos='https://cloud.r-project.org/')"

# Install Trajectory-specific packages
RUN R -e "BiocManager::install(c( \
    'slingshot', \
    'TrajectoryUtils', \
    'tradeSeq', \
    'SingleCellExperiment', \
    'scater', \
    'scran' \
), ask=FALSE)"

# Install monocle3 (from GitHub)
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org/')" && \
    R -e "remotes::install_github('cole-trapnell-lab/monocle3')" || \
    echo "Monocle3 installation may have failed, continuing..."

# Install visualization packages
RUN R -e "install.packages(c( \
    'patchwork', \
    'ggplot2', \
    'RColorBrewer', \
    'viridis', \
    'ggrepel' \
), repos='https://cloud.r-project.org/')"

# Install utility packages
RUN R -e "install.packages(c( \
    'tidyverse', \
    'rmarkdown', \
    'knitr', \
    'devtools' \
), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /home/rstudio/workspace

# Expose RStudio port
EXPOSE 8787

# Default command
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize", "0"]
