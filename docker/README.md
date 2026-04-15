# Docker Containers for BMBL Analysis Notebooks

This directory contains Docker configurations for creating reproducible, containerized environments for BMBL workflows.

## What is Docker?

Docker packages your entire computing environment (R version, packages, system tools) into a portable container. This ensures:
- **Reproducibility**: Same results on any computer
- **Portability**: Works on Mac, Windows, Linux
- **Isolation**: Doesn't interfere with your system

## Quick Start

### Prerequisites

Install Docker Desktop:
- **Mac/Windows**: [Download Docker Desktop](https://www.docker.com/products/docker-desktop)
- **Linux**: `sudo apt-get install docker.io docker-compose`

### Using Docker Compose (Recommended)

Start a specific workflow environment:

```bash
# scRNA-seq workflow
docker-compose up scrnaseq

# Trajectory analysis
docker-compose up trajectory

# scATAC-seq
docker-compose up scatacseq

# RNA-seq
docker-compose up rnaseq

# Spatial transcriptomics
docker-compose up spatial
```

Then open your browser to: **http://localhost:8788** (or port 8789-8792 for others)

**Login:**
- Username: `rstudio`
- Password: `bmbl2024`

### Building Images Manually

If you prefer not to use docker-compose:

```bash
# Build base image
docker build -t bmbl-base .

# Build specific workflow
docker build -t bmbl-scrnaseq -f docker/scRNAseq/Dockerfile .

# Run
docker run -p 8787:8787 -e PASSWORD=bmbl2024 bmbl-scrnaseq
```

## Available Containers

| Service | Port | Description |
|---------|------|-------------|
| `rstudio` | 8787 | Base environment with all packages |
| `scrnaseq` | 8788 | Single-cell RNA-seq (Seurat) |
| `trajectory` | 8789 | Pseudotime analysis (Slingshot) |
| `scatacseq` | 8790 | Single-cell ATAC-seq (Signac) |
| `rnaseq` | 8791 | Bulk RNA-seq (DESeq2/edgeR) |
| `spatial` | 8792 | Spatial transcriptomics |

## Data Management

### Mounting Your Data

Data is mounted from your local machine into the container:

```yaml
volumes:
  - ./data:/home/rstudio/data        # Your input data
  - ./outputs:/home/rstudio/outputs  # Results go here
```

Place your data in the `data/` folder before starting the container.

### Persisting R Packages

R packages are stored in a Docker volume so they persist between restarts:

```yaml
volumes:
  - r-packages:/usr/local/lib/R/site-library
```

## Troubleshooting

### Port Already in Use

If you get "port already allocated":
```bash
# Find what's using the port
lsof -i :8787

# Kill it or use a different port in docker-compose.yml
```

### Container Won't Start

Check logs:
```bash
docker-compose logs scrnaseq
```

### Rebuilding After Changes

```bash
# Force rebuild
docker-compose up --build scrnaseq
```

### Disk Space Issues

Docker images are large (~5-10 GB each). Clean up:
```bash
# Remove unused images
docker image prune

# Remove all stopped containers
docker container prune
```

## For Manuscripts

To create a fully reproducible environment for publication:

1. Build the container: `docker-compose up --build scrnaseq`
2. Run your analysis in RStudio
3. Document the image used: `docker images bmbl-scrnaseq`
4. Include in methods: "Analysis performed in Docker container bmbl-scrnaseq:v1.0"

## Customization

To add packages to a container:

1. Edit the relevant `Dockerfile`
2. Add `RUN R -e "install.packages('package_name')"`
3. Rebuild: `docker-compose up --build <service>`

## Questions?

Contact: Shaopeng Gu (shaopeng.gu@osumc.edu)
