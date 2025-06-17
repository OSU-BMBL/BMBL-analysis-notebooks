# SRA Data Fetcher Workflow

## Overview

- **What it does:** A parallel data fetching tool for downloading SRA (Sequence Read Archive) data using SRA Toolkit. The Sequence Read Archive (SRA) is a public repository of high-throughput sequencing data, primarily from next-generation sequencing platforms. It serves as a crucial resource for researchers to access and analyze genomic, transcriptomic, and metagenomic data. This tool facilitates efficient downloading of SRA datasets, which are essential for various bioinformatics analyses including gene expression studies, variant calling, metagenomic analysis, and comparative genomics. The parallel downloading capability significantly reduces the time required to retrieve multiple datasets, making it particularly useful for large-scale genomic studies and meta-analyses.
- **Who it's for:** Bioinformatics researchers and data scientists who need to efficiently download multiple SRA datasets.
- **Key Features:**
  - Parallel downloading of multiple SRA datasets
  - Configurable output and temporary directories
  - OSC (Ohio Supercomputer Center) integration support
  - JSON-based input configuration

---

## Getting Started

### Prerequisites

- **Software:** 
  - Python 3.12
  - SRA Toolkit 3.0.2
- **Knowledge:** Basic understanding of SRA data formats and command-line operations

### Instruction

1. Ensure you have Python 3.12 and SRA Toolkit 3.0.2 installed
2. Prepare a JSON file with SRA IDs to download
3. Create output and temporary directories if using custom paths

---

## Usage

### 1. Input Data

- **Required Format:** JSON file containing SRA IDs
- **Directory Structure:**
  ```
  SRA_Data_Fetcher/
  ├── sra_ids.json
  ├── fetch-with-sratoolkit.py
  └── osc_script.sh
  ```

### 2. Running the Workflow

Basic command:
```bash
python fetch-with-sratoolkit.py [-h] [-o OUTPUTDIR] [-t TEMPDIR] jsonfile
```

```
fetch-with-sratoolkit.py [-h] [-o OUTPUTDIR] [-t TEMPDIR] jsonfile

positional arguments:
  jsonfile              Json file specifying sra file ids to download

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUTDIR, --outputDir OUTPUTDIR
                        Output directory (optional), the default output directory is the current directory
  -t TEMPDIR, --tempDir TEMPDIR
                        Temporary directory (optional), the default temporary directory is the current
                        directory
```

#### Example:
```bash
python fetch-with-sratoolkit.py -o /output -t /tmp ./sra_ids.json
```
The `jsonfile` is to specify SRA file ids, its format is like:
```
{
    "sraIds": [
        "SRR2932830",
        "SRR2932831",
        "SRR2932832",
        ...
    ]
}
```
### 3. Testing with Sample Data

- **Command:** `python fetch-with-sratoolkit.py -o ./ -t ./ ./sra_ids.json`
- **Expected Output:** Downloaded SRA files in the specified output directory

### 4. Pipeline Output

The workflow downloads the specified SRA files to the output directory. Each SRA file will be downloaded and processed according to the SRA Toolkit specifications.

---

## OSC Integration

For running on the Ohio Supercomputer Center (OSC), use the provided `osc_script.sh`:

```bash
#!/usr/bin/bash
#SBATCH --job-name fetch_sra_with_sratoolkit
#SBATCH --account XXXXXXXX
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=128GB

date
module load python/3.12
module load sratoolkit/3.0.2

python ./fetch-with-sratoolkit.py -o ./ -t ./ ./sra_ids.json
date
```

---

**Author(s):** Shaohong Feng

**Tester(s):** Mirage Modi

**Contact Email:** shaohong.feng@osumc.edu
