# GEO Data Submission Workflow

## Introduction

This document outlines the workflow for submitting high-throughput functional genomic data to the Gene Expression Omnibus (GEO). The process involves preparing raw and processed data, filling out a metadata spreadsheet, and uploading the data and metadata to GEO.

## Input

The input to this pipeline includes:

- Raw data files: These should be FASTQ files containing raw data. All files should be moved to the same folder, not in separate subfolders for each sample.
- Processed data files: For RNA-seq and scRNA-seq data, these typically include a count matrix and metadata matrix in CSV format. The metadata should ideally include experimental data information and cell type annotation. You can use [export_data.rmd](./export_data.rmd) to export Seurat scRNA-seq object to CSV matrices.
- Metadata spreadsheet: This is a spreadsheet that contains descriptive information about the overall study, individual samples, all protocols, and references to processed and raw data file names. For collaborations, you should ask collaborators to fill in as much experimental information as possible. File names and MD5 checksums should be provided by the analyst. The example template can be found at [seq_template.xlsx](./seq_template.xlsx)

## Output

The output of this pipeline is a complete GEO submission. Initially, it is set to private upon upload. When the associated paper is published, you should set the submission to public and update the citation information.

## Steps

1. Visit the [GEO sequence data submission page](https://www.ncbi.nlm.nih.gov/geo/info/seq.html) and log in.

2. Download the metadata spreadsheet.

3. Transfer all your raw and processed data files to the GEO FTP server. Refer to the [GEO FTP submission page](https://www.ncbi.nlm.nih.gov/geo/info/submissionftp.html) for detailed instructions. The easiest way is to prepare your data in the same folder structure as used at the Ohio Supercomputer Center (OSC). You can use 'lftp' to accomplish this. The FTP server credentials can also be found on the same page. The password is temporary and changes every few days.

```bash
# Replace the ** part with your FTP password, and replace '/path/to/wd' with the path to your local directory containing the data to be uploaded.
# Also replace 'cankun.wang@orcid_6LOT6a04' with your personal GEO FTP folder.

cd /path/to/wd

lftp ftp://geoftp:********@ftp-private.ncbi.nlm.nih.gov
cd uploads/cankun.wang@orcid_6LOT6a04
mirror -R .
```

4. Once the data is uploaded, return to the GEO website to upload the metadata spreadsheet.

5. If all goes well, a data accession number should be created within a week. If there are any issues, you will receive an email from the GEO staff explaining what needs to be corrected.

## Directory structure


Here is an example of the folder structure for a GEO submission, as used at the Ohio Supercomputer Center (OSC):

![GEO Directory Structure](https://www.ncbi.nlm.nih.gov/geo/img/directory_structure_diagram.jpg)

This structure helps to organize the raw and processed data files along with the metadata spreadsheet. Following this structure can facilitate a smooth and efficient data submission process to GEO.

For RNA-seq data, the directory structure should include a "raw_data" folder to store all FASTQ files and a "processed_data" folder to store the count matrix and metadata matrix. These files are typically named counts.csv and metadata.csv, respectively. It is recommended to compress these files to save storage space and facilitate faster data transfer. The compressed files should be named counts.csv.gz and metadata.csv.gz.

## Contact

Author: Cankun Wang, cankun.wang@osumc.edu


