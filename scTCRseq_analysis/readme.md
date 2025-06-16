
# Metastasis_TCR

### 📌 Purpose

This pipeline integrates single-cell RNA sequencing (scRNA-seq) and T-cell receptor (TCR) sequencing data from metastatic and primary tumor samples. It is designed to analyze clonal expansion and TCR diversity across patients, enabling the investigation of immune dynamics in cancer metastasis.

### 📂 Input Files

- **TCR Annotation CSVs**  
  Output from 10x Genomics VDJ pipeline, e.g.,  
  `/.../SRR13737401_vdj/outs/filtered_contig_annotations.csv`

- **Gene Expression Data**  
  From processed scRNA-seq count matrices:
  - `RNA_matrix.mtx`: Matrix Market file with counts
  - `barcodes.csv`: Cell barcode file
  - `features.csv`: Gene feature names

These files are organized by patient and condition (e.g., primary tumor vs lymph node metastasis).

### 🧮 Output

- `MuData` object (`mdata`) containing:
  - `"gex"`: single-cell gene expression profiles
  - `"airr"`: associated TCR annotations
  - Cell-level labels for visualization and downstream analysis

### 🔧 Tools Used

- `scanpy`, `anndata`, `muon`, `scirpy`: For handling scRNA-seq and immune repertoire data
- `scipy`, `pandas`, `matplotlib`: For preprocessing and visualization

Author: Xiaoying Wang
Tested by: Shaopeng Gu
