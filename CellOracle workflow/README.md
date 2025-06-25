# Building GRN and in-silico gene perturbation using celloracle


## Pipeline input
#peak file from scATACseq and Cicero file to build your custom TF-gene linkage, however you can load prebuilt one using co.data

The Main input for cellorcle in scRNAseq in anndata object 

## Pipeline output
The output is dictionary where key is the celltype and value is the GRN (celltype specific GRN)


## Contact

Author: Ahmed Ghobashi
Tested by: Shaopeng Gu
## Methods for manuscript

Simulated gene perturbation was performed using CellOracle . We constructed our own gene-regulatory networks (GRN) from our scATAC-seq data using Cicero package. First, our scATC-seq data were converted to cicero object and cis-regulatory interactions were identified using run_cicero. Transcription start sites (TSS) were annotated using hg38 genome as reference using CellOracle. TSS information was integrated with cis-regulatory interactions using integrate_tss_peak_with_cicer. Peaks with weak coaccessibility scores using integrated.coaccess >= 0.8 command. Transcription factor motif scan was performed using tfi.scan. the final GRN and 3000 highly variable gene expression from ISX9 gene matrix were loaded into the CellOracle object using co.import_TF_data and co.import_anndata_as_raw_count, respectively. We constructed a
cluster-specific GRN for all clusters using oc.get_links and kept only network edges with p-value <=0.01. To simulate gene overexpression or knockout, we perturbed the gene expression to 1.5 or 0, respectively in the oc.simulate_shift function.


