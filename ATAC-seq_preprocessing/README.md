# ATAC-seq preprocessing using nfcore

The folder contains the ATAC-seq preprocessing workflow from nfcore/atacseq

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The Nextflow DSL2 implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. 


To run on your data, prepare a tab-separated samplesheet with your input data. Please follow the documentation on samplesheets for more details. An example samplesheet for running the pipeline looks as follows:

```csv
sample,fastq_1,fastq_2,replicate
CONTROL,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,1
CONTROL,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,2
CONTROL,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,3

```

Next, you can run the pipeline using:
```bash
nextflow run nf-core/atacseq --input samplesheet.csv --outdir <OUTDIR> --genome GRCh3/ --read_length <50|100|150|200> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
```

Note that the read length can be checked from head of fastq sequnces, and docker should be used in vp03 environment.