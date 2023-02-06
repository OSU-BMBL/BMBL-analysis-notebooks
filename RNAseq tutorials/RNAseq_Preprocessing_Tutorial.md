# RNAseq Preprocessing Tutorial

## Prepare Data and Tools

In order to use the following tutorial, users must have an OSC account and use PCON0022.
This tutorial will use example data found in the Example_Data folder under RNA_seq_tutorial_fastq, or it can be found here in OSC:
```
/fs/ess/PCON0022/bmbl_notebooks/RNAseq_tutorial_fastq
```


### 1. Tools

Assuming users are using PCON0022 in OSC, all tools needed for RNAseq preprocessing can be found here:

```
/fs/ess/PCON0022/tools
```

### 2. Reference Genome Data

- Identify the species that was used to generate the RNAseq data

- If using mouse data:
  
  ```
  /fs/ess/PCON0022/tools/genome/Mus_musculus.GRCm38.99
  ```

- If using human data:
  
  ```
  /fs/ess/PCON0022/tools/genome/Homo_sapiens.GRCh38.99
  ```

### 3. Prepare Data

1. **Example file directory structure.** Each folder is an RNA-seq sample and contains all fastq files. Put your fastq files into folders in the following format:
   
   ![](./fastqs.png)
   
   ```
   /working_directory/folder1/fastq1.fastq.gz 
   ```

2. **Prepare a fastq file list.** The file will be named fastq_list.txt and will contain three columns: the first two specify two pair-end files of a sample, and the third column is the sample name, separated by a tab.
   
  ![](./list.png#center)
   
  To generate a **fastq_list.txt** automat, run the following code:
   
   ```
   chmod +x *
   module load gnu mkl R/4.0.2
   Rscript build_fastq_list.R
   ```

---

## Alignment

The following code is used to align the data. This is in a file named **run_primary_alignment.sh**

Note: users must set correct working directory and select the correct relevant reference data.


After modifying the above code, submit the reads alignment jobs in the Pitzer Shell Acess Cluster is OSC. Be sure to set your working directory in the cluster before running the following code:

```
./submit_primary_alignment.sh
```
Running the primary alignment will generate the following directories:
- **alignment_out** : contains .bam, .sam, and .sorted.bam files
- **fastp_out** : contains R1 and R2 fastq.gz files
- **result** : contains a pre-alignment directory with .fastp.json and html files


## Quantification

After the alignment is finished, we will submit quantification to generate an RNA-seq count matrix.

To know if the alignment is finished, run this code:

```
squeue -u USERNAME
```

If there are no current jobs, continue with quantification.

The code for the file **run_quantification.sh** is included in this Github folder. Users will need to modify the working directory and reference genome.


In the Pitzer Shell Access Cluster in OSC, set your correct working directory and run the following code:

```
sbatch run_quantification.sh
```
Running the quanitification will generate a slurm.out file in the working directory, as well as out.txt, out.txt.summary, and out2.txt files.

To generate counts.csv and meta.csv files:

1. Download out2.txt from the newly generate result file in your working directory to your PC

2. Copy all content in the file and paste it into excel

3. Save the file as counts.csv, ensuring it is in csv format

4. Create a meta.csv file. Each sample name should correspond to column names in counts.csv

An example of what these outputs should look like can be found in the Example_Data folder.

## MultiQC (optional)

MultiQC can be used to create an automated quality check report

### 1. Install multiQC

```
ml python
conda create -n multiqx python=3.8
source activate mutliqc
pip install multiqc
```

### 2. Generate quality check report

Enter the results folder containing **out.txt.summary** and **out.txt**

![](./outs.png#center)

Then run the following code to generate a quality check report.

```
ml python
source activate multiqc
multiqc *
```

### 3. Access the report
In the same result folder, download the multiqc_report.html file to your PC. The result will resmeble the image below.
![](./multiqc.png#center)