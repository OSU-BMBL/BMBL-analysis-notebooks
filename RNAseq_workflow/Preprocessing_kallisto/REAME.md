This RNA-seq data preprocessing tutorial is based on [kallisto](https://pachterlab.github.io/kallisto/manual.html) and [txiport](https://bioconductor.org/packages/release/bioc/html/tximport.html).
* **kallisto**
    * Introduction: kallisto uses a pseudoalignment algorithm to determine the abundance of transcripts in a sample rapidly, without the need for alignment. The outputs kallisto consist of transcripts per million (TPM) and estimated counts (float number).
    * Installation: The most convenient way is to install in a conda envrioment.
    ``` bash
    conda install -c kallisto bioconda
    ```
* **tximport**
    * Introduction: tximport imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages.
    Installation:
    ``` R
    BiocManager::install("tximport")
    ```

**Workflow**
* Run `run_quantification.sh`, change some paths accordingly, you may need to download transcript files and gene annotation files from [GENCODE](https://www.gencodegenes.org/).
* Run `txiport.R`, change some paths accordingly, output the txiport object, which generate abundance(TPM), counts, length data slots.

For more usuage, please refer to the following link:
* [kallisto](https://pachterlab.github.io/kallisto/manual.html)
* [txiport](https://bioconductor.org/packages/release/bioc/html/tximport.html)

