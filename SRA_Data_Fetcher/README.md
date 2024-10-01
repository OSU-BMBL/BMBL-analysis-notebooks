# SRA data fetcher
This package is to download SRA data with SRA Toolkit in parallel.

### Usage

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
Note: please make sure the output directory and the temporary directory are exist.

Example:
```
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

### Usage in OSC
We also provide `osc_script.sh` script to download SRA data in OSC. Here is the script, please modify cluster configuration and the output directory and the temp directory as needed
```
#!/usr/bin/bash
#SBATCH --job-name fetch_sra_with_sratoolkit
#SBATCH --account PCON0100
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=128GB

date
module load python/3.9-2022.05
module load sratoolkit/2.11.2

cd $HOME
python ./fetch-with-sratoolkit.py -o ./ -t ./ ./sra_ids.json
date
```


