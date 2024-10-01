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


