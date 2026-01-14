# GEO Download Workflow

Download sample-level count matrices from NCBI GEO using HTTPS.

## Overview

GEO datasets often provide data in two ways:
1. **Series-level files**: Combined archives (RAW.tar), metadata sheets
2. **Sample-level files**: Individual files per GSM (count matrices, barcodes)

This workflow handles **sample-level downloads** - useful when:
- Combined archive doesn't include gene annotations
- You need specific samples, not the entire dataset
- Per-sample processing is required

## When to Use

**Use this workflow when:**
- Downloading scRNA-seq count matrices from GEO
- Dataset has per-sample supplementary files (TSV, MTX, CSV)
- You need individual sample files for processing
- FTP downloads are failing (GEO FTP is often unreliable)

**Don't use when:**
- You only need the combined RAW.tar archive
- Data is available through GEO2R
- You're downloading SRA/FASTQ files (use SRA toolkit instead)

## Quick Start

```bash
# 1. Install dependencies
python 0_install_packages.py

# 2. Create sample mapping file (samples.json)
# 3. Run download
python 1_download_samples.py --samples samples.json --output ./data/
```

## Step-by-Step Guide

### Step 1: Identify Your Samples

Find the GSM IDs and their associated filenames from:
- The GEO series page (GSExxxxx)
- The `geo_sample_metadata.csv` in your data folder
- NCBI's FTP listing

### Step 2: Create Sample Mapping

Create a JSON file mapping GSM IDs to filenames:

```json
{
  "GSM4716780": "GSM4716780_sample1_counts.tsv.gz",
  "GSM4716781": "GSM4716781_sample2_counts.tsv.gz",
  "GSM4716782": "GSM4716782_sample3_counts.tsv.gz"
}
```

Save as `samples.json`.

### Step 3: Run the Download

```bash
python 1_download_samples.py --samples samples.json --output ./raw_data/
```

The script will:
- Skip files that already exist
- Add delays between downloads (avoids rate limiting)
- Report progress and any failures

### Step 4: Verify Downloads

```bash
ls -lh ./raw_data/*.gz
```

## Usage Examples

### Basic Usage

```bash
python 1_download_samples.py --samples samples.json --output ./AD015/
```

**Expected output:**
```
[1/10] Downloading GSM4716780...
  Saved: GSM4716780_counts.tsv.gz
[2/10] Downloading GSM4716781...
  Saved: GSM4716781_counts.tsv.gz
[3/10] Skipping GSM4716782 (exists)
...
Done: 9 downloaded, 1 skipped
```

### With Custom Delay

For large datasets, increase delay to avoid rate limiting:

```bash
python 1_download_samples.py --samples samples.json --output ./data/ --delay 1.0
```

## How It Works

The script builds GEO HTTPS download URLs:

```
https://www.ncbi.nlm.nih.gov/geo/download/?acc={GSM}&format=file&file={filename}
```

This is more reliable than FTP because:
- HTTPS works through firewalls
- Better error handling
- No anonymous FTP issues

## Troubleshooting

### "Connection refused" or timeout errors

**Problem**: NCBI servers are busy or rate limiting you.

**Solution**:
- Increase `--delay` to 2.0 or higher
- Try again later (off-peak hours)
- Download in smaller batches

### File downloaded but empty or corrupted

**Problem**: Download interrupted or filename mismatch.

**Solution**:
- Delete the file and re-run (script skips existing files)
- Verify the filename exactly matches GEO listing

### "404 Not Found" errors

**Problem**: File doesn't exist at expected URL.

**Solution**:
- Verify GSM ID is correct
- Check if filename includes version suffix (e.g., `_v2`)
- Some datasets use different naming conventions

## Tips

1. **Start small**: Download one sample first to verify the filename format
2. **Check file sizes**: Compare with sizes listed on GEO
3. **Organize by dataset**: Use folders like `./AD015/`, `./AD016/`
4. **Resume safely**: Script skips existing files
5. **Next step**: Use `scRNAseq_H5AD_conversion_workflow` to convert to H5AD format

## Files

| File | Description |
|------|-------------|
| `0_install_packages.py` | Install dependencies |
| `1_download_samples.py` | Main download script |

## Contact

BMBL Lab - https://u.osu.edu/bmbl/
