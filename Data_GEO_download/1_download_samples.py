#!/usr/bin/env python3
"""
Download sample-level files from GEO using HTTPS.

Usage:
    python download.py --gse GSE123456 --output /path/to/output
    python download.py --samples samples.json --output /path/to/output
"""

import argparse
import json
import os
import time
import urllib.parse
from pathlib import Path

import requests


def get_download_url(gsm: str, filename: str) -> str:
    """Build GEO HTTPS download URL."""
    encoded = urllib.parse.quote(filename)
    return f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gsm}&format=file&file={encoded}"


def download_file(url: str, filepath: Path, timeout: int = 120) -> bool:
    """Download a file with streaming."""
    try:
        r = requests.get(url, stream=True, timeout=timeout)
        r.raise_for_status()
        with open(filepath, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except Exception as e:
        print(f"  Error: {e}")
        return False


def download_samples(
    samples: dict[str, str], output_dir: Path, delay: float = 0.5
) -> tuple[int, int]:
    """
    Download all samples.

    Args:
        samples: Dict mapping GSM IDs to filenames
        output_dir: Output directory
        delay: Delay between downloads (seconds)

    Returns:
        Tuple of (success_count, skip_count)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    success = 0
    skipped = 0

    for i, (gsm, filename) in enumerate(samples.items(), 1):
        filepath = output_dir / filename

        if filepath.exists():
            print(f"[{i}/{len(samples)}] Skipping {gsm} (exists)")
            skipped += 1
            continue

        print(f"[{i}/{len(samples)}] Downloading {gsm}...")
        url = get_download_url(gsm, filename)

        if download_file(url, filepath):
            success += 1
            print(f"  Saved: {filepath.name}")
        else:
            print(f"  Failed: {gsm}")

        time.sleep(delay)

    return success, skipped


def main():
    parser = argparse.ArgumentParser(description="Download GEO sample files")
    parser.add_argument("--gse", help="GEO Series ID (e.g., GSE123456)")
    parser.add_argument(
        "--samples", help="JSON file with GSM->filename mapping"
    )
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument(
        "--delay", type=float, default=0.5, help="Delay between downloads"
    )
    args = parser.parse_args()

    if args.samples:
        with open(args.samples) as f:
            samples = json.load(f)
    elif args.gse:
        print(f"Fetching sample list for {args.gse}...")
        print("Note: You may need to manually create a samples.json file")
        print("with GSM IDs mapped to filenames.")
        return
    else:
        parser.error("Either --gse or --samples is required")

    output_dir = Path(args.output)
    success, skipped = download_samples(samples, output_dir, args.delay)

    print(f"\nDone: {success} downloaded, {skipped} skipped")


if __name__ == "__main__":
    main()
