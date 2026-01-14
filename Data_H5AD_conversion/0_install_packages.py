#!/usr/bin/env python3
"""
Install required packages for scRNA-seq H5AD conversion workflow.

Run this script first:
    python 0_install_packages.py
"""

import subprocess
import sys

packages = [
    "anndata",
    "scipy",
    "pandas",
    "numpy",
    "mygene",
]

if __name__ == "__main__":
    print("Installing required packages...")
    for pkg in packages:
        print(f"  Installing {pkg}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
    print("Done!")
