#!/usr/bin/env python3
"""
Convert scRNA-seq count matrices to H5AD format.

Usage:
    python convert.py --input sample.tsv.gz --output sample.h5ad --species mouse
    python convert.py --input-dir ./samples --output-dir ./h5ad --species mouse
"""

import argparse
import gzip
from pathlib import Path

import anndata as ad
import mygene
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix


def convert_ensembl_to_symbols(
    ensembl_ids: list[str], species: str = "mouse"
) -> tuple[list[str], dict[str, str]]:
    """
    Convert Ensembl IDs to gene symbols using mygene.

    Returns:
        Tuple of (gene_symbols, ensembl_to_symbol_mapping)
    """
    mg = mygene.MyGeneInfo()
    results = mg.querymany(
        ensembl_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species=species,
        verbose=False,
    )

    # Build mapping
    ensembl_to_symbol = {}
    for r in results:
        query = r.get("query", "")
        symbol = r.get("symbol", query)  # fallback to Ensembl ID
        ensembl_to_symbol[query] = symbol

    # Convert with fallback
    gene_symbols = [ensembl_to_symbol.get(eid, eid) for eid in ensembl_ids]

    return gene_symbols, ensembl_to_symbol


def make_unique_symbols(symbols: list[str]) -> list[str]:
    """Make gene symbols unique by appending suffixes to duplicates."""
    seen = {}
    unique = []
    for sym in symbols:
        if sym in seen:
            seen[sym] += 1
            unique.append(f"{sym}_{seen[sym]}")
        else:
            seen[sym] = 0
            unique.append(sym)
    return unique


def load_count_matrix(filepath: Path) -> pd.DataFrame:
    """Load count matrix from TSV/CSV (optionally gzipped)."""
    suffix = "".join(filepath.suffixes)

    if ".gz" in suffix:
        compression = "gzip"
    else:
        compression = None

    if ".tsv" in suffix or ".txt" in suffix:
        sep = "\t"
    else:
        sep = ","

    return pd.read_csv(filepath, sep=sep, compression=compression, index_col=0)


def create_anndata(
    counts: np.ndarray,
    cell_ids: list[str],
    gene_symbols: list[str],
    ensembl_ids: list[str],
    sample_name: str,
    metadata: dict | None = None,
) -> ad.AnnData:
    """Create AnnData object from count matrix."""
    # Ensure sparse matrix (cells x genes)
    if not isinstance(counts, csr_matrix):
        counts = csr_matrix(counts)

    # Build obs (cell metadata)
    obs = pd.DataFrame(index=cell_ids)
    obs["sample"] = sample_name

    # Build var (gene metadata)
    var = pd.DataFrame(index=gene_symbols)
    var["ensembl_id"] = ensembl_ids

    # Create AnnData
    adata = ad.AnnData(X=counts, obs=obs, var=var)

    # Add metadata to uns
    if metadata:
        for key, value in metadata.items():
            adata.uns[key] = value

    return adata


def save_preview(adata: ad.AnnData, filepath: Path, n_cells: int = 1000, n_genes: int = 100):
    """Save a preview CSV with subset of data."""
    n_cells = min(n_cells, adata.n_obs)
    n_genes = min(n_genes, adata.n_vars)

    preview = pd.DataFrame(
        adata.X[:n_cells, :n_genes].toarray(),
        index=adata.obs_names[:n_cells],
        columns=adata.var_names[:n_genes],
    )
    preview.to_csv(filepath)


def convert_sample(
    input_path: Path,
    output_path: Path,
    species: str = "mouse",
    sample_name: str | None = None,
    metadata: dict | None = None,
    save_preview_csv: bool = True,
) -> ad.AnnData:
    """Convert a single sample to H5AD."""
    if sample_name is None:
        sample_name = input_path.stem.replace(".tsv", "").replace(".csv", "")

    print(f"Loading {input_path.name}...")
    df = load_count_matrix(input_path)

    # Determine orientation (genes x cells or cells x genes)
    # Assume genes are rows if row count > col count
    if df.shape[0] > df.shape[1]:
        # genes x cells -> transpose to cells x genes
        ensembl_ids = list(df.index)
        cell_ids = list(df.columns)
        counts = df.values.T
    else:
        # cells x genes
        cell_ids = list(df.index)
        ensembl_ids = list(df.columns)
        counts = df.values

    print(f"  Shape: {len(cell_ids)} cells x {len(ensembl_ids)} genes")

    # Convert gene IDs
    print(f"  Converting Ensembl IDs to symbols ({species})...")
    gene_symbols, _ = convert_ensembl_to_symbols(ensembl_ids, species)
    gene_symbols = make_unique_symbols(gene_symbols)

    # Create AnnData
    print("  Creating AnnData...")
    adata = create_anndata(
        counts=csr_matrix(counts),
        cell_ids=cell_ids,
        gene_symbols=gene_symbols,
        ensembl_ids=ensembl_ids,
        sample_name=sample_name,
        metadata=metadata,
    )

    # Save
    print(f"  Saving {output_path.name}...")
    adata.write_h5ad(output_path, compression="gzip")

    if save_preview_csv:
        preview_path = output_path.with_suffix(".preview.csv")
        save_preview(adata, preview_path)
        print(f"  Saved preview: {preview_path.name}")

    return adata


def main():
    parser = argparse.ArgumentParser(description="Convert scRNA-seq to H5AD")
    parser.add_argument("--input", help="Input count matrix file")
    parser.add_argument("--input-dir", help="Directory with input files")
    parser.add_argument("--output", help="Output H5AD file")
    parser.add_argument("--output-dir", help="Output directory for batch")
    parser.add_argument(
        "--species", default="mouse", choices=["mouse", "human"],
        help="Species for gene symbol conversion"
    )
    parser.add_argument(
        "--no-preview", action="store_true", help="Skip preview CSV"
    )
    args = parser.parse_args()

    if args.input:
        input_path = Path(args.input)
        output_path = Path(args.output) if args.output else input_path.with_suffix(".h5ad")
        convert_sample(input_path, output_path, args.species, save_preview_csv=not args.no_preview)

    elif args.input_dir:
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir) if args.output_dir else input_dir / "h5ad"
        output_dir.mkdir(parents=True, exist_ok=True)

        patterns = ["*.tsv.gz", "*.tsv", "*.csv.gz", "*.csv"]
        files = []
        for pattern in patterns:
            files.extend(input_dir.glob(pattern))

        print(f"Found {len(files)} files to process")
        for f in files:
            output_path = output_dir / f"{f.stem.replace('.tsv', '').replace('.csv', '')}.h5ad"
            convert_sample(f, output_path, args.species, save_preview_csv=not args.no_preview)
    else:
        parser.error("Either --input or --input-dir is required")


if __name__ == "__main__":
    main()
