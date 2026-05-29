"""
utils/io.py
-----------
Reads count matrices in multiple formats and builds AnnData objects.

WHY ANNDATA AS THE CONTRACT:
    AnnData is the standard container for single-cell data.
    - adata.X   → count matrix (cells × genes)
    - adata.obs → cell-level metadata (condition, timepoint, etc.)
    - adata.var → gene-level metadata (gene names, is_mito, etc.)

    Every downstream module — RNA velocity, GRN, trajectory —
    reads the same adata.obs without needing to re-load metadata.
    This is the core of the future-proof design.

SUPPORTED INPUT FORMATS:
    - 10x HDF5   (.h5)
    - 10x MEX    (directory with matrix.mtx, barcodes.tsv, features.tsv)
    - Loom        (.loom)
    - CSV/TSV     (cells × genes, cell IDs as index)
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path
from typing import Optional

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse as sp

from utils.config import METADATA_CONFIG, OPTIONAL_COLS
from utils.metadata import (
    load_metadata,
    validate_metadata,
    align_metadata,
    fill_missing_optional,
    ValidationResult,
)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------
def build_anndata(
    counts_path: str,
    metadata_path: str,
) -> tuple[ad.AnnData, ValidationResult]:
    """
    Build a validated AnnData object from a counts matrix and metadata CSV.

    This is the single function called by the API and dashboard.
    It handles format detection, metadata validation, and alignment.

    Args:
        counts_path   : path to count matrix (auto-detects format)
        metadata_path : path to metadata.csv

    Returns:
        (adata, validation_result)
        adata             — fully constructed AnnData with metadata in adata.obs
        validation_result — contains any warnings / errors from metadata checks

    Raises:
        ValueError for unrecoverable errors (missing required columns, no overlap)
    """
    # 1. Load counts ───────────────────────────────────────────────────────
    adata = _read_counts(counts_path)

    # 2. Load metadata ─────────────────────────────────────────────────────
    meta_df = load_metadata(metadata_path)

    # 3. Validate metadata ─────────────────────────────────────────────────
    val_result = validate_metadata(meta_df)

    # Hard stop on errors — don't build a corrupt AnnData
    if not val_result.is_valid:
        raise ValueError(
            "Metadata validation failed:\n" + "\n".join(val_result.errors)
        )

    # 4. Align metadata to counts matrix cell order ────────────────────────
    # This ensures adata.X[i] and adata.obs.iloc[i] always refer to the same cell
    aligned_meta = align_metadata(list(adata.obs_names), meta_df)

    # Subset adata to only cells that exist in both counts and metadata
    adata = adata[aligned_meta.index, :].copy()

    # 5. Fill missing optional columns with 'unknown' ──────────────────────
    aligned_meta = fill_missing_optional(aligned_meta)

    # 6. Attach metadata to adata.obs ──────────────────────────────────────
    #
    # WHY adata.obs:
    #   This is the AnnData-standard location for cell metadata.
    #   All downstream tools (scVelo, pySCENIC, CellOracle) read from here.
    #   RNA velocity will read adata.obs["timepoint"] for pseudotime ordering.
    #   GRN analysis will read adata.obs["condition"] to split regulon networks.
    #
    adata.obs = aligned_meta

    # 7. Annotate mitochondrial genes ──────────────────────────────────────
    _flag_mito_genes(adata)

    # 8. Store pipeline provenance in adata.uns ────────────────────────────
    #
    # adata.uns = unstructured metadata (pipeline params, software versions, etc.)
    # Future modules append their own entries here for full reproducibility.
    #
    adata.uns["pipeline"] = {
        "counts_source": str(counts_path),
        "metadata_source": str(metadata_path),
        "metadata_columns": list(aligned_meta.columns),
    }

    return adata, val_result


# ---------------------------------------------------------------------------
# Format detection + readers
# ---------------------------------------------------------------------------
def _read_counts(path: str) -> ad.AnnData:
    """
    Auto-detect count matrix format and load into AnnData.
    Returns AnnData with only adata.X, adata.obs_names, adata.var_names.
    Metadata is attached later in build_anndata().
    """
    path = Path(path)

    if path.is_dir():
        return _read_10x_mex(path)

    ext = "".join(path.suffixes).lower()

    if ext in (".h5",):
        return _read_10x_h5(path)
    elif ext in (".loom",):
        return _read_loom(path)
    elif ext in (".csv", ".tsv", ".txt"):
        return _read_dense_csv(path)
    elif ext in (".h5ad",):
        # Already an AnnData — load as-is, metadata will be overlaid
        return ad.read_h5ad(path)
    else:
        raise ValueError(
            f"Unrecognized counts matrix format: '{ext}'. "
            f"Supported: 10x HDF5 (.h5), 10x MEX (directory), "
            f"Loom (.loom), CSV/TSV (.csv/.tsv), AnnData (.h5ad)"
        )


def _read_10x_h5(path: Path) -> ad.AnnData:
    """Read 10x Genomics HDF5 format."""
    import scanpy as sc
    adata = sc.read_10x_h5(str(path))
    adata.var_names_make_unique()
    return adata


def _read_10x_mex(directory: Path) -> ad.AnnData:
    """Read 10x Genomics MEX format (matrix.mtx + barcodes.tsv + features.tsv)."""
    import scanpy as sc
    adata = sc.read_10x_mtx(str(directory), var_names="gene_symbols")
    adata.var_names_make_unique()
    return adata


def _read_loom(path: Path) -> ad.AnnData:
    """Read Loom format."""
    adata = ad.read_loom(str(path))
    adata.var_names_make_unique()
    return adata


def _read_dense_csv(path: Path) -> ad.AnnData:
    """
    Read dense CSV/TSV where rows = cells, columns = genes.
    First column must be cell IDs.
    """
    sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(path, index_col=0, sep=sep)

    # Validate it looks like a counts matrix
    if df.shape[0] == 0 or df.shape[1] == 0:
        raise ValueError("Counts matrix CSV is empty.")

    # Check for non-numeric values
    non_numeric = df.select_dtypes(exclude="number").columns.tolist()
    if non_numeric:
        raise ValueError(
            f"Counts matrix has non-numeric columns: {non_numeric[:5]}. "
            f"Ensure only gene count columns are present (no metadata columns)."
        )

    X = sp.csr_matrix(df.values)
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=df.index),
        var=pd.DataFrame(index=df.columns),
    )
    return adata


# ---------------------------------------------------------------------------
# Gene annotation helpers
# ---------------------------------------------------------------------------
def _flag_mito_genes(adata: ad.AnnData) -> None:
    """
    Identify mitochondrial genes and store flag in adata.var.

    WHY:
        MT gene % is a key QC metric.
        Flagging here once means QC module doesn't need to know
        how to detect MT genes — it just reads adata.var["is_mito"].
    """
    # Handles human (MT-) and mouse (mt-) naming conventions
    mito_mask = (
        adata.var_names.str.startswith("MT-") |
        adata.var_names.str.startswith("mt-")
    )
    adata.var["is_mito"] = mito_mask

    if mito_mask.sum() == 0:
        warnings.warn(
            "No mitochondrial genes detected (expected names starting with 'MT-' or 'mt-'). "
            "MT% QC metric will be zero for all cells.",
            UserWarning,
            stacklevel=2,
        )


# ---------------------------------------------------------------------------
# Persist / reload
# ---------------------------------------------------------------------------
def save_anndata(adata: ad.AnnData, output_path: str) -> None:
    """Save AnnData to h5ad. All metadata in adata.obs is preserved."""
    adata.write_h5ad(output_path)


def load_anndata(path: str) -> ad.AnnData:
    """Load a previously saved AnnData. Metadata in adata.obs is restored."""
    return ad.read_h5ad(path)
