"""
utils/metadata.py
-----------------
All metadata operations live here.

WHY CENTRALIZED:
    Keeps validation, alignment, and loading logic in one place.
    Every other module imports from here — no duplication.

CONTRACT:
    Input  → raw CSV path (str) or DataFrame
    Output → clean, validated, cell-aligned DataFrame
             ready to be stored in adata.obs

FUTURE MODULE USE:
    RNA velocity : adata.obs[config["timepoint_col"]] → cell ordering
    GRN          : adata.obs[config["condition_col"]] → per-condition networks
    Trajectory   : adata.obs[config["batch_col"]]     → batch-aware analysis
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
import pandas as pd

from utils.config import METADATA_CONFIG, REQUIRED_COLS, OPTIONAL_COLS


# ---------------------------------------------------------------------------
# Validation result — carries warnings and errors cleanly
# ---------------------------------------------------------------------------
@dataclass
class ValidationResult:
    is_valid: bool = True
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    def add_error(self, msg: str):
        self.errors.append(msg)
        self.is_valid = False

    def add_warning(self, msg: str):
        self.warnings.append(msg)

    def summary(self) -> str:
        lines = []
        for e in self.errors:
            lines.append(f"ERROR: {e}")
        for w in self.warnings:
            lines.append(f"WARNING: {w}")
        return "\n".join(lines) if lines else "All checks passed."


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
def load_metadata(filepath: str) -> pd.DataFrame:
    """
    Load metadata CSV into a DataFrame.

    Sets cell_id as the index so it aligns with AnnData convention
    (adata.obs is indexed by cell barcodes / IDs).

    Args:
        filepath: path to metadata.csv

    Returns:
        DataFrame with cell_id as index.

    Raises:
        ValueError if file cannot be read or cell_id column is missing.
    """
    try:
        meta_df = pd.read_csv(filepath, dtype=str)  # read all as str first; parse types later
    except Exception as e:
        raise ValueError(f"Could not read metadata file: {e}")

    cell_id_col = METADATA_CONFIG["cell_id_col"]

    if cell_id_col not in meta_df.columns:
        raise ValueError(
            f"Metadata CSV must have a '{cell_id_col}' column. "
            f"Found columns: {list(meta_df.columns)}"
        )

    # Set cell_id as index — this is how adata.obs is structured in AnnData
    meta_df = meta_df.set_index(cell_id_col)

    return meta_df


# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------
def validate_metadata(meta_df: pd.DataFrame) -> ValidationResult:
    """
    Run all validation checks on metadata DataFrame.

    Checks:
        1. Required columns exist
        2. No duplicate cell IDs
        3. Missing values in required columns
        4. Missing values in optional columns (warnings only)

    Args:
        meta_df: DataFrame with cell_id as index (output of load_metadata)

    Returns:
        ValidationResult with .is_valid, .errors, .warnings
    """
    result = ValidationResult()
    condition_col = METADATA_CONFIG["condition_col"]

    # 1. Required columns ─────────────────────────────────────────────────
    for col in REQUIRED_COLS:
        if col == METADATA_CONFIG["cell_id_col"]:
            continue  # already enforced as index in load_metadata
        if col not in meta_df.columns:
            result.add_error(
                f"Required column '{col}' is missing from metadata.csv. "
                f"Found: {list(meta_df.columns)}"
            )

    # 2. Duplicate cell IDs ───────────────────────────────────────────────
    if meta_df.index.duplicated().any():
        dupes = meta_df.index[meta_df.index.duplicated()].tolist()
        result.add_error(
            f"Duplicate cell IDs found in metadata: {dupes[:5]}"
            + (" ..." if len(dupes) > 5 else "")
        )

    # Stop here if required columns are missing — further checks would error
    if not result.is_valid:
        return result

    # 3. Missing values in required columns ───────────────────────────────
    if condition_col in meta_df.columns:
        n_missing = meta_df[condition_col].isna().sum()
        if n_missing > 0:
            result.add_error(
                f"Required column '{condition_col}' has {n_missing} missing values. "
                f"All cells must have a condition label."
            )

    # 4. Missing values in optional columns (warn, don't fail) ────────────
    for col in OPTIONAL_COLS:
        if col in meta_df.columns:
            n_missing = meta_df[col].isna().sum()
            if n_missing > 0:
                pct = 100 * n_missing / len(meta_df)
                result.add_warning(
                    f"Optional column '{col}' has {n_missing} missing values ({pct:.1f}%). "
                    f"Cells without a value will be labeled 'unknown'."
                )

    return result


# ---------------------------------------------------------------------------
# Align
# ---------------------------------------------------------------------------
def align_metadata(counts_cell_ids: list, meta_df: pd.DataFrame) -> pd.DataFrame:
    """
    Align metadata rows to the cell order in the counts matrix.

    WHY THIS MATTERS:
        AnnData requires adata.obs to have the exact same row order
        as adata.X (the counts matrix). If they don't match, gene
        expression and metadata get silently misassigned to wrong cells.

    Args:
        counts_cell_ids : list of cell IDs from the counts matrix (in order)
        meta_df         : DataFrame with cell_id as index

    Returns:
        Aligned DataFrame matching counts matrix cell order.

    Raises:
        ValueError if cell IDs don't overlap sufficiently.
    """
    counts_set = set(counts_cell_ids)
    meta_set   = set(meta_df.index)

    in_counts_not_meta = counts_set - meta_set
    in_meta_not_counts = meta_set - counts_set
    overlap            = counts_set & meta_set

    if len(overlap) == 0:
        raise ValueError(
            "No cell IDs overlap between counts matrix and metadata.csv. "
            "Check that cell_id values match exactly (case-sensitive, no extra spaces). "
            f"Example counts IDs: {list(counts_cell_ids)[:3]} | "
            f"Example metadata IDs: {list(meta_df.index)[:3]}"
        )

    if in_counts_not_meta:
        warnings.warn(
            f"{len(in_counts_not_meta)} cells in counts matrix have no metadata entry. "
            f"These cells will be dropped. Example: {list(in_counts_not_meta)[:3]}",
            UserWarning,
            stacklevel=2,
        )

    if in_meta_not_counts:
        warnings.warn(
            f"{len(in_meta_not_counts)} metadata rows have no matching cell in counts. "
            f"These rows will be ignored. Example: {list(in_meta_not_counts)[:3]}",
            UserWarning,
            stacklevel=2,
        )

    # Keep only cells that appear in both, in counts matrix order
    aligned_ids = [cid for cid in counts_cell_ids if cid in meta_set]
    return meta_df.loc[aligned_ids]


# ---------------------------------------------------------------------------
# Fill missing optional columns
# ---------------------------------------------------------------------------
def fill_missing_optional(meta_df: pd.DataFrame) -> pd.DataFrame:
    """
    Fill NaN values in optional columns with 'unknown'.

    Keeps downstream code clean — no NaN checks scattered everywhere.
    Modules can filter on 'unknown' to exclude cells without that label.
    """
    meta_df = meta_df.copy()
    for col in OPTIONAL_COLS:
        if col in meta_df.columns:
            meta_df[col] = meta_df[col].fillna("unknown")
    return meta_df


# ---------------------------------------------------------------------------
# Convenience: get present optional columns
# ---------------------------------------------------------------------------
def get_available_filter_cols(meta_df: pd.DataFrame) -> List[str]:
    """
    Return list of columns that can be used as dashboard filters.
    Only returns columns that are actually present in metadata.

    Used by the dashboard to dynamically build filter widgets.
    """
    filter_cols = []
    for key in ["condition_col", "timepoint_col", "replicate_col",
                "batch_col", "species_col", "treatment_col"]:
        col = METADATA_CONFIG[key]
        if col in meta_df.columns:
            filter_cols.append(col)
    return filter_cols
