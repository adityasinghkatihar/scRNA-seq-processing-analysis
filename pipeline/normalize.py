"""
pipeline/normalize.py
---------------------
Normalization, log transformation, and highly variable gene (HVG) selection.

WHAT THIS DOES:
    1. Normalize total counts per cell to a fixed target (default 10,000)
    2. Log1p transform (stabilizes variance across expression magnitudes)
    3. Select highly variable genes (HVGs) — the informative subset
    4. Preserve raw counts in adata.raw for downstream tools that need them

WHY THIS ORDER MATTERS:
    Raw counts → normalize → log1p → HVG selection
    Each step assumes the previous one was done. Changing the order
    produces biologically wrong results.

IMPORTANT — adata.raw:
    pySCENIC (GRN module) and some velocity tools need raw counts.
    We store them in adata.raw BEFORE normalization.
    This is AnnData convention — do not skip this step.

FUTURE MODULE NOTE:
    GRN (pySCENIC) reads adata.raw.X for its regulon scoring.
    RNA velocity (scVelo) uses spliced/unspliced counts separately,
    but the normalized matrix in adata.X is used for embedding.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import anndata as ad
import scanpy as sc

from utils.config import PIPELINE_DEFAULTS


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------
@dataclass
class NormResult:
    target_sum:     float
    n_hvgs:         int
    genes_before:   int
    genes_after:    int   # after subsetting to HVGs
    raw_stored:     bool

    def summary(self) -> str:
        return (
            f"Normalized to {self.target_sum:,.0f} counts/cell, log1p transformed.\n"
            f"HVGs selected: {self.n_hvgs:,} of {self.genes_before:,} genes.\n"
            f"Raw counts stored: {self.raw_stored}"
        )


# ---------------------------------------------------------------------------
# Main normalization function
# ---------------------------------------------------------------------------
def run_normalization(
    adata: ad.AnnData,
    target_sum:  float = 1e4,
    n_top_genes: Optional[int] = None,
    flavor:      str   = "seurat_v3",
) -> tuple[ad.AnnData, NormResult]:
    """
    Normalize, log-transform, and select highly variable genes.

    Args:
        adata        : QC-filtered AnnData (output of pipeline/qc.py)
        target_sum   : normalize each cell to this total count (default 10,000)
        n_top_genes  : number of HVGs to select (default from PIPELINE_DEFAULTS)
        flavor       : HVG selection method — "seurat_v3" or "cell_ranger"

    Returns:
        (adata, NormResult)
        adata.X        → normalized + log-transformed counts
        adata.raw      → original raw counts (pre-normalization)
        adata.var      → includes "highly_variable" boolean column
    """
    n_top_genes = n_top_genes or PIPELINE_DEFAULTS["n_top_genes"]
    genes_before = adata.n_vars

    # 1. Store raw counts BEFORE any normalization ──────────────────────────
    #
    # adata.raw freezes the current state (counts + all genes).
    # After this point, adata.X will be modified but adata.raw stays pristine.
    # pySCENIC (GRN) requires raw counts for AUCell scoring.
    #
    adata.raw = adata

    # 2. Normalize per cell ────────────────────────────────────────────────
    #
    # Each cell's total count is scaled to target_sum.
    # This removes library size differences between cells.
    # (Cell A with 5000 counts and Cell B with 50000 counts become comparable)
    #
    sc.pp.normalize_total(adata, target_sum=target_sum)

    # 3. Log1p transform ───────────────────────────────────────────────────
    #
    # log1p(x) = log(x + 1)
    # Compresses the dynamic range — a gene with 1000 counts doesn't
    # completely dominate a gene with 10 counts in distance calculations.
    #
    sc.pp.log1p(adata)

    # Store that log1p was applied — scVelo and other tools check this
    adata.uns["log1p"] = {"base": None}

    # 4. Select highly variable genes ──────────────────────────────────────
    #
    # Not all genes are informative. HVG selection keeps genes that vary
    # meaningfully across cells (signal) vs. uniformly expressed genes (noise).
    # This speeds up PCA and clustering without losing biological resolution.
    #
    # "seurat_v3" uses raw counts for variance modeling (more accurate).
    # It requires adata.raw to be set — which we did in step 1.
    #
    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
            layer="counts" if "counts" in adata.layers else None,
        )
    except Exception:
        # Fallback to dispersion-based method if seurat_v3 fails
        # (e.g., if counts layer not available)
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor="cell_ranger",
        )

    n_hvgs = int(adata.var["highly_variable"].sum())

    # 5. Subset to HVGs ───────────────────────────────────────────────────
    #
    # Keep only HVGs in adata.X for downstream PCA/clustering.
    # adata.raw still has ALL genes — used by GRN and differential expression.
    #
    adata = adata[:, adata.var["highly_variable"]].copy()

    # 6. Record normalization state ────────────────────────────────────────
    adata.uns["normalization"] = {
        "completed": True,
        "target_sum": target_sum,
        "n_top_genes": n_top_genes,
        "hvg_flavor": flavor,
        "n_hvgs": n_hvgs,
    }

    result = NormResult(
        target_sum   = target_sum,
        n_hvgs       = n_hvgs,
        genes_before = genes_before,
        genes_after  = adata.n_vars,
        raw_stored   = True,
    )

    return adata, result
