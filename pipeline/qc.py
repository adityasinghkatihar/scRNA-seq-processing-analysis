"""
pipeline/qc.py
--------------
Quality control: filter low-quality cells and uninformative genes.

WHAT THIS DOES:
    1. Calculate per-cell QC metrics (gene count, total counts, MT%)
    2. Filter cells that fall outside acceptable thresholds
    3. Filter genes expressed in too few cells
    4. Store QC metrics in adata.obs so the dashboard can visualize them

WHY QC MATTERS:
    Raw scRNA-seq data contains:
    - Dead/damaged cells (high MT%, low gene count)
    - Empty droplets (very low counts)
    - Doublets (unusually high gene/count numbers)
    Keeping these in distorts clustering and downstream analysis.

FUTURE MODULE NOTE:
    adata.obs will carry QC columns (n_genes_by_counts, pct_counts_mt, etc.)
    RNA velocity and GRN modules can use these to further filter if needed.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional

import anndata as ad
import scanpy as sc

from utils.config import PIPELINE_DEFAULTS


# ---------------------------------------------------------------------------
# QC result — structured summary returned to API / dashboard
# ---------------------------------------------------------------------------
@dataclass
class QCResult:
    cells_before: int
    cells_after:  int
    genes_before: int
    genes_after:  int
    cells_removed: int
    genes_removed: int
    removal_reasons: dict  # {"low_genes": N, "high_mt": N, ...}
    thresholds_used: dict

    def summary(self) -> str:
        return (
            f"Cells: {self.cells_before:,} → {self.cells_after:,} "
            f"({self.cells_removed:,} removed)\n"
            f"Genes: {self.genes_before:,} → {self.genes_after:,} "
            f"({self.genes_removed:,} removed)"
        )


# ---------------------------------------------------------------------------
# Main QC function
# ---------------------------------------------------------------------------
def run_qc(
    adata: ad.AnnData,
    min_genes:   Optional[int]   = None,
    max_genes:   Optional[int]   = None,
    max_mt_pct:  Optional[float] = None,
    min_cells:   int             = 3,
) -> tuple[ad.AnnData, QCResult]:
    """
    Run quality control filtering on an AnnData object.

    All thresholds default to PIPELINE_DEFAULTS if not provided.
    This allows per-run overrides via API request body.

    Args:
        adata      : AnnData with adata.obs populated (output of io.build_anndata)
        min_genes  : minimum genes per cell (cells below this are filtered)
        max_genes  : maximum genes per cell (cells above may be doublets)
        max_mt_pct : maximum mitochondrial gene percentage per cell
        min_cells  : minimum cells a gene must appear in (genes below are filtered)

    Returns:
        (filtered_adata, QCResult)
    """
    # Use defaults if not overridden
    min_genes  = min_genes  if min_genes  is not None else PIPELINE_DEFAULTS["min_genes"]
    max_genes  = max_genes  if max_genes  is not None else PIPELINE_DEFAULTS["max_genes"]
    max_mt_pct = max_mt_pct if max_mt_pct is not None else PIPELINE_DEFAULTS["max_mt_pct"]

    cells_before = adata.n_obs
    genes_before = adata.n_vars

    # 1. Calculate QC metrics ──────────────────────────────────────────────
    #
    # Adds to adata.obs:
    #   n_genes_by_counts  — number of genes with > 0 counts
    #   total_counts       — total UMI count per cell
    #   pct_counts_mt      — % of counts from mitochondrial genes
    #
    # Adds to adata.var:
    #   n_cells_by_counts  — number of cells expressing this gene
    #   mean_counts        — mean expression across cells
    #
    qc_vars = ["is_mito"] if "is_mito" in adata.var.columns else []
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, inplace=True)

    # 2. Track removal reasons ─────────────────────────────────────────────
    removal_reasons = {}

    low_gene_mask  = adata.obs["n_genes_by_counts"] < min_genes
    high_gene_mask = adata.obs["n_genes_by_counts"] > max_genes
    removal_reasons["low_genes"]  = int(low_gene_mask.sum())
    removal_reasons["high_genes"] = int(high_gene_mask.sum())

    if "pct_counts_is_mito" in adata.obs.columns:
        high_mt_mask = adata.obs["pct_counts_is_mito"] > max_mt_pct
        removal_reasons["high_mt_pct"] = int(high_mt_mask.sum())
    else:
        high_mt_mask = None
        removal_reasons["high_mt_pct"] = 0

    # 3. Filter cells ──────────────────────────────────────────────────────
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_cells(adata, max_genes=max_genes)

    if high_mt_mask is not None:
        # Re-evaluate after cell filtering above (obs index may have changed)
        if "pct_counts_is_mito" in adata.obs.columns:
            adata = adata[adata.obs["pct_counts_is_mito"] <= max_mt_pct, :].copy()

    # 4. Filter genes ──────────────────────────────────────────────────────
    sc.pp.filter_genes(adata, min_cells=min_cells)

    cells_after = adata.n_obs
    genes_after = adata.n_vars

    # 5. Store QC status in adata.uns for downstream access ────────────────
    adata.uns["qc"] = {
        "completed": True,
        "thresholds": {
            "min_genes":   min_genes,
            "max_genes":   max_genes,
            "max_mt_pct":  max_mt_pct,
            "min_cells":   min_cells,
        },
        "cells_before": cells_before,
        "cells_after":  cells_after,
        "genes_before": genes_before,
        "genes_after":  genes_after,
    }

    result = QCResult(
        cells_before   = cells_before,
        cells_after    = cells_after,
        genes_before   = genes_before,
        genes_after    = genes_after,
        cells_removed  = cells_before - cells_after,
        genes_removed  = genes_before - genes_after,
        removal_reasons= removal_reasons,
        thresholds_used= adata.uns["qc"]["thresholds"],
    )

    return adata, result
