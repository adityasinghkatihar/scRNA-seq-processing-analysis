"""
pipeline/annotate.py
--------------------
Cell type annotation using celltypist (pre-trained models) with a
marker-gene score fallback for custom or unsupported tissue types.

STRATEGY:
    1. Try celltypist with an appropriate pre-trained model
    2. If model unavailable or tissue not supported → fall back to
       Scanpy's score_genes() with user-provided marker dict
    3. Store results in adata.obs["cell_type"] in both cases

WHY celltypist FIRST:
    Pre-trained models cover 40+ tissue types and require no manual
    marker curation. One line of code vs. curating marker lists.
    Download happens once and models are cached locally.

WHY KEEP THE FALLBACK:
    Custom organisms, rare tissues, or in-house marker lists need
    the manual route. The fallback ensures the module never dead-ends.

FUTURE MODULE NOTE:
    adata.obs["cell_type"]   — used by GRN to build per-cell-type regulons
    adata.obs["leiden"]      — Leiden labels are kept alongside cell types
                               so users can compare clustering vs annotation
    Cell type labels feed directly into the biological interpretation layer.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


# ---------------------------------------------------------------------------
# Built-in marker gene sets for common cell types
# Used as fallback when celltypist model is unavailable
# Users can override this with their own dict via the API
# ---------------------------------------------------------------------------
DEFAULT_MARKERS: Dict[str, List[str]] = {
    "T_cell":        ["CD3D", "CD3E", "CD3G", "CD2"],
    "B_cell":        ["CD19", "MS4A1", "CD79A", "CD79B"],
    "NK_cell":       ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "Monocyte":      ["CD14", "LYZ", "CST3", "FCGR3A"],
    "Dendritic":     ["FCER1A", "CST3", "IL3RA", "CLEC4C"],
    "Macrophage":    ["CD68", "CSF1R", "MRC1", "MARCO"],
    "Neutrophil":    ["ELANE", "MPO", "AZU1", "PRTN3"],
    "Platelet":      ["PPBP", "PF4", "GP9", "ITGA2B"],
    "Erythrocyte":   ["HBB", "HBA1", "HBA2", "GYPA"],
    "Epithelial":    ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Fibroblast":    ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial":   ["PECAM1", "VWF", "CDH5", "ENG"],
}


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------
@dataclass
class AnnotationResult:
    method:             str            # "celltypist" or "marker_score"
    model_used:         Optional[str]  # celltypist model name if applicable
    cell_type_counts:   dict           # {"T_cell": 345, "B_cell": 210, ...}
    n_unassigned:       int

    def summary(self) -> str:
        top = sorted(self.cell_type_counts.items(), key=lambda x: -x[1])[:5]
        top_str = ", ".join(f"{k}: {v}" for k, v in top)
        return (
            f"Method: {self.method} ({self.model_used or 'no model'})\n"
            f"Top cell types: {top_str}\n"
            f"Unassigned: {self.n_unassigned}"
        )


# ---------------------------------------------------------------------------
# Main annotation function
# ---------------------------------------------------------------------------
def run_annotation(
    adata: ad.AnnData,
    celltypist_model: Optional[str]          = "Immune_All_Low.pkl",
    marker_dict:      Optional[Dict[str, List[str]]] = None,
    min_score:        float                  = 0.5,
) -> tuple[ad.AnnData, AnnotationResult]:
    """
    Annotate cells using celltypist or marker gene scoring.

    Args:
        adata             : clustered AnnData (output of pipeline/cluster.py)
        celltypist_model  : name of pre-trained celltypist model.
                            Full list: celltypist.models.get_all_models()
                            Set to None to skip and use marker scoring.
        marker_dict       : custom marker dict {cell_type: [gene1, gene2, ...]}.
                            If None, DEFAULT_MARKERS is used as fallback.
        min_score         : minimum celltypist probability to accept a label.
                            Cells below this are labeled "Unassigned".

    Returns:
        (adata, AnnotationResult)
        adata.obs["cell_type"]              → final cell type label
        adata.obs["cell_type_confidence"]   → confidence score (celltypist only)
    """
    method_used = None
    model_used  = None

    # Try celltypist first ─────────────────────────────────────────────────
    if celltypist_model is not None:
        try:
            adata, method_used, model_used = _annotate_celltypist(
                adata, celltypist_model, min_score
            )
        except ImportError:
            warnings.warn(
                "celltypist not installed. Falling back to marker gene scoring. "
                "Install with: pip install celltypist",
                UserWarning,
                stacklevel=2,
            )
        except Exception as e:
            warnings.warn(
                f"celltypist annotation failed ({e}). Falling back to marker gene scoring.",
                UserWarning,
                stacklevel=2,
            )

    # Fallback: marker gene scoring ────────────────────────────────────────
    if method_used is None:
        markers = marker_dict or DEFAULT_MARKERS
        adata, method_used = _annotate_marker_score(adata, markers)

    # Build result summary ─────────────────────────────────────────────────
    cell_type_counts = adata.obs["cell_type"].value_counts().to_dict()
    n_unassigned     = int(cell_type_counts.get("Unassigned", 0))

    # Store annotation state
    adata.uns["annotation"] = {
        "completed":   True,
        "method":      method_used,
        "model":       model_used,
        "n_cell_types": len(cell_type_counts),
    }

    result = AnnotationResult(
        method           = method_used,
        model_used       = model_used,
        cell_type_counts = {k: int(v) for k, v in cell_type_counts.items()},
        n_unassigned     = n_unassigned,
    )

    return adata, result


# ---------------------------------------------------------------------------
# celltypist backend
# ---------------------------------------------------------------------------
def _annotate_celltypist(
    adata: ad.AnnData,
    model_name: str,
    min_score: float,
) -> tuple[ad.AnnData, str, str]:
    """
    Run celltypist annotation.
    Returns (adata, method_str, model_name)
    """
    import celltypist
    from celltypist import models

    # Download model if not cached (happens once)
    model = models.Model.load(model=model_name)

    # celltypist expects log1p-normalized data with 10,000 counts/cell
    # which is exactly what pipeline/normalize.py produces
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=True,   # smooth predictions using KNN graph
    )

    # Transfer predictions to adata.obs
    pred_adata = predictions.to_adata()
    adata.obs["cell_type"]            = pred_adata.obs["majority_voting"]
    adata.obs["cell_type_confidence"] = pred_adata.obs["conf_score"]

    # Apply minimum confidence threshold
    low_conf = adata.obs["cell_type_confidence"] < min_score
    adata.obs.loc[low_conf, "cell_type"] = "Unassigned"

    return adata, "celltypist", model_name


# ---------------------------------------------------------------------------
# Marker gene score backend
# ---------------------------------------------------------------------------
def _annotate_marker_score(
    adata: ad.AnnData,
    marker_dict: Dict[str, List[str]],
) -> tuple[ad.AnnData, str]:
    """
    Score each cell type using Scanpy's score_genes(), then assign
    the highest-scoring label to each cell.

    score_genes() adds a column to adata.obs for each cell type.
    We then take the argmax across all score columns.
    """
    score_cols = []

    for cell_type, markers in marker_dict.items():
        # Only score with genes that exist in the dataset
        available = [g for g in markers if g in adata.var_names]
        if len(available) < 2:
            # Too few markers to produce a meaningful score — skip
            continue

        col = f"score_{cell_type}"
        sc.tl.score_genes(adata, gene_list=available, score_name=col)
        score_cols.append((cell_type, col))

    if not score_cols:
        warnings.warn(
            "No marker genes found in dataset. All cells labeled 'Unassigned'.",
            UserWarning,
            stacklevel=2,
        )
        adata.obs["cell_type"] = "Unassigned"
        return adata, "marker_score"

    # Build score matrix and find best label per cell
    score_df = pd.DataFrame(
        {ct: adata.obs[col] for ct, col in score_cols},
        index=adata.obs_names,
    )

    # Assign cell type with highest score
    adata.obs["cell_type"] = score_df.idxmax(axis=1)

    # Cells where all scores are negative get "Unassigned"
    all_negative = (score_df < 0).all(axis=1)
    adata.obs.loc[all_negative, "cell_type"] = "Unassigned"

    # Clean up temporary score columns from adata.obs
    temp_cols = [col for _, col in score_cols]
    adata.obs.drop(columns=temp_cols, inplace=True, errors="ignore")

    return adata, "marker_score"
