"""
pipeline/pathway.py
-------------------
Per-cluster pathway enrichment analysis using gseapy.

WHAT THIS DOES:
    For each cluster (or cell type):
    1. Find differentially expressed genes (DEGs) vs. all other cells
    2. Run enrichment against curated gene sets (GO, KEGG, Reactome, MSigDB)
    3. Return structured JSON — one entry per cluster with top pathways

WHY STRUCTURED JSON OUTPUT:
    This JSON is the input to the biological interpretation layer.
    The LLM agent (future) reads it to generate narrative summaries.
    Keeping it structured and consistent now saves a rewrite later.

OUTPUT SCHEMA (per cluster):
    {
        "cluster_id": "0",
        "cell_type": "T_cell",          # from annotation step
        "n_cells": 345,
        "top_degs": ["gene1", ...],
        "pathways": [
            {
                "term": "T cell activation",
                "gene_set": "GO_Biological_Process_2021",
                "p_value": 0.0001,
                "adjusted_p": 0.002,
                "overlap_genes": ["CD3D", "CD69", ...]
            },
            ...
        ]
    }

FUTURE MODULE NOTE:
    This JSON is consumed by:
    - dashboard/app.py        → pathway cards per cluster
    - interpretation layer    → LLM narrative generation
    - GRN module              → cross-reference regulons with enriched pathways
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import anndata as ad
import pandas as pd
import scanpy as sc


# Gene sets to query — ordered by priority
# Reduced list avoids rate-limiting the Enrichr API
DEFAULT_GENE_SETS = [
    "GO_Biological_Process_2021",
    "KEGG_2021_Human",
    "Reactome_2022",
]

# Mouse-specific alternative when species = mouse
MOUSE_GENE_SETS = [
    "GO_Biological_Process_2021",
    "KEGG_2019_Mouse",
]


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------
@dataclass
class PathwayResult:
    n_clusters:   int
    gene_sets:    List[str]
    results:      List[dict]  # one entry per cluster (schema above)
    failed_clusters: List[str] = field(default_factory=list)

    def to_json_ready(self) -> List[dict]:
        """Return results list ready for json.dumps()."""
        return self.results

    def summary(self) -> str:
        successful = self.n_clusters - len(self.failed_clusters)
        return (
            f"Pathway analysis complete: {successful}/{self.n_clusters} clusters.\n"
            f"Gene sets used: {', '.join(self.gene_sets)}"
        )


# ---------------------------------------------------------------------------
# Main pathway function
# ---------------------------------------------------------------------------
def run_pathway_analysis(
    adata: ad.AnnData,
    cluster_col:   str          = "leiden",
    gene_sets:     Optional[List[str]] = None,
    n_top_degs:    int          = 100,
    p_cutoff:      float        = 0.05,
    species:       str          = "human",
) -> tuple[ad.AnnData, PathwayResult]:
    """
    Run per-cluster differential expression + pathway enrichment.

    Args:
        adata        : annotated AnnData (output of pipeline/annotate.py)
        cluster_col  : column in adata.obs to group cells by (default "leiden")
        gene_sets    : list of Enrichr gene set databases to query
        n_top_degs   : number of top DEGs per cluster to use for enrichment
        p_cutoff     : adjusted p-value cutoff for significant pathways
        species      : "human" or "mouse" — affects gene set selection

    Returns:
        (adata, PathwayResult)
        adata.uns["pathway_results"] → same as PathwayResult.results (for persistence)
    """
    try:
        import gseapy as gp
    except ImportError:
        raise ImportError(
            "gseapy is required for pathway analysis. "
            "Install with: pip install gseapy"
        )

    if cluster_col not in adata.obs.columns:
        raise ValueError(
            f"Cluster column '{cluster_col}' not found in adata.obs. "
            f"Run pipeline/cluster.py first."
        )

    # Select gene sets based on species
    if gene_sets is None:
        gene_sets = MOUSE_GENE_SETS if species.lower() == "mouse" else DEFAULT_GENE_SETS

    clusters = sorted(adata.obs[cluster_col].unique().tolist())
    all_results   = []
    failed_clusters = []

    # 1. Compute DEGs for all clusters at once (more efficient than per-cluster)
    _compute_degs(adata, cluster_col)

    # 2. Per-cluster enrichment ────────────────────────────────────────────
    for cluster_id in clusters:
        try:
            cluster_result = _enrich_cluster(
                adata       = adata,
                cluster_id  = cluster_id,
                cluster_col = cluster_col,
                gene_sets   = gene_sets,
                n_top_degs  = n_top_degs,
                p_cutoff    = p_cutoff,
                gp          = gp,
            )
            all_results.append(cluster_result)

        except Exception as e:
            warnings.warn(
                f"Pathway enrichment failed for cluster {cluster_id}: {e}",
                UserWarning,
                stacklevel=2,
            )
            failed_clusters.append(str(cluster_id))
            # Add empty entry so downstream code doesn't break
            all_results.append(_empty_cluster_entry(cluster_id, adata, cluster_col))

    # 3. Store in adata.uns for persistence ────────────────────────────────
    # AnnData/HDF5 can't reliably serialize list[dict] with mixed nested values.
    # Persist pathway payload as JSON text so writes are stable across versions.
    adata.uns["pathway_results"] = json.dumps(all_results)
    adata.uns["pathway_config"] = {
        "gene_sets": gene_sets,
        "n_top_degs": n_top_degs,
        "p_cutoff": p_cutoff,
        "species": species,
    }

    result = PathwayResult(
        n_clusters      = len(clusters),
        gene_sets       = gene_sets,
        results         = all_results,
        failed_clusters = failed_clusters,
    )

    return adata, result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _compute_degs(adata: ad.AnnData, cluster_col: str) -> None:
    """
    Run Wilcoxon rank-sum test for all clusters vs. rest.
    Results stored in adata.uns["rank_genes_groups"].
    Uses adata.raw if available (raw counts preferred for DEG testing).
    """
    # Use raw counts for DEG testing if available
    # Raw counts better reflect true expression differences
    test_adata = adata.raw.to_adata() if adata.raw is not None else adata

    # Re-attach obs so cluster labels are available
    test_adata.obs = adata.obs

    sc.tl.rank_genes_groups(
        test_adata,
        groupby=cluster_col,
        method="wilcoxon",
        use_raw=False,
        pts=True,      # compute fraction of cells expressing each gene
    )
    # Copy DEG results back to main adata
    adata.uns["rank_genes_groups"] = test_adata.uns["rank_genes_groups"]


def _enrich_cluster(
    adata: ad.AnnData,
    cluster_id: str,
    cluster_col: str,
    gene_sets: List[str],
    n_top_degs: int,
    p_cutoff: float,
    gp,
) -> dict:
    """Run enrichment for a single cluster and return structured result dict."""

    # Extract top DEGs for this cluster
    degs = sc.get.rank_genes_groups_df(
        adata,
        group=str(cluster_id),
        key="rank_genes_groups",
    )

    # Filter to upregulated genes only (logfoldchange > 0)
    degs = degs[degs["logfoldchanges"] > 0].head(n_top_degs)
    top_genes = degs["names"].tolist()

    if len(top_genes) < 5:
        return _empty_cluster_entry(cluster_id, adata, cluster_col,
                                     reason="insufficient DEGs")

    # Get cell type label for this cluster (if annotation was run)
    cell_type = _get_cluster_cell_type(adata, cluster_id, cluster_col)
    n_cells   = int((adata.obs[cluster_col] == cluster_id).sum())

    # Run Enrichr via gseapy
    enr = gp.enrichr(
        gene_list   = top_genes,
        gene_sets   = gene_sets,
        outdir      = None,     # no file output — results in memory only
        verbose     = False,
        cutoff      = p_cutoff,
    )

    # Parse results into clean list
    pathways = _parse_enrichr_results(enr.results, p_cutoff)

    return {
        "cluster_id":    str(cluster_id),
        "cell_type":     cell_type,
        "n_cells":       n_cells,
        "top_degs":      top_genes[:20],   # top 20 for display
        "n_degs_tested": len(top_genes),
        "pathways":      pathways,
    }


def _parse_enrichr_results(results_df: pd.DataFrame, p_cutoff: float) -> List[dict]:
    """Convert gseapy results DataFrame into clean list of dicts."""
    if results_df is None or results_df.empty:
        return []

    # gseapy column names vary slightly by version — handle both
    pval_col = "Adjusted P-value" if "Adjusted P-value" in results_df.columns else "P-value"
    sig = results_df[results_df[pval_col] <= p_cutoff].copy()
    sig = sig.sort_values(pval_col).head(10)  # top 10 per cluster

    pathways = []
    for _, row in sig.iterrows():
        # Parse overlap genes from "gene1;gene2;..." format
        genes_str = row.get("Genes", "")
        overlap   = [g.strip() for g in genes_str.split(";") if g.strip()]

        pathways.append({
            "term":          row.get("Term", "Unknown"),
            "gene_set":      row.get("Gene_set", ""),
            "p_value":       float(row.get("P-value", 1.0)),
            "adjusted_p":    float(row.get("Adjusted P-value", 1.0)),
            "overlap_genes": overlap[:10],   # cap at 10 for readability
            "odds_ratio":    float(row.get("Odds Ratio", 0.0)),
        })

    return pathways


def _get_cluster_cell_type(
    adata: ad.AnnData,
    cluster_id: str,
    cluster_col: str,
) -> str:
    """Get majority cell type label for cells in this cluster."""
    if "cell_type" not in adata.obs.columns:
        return "Unknown"

    mask = adata.obs[cluster_col] == cluster_id
    if mask.sum() == 0:
        return "Unknown"

    return str(adata.obs.loc[mask, "cell_type"].mode()[0])


def _empty_cluster_entry(
    cluster_id: str,
    adata: ad.AnnData,
    cluster_col: str,
    reason: str = "analysis failed",
) -> dict:
    """Return a safe empty entry for clusters where enrichment failed."""
    n_cells   = int((adata.obs[cluster_col] == cluster_id).sum())
    cell_type = _get_cluster_cell_type(adata, cluster_id, cluster_col)

    return {
        "cluster_id":    str(cluster_id),
        "cell_type":     cell_type,
        "n_cells":       n_cells,
        "top_degs":      [],
        "n_degs_tested": 0,
        "pathways":      [],
        "note":          reason,
    }
