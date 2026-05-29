"""
pipeline/cluster.py
-------------------
Dimensionality reduction, nearest-neighbor graph, and clustering.

PIPELINE ORDER:
    HVG-normalized data → PCA → KNN graph → UMAP + Leiden clustering

WHAT EACH STEP DOES:
    PCA       : reduces ~2000 HVGs to 50 principal components.
                Linear compression. Fast. Removes technical noise.

    KNN graph : builds a graph where each cell connects to its
                k nearest neighbors in PCA space. This graph is
                the foundation for both UMAP and Leiden.

    UMAP      : non-linear 2D embedding for visualization only.
                Do NOT use UMAP coordinates for analysis — only for plots.

    Leiden    : community detection on the KNN graph → cluster labels.
                Resolution controls granularity: higher = more clusters.

FUTURE MODULE NOTE:
    adata.obsm["X_pca"]   — used by scVelo for velocity graph
    adata.obsm["X_umap"]  — used by dashboard and CellOracle visualization
    adata.obs["leiden"]   — used by pySCENIC to split per-cluster regulons
    adata.obsp["connectivities"] — used by RNA velocity transition matrix
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
class ClusterResult:
    n_clusters:   int
    resolution:   float
    n_pcs:        int
    n_neighbors:  int
    cluster_sizes: dict   # {"0": 145, "1": 302, ...}

    def summary(self) -> str:
        return (
            f"Found {self.n_clusters} clusters at resolution {self.resolution}.\n"
            f"PCA components: {self.n_pcs}, Neighbors: {self.n_neighbors}.\n"
            f"Cluster sizes: { {k: v for k, v in list(self.cluster_sizes.items())[:5]} }"
            + (" ..." if len(self.cluster_sizes) > 5 else "")
        )


# ---------------------------------------------------------------------------
# Main clustering function
# ---------------------------------------------------------------------------
def run_clustering(
    adata: ad.AnnData,
    n_pcs:       Optional[int]   = None,
    n_neighbors: Optional[int]   = None,
    resolution:  Optional[float] = None,
    random_state: int            = 42,
) -> tuple[ad.AnnData, ClusterResult]:
    """
    Run PCA → KNN → UMAP → Leiden clustering on normalized AnnData.

    Args:
        adata        : normalized AnnData (output of pipeline/normalize.py)
        n_pcs        : number of PCA components (default from PIPELINE_DEFAULTS)
        n_neighbors  : number of nearest neighbors for KNN graph
        resolution   : Leiden clustering resolution (higher = more clusters)
        random_state : random seed for reproducibility

    Returns:
        (adata, ClusterResult)
        adata.obsm["X_pca"]  → PCA embedding
        adata.obsm["X_umap"] → UMAP embedding
        adata.obs["leiden"]  → cluster labels (string: "0", "1", ...)
    """
    n_pcs       = n_pcs       or PIPELINE_DEFAULTS["n_pcs"]
    n_neighbors = n_neighbors or PIPELINE_DEFAULTS["n_neighbors"]
    resolution  = resolution  or PIPELINE_DEFAULTS["leiden_resolution"]

    # 1. Scale gene expression ─────────────────────────────────────────────
    #
    # Zero-centers and scales each gene to unit variance.
    # This prevents high-expression genes from dominating PCA.
    # We clip at max_value=10 to reduce the influence of extreme outliers.
    #
    sc.pp.scale(adata, max_value=10)

    # 2. PCA ───────────────────────────────────────────────────────────────
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=random_state)

    # 3. KNN graph ─────────────────────────────────────────────────────────
    #
    # Builds adata.obsp["connectivities"] and adata.obsp["distances"].
    # n_pcs controls how many PC dimensions to use for distance calculation.
    #
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)

    # 4. UMAP ──────────────────────────────────────────────────────────────
    sc.tl.umap(adata, random_state=random_state)

    # 5. Leiden clustering ─────────────────────────────────────────────────
    #
    # Stores cluster labels in adata.obs["leiden"]
    # Labels are strings: "0", "1", "2", ...
    #
    sc.tl.leiden(adata, resolution=resolution, random_state=random_state)

    # 6. Record clustering state ───────────────────────────────────────────
    cluster_sizes = adata.obs["leiden"].value_counts().to_dict()
    n_clusters    = len(cluster_sizes)

    adata.uns["clustering"] = {
        "completed":   True,
        "method":      "leiden",
        "resolution":  resolution,
        "n_clusters":  n_clusters,
        "n_pcs":       n_pcs,
        "n_neighbors": n_neighbors,
    }

    result = ClusterResult(
        n_clusters    = n_clusters,
        resolution    = resolution,
        n_pcs         = n_pcs,
        n_neighbors   = n_neighbors,
        cluster_sizes = {k: int(v) for k, v in cluster_sizes.items()},
    )

    return adata, result
