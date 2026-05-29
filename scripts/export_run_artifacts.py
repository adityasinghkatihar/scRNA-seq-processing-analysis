"""
Export full pipeline run artifacts: plots, tables, summaries, logs.
Usage: python scripts/export_run_artifacts.py <job_id> <output_dir>
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

import anndata as ad


def _ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def _write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=str), encoding="utf-8")


def _export_qc_summary(adata: ad.AnnData, out_dir: Path) -> dict:
    qc = adata.uns.get("qc", {})
    summary = {
        "cells_before": qc.get("cells_before"),
        "cells_after": qc.get("cells_after"),
        "genes_before": qc.get("genes_before"),
        "genes_after": qc.get("genes_after"),
        "thresholds": qc.get("thresholds", {}),
    }
    if summary["cells_before"] and summary["cells_after"] is not None:
        summary["cells_removed"] = summary["cells_before"] - summary["cells_after"]
    if summary["genes_before"] and summary["genes_after"] is not None:
        summary["genes_removed"] = summary["genes_before"] - summary["genes_after"]

    _write_json(out_dir / "qc_summary.json", summary)

    obs_cols = [c for c in ["n_genes_by_counts", "total_counts", "pct_counts_is_mito"] if c in adata.obs.columns]
    if obs_cols:
        adata.obs[obs_cols].describe().to_csv(out_dir / "qc_metrics_summary.csv")

    return summary


def _export_cluster_summary(adata: ad.AnnData, out_dir: Path) -> dict:
    clustering = adata.uns.get("clustering", {})
    sizes = adata.obs["leiden"].value_counts().sort_index().to_dict() if "leiden" in adata.obs else {}
    summary = {
        **clustering,
        "cluster_sizes": {str(k): int(v) for k, v in sizes.items()},
    }
    _write_json(out_dir / "cluster_summary.json", summary)
    if sizes:
        pd.DataFrame({"cluster": list(sizes.keys()), "n_cells": list(sizes.values())}).to_csv(
            out_dir / "cluster_sizes.csv", index=False
        )
    return summary


def _export_annotation_summary(adata: ad.AnnData, out_dir: Path) -> dict:
    ann = adata.uns.get("annotation", {})
    counts = adata.obs["cell_type"].value_counts().to_dict() if "cell_type" in adata.obs else {}
    summary = {
        **ann,
        "cell_type_counts": {str(k): int(v) for k, v in counts.items()},
        "n_unassigned": int(counts.get("Unassigned", 0)),
    }
    _write_json(out_dir / "annotation_summary.json", summary)
    if counts:
        pd.DataFrame({"cell_type": list(counts.keys()), "n_cells": list(counts.values())}).to_csv(
            out_dir / "cell_type_counts.csv", index=False
        )
    return summary


def _export_marker_tables(adata: ad.AnnData, out_dir: Path) -> None:
    markers_dir = _ensure_dir(out_dir / "marker_tables")
    if "rank_genes_groups" not in adata.uns:
        return

    groups = adata.uns["rank_genes_groups"]["names"].dtype.names
    all_rows = []
    for group in groups:
        df = sc.get.rank_genes_groups_df(adata, group=str(group))
        df.insert(0, "cluster", group)
        df.to_csv(markers_dir / f"markers_cluster_{group}.csv", index=False)
        all_rows.append(df.head(50))

    if all_rows:
        pd.concat(all_rows, ignore_index=True).to_csv(markers_dir / "markers_top50_all_clusters.csv", index=False)


def _export_pathway_tables(adata: ad.AnnData, out_dir: Path) -> None:
    pathway_dir = _ensure_dir(out_dir / "pathway")
    raw = adata.uns.get("pathway_results")
    if raw is None:
        return

    results = json.loads(raw) if isinstance(raw, str) else raw
    _write_json(pathway_dir / "pathway_results.json", results)

    rows = []
    for cluster in results:
        for pw in cluster.get("pathways", []):
            rows.append({
                "cluster_id": cluster.get("cluster_id"),
                "cell_type": cluster.get("cell_type"),
                "term": pw.get("term"),
                "gene_set": pw.get("gene_set"),
                "p_value": pw.get("p_value"),
                "adjusted_p": pw.get("adjusted_p"),
                "overlap_genes": ";".join(pw.get("overlap_genes", [])),
            })
    if rows:
        pd.DataFrame(rows).to_csv(pathway_dir / "pathway_enrichment.csv", index=False)


def _export_plots(adata: ad.AnnData, plots_dir: Path) -> None:
    sc.settings.figdir = str(plots_dir)
    sc.settings.set_figure_params(dpi=120, facecolor="white", frameon=False)

    # QC plots
    qc_cols = [c for c in ["n_genes_by_counts", "total_counts", "pct_counts_is_mito"] if c in adata.obs.columns]
    if qc_cols:
        sc.pl.violin(adata, keys=qc_cols, jitter=0.4, multi_panel=True, show=False, save="_qc_violin.png")
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_is_mito" if "pct_counts_is_mito" in adata.obs else None, show=False, save="_qc_scatter_counts_genes.png")

    # Embeddings
    if "X_umap" in adata.obsm:
        color_keys = [k for k in ["leiden", "cell_type", "condition"] if k in adata.obs.columns]
        if color_keys:
            sc.pl.umap(adata, color=color_keys, wspace=0.4, show=False, save="_umap_overview.png")
        if "leiden" in adata.obs.columns:
            sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False, save="_umap_leiden.png")
        if "cell_type" in adata.obs.columns:
            sc.pl.umap(adata, color="cell_type", legend_loc="right margin", show=False, save="_umap_cell_type.png")

    # Marker dotplot for top clusters
    if "rank_genes_groups" in adata.uns and "leiden" in adata.obs.columns:
        try:
            sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False, save="_markers_dotplot.png")
            sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, show=False, save="_markers_heatmap.png")
        except Exception:
            pass


def _export_metadata_tables(adata: ad.AnnData, out_dir: Path) -> None:
    meta_dir = _ensure_dir(out_dir / "metadata")
    adata.obs.to_csv(meta_dir / "cell_metadata.csv")
    adata.var.to_csv(meta_dir / "gene_metadata.csv")


def _write_run_log(out_dir: Path, job_id: str, job_dir: Path, adata: ad.AnnData, summaries: dict) -> None:
    log_dir = _ensure_dir(out_dir / "logs")
    lines = [
        "Pipeline Run Export Log",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        f"Job ID: {job_id}",
        f"Job directory: {job_dir}",
        "",
        "=== QC ===",
        json.dumps(summaries.get("qc", {}), indent=2),
        "",
        "=== Clustering ===",
        json.dumps(summaries.get("cluster", {}), indent=2),
        "",
        "=== Annotation ===",
        json.dumps(summaries.get("annotation", {}), indent=2),
        "",
        "=== Normalization ===",
        json.dumps(adata.uns.get("normalization", {}), indent=2),
    ]
    (log_dir / "pipeline_run.log").write_text("\n".join(lines), encoding="utf-8")


def export(job_id: str, output_dir: str, base_outputs: str = "/tmp/scrna_outputs") -> Path:
    job_dir = Path(base_outputs) / job_id
    h5ad_path = job_dir / "output.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Missing processed h5ad: {h5ad_path}")

    out = _ensure_dir(Path(output_dir))
    for sub in ["plots", "marker_tables", "pathway", "metadata", "logs", "summaries", "data"]:
        _ensure_dir(out / sub)

    adata = ad.read_h5ad(h5ad_path)

    summaries = {
        "qc": _export_qc_summary(adata, out / "summaries"),
        "cluster": _export_cluster_summary(adata, out / "summaries"),
        "annotation": _export_annotation_summary(adata, out / "summaries"),
    }

    norm = adata.uns.get("normalization", {})
    _write_json(out / "summaries" / "normalization_summary.json", norm)

    _export_marker_tables(adata, out)
    _export_pathway_tables(adata, out)
    _export_plots(adata, out / "plots")
    _export_metadata_tables(adata, out)

    # Copy core data files
    import shutil
    for name in ["output.h5ad", "input.h5ad", "metadata.csv"]:
        src = job_dir / name
        if src.exists():
            shutil.copy2(src, out / "data" / name)

    # API results snapshot if available
    try:
        import requests
        resp = requests.get(f"http://localhost:8000/results/results/{job_id}", timeout=60)
        if resp.ok:
            _write_json(out / "summaries" / "api_results.json", resp.json())
    except Exception as e:
        (out / "logs" / "api_fetch_warning.txt").write_text(str(e), encoding="utf-8")

    _write_run_log(out, job_id, job_dir, adata, summaries)

    manifest = {
        "job_id": job_id,
        "exported_at": datetime.now(timezone.utc).isoformat(),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "artifacts": {
            "plots": sorted(p.name for p in (out / "plots").glob("*") if p.is_file()),
            "summaries": sorted(p.name for p in (out / "summaries").glob("*")),
            "marker_tables": sorted(p.name for p in (out / "marker_tables").glob("*")),
            "pathway": sorted(p.name for p in (out / "pathway").glob("*")),
            "metadata": sorted(p.name for p in (out / "metadata").glob("*")),
            "data": sorted(p.name for p in (out / "data").glob("*")),
            "logs": sorted(p.name for p in (out / "logs").glob("*")),
        },
    }
    _write_json(out / "MANIFEST.json", manifest)
    return out


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python scripts/export_run_artifacts.py <job_id> <output_dir>")
        sys.exit(1)
    result = export(sys.argv[1], sys.argv[2])
    print(f"Exported to {result}")
