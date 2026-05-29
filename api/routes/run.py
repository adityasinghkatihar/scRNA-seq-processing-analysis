"""
api/routes/run.py
-----------------
Triggers pipeline execution for an uploaded dataset.

POST /run/{job_id}
    - Accepts RunConfig (all optional — has defaults)
    - Launches pipeline as async background task
    - Returns immediately with job_id
    - Client polls GET /status/{job_id} for progress
"""

import asyncio
from pathlib import Path

from fastapi import APIRouter, BackgroundTasks, HTTPException, Request

from api.schemas import JobResponse, JobStatus, RunConfig
from utils.io import load_anndata, save_anndata

router = APIRouter()


# ---------------------------------------------------------------------------
# Route
# ---------------------------------------------------------------------------
@router.post("/{job_id}", response_model=JobResponse)
async def run_pipeline(
    job_id:           str,
    request:          Request,
    background_tasks: BackgroundTasks,
    config:           RunConfig = None,
):
    """
    Start the analysis pipeline for an uploaded dataset.

    Body (all optional):
        qc, norm, cluster, annotation, pathway configs.
        Set run_* flags to False to skip individual steps.
    """
    jobs = request.app.state.jobs

    if job_id not in jobs:
        raise HTTPException(status_code=404, detail=f"Job '{job_id}' not found.")

    job = jobs[job_id]

    if job["status"] == JobStatus.RUNNING:
        raise HTTPException(status_code=409, detail="Pipeline already running for this job.")

    if job["status"] == JobStatus.COMPLETED:
        raise HTTPException(status_code=409, detail="Pipeline already completed. Fetch results.")

    config = config or RunConfig()

    # Mark as running
    job["status"]   = JobStatus.RUNNING
    job["progress"] = 0
    job["step"]     = "starting"
    job["error"]    = None

    # Launch as background task — returns immediately to caller
    background_tasks.add_task(_run_pipeline_task, job_id, jobs, config)

    return JobResponse(
        job_id  = job_id,
        status  = JobStatus.RUNNING,
        message = f"Pipeline started. Poll GET /status/{job_id} for progress.",
    )


# ---------------------------------------------------------------------------
# Background pipeline task
# ---------------------------------------------------------------------------
async def _run_pipeline_task(job_id: str, jobs: dict, config: RunConfig):
    """
    Runs the full pipeline sequentially in a background asyncio task.
    Updates job["progress"] and job["step"] at each stage.
    Errors are caught and stored in job["error"] — never crash the server.
    """
    job = jobs[job_id]

    # Import pipeline modules here to keep startup fast
    from pipeline.qc         import run_qc
    from pipeline.normalize  import run_normalization
    from pipeline.cluster    import run_clustering
    from pipeline.annotate   import run_annotation
    from pipeline.pathway    import run_pathway_analysis

    try:
        adata = load_anndata(job["adata_path"])

        # Step 1: QC ───────────────────────────────────────────────────────
        if config.run_qc:
            _update(job, step="qc", progress=10)
            # Run in thread pool to avoid blocking the event loop
            adata, qc_result = await asyncio.to_thread(
                run_qc, adata,
                min_genes  = config.qc.min_genes,
                max_genes  = config.qc.max_genes,
                max_mt_pct = config.qc.max_mt_pct,
                min_cells  = config.qc.min_cells,
            )
            job["results"]["qc"] = qc_result
            _update(job, step="qc_done", progress=25)

        # Step 2: Normalize ────────────────────────────────────────────────
        if config.run_norm:
            _update(job, step="normalization", progress=30)
            adata, norm_result = await asyncio.to_thread(
                run_normalization, adata,
                target_sum  = config.norm.target_sum,
                n_top_genes = config.norm.n_top_genes,
            )
            job["results"]["norm"] = norm_result
            _update(job, step="norm_done", progress=45)

        # Step 3: Cluster ──────────────────────────────────────────────────
        if config.run_cluster:
            _update(job, step="clustering", progress=50)
            adata, cluster_result = await asyncio.to_thread(
                run_clustering, adata,
                n_pcs       = config.cluster.n_pcs,
                n_neighbors = config.cluster.n_neighbors,
                resolution  = config.cluster.resolution,
            )
            job["results"]["cluster"] = cluster_result
            _update(job, step="cluster_done", progress=70)

        # Step 4: Annotate ─────────────────────────────────────────────────
        if config.run_annotation:
            _update(job, step="annotation", progress=75)
            ct_model = (
                config.annotation.celltypist_model
                if config.annotation.method.value == "celltypist"
                else None
            )
            adata, ann_result = await asyncio.to_thread(
                run_annotation, adata,
                celltypist_model = ct_model,
                marker_dict      = config.annotation.custom_markers,
                min_score        = config.annotation.min_score,
            )
            job["results"]["annotation"] = ann_result
            _update(job, step="annotation_done", progress=85)

        # Step 5: Pathway ──────────────────────────────────────────────────
        if config.run_pathway:
            _update(job, step="pathway_analysis", progress=88)
            adata, pw_result = await asyncio.to_thread(
                run_pathway_analysis, adata,
                gene_sets  = config.pathway.gene_sets,
                n_top_degs = config.pathway.n_top_degs,
                p_cutoff   = config.pathway.p_cutoff,
                species    = config.pathway.species.value,
            )
            job["results"]["pathway"] = pw_result
            _update(job, step="pathway_done", progress=95)

        # Save final AnnData ───────────────────────────────────────────────
        output_path = str(Path(job["adata_path"]).parent / "output.h5ad")
        await asyncio.to_thread(save_anndata, adata, output_path)
        job["output_path"] = output_path

        _update(job, step="completed", progress=100, status=JobStatus.COMPLETED)

    except Exception as e:
        job["status"]   = JobStatus.FAILED
        job["error"]    = str(e)
        job["step"]     = "failed"
        job["progress"] = 0


def _update(job: dict, step: str, progress: int, status: JobStatus = JobStatus.RUNNING):
    job["step"]     = step
    job["progress"] = progress
    job["status"]   = status
