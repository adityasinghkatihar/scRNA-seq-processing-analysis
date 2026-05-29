"""
api/routes/results.py
---------------------
Status polling and result retrieval endpoints.

GET /status/{job_id}  → lightweight poll (progress %, current step)
GET /results/{job_id} → full results once completed
GET /download/{job_id}→ download the processed .h5ad file
"""

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import FileResponse

from api.schemas import (
    JobStatus,
    JobStatusResponse,
    PipelineResultResponse,
    QCResultResponse,
    ClusterResultResponse,
    AnnotationResultResponse,
)

router = APIRouter()


# ---------------------------------------------------------------------------
# Status poll — lightweight, called frequently
# ---------------------------------------------------------------------------
@router.get("/status/{job_id}", response_model=JobStatusResponse)
def get_status(job_id: str, request: Request):
    """
    Poll pipeline progress.
    Call this every 2–5 seconds while status is 'running'.
    """
    jobs = request.app.state.jobs
    job  = _get_job_or_404(jobs, job_id)

    return JobStatusResponse(
        job_id   = job_id,
        status   = job["status"],
        progress = job["progress"],
        step     = job["step"],
        error    = job.get("error"),
    )


# ---------------------------------------------------------------------------
# Full results — called once when status == "completed"
# ---------------------------------------------------------------------------
@router.get("/results/{job_id}", response_model=PipelineResultResponse)
def get_results(job_id: str, request: Request):
    """
    Retrieve full pipeline results.
    Only meaningful after status == 'completed'.
    """
    jobs = request.app.state.jobs
    job  = _get_job_or_404(jobs, job_id)

    if job["status"] == JobStatus.RUNNING:
        raise HTTPException(
            status_code=202,
            detail=f"Pipeline still running ({job['progress']}% — {job['step']}). Poll /status/{job_id}."
        )

    if job["status"] == JobStatus.FAILED:
        raise HTTPException(
            status_code=500,
            detail=f"Pipeline failed: {job.get('error', 'unknown error')}"
        )

    if job["status"] == JobStatus.PENDING:
        raise HTTPException(
            status_code=400,
            detail=f"Pipeline not started. Call POST /run/{job_id} first."
        )

    # Serialize result objects
    pipeline_results = job.get("results", {})

    qc_resp = None
    if "qc" in pipeline_results:
        qc = pipeline_results["qc"]
        qc_resp = QCResultResponse(
            cells_before    = qc.cells_before,
            cells_after     = qc.cells_after,
            genes_before    = qc.genes_before,
            genes_after     = qc.genes_after,
            cells_removed   = qc.cells_removed,
            genes_removed   = qc.genes_removed,
            removal_reasons = qc.removal_reasons,
            thresholds_used = qc.thresholds_used,
        )

    cluster_resp = None
    if "cluster" in pipeline_results:
        cl = pipeline_results["cluster"]
        cluster_resp = ClusterResultResponse(
            n_clusters    = cl.n_clusters,
            resolution    = cl.resolution,
            n_pcs         = cl.n_pcs,
            n_neighbors   = cl.n_neighbors,
            cluster_sizes = cl.cluster_sizes,
        )

    ann_resp = None
    if "annotation" in pipeline_results:
        an = pipeline_results["annotation"]
        ann_resp = AnnotationResultResponse(
            method           = an.method,
            model_used       = an.model_used,
            cell_type_counts = an.cell_type_counts,
            n_unassigned     = an.n_unassigned,
        )

    pathway_resp = None
    if "pathway" in pipeline_results:
        pathway_resp = pipeline_results["pathway"].to_json_ready()

    # Build plain-text summary
    summary_parts = []
    if qc_resp:
        summary_parts.append(pipeline_results["qc"].summary())
    if cluster_resp:
        summary_parts.append(pipeline_results["cluster"].summary())
    if ann_resp:
        summary_parts.append(pipeline_results["annotation"].summary())
    if "pathway" in pipeline_results:
        summary_parts.append(pipeline_results["pathway"].summary())

    return PipelineResultResponse(
        job_id      = job_id,
        status      = job["status"],
        qc          = qc_resp,
        clustering  = cluster_resp,
        annotation  = ann_resp,
        pathway     = pathway_resp,
        output_path = job.get("output_path"),
        summary     = "\n\n".join(summary_parts),
    )


# ---------------------------------------------------------------------------
# Download the processed .h5ad
# ---------------------------------------------------------------------------
@router.get("/download/{job_id}")
def download_h5ad(job_id: str, request: Request):
    """
    Download the fully processed AnnData (.h5ad) file.
    This file can be passed directly to the velocity or GRN repos.
    """
    jobs = request.app.state.jobs
    job  = _get_job_or_404(jobs, job_id)

    if job["status"] != JobStatus.COMPLETED:
        raise HTTPException(
            status_code=400,
            detail="Pipeline not completed. Check /status/{job_id}."
        )

    output_path = job.get("output_path")
    if not output_path:
        raise HTTPException(status_code=404, detail="Output file not found.")

    return FileResponse(
        path             = output_path,
        media_type       = "application/octet-stream",
        filename         = f"scrna_{job_id[:8]}_processed.h5ad",
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _get_job_or_404(jobs: dict, job_id: str) -> dict:
    if job_id not in jobs:
        raise HTTPException(
            status_code=404,
            detail=f"Job '{job_id}' not found. Upload files first via POST /upload."
        )
    return jobs[job_id]
