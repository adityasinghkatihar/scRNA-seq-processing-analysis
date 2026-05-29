"""
api/routes/upload.py
--------------------
Handles file uploads and initial AnnData construction.

POST /upload
    - Accepts counts matrix + metadata.csv
    - Validates metadata immediately (fast, synchronous)
    - Returns job_id for use in subsequent /run call
"""

import os
import uuid
from pathlib import Path

from fastapi import APIRouter, File, HTTPException, Request, UploadFile

from api.schemas import JobResponse, JobStatus
from utils.io import build_anndata
from utils.io import save_anndata

router = APIRouter()

# Allowed file extensions for counts matrix
ALLOWED_COUNTS_EXTENSIONS = {".h5", ".csv", ".tsv", ".h5ad", ".loom"}


@router.post("", response_model=JobResponse)
async def upload_files(
    request:      Request,
    counts_file:  UploadFile = File(..., description="Count matrix (.h5, .csv, .tsv, .h5ad)"),
    metadata_file: UploadFile = File(..., description="Metadata CSV (cell_id, condition, ...)"),
):
    """
    Upload counts matrix and metadata CSV.

    Creates a job entry and saves validated AnnData to disk.
    Returns job_id for use with POST /run/{job_id}.
    """
    jobs       = request.app.state.jobs
    output_dir = request.app.state.output_dir

    # Validate file extensions
    counts_ext = Path(counts_file.filename).suffix.lower()
    if counts_ext not in ALLOWED_COUNTS_EXTENSIONS:
        raise HTTPException(
            status_code=422,
            detail=f"Unsupported counts format '{counts_ext}'. "
                   f"Allowed: {sorted(ALLOWED_COUNTS_EXTENSIONS)}"
        )

    if not metadata_file.filename.endswith(".csv"):
        raise HTTPException(
            status_code=422,
            detail="metadata_file must be a .csv file."
        )

    job_id = str(uuid.uuid4())

    # Save uploaded files to a temp location
    job_dir = output_dir / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    counts_path = job_dir / f"counts{counts_ext}"
    meta_path   = job_dir / "metadata.csv"

    counts_path.write_bytes(await counts_file.read())
    meta_path.write_bytes(await metadata_file.read())

    # Build and validate AnnData immediately
    # This is fast (just IO + metadata validation) — fine to do synchronously
    try:
        adata, val_result = build_anndata(str(counts_path), str(meta_path))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to load data: {e}")

    # Save validated AnnData for the pipeline to pick up
    adata_path = job_dir / "input.h5ad"
    save_anndata(adata, str(adata_path))

    # Initialize job record
    jobs[job_id] = {
        "status":     JobStatus.PENDING,
        "progress":   0,
        "step":       "uploaded",
        "error":      None,
        "results":    {},
        "adata_path": str(adata_path),
        "warnings":   val_result.warnings,
        "cells":      adata.n_obs,
        "genes":      adata.n_vars,
    }

    return JobResponse(
        job_id  = job_id,
        status  = JobStatus.PENDING,
        message = (
            f"Uploaded {adata.n_obs:,} cells × {adata.n_vars:,} genes. "
            f"Warnings: {len(val_result.warnings)}. "
            f"Call POST /run/{job_id} to start pipeline."
        ),
    )
