"""
api/main.py
-----------
FastAPI application entry point.

JOB MODEL:
    POST /upload         → upload counts + metadata, get job_id
    POST /run/{job_id}   → start pipeline with config, runs async
    GET  /status/{job_id}→ poll progress
    GET  /results/{job_id}→ fetch full results when complete

WHY ASYNC:
    Large datasets (100k+ cells) take minutes to process.
    Synchronous requests would timeout. The client polls /status
    until "completed", then fetches /results once.

JOB STORE:
    In-memory dict for now (single-process). Swap for Redis
    when scaling to multi-worker deployment.

FUTURE:
    This API is the integration point for the velocity and GRN repos.
    Those services will POST to /results/{job_id} to fetch the .h5ad
    and start their own pipelines.
"""

import asyncio
import os
import tempfile
import uuid
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware

from api.routes import upload, run, results
from api.schemas import JobStatus


# ---------------------------------------------------------------------------
# Output directory for saved .h5ad files
# ---------------------------------------------------------------------------
OUTPUT_DIR = Path(os.getenv("OUTPUT_DIR", "/tmp/scrna_outputs"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# In-memory job store
# Structure: { job_id: { "status": ..., "progress": ..., "step": ...,
#                         "error": ..., "results": ..., "adata_path": ... } }
# ---------------------------------------------------------------------------
JOB_STORE: dict = {}


# ---------------------------------------------------------------------------
# App lifecycle
# ---------------------------------------------------------------------------
@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup / shutdown hooks."""
    # Expose job store and output dir via app state
    app.state.jobs       = JOB_STORE
    app.state.output_dir = OUTPUT_DIR
    yield
    # Cleanup temp files on shutdown (optional — comment out to persist)
    # shutil.rmtree(OUTPUT_DIR, ignore_errors=True)


# ---------------------------------------------------------------------------
# FastAPI app
# ---------------------------------------------------------------------------
app = FastAPI(
    title       = "scRNA-seq Pipeline API",
    description = "QC, normalization, clustering, annotation, and pathway analysis for scRNA-seq data.",
    version     = "0.1.0",
    lifespan    = lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins  = ["*"],   # tighten in production
    allow_methods  = ["*"],
    allow_headers  = ["*"],
)

# Register routers
app.include_router(upload.router,  prefix="/upload",  tags=["Upload"])
app.include_router(run.router,     prefix="/run",     tags=["Pipeline"])
app.include_router(results.router, prefix="/results", tags=["Results"])


# ---------------------------------------------------------------------------
# Health check
# ---------------------------------------------------------------------------
@app.get("/health", tags=["Health"])
def health():
    return {
        "status":    "ok",
        "jobs_active": sum(
            1 for j in JOB_STORE.values()
            if j.get("status") == JobStatus.RUNNING
        ),
    }
