"""
api/schemas.py
--------------
Pydantic schemas for all API request and response bodies.

WHY PYDANTIC:
    FastAPI uses these for automatic validation and OpenAPI docs.
    All pipeline parameters have defaults from PIPELINE_DEFAULTS —
    callers only need to supply overrides.
"""

from __future__ import annotations

from enum import Enum
from typing import Dict, List, Optional, Any

from pydantic import BaseModel, Field

from utils.config import PIPELINE_DEFAULTS


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------
class JobStatus(str, Enum):
    PENDING    = "pending"
    RUNNING    = "running"
    COMPLETED  = "completed"
    FAILED     = "failed"


class AnnotationMethod(str, Enum):
    CELLTYPIST   = "celltypist"
    MARKER_SCORE = "marker_score"


class Species(str, Enum):
    HUMAN = "human"
    MOUSE = "mouse"


# ---------------------------------------------------------------------------
# Pipeline run config — sent with every /run request
# All fields optional — defaults come from PIPELINE_DEFAULTS / pipeline modules
# ---------------------------------------------------------------------------
class QCConfig(BaseModel):
    min_genes:  int   = PIPELINE_DEFAULTS["min_genes"]
    max_genes:  int   = PIPELINE_DEFAULTS["max_genes"]
    max_mt_pct: float = PIPELINE_DEFAULTS["max_mt_pct"]
    min_cells:  int   = 3


class NormConfig(BaseModel):
    target_sum:  float = 1e4
    n_top_genes: int   = PIPELINE_DEFAULTS["n_top_genes"]


class ClusterConfig(BaseModel):
    n_pcs:       int   = PIPELINE_DEFAULTS["n_pcs"]
    n_neighbors: int   = PIPELINE_DEFAULTS["n_neighbors"]
    resolution:  float = PIPELINE_DEFAULTS["leiden_resolution"]


class AnnotationConfig(BaseModel):
    method:           AnnotationMethod = AnnotationMethod.CELLTYPIST
    celltypist_model: str              = "Immune_All_Low.pkl"
    min_score:        float            = 0.5
    # Optional custom marker dict: {"T_cell": ["CD3D", "CD3E"], ...}
    custom_markers:   Optional[Dict[str, List[str]]] = None


class PathwayConfig(BaseModel):
    gene_sets:  Optional[List[str]] = None  # None = use defaults for species
    n_top_degs: int                 = 100
    p_cutoff:   float               = 0.05
    species:    Species             = Species.HUMAN


class RunConfig(BaseModel):
    """Full pipeline run configuration. All sub-configs are optional."""
    qc:         QCConfig         = Field(default_factory=QCConfig)
    norm:       NormConfig       = Field(default_factory=NormConfig)
    cluster:    ClusterConfig    = Field(default_factory=ClusterConfig)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)
    pathway:    PathwayConfig    = Field(default_factory=PathwayConfig)

    # Which steps to run (skip steps by setting to False)
    run_qc:         bool = True
    run_norm:       bool = True
    run_cluster:    bool = True
    run_annotation: bool = True
    run_pathway:    bool = True


# ---------------------------------------------------------------------------
# Responses
# ---------------------------------------------------------------------------
class JobResponse(BaseModel):
    job_id:  str
    status:  JobStatus
    message: str = ""


class JobStatusResponse(BaseModel):
    job_id:   str
    status:   JobStatus
    progress: int              # 0–100
    step:     str              # current pipeline step name
    message:  str = ""
    error:    Optional[str] = None


class QCResultResponse(BaseModel):
    cells_before:    int
    cells_after:     int
    genes_before:    int
    genes_after:     int
    cells_removed:   int
    genes_removed:   int
    removal_reasons: Dict[str, int]
    thresholds_used: Dict[str, Any]


class ClusterResultResponse(BaseModel):
    n_clusters:    int
    resolution:    float
    n_pcs:         int
    n_neighbors:   int
    cluster_sizes: Dict[str, int]


class AnnotationResultResponse(BaseModel):
    method:           str
    model_used:       Optional[str]
    cell_type_counts: Dict[str, int]
    n_unassigned:     int


class PipelineResultResponse(BaseModel):
    """Full results returned by /results/{job_id}"""
    job_id:      str
    status:      JobStatus
    qc:          Optional[QCResultResponse]          = None
    clustering:  Optional[ClusterResultResponse]     = None
    annotation:  Optional[AnnotationResultResponse]  = None
    pathway:     Optional[List[dict]]                = None
    output_path: Optional[str]                       = None   # path to saved .h5ad
    summary:     Optional[str]                       = None
