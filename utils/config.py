"""
utils/config.py
---------------
Single source of truth for metadata column names.

WHY THIS EXISTS:
    Without this, every module hardcodes "condition", "timepoint", etc.
    If a user's CSV has "Condition" or "time_point", everything breaks.
    Change column names here ONCE — all modules adapt automatically.

HOW FUTURE MODULES USE THIS:
    RNA velocity  → reads config["timepoint_col"] to order pseudotime
    GRN analysis  → reads config["condition_col"] to split regulon networks
    Trajectory    → reads config["replicate_col"] to average across replicates
"""

# ---------------------------------------------------------------------------
# Metadata column name mapping
# Key   = internal name used throughout the codebase
# Value = actual column name expected in the user's metadata.csv
# ---------------------------------------------------------------------------
METADATA_CONFIG = {
    "cell_id_col":    "cell_id",     # required — must match count matrix index
    "condition_col":  "condition",   # required — e.g. "control", "treated"
    "timepoint_col":  "timepoint",   # optional — e.g. "0h", "24h", "48h"
    "replicate_col":  "replicate",   # optional — e.g. "rep1", "rep2"
    "batch_col":      "batch",       # optional — used for batch correction
    "species_col":    "species",     # optional — e.g. "human", "mouse"
    "treatment_col":  "treatment",   # optional — e.g. "drug_A", "vehicle"
}

# Columns that MUST be present in every metadata.csv upload
REQUIRED_COLS = ["cell_id", "condition"]

# Columns that are optional but will be used if present
OPTIONAL_COLS = ["timepoint", "replicate", "batch", "species", "treatment"]

# All recognized columns (required + optional)
ALL_EXPECTED_COLS = REQUIRED_COLS + OPTIONAL_COLS

# ---------------------------------------------------------------------------
# Pipeline defaults — QC thresholds, clustering params etc.
# These are separate from metadata config intentionally.
# ---------------------------------------------------------------------------
PIPELINE_DEFAULTS = {
    "min_genes":       200,    # cells with fewer genes are filtered
    "max_genes":       6000,   # cells with more genes may be doublets
    "max_mt_pct":      20.0,   # mitochondrial gene % cutoff
    "leiden_resolution": 0.5,  # clustering granularity
    "n_top_genes":     2000,   # highly variable genes
    "n_pcs":           50,
    "n_neighbors":     15,
}
