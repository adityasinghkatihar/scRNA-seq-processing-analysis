"""
dashboard/app.py
----------------
Streamlit dashboard for scRNA-seq analysis.

METADATA INTEGRATION:
    - User uploads counts matrix + metadata.csv
    - Metadata is validated and aligned via utils/metadata.py
    - All filters are built dynamically from adata.obs columns
    - No column names are hardcoded in the UI layer

FUTURE READINESS:
    - This dashboard passes adata to analysis modules (QC, clustering)
    - Downstream tabs (velocity, GRN) will read adata.obs directly
    - Session state holds one adata object shared across all tabs
"""

import warnings
from io import StringIO
from pathlib import Path
import tempfile
import os

import numpy as np
import pandas as pd
import streamlit as st
import anndata as ad

# Internal modules
from utils.config import METADATA_CONFIG, PIPELINE_DEFAULTS
from utils.io import build_anndata, save_anndata
from utils.metadata import get_available_filter_cols, ValidationResult


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="scRNA-seq Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ---------------------------------------------------------------------------
# Session state initialization
# ---------------------------------------------------------------------------
def init_session():
    """Initialize session state keys on first load."""
    defaults = {
        "adata": None,              # the central AnnData object
        "val_result": None,         # ValidationResult from last upload
        "active_filters": {},       # {col: [selected_values]}
        "upload_error": None,
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


# ---------------------------------------------------------------------------
# Upload + build AnnData
# ---------------------------------------------------------------------------
def render_upload_panel():
    """Sidebar panel for file uploads and AnnData construction."""
    st.sidebar.header("Upload Data")

    counts_file = st.sidebar.file_uploader(
        "Counts Matrix",
        type=["h5", "csv", "tsv", "h5ad"],
        help="10x HDF5 (.h5), dense CSV, or AnnData (.h5ad)",
    )
    meta_file = st.sidebar.file_uploader(
        "metadata.csv",
        type=["csv"],
        help="Must contain: cell_id, condition. Optional: timepoint, replicate, batch, species, treatment",
    )

    # Show metadata format hint
    with st.sidebar.expander("metadata.csv format"):
        st.markdown("""
**Required columns:**
- `cell_id` — must match cell barcodes in counts matrix
- `condition` — e.g. control, treated

**Optional columns:**
- `timepoint` — e.g. 0h, 24h, 48h
- `replicate` — e.g. rep1, rep2
- `batch` — e.g. batch_1, batch_2
- `species` — e.g. human, mouse
- `treatment` — e.g. drug_A, vehicle
        """)

    if st.sidebar.button("Load & Validate", type="primary", disabled=(not counts_file or not meta_file)):
        _load_files(counts_file, meta_file)

    # Show validation messages
    if st.session_state.upload_error:
        st.sidebar.error(st.session_state.upload_error)

    if st.session_state.val_result:
        val: ValidationResult = st.session_state.val_result
        if val.warnings:
            for w in val.warnings:
                st.sidebar.warning(w)
        if val.is_valid:
            st.sidebar.success(f"Loaded {st.session_state.adata.n_obs:,} cells × {st.session_state.adata.n_vars:,} genes")


def _load_files(counts_file, meta_file):
    """Write uploads to temp files and call build_anndata."""
    st.session_state.upload_error = None
    st.session_state.val_result   = None

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write counts file
        counts_suffix = Path(counts_file.name).suffix
        counts_path   = os.path.join(tmpdir, f"counts{counts_suffix}")
        with open(counts_path, "wb") as f:
            f.write(counts_file.getvalue())

        # Write metadata file
        meta_path = os.path.join(tmpdir, "metadata.csv")
        with open(meta_path, "wb") as f:
            f.write(meta_file.getvalue())

        try:
            adata, val_result = build_anndata(counts_path, meta_path)
            st.session_state.adata      = adata
            st.session_state.val_result = val_result
            # Reset filters when new data is loaded
            st.session_state.active_filters = {}

        except ValueError as e:
            st.session_state.upload_error = str(e)
        except Exception as e:
            st.session_state.upload_error = f"Unexpected error: {e}"


# ---------------------------------------------------------------------------
# Dynamic metadata filters
# ---------------------------------------------------------------------------
def render_metadata_filters(adata: ad.AnnData) -> ad.AnnData:
    """
    Build filter widgets dynamically from adata.obs columns.

    WHY DYNAMIC:
        We don't know which optional columns the user provided.
        Instead of hardcoding a filter for every possible column,
        we inspect adata.obs and only show filters for present columns.

    Returns filtered adata subset.
    """
    st.sidebar.header("Filter Cells")

    # Get columns that can be used as filters
    filter_cols = get_available_filter_cols(adata.obs)

    if not filter_cols:
        st.sidebar.info("No filterable metadata columns found.")
        return adata

    active_filters = {}

    for col in filter_cols:
        # Get unique values; sort for consistent ordering
        unique_vals = sorted(adata.obs[col].dropna().unique().tolist())

        if len(unique_vals) <= 1:
            continue  # No point showing a filter with one option

        # Use multiselect — all selected by default
        selected = st.sidebar.multiselect(
            label=col.replace("_", " ").title(),
            options=unique_vals,
            default=unique_vals,
            key=f"filter_{col}",
        )

        if selected and set(selected) != set(unique_vals):
            active_filters[col] = selected

    st.session_state.active_filters = active_filters

    # Apply all active filters
    mask = pd.Series([True] * adata.n_obs, index=adata.obs_names)
    for col, selected_vals in active_filters.items():
        mask &= adata.obs[col].isin(selected_vals)

    filtered = adata[mask.values, :].copy()

    # Show filter summary
    if active_filters:
        st.sidebar.caption(
            f"Showing {filtered.n_obs:,} of {adata.n_obs:,} cells after filtering"
        )

    return filtered


# ---------------------------------------------------------------------------
# Main dashboard tabs
# ---------------------------------------------------------------------------
def render_overview_tab(adata: ad.AnnData):
    """Dataset overview — metadata summary table and cell counts."""
    st.header("Dataset Overview")

    col1, col2, col3 = st.columns(3)
    col1.metric("Total Cells", f"{adata.n_obs:,}")
    col2.metric("Genes", f"{adata.n_vars:,}")
    col3.metric("Metadata Columns", len(adata.obs.columns))

    st.subheader("Metadata Summary")

    # Show value counts for each categorical metadata column
    filter_cols = get_available_filter_cols(adata.obs)
    if filter_cols:
        tabs = st.tabs([c.replace("_", " ").title() for c in filter_cols])
        for tab, col in zip(tabs, filter_cols):
            with tab:
                counts = adata.obs[col].value_counts().reset_index()
                counts.columns = [col, "Cell Count"]
                counts["Percentage"] = (counts["Cell Count"] / adata.n_obs * 100).round(1)
                st.dataframe(counts, use_container_width=True)
    else:
        st.info("Upload data with metadata to see summary.")

    # Raw metadata preview
    with st.expander("Raw metadata (first 50 rows)"):
        st.dataframe(adata.obs.head(50), use_container_width=True)


def render_qc_tab(adata: ad.AnnData):
    """
    QC metrics tab.
    
    NOTE: Full QC pipeline (scanpy.pp.calculate_qc_metrics) will be
    called here in pipeline/qc.py. This is the placeholder structure.
    """
    st.header("Quality Control")

    try:
        import scanpy as sc

        # Calculate QC metrics if not already done
        if "n_genes_by_counts" not in adata.obs.columns:
            with st.spinner("Calculating QC metrics..."):
                sc.pp.calculate_qc_metrics(
                    adata,
                    qc_vars=["is_mito"] if "is_mito" in adata.var.columns else [],
                    inplace=True,
                )

        col_a, col_b = st.columns(2)

        with col_a:
            st.subheader("Genes per Cell")
            if "n_genes_by_counts" in adata.obs.columns:
                st.bar_chart(
                    adata.obs["n_genes_by_counts"]
                    .value_counts()
                    .sort_index()
                    .head(100)
                )

        with col_b:
            st.subheader("Counts per Cell")
            if "total_counts" in adata.obs.columns:
                st.bar_chart(
                    adata.obs["total_counts"]
                    .value_counts()
                    .sort_index()
                    .head(100)
                )

        # QC thresholds — pulled from PIPELINE_DEFAULTS (not hardcoded)
        st.subheader("QC Thresholds")
        c1, c2, c3 = st.columns(3)
        c1.metric("Min Genes", PIPELINE_DEFAULTS["min_genes"])
        c2.metric("Max Genes", PIPELINE_DEFAULTS["max_genes"])
        c3.metric("Max MT%",  f"{PIPELINE_DEFAULTS['max_mt_pct']}%")

    except ImportError:
        st.error("scanpy is required for QC. Install: pip install scanpy")
    except Exception as e:
        st.error(f"QC calculation failed: {e}")


def render_metadata_tab(adata: ad.AnnData):
    """
    Detailed metadata inspection tab.
    Shows which columns are present and flags any data quality issues.
    """
    st.header("Metadata Inspector")

    condition_col = METADATA_CONFIG["condition_col"]
    timepoint_col = METADATA_CONFIG["timepoint_col"]

    # Column presence report
    st.subheader("Column Status")

    all_tracked = {
        "cell_id":   "Required — cell identifiers",
        "condition": "Required — experimental condition",
        "timepoint": "Optional — used by RNA velocity",
        "replicate": "Optional — used for batch-aware analysis",
        "batch":     "Optional — used for batch correction",
        "species":   "Optional — organism source",
        "treatment": "Optional — perturbation label",
    }

    rows = []
    for col, description in all_tracked.items():
        if col == "cell_id":
            present = True  # always present as index
        else:
            present = col in adata.obs.columns
        n_missing = adata.obs[col].isna().sum() if present and col != "cell_id" else 0
        rows.append({
            "Column": col,
            "Status": "Present" if present else "Not uploaded",
            "Missing Values": n_missing if present and col != "cell_id" else "—",
            "Description": description,
        })

    status_df = pd.DataFrame(rows)
    st.dataframe(status_df, use_container_width=True, hide_index=True)

    # Cell counts per condition
    if condition_col in adata.obs.columns:
        st.subheader(f"Cells per {condition_col.title()}")
        st.bar_chart(adata.obs[condition_col].value_counts())


# ---------------------------------------------------------------------------
# App entry point
# ---------------------------------------------------------------------------
def main():
    init_session()

    st.title("scRNA-seq Analysis Pipeline")
    st.caption("Upload your counts matrix and metadata to begin.")

    # Sidebar: upload + filters
    render_upload_panel()

    adata = st.session_state.adata

    if adata is None:
        st.info("Upload a counts matrix and metadata.csv to get started.")

        # Show example metadata format
        st.subheader("Expected metadata.csv format")
        example = pd.DataFrame({
            "cell_id":    ["AAACCT", "AAGCCT", "ATGCCT"],
            "condition":  ["control", "treated", "control"],
            "timepoint":  ["0h", "24h", "48h"],
            "replicate":  ["rep1", "rep1", "rep2"],
            "batch":      ["batch_1", "batch_1", "batch_2"],
            "species":    ["human", "human", "human"],
            "treatment":  ["vehicle", "drug_A", "vehicle"],
        })
        st.dataframe(example, use_container_width=True, hide_index=True)
        return

    # Apply metadata filters — returns a filtered adata view
    filtered_adata = render_metadata_filters(adata)

    # Main content tabs
    tab_overview, tab_qc, tab_metadata = st.tabs([
        "Overview",
        "QC Metrics",
        "Metadata Inspector",
    ])

    with tab_overview:
        render_overview_tab(filtered_adata)

    with tab_qc:
        render_qc_tab(filtered_adata)

    with tab_metadata:
        render_metadata_tab(filtered_adata)


if __name__ == "__main__":
    main()
