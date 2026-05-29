"""
tests/test_metadata.py
-----------------------
Unit tests for metadata validation and alignment logic.
Run with: pytest tests/test_metadata.py -v
"""

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import anndata as ad

from utils.metadata import (
    validate_metadata,
    align_metadata,
    fill_missing_optional,
    get_available_filter_cols,
)
from utils.config import METADATA_CONFIG


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def make_meta(include_optional=True, extra_cols=None):
    """Build a valid metadata DataFrame for testing."""
    data = {
        "condition": ["control", "treated", "control", "treated"],
    }
    if include_optional:
        data.update({
            "timepoint":  ["0h", "24h", "0h", "24h"],
            "replicate":  ["rep1", "rep1", "rep2", "rep2"],
            "species":    ["human", "human", "human", "human"],
        })
    if extra_cols:
        data.update(extra_cols)

    df = pd.DataFrame(data, index=["cell_A", "cell_B", "cell_C", "cell_D"])
    df.index.name = "cell_id"
    return df


# ---------------------------------------------------------------------------
# validate_metadata
# ---------------------------------------------------------------------------
class TestValidateMetadata:

    def test_valid_metadata_passes(self):
        meta = make_meta()
        result = validate_metadata(meta)
        assert result.is_valid
        assert len(result.errors) == 0

    def test_missing_required_condition_fails(self):
        meta = make_meta()
        meta = meta.drop(columns=["condition"])
        result = validate_metadata(meta)
        assert not result.is_valid
        assert any("condition" in e for e in result.errors)

    def test_duplicate_cell_ids_fails(self):
        meta = make_meta()
        meta = pd.concat([meta, meta.iloc[[0]]])  # add duplicate row
        result = validate_metadata(meta)
        assert not result.is_valid
        assert any("Duplicate" in e for e in result.errors)

    def test_missing_condition_values_fails(self):
        meta = make_meta()
        meta.loc["cell_A", "condition"] = None
        result = validate_metadata(meta)
        assert not result.is_valid
        assert any("missing values" in e for e in result.errors)

    def test_missing_optional_values_warns(self):
        meta = make_meta()
        meta.loc["cell_A", "timepoint"] = None
        result = validate_metadata(meta)
        assert result.is_valid                            # still valid
        assert any("timepoint" in w for w in result.warnings)  # but warned

    def test_minimal_metadata_valid(self):
        # Only required columns — should still pass
        meta = make_meta(include_optional=False)
        result = validate_metadata(meta)
        assert result.is_valid


# ---------------------------------------------------------------------------
# align_metadata
# ---------------------------------------------------------------------------
class TestAlignMetadata:

    def test_perfect_overlap_preserves_order(self):
        meta = make_meta()
        counts_ids = ["cell_D", "cell_B", "cell_A", "cell_C"]  # different order
        aligned = align_metadata(counts_ids, meta)
        assert list(aligned.index) == counts_ids

    def test_extra_metadata_rows_ignored(self):
        meta = make_meta()
        counts_ids = ["cell_A", "cell_B"]  # subset of metadata
        with pytest.warns(UserWarning, match="no matching cell"):
            aligned = align_metadata(counts_ids, meta)
        assert list(aligned.index) == counts_ids

    def test_extra_counts_cells_dropped_with_warning(self):
        meta = make_meta()
        counts_ids = ["cell_A", "cell_B", "cell_UNKNOWN"]  # unknown cell
        with pytest.warns(UserWarning, match="no metadata entry"):
            aligned = align_metadata(counts_ids, meta)
        assert "cell_UNKNOWN" not in aligned.index

    def test_no_overlap_raises(self):
        meta = make_meta()
        counts_ids = ["barcode_1", "barcode_2"]  # completely different
        with pytest.raises(ValueError, match="No cell IDs overlap"):
            align_metadata(counts_ids, meta)


# ---------------------------------------------------------------------------
# fill_missing_optional
# ---------------------------------------------------------------------------
class TestFillMissingOptional:

    def test_fills_nan_with_unknown(self):
        meta = make_meta()
        meta.loc["cell_A", "timepoint"] = None
        filled = fill_missing_optional(meta)
        assert filled.loc["cell_A", "timepoint"] == "unknown"

    def test_does_not_modify_existing_values(self):
        meta = make_meta()
        filled = fill_missing_optional(meta)
        assert filled.loc["cell_A", "timepoint"] == "0h"

    def test_does_not_mutate_original(self):
        meta = make_meta()
        meta.loc["cell_A", "timepoint"] = None
        fill_missing_optional(meta)           # should not modify meta
        assert pd.isna(meta.loc["cell_A", "timepoint"])


# ---------------------------------------------------------------------------
# get_available_filter_cols
# ---------------------------------------------------------------------------
class TestGetAvailableFilterCols:

    def test_returns_present_cols_only(self):
        meta = make_meta(include_optional=False)  # only condition
        cols = get_available_filter_cols(meta)
        assert "condition" in cols
        assert "timepoint" not in cols

    def test_returns_all_when_all_present(self):
        meta = make_meta(include_optional=True)
        cols = get_available_filter_cols(meta)
        assert "condition" in cols
        assert "timepoint" in cols
        assert "replicate" in cols
