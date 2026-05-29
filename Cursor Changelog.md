# Cursor Changelog

### Fixed

- **Repository layout**: Restored expected package structure (`api/`, `utils/`, `pipeline/`, `dashboard/`, `tests/`) so imports like `api.main`, `utils.io`, and `pipeline.qc` resolve correctly in Docker and locally.
- **Dockerfile**: Removed accidental Markdown code fences (` ```dockerfile ` / ` ``` `) that caused `docker build` parse errors.
- **docker-compose.yml**: Removed obsolete top-level `version` key (Compose v2 warning).
- **Pathway persistence** (`pipeline/pathway.py`): Store `adata.uns["pathway_results"]` as JSON string instead of `list[dict]` to fix H5AD save error: `Can't implicitly convert non-string objects to strings`.

### Added

- **`scripts/export_run_artifacts.py`**: Export plots, marker tables, pathway tables, QC/cluster/annotation summaries, metadata CSVs, and run logs after a completed job.
- **Package `__init__.py` files** for `api`, `api/routes`, `utils`, `pipeline`, `dashboard`, and `tests`.

### Verified

- `docker compose up --build` — API (8000) and dashboard (8501) start successfully.
- `pytest tests/test_metadata.py` — 15 passed.
- End-to-end PBMC3k run — upload → run → completed; processed `.h5ad` downloadable.
