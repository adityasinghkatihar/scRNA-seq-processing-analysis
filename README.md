# scrna-pipeline

scRNA-seq preprocessing and analysis service.
Standalone repo. Outputs `.h5ad` consumed by velocity and GRN repos.

## Repository layout

```
scrna-pipeline/
├── utils/
│   ├── config.py        ← metadata column names + pipeline defaults (edit here)
│   ├── metadata.py      ← load, validate, align metadata
│   └── io.py            ← AnnData builder (counts + metadata → .h5ad)
├── pipeline/
│   ├── qc.py            ← filter cells/genes by QC metrics
│   ├── normalize.py     ← normalize, log1p, HVG selection
│   ├── cluster.py       ← PCA → KNN → UMAP → Leiden
│   ├── annotate.py      ← celltypist / marker gene scoring
│   └── pathway.py       ← per-cluster DEG + Enrichr enrichment → JSON
├── api/
│   ├── main.py          ← FastAPI app
│   ├── schemas.py       ← Pydantic request/response models
│   └── routes/
│       ├── upload.py    ← POST /upload
│       ├── run.py       ← POST /run/{job_id}
│       └── results.py   ← GET /status, /results, /download
├── dashboard/
│   └── app.py           ← Streamlit UI with dynamic metadata filters
├── tests/
│   └── test_metadata.py ← metadata validation + alignment tests
├── scripts/
│   └── export_run_artifacts.py ← export plots, tables, summaries after a run
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```

## Quickstart (local, no Docker)

```bash
pip install -r requirements.txt

# Run API
uvicorn api.main:app --reload --port 8000

# Run dashboard (separate terminal)
streamlit run dashboard/app.py
```

## Quickstart (Docker)

```bash
docker-compose up --build
# API:       http://localhost:8000
# Dashboard: http://localhost:8501
# API docs:  http://localhost:8000/docs
```

## API workflow

```bash
# 1. Upload files
curl -X POST http://localhost:8000/upload \
  -F "counts_file=@counts.h5" \
  -F "metadata_file=@metadata.csv"
# → {"job_id": "abc-123", "status": "pending"}

# 2. Start pipeline (all defaults)
curl -X POST http://localhost:8000/run/abc-123

# 3. Poll status
curl http://localhost:8000/results/status/abc-123
# → {"progress": 70, "step": "clustering"}

# 4. Fetch results
curl http://localhost:8000/results/results/abc-123

# 5. Download .h5ad (input for velocity / GRN repos)
curl -O http://localhost:8000/results/download/abc-123
```

## metadata.csv format

| cell_id   | condition | timepoint | replicate | batch   | species | treatment |
|-----------|-----------|-----------|-----------|---------|---------|-----------|
| AAACCT... | control   | 0h        | rep1      | batch_1 | human   | vehicle   |
| AAGCCT... | treated   | 24h       | rep1      | batch_1 | human   | drug_A    |

- `cell_id` and `condition` are **required**
- All other columns are **optional**

## Run tests

```bash
pytest tests/ -v
```

## AnnData contract for downstream repos

The output `.h5ad` contains:

| Location | Content | Used by |
|---|---|---|
| `adata.X` | normalized + log1p counts (HVGs) | clustering, visualization |
| `adata.raw` | raw counts (all genes) | pySCENIC, DEG testing |
| `adata.obs` | all metadata columns + QC metrics + leiden + cell_type | velocity, GRN, trajectory |
| `adata.obsm["X_pca"]` | PCA embedding | scVelo |
| `adata.obsm["X_umap"]` | UMAP embedding | visualization |
| `adata.obsp["connectivities"]` | KNN graph | RNA velocity transition matrix |
| `adata.uns["pathway_results"]` | pathway JSON | interpretation layer |
