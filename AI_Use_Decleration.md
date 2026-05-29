# AI Use Decleration.md

Core principle: Transparency handled confidently is a professional strength. Discovered
concealment is disqualifying. These materials are designed to present honest disclosure as evidence
of scientific maturity, not as an apology.


## What AI assisted with

Code implementation and boilerplate generation - including FastAPI route scaffolding, Pydantic schema structure, and pytest fixture setup. Docstring drafting and clarification of complex passages. It was built with AI coding assistance of (Anthropic) Claude and validated end-to-end on the 10xPBMC3k dataset. 

Debugging assistance-
Cursor and (OpenAI) Codex were used to accelerate implementation of boilerplate and to assist
in debugging specific issues. 
For specific technical issues, 
1. Repository layout: Restored expected package structure (api/, utils/, pipeline/, dashboard/, tests/) so imports like api.main, utils.io,
   and pipeline.qc resolve correctly in Docker and locally.
2. Pathway persistence (pipeline/pathway.py): Store adata.uns["pathway_results"] as JSON string
   instead of list[dict] to fix H5AD save error: Can't implicitly convert non-string objects to strings.
A Changelog for Cursor has been provided for further disclosure. File name - Cursor Changelog.md

## What I did independently

All biological and architectural decisions. This includes: the choice to store adata.raw before normalisation to
preserve raw counts for GRN and DEG testing; the selection of Wilcoxon rank-sum for differential expression
(non-parametric, appropriate for sparse overdispersed count data); the celltypist majority-voting configuration and
confidence thresholding logic; the async job model to handle datasets exceeding HTTP timeout limits; and the
JSON serialisation workaround for pathway results. The full pipeline order (QC - normalise - cluster - annotate
- pathway), the AnnData data contract connecting this pipeline to downstream RNA velocity and GRN modules,
and all parameter defaults and their biological justification were designed independently.

## Validation

The complete pipeline was run on the 10x PBMC3k public dataset. All outputs - UMAP embeddings, cluster
assignments, cell type annotations, pathway enrichment results - were reviewed. A pathway serialisation bug
encountered during this run was diagnosed and resolved using cursor. All 15 unit tests in tests/test_metadata.py pass. The
processed .h5ad is downloadable and is compatible with scVelo and pySCENIC.


### Pipeline Audit Report

Codex was used as an independent third-party reviewer in the preparation of the audit report for the single-cell RNA-seq pipeline. 
The audit was conducted using the public 10x Genomics PBMC3k dataset as a representative input and included a review of the pipeline's generated outputs, workflow structure, and analytical results. 
The findings and observations were compiled into **Final_scRNA_Pipeline_Audit_Report.docx**. 
All conclusions and recommendations were reviewed and validated by the author.


In addition, the processed results were manually compared against the published PBMC3k reference outputs
to verify the correctness and consistency of the pipeline's analytical results. 
Link to Preprocessing and clustering 3k PBMCs (legacy workflow) - https://scanpy.readthedocs.io/en/latest/tutorials/basics/clustering-2017.html


## Statement

I consider effective and transparent use of AI tools a professional competency in modern research. AI accelerated
implementation; it did not substitute for biological reasoning, scientific judgement, troubleshooting, or interpretation.
I can explain and defend every design decision in this codebase.
