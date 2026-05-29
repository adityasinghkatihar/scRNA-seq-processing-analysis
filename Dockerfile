# scrna-pipeline/Dockerfile
# Multi-stage build: keeps final image lean by separating build deps

# ── Stage 1: dependency builder ──────────────────────────────────────────
FROM python:3.11-slim AS builder

WORKDIR /build

# System deps needed to compile scipy / anndata wheels
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    libhdf5-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

# Install Python dependencies with extended timeout/retries
RUN pip install \
    --default-timeout=1000 \
    --retries=10 \
    --prefix=/install \
    --no-cache-dir \
    -r requirements.txt


# ── Stage 2: runtime image ────────────────────────────────────────────────
FROM python:3.11-slim AS runtime

WORKDIR /app

# Only runtime system libs (HDF5 for anndata)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    hdf5-tools \
    && rm -rf /var/lib/apt/lists/*

# Copy installed packages from builder
COPY --from=builder /install /usr/local

# Copy application code
COPY . .

# Create output directory (bind-mount in production)
RUN mkdir -p /tmp/scrna_outputs

# ── Environment ───────────────────────────────────────────────────────────
ENV OUTPUT_DIR=/tmp/scrna_outputs
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Expose both FastAPI and Streamlit ports
EXPOSE 8000
EXPOSE 8501

# ── Default: run API ──────────────────────────────────────────────────────
# Override CMD to run the dashboard instead:
# docker run ... streamlit run dashboard/app.py --server.port 8501
CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]
