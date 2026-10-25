#!/usr/bin/env bash
# Full pipeline, from download to replication table and figures.
# Needs: python (requirements.txt) as $PY or python on PATH; Rscript with edgeR, limma, sva, isa2
# (R/install.R) as $RSCRIPT or on PATH. Optional: COS_DATA_DIR / COS_WORK_DIR for large files.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PY:-python}"
"$PY" scripts/fetch_data.py
"$PY" scripts/01_prepare.py
"$PY" scripts/02_published_signatures.py
"$PY" scripts/03_denovo_dog.py          # canine only; freezes results/frozen_signatures.json
"$PY" scripts/04_transfer_target.py     # verifies the frozen hash before reading TARGET
"$PY" scripts/05_immune_pathways.py
"$PY" scripts/06_sensitivity.py
"$PY" scripts/08_replication_table.py
"$PY" scripts/07_figures.py
"$PY" -m pytest -q
