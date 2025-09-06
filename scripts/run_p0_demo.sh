#!/usr/bin/env bash
# -*- coding: ascii -*-
set -euo pipefail
IFS=$'\n\t'

# Ensure conda is available
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found"; exit 1
fi

# Activate conda environment with proper shell hook
eval "$(conda shell.bash hook)"
conda activate halo-p0

# Key dependency preflight check (hard failure)
python - <<'PY'
import sys
missing = []
for m in ("rdkit", "pandas", "pyarrow", "yaml"):
    try:
        __import__(m)
    except Exception:
        missing.append(m)
if missing:
    print("ERROR: Missing Python deps:", ", ".join(missing))
    sys.exit(1)
print("Deps preflight OK")
PY

echo "Starting P0 halogenation demo..."

# 1) Parse names -> parents.smi
echo "Step 1: Ingesting flavonoid names and standardizing..."
python -m halogenator.cli ingest -c configs/p0.yaml
[ -f "data/output/p0/parents.smi" ] || { echo "ERROR: parents.smi not generated"; exit 1; }
echo "  -> parents.smi generated successfully"

# 2) k=1 single substitution (R1-R5 + QC + dedupe)
echo "Step 2: Generating k=1 halogenated products (flavonoids)..."
python -m halogenator.cli k1 -c configs/p0.yaml --subset flavonoids
[ -f "data/output/p0/products_k1.parquet" ] || { echo "ERROR: flavonoid products not generated"; exit 1; }
echo "  -> flavonoid products generated successfully"

echo "Step 2b: Generating k=1 halogenated products (probes)..."
python -m halogenator.cli k1 -c configs/p0.yaml --subset probes
[ -f "data/output/p0/products_k1_probes.parquet" ] || { echo "ERROR: probe products not generated"; exit 1; }
echo "  -> probe products generated successfully"

echo "Step 2c: Generating k=1 halogenated products (all)..."
python -m halogenator.cli k1 -c configs/p0.yaml --subset all
[ -f "data/output/p0/products_k1_all.parquet" ] || { echo "ERROR: combined products not generated"; exit 1; }
echo "  -> combined products generated successfully"

# 3) Generate summary reports
echo "Step 3: Generating summary report (flavonoids)..."
python -m halogenator.cli report -c configs/p0.yaml --subset flavonoids
[ -f "data/output/p0/summary_k1.csv" ] || { echo "ERROR: flavonoid summary not generated"; exit 1; }
[ -f "data/output/p0/summary_k1_pivot.csv" ] || { echo "ERROR: flavonoid pivot not generated"; exit 1; }
echo "  -> flavonoid report generated successfully"

echo "Step 3b: Generating summary report (probes)..."
python -m halogenator.cli report -c configs/p0.yaml --subset probes
[ -f "data/output/p0/summary_k1_probes.csv" ] || { echo "ERROR: probe summary not generated"; exit 1; }
[ -f "data/output/p0/summary_k1_probes_pivot.csv" ] || { echo "ERROR: probe pivot not generated"; exit 1; }
echo "  -> probe report generated successfully"

echo "Step 3c: Generating summary report (all)..."
python -m halogenator.cli report -c configs/p0.yaml --subset all
[ -f "data/output/p0/summary_k1_all.csv" ] || { echo "ERROR: combined summary not generated"; exit 1; }
[ -f "data/output/p0/summary_k1_all_pivot.csv" ] || { echo "ERROR: combined pivot not generated"; exit 1; }
echo "  -> combined report generated successfully"

echo ""
echo "P0 finished. Outputs in data/output/p0/"
echo ""
echo "Compare pivot tables:"
echo "  Flavonoids: data/output/p0/summary_k1_pivot.csv"
echo "  Probes:     data/output/p0/summary_k1_probes_pivot.csv"
echo "  All:        data/output/p0/summary_k1_all_pivot.csv"
echo "              data/output/p0/summary_k1_all_flavonoids_pivot.csv"
echo "              data/output/p0/summary_k1_all_probes_pivot.csv"