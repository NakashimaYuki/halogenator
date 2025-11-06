#!/bin/bash
# ETCM 2000 Flavonoids Batch Processing Script
# Complete workflow from Excel ingestion to halogenated product generation

set -e  # Exit on any error

# Configuration - ADJUST THESE FOR YOUR SETUP
ETCM_INPUT="flavonoid_smiles_data_ETCM.xlsx"
K2_CONFIG="configs/etcm_flavonoids.yaml"
K3_CONFIG="configs/etcm_flavonoids_k3.yaml"
OUTPUT_BASE="out/etcm2000"
PARENTS_DIR="out/etcm2000/parents"
P1_K2_DIR="out/etcm2000/p1_k2"
P1_K3_DIR="out/etcm2000/p1_k3"
AUDIT_DIR="out/etcm2000/audit"
SAMPLE_SIZE=0  # Set >0 for testing (e.g., 50), 0 for full batch

# QA conflict tolerance settings
K2_TOL=1  # k=2 conflict tolerance (recommended: 1)
K3_TOL=2  # k=3 conflict tolerance (recommended: 2, more relaxed)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

echo_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

echo_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check prerequisites
echo_info "Checking prerequisites..."

if [ ! -f "$ETCM_INPUT" ]; then
    echo_error "ETCM input file not found: $ETCM_INPUT"
    echo "Please ensure the ETCM Excel file is in the project root directory."
    exit 1
fi

# Removed old CONFIG check - now using K2_CONFIG and K3_CONFIG

# Validate prerequisites
echo_info "Validating prerequisites..."

if [ ! -f "$ETCM_INPUT" ]; then
    echo_error "ETCM input file not found: $ETCM_INPUT"
    echo "Please ensure the ETCM Excel file is in the project root directory."
    exit 1
fi

if [ ! -f "$K2_CONFIG" ]; then
    echo_error "k=2 configuration not found: $K2_CONFIG"
    exit 1
fi

# Check if halogenator CLI is available
if ! python -m halogenator.cli --help &> /dev/null; then
    echo_error "Halogenator CLI not available. Please install or run from project directory."
    exit 1
fi

# Check pandas/pyarrow availability for statistics
PANDAS_AVAILABLE=0
if python scripts/batch_stats.py --check-deps &> /dev/null; then
    PANDAS_AVAILABLE=1
    echo_info "Pandas/PyArrow available for statistics"
else
    echo_warn "Pandas/PyArrow not available - statistics will be skipped or simplified"
fi

# Create directories
mkdir -p "$OUTPUT_BASE"
mkdir -p "$PARENTS_DIR"
mkdir -p "$P1_K2_DIR"
mkdir -p "$P1_K3_DIR"
mkdir -p "$AUDIT_DIR"

echo_info "Starting ETCM 2000 flavonoids batch processing workflow..."

# Step 1: ETCM Data Ingestion
echo_info "Step 1: Converting ETCM Excel to SMILES..."

SAMPLE_ARGS=""
if [ "$SAMPLE_SIZE" -gt 0 ]; then
    SAMPLE_ARGS="--sample $SAMPLE_SIZE"
    echo_warn "Running in SAMPLE MODE with $SAMPLE_SIZE molecules"
fi

python -m halogenator.cli etcm-ingest \
    -i "$ETCM_INPUT" \
    -o "$PARENTS_DIR/parents.smi" \
    --out-meta "$PARENTS_DIR/parents_meta.csv" \
    $SAMPLE_ARGS

if [ $? -ne 0 ]; then
    echo_error "ETCM ingestion failed"
    exit 1
fi

MOLECULE_COUNT=$(wc -l < "$PARENTS_DIR/parents.smi")
echo_info "Initial ingestion: $MOLECULE_COUNT flavonoid parent molecules"

# Sampling validation and correction
if [ "${SAMPLE_SIZE:-0}" -gt 0 ]; then
    if [ "$MOLECULE_COUNT" -gt "$SAMPLE_SIZE" ]; then
        # Fallback: Force trim to SAMPLE_SIZE lines
        echo_warn "More molecules ($MOLECULE_COUNT) than requested ($SAMPLE_SIZE), trimming to sample size"
        head -n "$SAMPLE_SIZE" "$PARENTS_DIR/parents.smi" > "$PARENTS_DIR/parents.sampled.smi"
        mv "$PARENTS_DIR/parents.sampled.smi" "$PARENTS_DIR/parents.smi"
        MOLECULE_COUNT=$(wc -l < "$PARENTS_DIR/parents.smi")
        echo_info "Forced parents.smi to $MOLECULE_COUNT lines for POC sampling"
    elif [ "$MOLECULE_COUNT" -lt "$SAMPLE_SIZE" ]; then
        echo_warn "After dedupe got $MOLECULE_COUNT (< $SAMPLE_SIZE). Proceeding with available parents."
    else
        echo_info "Sample size $MOLECULE_COUNT matches requested $SAMPLE_SIZE"
    fi
fi

echo_info "Final parent molecules count: $MOLECULE_COUNT"

# Step 2: k=2 Halogenation Enumeration
echo_info "Step 2: Running k=2 halogenation enumeration..."
echo_info "Using parent molecules from: $PARENTS_DIR/parents.smi"

python -m halogenator.cli enum \
    -c "$K2_CONFIG" \
    --subset flavonoids \
    --outdir "$P1_K2_DIR" \
    --qa-completion-mode distribute \
    --qa-conflict-tolerance $K2_TOL \
    --qa-max-warnings 500

if [ $? -ne 0 ]; then
    echo_error "k=2 enumeration failed"
    exit 1
fi

# Step 3: Validate k=2 Results
echo_info "Step 3: Validating k=2 results..."

K2_PRODUCTS="$P1_K2_DIR/products_k2.parquet"
K2_QA="$P1_K2_DIR/qa_summary.json"

if [ ! -f "$K2_PRODUCTS" ]; then
    echo_error "k=2 products not found: $K2_PRODUCTS"
    exit 1
fi

if [ ! -f "$K2_QA" ]; then
    echo_warn "k=2 QA summary not found: $K2_QA"
fi

# Create standardized file names for k=2 products
echo_info "Creating standardized product file names..."
cp "$K2_PRODUCTS" "$P1_K2_DIR/products.parquet"

# Generate CSV version using Python
if [ "$PANDAS_AVAILABLE" -eq 1 ]; then
    python -c "
import pandas as pd
df = pd.read_parquet('$K2_PRODUCTS')
df.to_csv('$P1_K2_DIR/products.csv', index=False)
print('Created $P1_K2_DIR/products.csv')
"
else
    echo_warn "Pandas not available - skipping CSV generation"
fi

# Display k=2 statistics
python scripts/batch_stats.py --products "$K2_PRODUCTS" --qa "$K2_QA" --label "k=2"

# Generate audit sample
echo_info "Generating audit sample..."
python scripts/batch_stats.py --products "$K2_PRODUCTS" --audit "$AUDIT_DIR" --audit-size 100

# Step 4: Optional k=3 Enumeration
echo_info "Step 4: k=3 enumeration (optional - resource intensive)..."

# Check for non-interactive mode via environment variable
if [ -n "$RUN_K3" ]; then
    echo_info "RUN_K3 environment variable set to: $RUN_K3"
else
    echo "Do you want to run k=3 enumeration? (y/N): "
    read -t 30 -r RUN_K3 || RUN_K3="n"
fi

if [[ "$RUN_K3" =~ ^[Yy]$ ]]; then
    echo_info "Running k=3 enumeration (tri-substituted products)..."
    echo_info "Using parent molecules from: $PARENTS_DIR/parents.smi"

    python -m halogenator.cli enum \
        -c "$K3_CONFIG" \
        --subset flavonoids \
        --outdir "$P1_K3_DIR" \
        --qa-completion-mode distribute \
        --qa-conflict-tolerance $K3_TOL \
        --qa-max-warnings 1000

    if [ $? -eq 0 ]; then
        echo_info "k=3 enumeration completed successfully"

        # Create standardized file names for k=3 products
        K3_PRODUCTS="$P1_K3_DIR/products_k3.parquet"
        if [ -f "$K3_PRODUCTS" ]; then
            echo_info "Creating standardized k=3 product file names..."
            cp "$K3_PRODUCTS" "$P1_K3_DIR/products.parquet"

            # Generate CSV version using Python
            if [ "$PANDAS_AVAILABLE" -eq 1 ]; then
                python -c "
import pandas as pd
df = pd.read_parquet('$K3_PRODUCTS')
df.to_csv('$P1_K3_DIR/products.csv', index=False)
print('Created $P1_K3_DIR/products.csv')
"
            else
                echo_warn "Pandas not available - skipping k=3 CSV generation"
            fi
        fi

        # Display k=3 statistics
        K3_QA="$P1_K3_DIR/qa_summary.json"
        if [ -f "$K3_PRODUCTS" ]; then
            echo_info "k=3 statistics:"
            python scripts/batch_stats.py --products "$K3_PRODUCTS" --qa "$K3_QA" --label "k=3"
        fi
    else
        echo_warn "k=3 enumeration failed but continuing with k=2 results"
    fi
else
    echo_info "Skipping k=3 enumeration"
fi

# Final Summary
echo ""
echo "=================================================================="
echo "            ETCM 2000 FLAVONOIDS BATCH COMPLETE"
echo "=================================================================="
echo ""
echo "Processing Summary:"
echo "  Input file: $ETCM_INPUT"
echo "  Parent molecules: $MOLECULE_COUNT"
echo "  k=2 configuration: $K2_CONFIG"

if [ "$SAMPLE_SIZE" -gt 0 ]; then
    echo_warn "  MODE: SAMPLE RUN ($SAMPLE_SIZE molecules)"
    echo "  To run full batch: set SAMPLE_SIZE=0 and rerun"
else
    echo "  MODE: FULL BATCH"
fi

echo ""
echo "Generated Outputs:"
echo "  k=2 Products:  $OUTPUT_BASE/products_k2.parquet"
echo "  k=2 CSV:       $OUTPUT_BASE/summary_k2.csv"
echo "  k=2 QA:        $OUTPUT_BASE/qa_summary.json"

if [ -f "$OUTPUT_BASE/k3/products_k3.parquet" ]; then
    echo "  k=3 Products:  $OUTPUT_BASE/k3/products_k3.parquet"
    echo "  k=3 CSV:       $OUTPUT_BASE/k3/summary_k3.csv"
    echo "  k=3 QA:        $OUTPUT_BASE/k3/qa_summary.json"
fi

echo ""
echo "Next Steps:"
echo "1. Review QA summaries for data quality issues"
echo "2. Load products into your analysis pipeline"
echo "3. Consider ADMET/property prediction on the generated library"
echo ""
echo_info "Batch processing workflow complete!"