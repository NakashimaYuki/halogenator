#!/usr/bin/env bash
# ==============================================================================
# Parallel Flavonoid Shard Enumeration (Linux/macOS)
# ==============================================================================
#
# Purpose: Run halogenation enumeration on all shard files in parallel
#
# Features:
# - Controlled parallelism (avoid CPU/memory oversubscription)
# - Individual log files per shard
# - Progress tracking
# - Error handling and resume capability
#
# Usage:
#   bash scripts/03_enum_shards.bash [max_parallel_jobs]
#
# Example:
#   bash scripts/03_enum_shards.bash 4  # Run 4 shards in parallel
#
# ==============================================================================

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# ==============================================================================
# Configuration
# ==============================================================================

# Paths
SHARD_DIR="data/work/shards"
OUTPUT_ROOT="data/output/cnpd_flav_k2"
CONFIG="configs/flavonoids_k2_prod.yaml"
LOG_DIR="data/output/cnpd_flav_k2_logs"

# Parallelism (default: 4, override with command-line argument)
MAX_PARALLEL="${1:-4}"

# RDKit threads per job (adjust based on total CPU cores)
# Formula: total_cores / max_parallel (with some headroom)
# Example: 16 cores, 4 parallel jobs → 3-4 threads each
RDKIT_THREADS="${RDKIT_THREADS:-4}"

# ==============================================================================
# Validation
# ==============================================================================

echo "================================================================================"
echo "Flavonoid Shard Enumeration - Parallel Execution"
echo "================================================================================"
echo ""

# Check if shard directory exists
if [[ ! -d "$SHARD_DIR" ]]; then
    echo "ERROR: Shard directory not found: $SHARD_DIR"
    echo "Please run scripts/02_make_shards.py first"
    exit 1
fi

# Check if config exists
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Configuration file not found: $CONFIG"
    exit 1
fi

# Count shards
SHARD_FILES=("$SHARD_DIR"/flav_shard_*.smi)
NUM_SHARDS=${#SHARD_FILES[@]}

if [[ $NUM_SHARDS -eq 0 ]]; then
    echo "ERROR: No shard files found in $SHARD_DIR"
    exit 1
fi

echo "Configuration:"
echo "  Shard directory: $SHARD_DIR"
echo "  Number of shards: $NUM_SHARDS"
echo "  Output root: $OUTPUT_ROOT"
echo "  Config file: $CONFIG"
echo "  Max parallel jobs: $MAX_PARALLEL"
echo "  RDKit threads per job: $RDKIT_THREADS"
echo ""

# Create output directories
mkdir -p "$OUTPUT_ROOT"
mkdir -p "$LOG_DIR"

# ==============================================================================
# Progress Tracking
# ==============================================================================

PROGRESS_FILE="$LOG_DIR/progress.txt"
COMPLETED_FILE="$LOG_DIR/completed_shards.txt"
FAILED_FILE="$LOG_DIR/failed_shards.txt"

# Initialize progress files
: > "$PROGRESS_FILE"
[[ ! -f "$COMPLETED_FILE" ]] && : > "$COMPLETED_FILE"
[[ ! -f "$FAILED_FILE" ]] && : > "$FAILED_FILE"

# Load completed shards (for resume capability)
declare -A COMPLETED
while IFS= read -r shard; do
    COMPLETED["$shard"]=1
done < "$COMPLETED_FILE"

# ==============================================================================
# Enumeration Function
# ==============================================================================

run_shard() {
    local smi_file="$1"
    local shard_name=$(basename "$smi_file" .smi)

    # Check if already completed
    if [[ -n "${COMPLETED[$shard_name]:-}" ]]; then
        echo "[SKIP] $shard_name (already completed)"
        return 0
    fi

    local outdir="$OUTPUT_ROOT/${shard_name}"
    local logfile="$LOG_DIR/${shard_name}.log"

    echo "[START] $shard_name → $outdir"

    # Create output directory
    mkdir -p "$outdir"

    # Run enumeration
    # Note: No --out-structure hierarchical (uses YAML config: structure: flat)
    if python -m halogenator.cli enum \
        -c "$CONFIG" \
        --rdkit-threads "$RDKIT_THREADS" \
        --outdir "$outdir" \
        -i "$smi_file" \
        > "$logfile" 2>&1; then

        # Success
        echo "[DONE] $shard_name"
        echo "$shard_name" >> "$COMPLETED_FILE"
        echo "$(date '+%Y-%m-%d %H:%M:%S') SUCCESS $shard_name" >> "$PROGRESS_FILE"
        return 0
    else
        # Failure
        echo "[FAIL] $shard_name (see $logfile)"
        echo "$shard_name" >> "$FAILED_FILE"
        echo "$(date '+%Y-%m-%d %H:%M:%S') FAILED $shard_name" >> "$PROGRESS_FILE"
        return 1
    fi
}

export -f run_shard
export OUTPUT_ROOT CONFIG LOG_DIR RDKIT_THREADS PROGRESS_FILE COMPLETED_FILE FAILED_FILE
export -A COMPLETED

# ==============================================================================
# Parallel Execution
# ==============================================================================

echo "================================================================================"
echo "Starting parallel enumeration..."
echo "================================================================================"
echo ""

START_TIME=$(date +%s)

# Use GNU parallel if available, otherwise fall back to manual job control
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel (preferred method)"
    printf '%s\n' "${SHARD_FILES[@]}" | parallel -j "$MAX_PARALLEL" --eta run_shard {}
else
    echo "Using bash job control (GNU parallel not found)"
    echo "Tip: Install GNU parallel for better progress tracking"
    echo ""

    # Manual parallel execution with job control
    for smi_file in "${SHARD_FILES[@]}"; do
        # Wait if max jobs reached
        while [[ $(jobs -r -p | wc -l) -ge $MAX_PARALLEL ]]; do
            sleep 1
        done

        # Launch job in background
        run_shard "$smi_file" &
    done

    # Wait for all remaining jobs
    wait
fi

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

# ==============================================================================
# Summary
# ==============================================================================

echo ""
echo "================================================================================"
echo "Enumeration Complete"
echo "================================================================================"
echo ""

# Count results
COMPLETED_COUNT=$(wc -l < "$COMPLETED_FILE" | tr -d ' ')
FAILED_COUNT=$(wc -l < "$FAILED_FILE" | tr -d ' ')

echo "Results:"
echo "  Total shards: $NUM_SHARDS"
echo "  Completed: $COMPLETED_COUNT"
echo "  Failed: $FAILED_COUNT"
echo ""
echo "Elapsed time: ${ELAPSED}s ($(printf '%02d:%02d:%02d' $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))))"
echo ""

if [[ $FAILED_COUNT -gt 0 ]]; then
    echo "Failed shards:"
    cat "$FAILED_FILE"
    echo ""
    echo "Review log files in $LOG_DIR for error details"
    exit 1
fi

echo "All shards completed successfully!"
echo ""
echo "Next step: Merge shards and run QC"
echo "  python scripts/04_merge_and_qc.py"
echo ""
