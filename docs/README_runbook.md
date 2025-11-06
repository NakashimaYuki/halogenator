# CNPD-ETCM Halogenated Flavonoid Library Pipeline - Operational Runbook

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Pipeline Architecture](#pipeline-architecture)
4. [Step-by-Step Execution Guide](#step-by-step-execution-guide)
5. [Configuration Reference](#configuration-reference)
6. [Troubleshooting](#troubleshooting)
7. [Performance Tuning](#performance-tuning)
8. [Output Files Reference](#output-files-reference)

---

## Overview

This pipeline generates a comprehensive **halogenated flavonoid library** (k≤2) from the **CNPD-ETCM natural product database** (~70,000 compounds).

### Key Features

- **Robust flavonoid extraction** using SMARTS patterns + morphological heuristics
- **Two-stage filtering** (loose → strict) to minimize false negatives
- **RDKit standardization** (largest fragment, metal disconnect, uncharge, tautomer canonical)
- **Parallel shard-based enumeration** for scalability
- **Production configuration**:
  - k_max = 2 (up to double substitution)
  - Deduplication enabled (clean structural library)
  - Sugar mask enabled (protect glycosidic groups)
  - R2 (prenylation) enabled with fallback
  - R6_methyl enabled (step + macro modes)
  - All halogens: F, Cl, Br, I
- **Comprehensive QC validation** suite
- **Statistical analysis** and reporting

### Pipeline Stages

```
CNPD-ETCM Excel
      ↓
[1] Extract Flavonoids (SMARTS + morphology + standardization)
      ↓
[2] Generate Shards (2-5k compounds per shard)
      ↓
[3] Parallel Enumeration (k≤2 halogenation)
      ↓
[4] Merge & QC (global dedup + validation)
      ↓
[5] Statistics & Reports
      ↓
Final Library: products_k2.parquet
```

---

## Prerequisites

### Software Requirements

- **Python 3.8+**
- **RDKit** (`conda install -c conda-forge rdkit`)
- **pandas** (`pip install pandas`)
- **pyarrow** (`pip install pyarrow`) for parquet support
- **openpyxl** (`pip install openpyxl`) for Excel reading
- **halogenator** package (this project)

### System Requirements

**Minimum:**
- 16 GB RAM
- 4 CPU cores
- 50 GB free disk space

**Recommended:**
- 32+ GB RAM
- 8+ CPU cores
- 100+ GB free disk space (for temporary files + outputs)

### Input Data

- `data/input/CNPD-ETCM-合并去重.xlsx`
  - Should contain columns: `Smiles`, `Name`, `ID`, `Pubchem-ID` (or similar)
  - ~70,000 natural product compounds

---

## Pipeline Architecture

### Directory Structure

```
project/
├─ data/
│  ├─ input/
│  │  └─ CNPD-ETCM-合并去重.xlsx          # Input Excel file
│  ├─ work/
│  │  ├─ cnpd_raw.parquet                 # Standardized database
│  │  ├─ flavonoids_candidates.parquet    # Stage I (loose filtering)
│  │  ├─ flavonoids_final.parquet         # Stage II (strict filtering)
│  │  └─ shards/
│  │     ├─ flav_shard_0001.smi           # Shard files
│  │     └─ ...
│  └─ output/
│     ├─ cnpd_flav_k2/                    # Shard enumeration outputs
│     │  ├─ flav_shard_0001/
│     │  │  └─ products_k2.parquet
│     │  └─ ...
│     ├─ cnpd_flav_k2_logs/               # Enumeration logs
│     └─ haloflav_k2/                     # Final merged library
│        ├─ products_k2.parquet           # Main output file
│        ├─ by_rule.csv
│        ├─ parents_coverage.csv
│        ├─ macro_summary.csv
│        ├─ overall_stats.json
│        ├─ SUMMARY_REPORT.txt
│        └─ RUN_METADATA.json
├─ configs/
│  └─ flavonoids_k2_prod.yaml             # Production configuration
├─ scripts/
│  ├─ 01_extract_flavonoids.py
│  ├─ 02_make_shards.py
│  ├─ 03_enum_shards.bash
│  ├─ 03_enum_shards.ps1
│  ├─ 04_merge_and_qc.py
│  └─ 05_summaries.py
└─ README_runbook.md                      # This file
```

---

## Step-by-Step Execution Guide

### Step 0: Environment Setup

```bash
# Activate conda environment (if using conda)
conda activate halogenator

# Verify RDKit installation
python -c "from rdkit import Chem; print('RDKit OK')"

# Verify halogenator package
python -c "from halogenator import cli; print('Halogenator OK')"
```

---

### Step 1: Extract Flavonoids from CNPD-ETCM

**Purpose:** Identify and extract flavonoid compounds using robust filtering

**Command:**

```bash
python scripts/01_extract_flavonoids.py
```

**What it does:**
1. Loads `CNPD-ETCM-合并去重.xlsx`
2. Standardizes molecules (RDKit: largest fragment, metal disconnect, uncharge, tautomer canonical)
3. Deduplicates by InChIKey at source data level
4. Stage I filtering: SMARTS OR morphology OR name keywords → candidates
5. Stage II filtering: SMARTS AND morphology≥2 → final flavonoids

**Output:**
- `data/work/cnpd_raw.parquet` - Standardized full database
- `data/work/flavonoids_candidates.parquet` - Stage I candidates
- `data/work/flavonoids_final.parquet` - Final flavonoid set

**Expected Results:**
- ~70,000 total compounds
- ~5,000-15,000 flavonoid candidates (Stage I)
- ~3,000-10,000 final flavonoids (Stage II)

**Time:** 10-30 minutes (depends on dataset size)

**Validation:**
```bash
# Check final flavonoid count
python -c "import pandas as pd; df = pd.read_parquet('data/work/flavonoids_final.parquet'); print(f'Flavonoids: {len(df)}')"
```

---

### Step 2: Generate Shard Files

**Purpose:** Split flavonoids into manageable chunks for parallel processing

**Command:**

```bash
python scripts/02_make_shards.py \
  --input data/work/flavonoids_final.parquet \
  --output data/work/shards \
  --shard-size 3000
```

**Parameters:**
- `--shard-size`: Number of compounds per shard
  - Recommended: 2000-5000
  - Smaller shards: more parallelism, more overhead
  - Larger shards: less overhead, higher memory usage

**Output:**
- `data/work/shards/flav_shard_0001.smi`
- `data/work/shards/flav_shard_0002.smi`
- ...
- `data/work/shards/shard_manifest.csv` - Metadata

**Format of .smi files:**
```
<SMILES> <PARENT_NAME>
COc1cc... M1_000123
...
```

**Time:** 1-5 minutes

**Validation:**
```bash
# Check shard count
ls data/work/shards/flav_shard_*.smi | wc -l

# Check manifest
cat data/work/shards/shard_manifest.csv
```

---

### Step 3: Parallel Enumeration

**Purpose:** Run halogenation enumeration on all shards in parallel

#### Linux/macOS

**Command:**

```bash
bash scripts/03_enum_shards.bash 4
```

Where `4` is the number of parallel jobs.

#### Windows

**Command:**

```powershell
powershell -ExecutionPolicy Bypass -File scripts/03_enum_shards.ps1 -MaxParallel 4 -RdkitThreads 4
```

**Parameters:**
- `MaxParallel`: Number of shards to process simultaneously
  - Formula: `total_cpu_cores / rdkit_threads`
  - Example: 16 cores, 4 threads per job → 4 parallel jobs
- `RdkitThreads`: RDKit parallelism per job
  - Typical: 2-8 threads
  - Higher is not always better (diminishing returns)

**Resume capability:**

If enumeration is interrupted, you can resume:

```bash
# Linux/macOS - automatically resumes
bash scripts/03_enum_shards.bash 4

# Windows - use -Resume flag
powershell -File scripts/03_enum_shards.ps1 -MaxParallel 4 -Resume
```

**What it does:**
1. Reads configuration from `configs/flavonoids_k2_prod.yaml`
2. For each shard:
   - Runs `python -m halogenator.cli enum -c <config> -i <shard.smi> --outdir <output>`
   - Generates `products_k2.parquet`, `by_rule.csv`, QA reports
3. Tracks progress in `data/output/cnpd_flav_k2_logs/`
4. Logs successes and failures

**Output:**
- `data/output/cnpd_flav_k2/flav_shard_0001/products_k2.parquet`
- `data/output/cnpd_flav_k2/flav_shard_0002/products_k2.parquet`
- ...
- `data/output/cnpd_flav_k2_logs/flav_shard_0001.log`
- `data/output/cnpd_flav_k2_logs/progress.txt`
- `data/output/cnpd_flav_k2_logs/completed_shards.txt`
- `data/output/cnpd_flav_k2_logs/failed_shards.txt`

**Time:** Highly variable
- Small library (1000 flavonoids): 30 minutes - 2 hours
- Medium library (5000 flavonoids): 3-8 hours
- Large library (10000 flavonoids): 8-24 hours

**Progress monitoring:**

```bash
# Linux/macOS
tail -f data/output/cnpd_flav_k2_logs/progress.txt

# Check completed count
wc -l data/output/cnpd_flav_k2_logs/completed_shards.txt

# Check for failures
cat data/output/cnpd_flav_k2_logs/failed_shards.txt
```

**Validation:**

```bash
# Check that all shards completed
diff <(ls data/work/shards/flav_shard_*.smi | wc -l) \
     <(ls data/output/cnpd_flav_k2/flav_shard_*/products_k2.parquet | wc -l)
```

---

### Step 4: Merge Shards and Run QC

**Purpose:** Merge all shard outputs, perform global deduplication, run QC validation

**Command:**

```bash
python scripts/04_merge_and_qc.py \
  --shard-output-root data/output/cnpd_flav_k2 \
  --output-dir data/output/haloflav_k2 \
  --config configs/flavonoids_k2_prod.yaml
```

**Options:**
- `--skip-qc`: Skip quality control validation (not recommended)
- `--skip-dedup`: Skip global deduplication (not recommended)

**What it does:**
1. Finds all `products_k2.parquet` files in shard directories
2. Concatenates into single dataframe
3. Performs global deduplication by InChIKey
4. Saves merged file
5. Runs QC validation suite:
   - Field consistency check
   - Structural consistency check
   - Parent chain integrity check
   - CI validation gate
6. Generates metadata file

**Output:**
- `data/output/haloflav_k2/products_k2.parquet` - **Main output file**
- `data/output/haloflav_k2/products_k2_pre_dedup.parquet` - Pre-dedup version (debug)
- `data/output/haloflav_k2/RUN_METADATA.json` - Run metadata

**QC Validation:**

All QC checks must pass for production use. If any fail:

1. Review error messages in console output
2. Check individual shard logs for problematic shards
3. Investigate data quality issues
4. Re-run failed shards if needed

**Time:** 5-20 minutes (depends on library size)

**Validation:**

```bash
# Check final record count
python -c "import pandas as pd; df = pd.read_parquet('data/output/haloflav_k2/products_k2.parquet'); print(f'Final library: {len(df):,} products')"

# Check metadata
cat data/output/haloflav_k2/RUN_METADATA.json
```

---

### Step 5: Generate Statistics and Reports

**Purpose:** Analyze the final library and generate comprehensive statistics

**Command:**

```bash
python scripts/05_summaries.py \
  --input data/output/haloflav_k2/products_k2.parquet \
  --output-dir data/output/haloflav_k2
```

**What it does:**
1. Loads merged library
2. Generates statistical summaries:
   - Product counts by rule × halogen × k-level
   - Per-parent productivity analysis
   - Macro substitution statistics
   - Overall library metrics
3. Creates human-readable report

**Output:**
- `data/output/haloflav_k2/by_rule.csv`
- `data/output/haloflav_k2/parents_coverage.csv`
- `data/output/haloflav_k2/macro_summary.csv`
- `data/output/haloflav_k2/overall_stats.json`
- `data/output/haloflav_k2/SUMMARY_REPORT.txt`

**Time:** 2-10 minutes

**Review:**

```bash
# View summary report
cat data/output/haloflav_k2/SUMMARY_REPORT.txt

# Explore CSV files
head data/output/haloflav_k2/by_rule.csv
head data/output/haloflav_k2/parents_coverage.csv
```

---

## Configuration Reference

### Production Configuration: `configs/flavonoids_k2_prod.yaml`

```yaml
subset: cnpd_etcm_flavonoids
k_max: 2                      # Maximum substitution level
halogens: [F, Cl, Br, I]      # All halogens
rules: [R1, R2, R3, R6_methyl]

engine_cfg:
  budget_mode: ops            # Count by operations
  rdkit_threads: 8            # Adjust based on CPU cores

dedup:
  enable: true                # Enable deduplication

sugar_cfg:
  mode: heuristic             # Protect glycosidic groups

constraints:
  enable: true                # Chemical feasibility constraints

rules_cfg:
  r2:
    enable: true
    fallback: true            # Broader prenylation coverage

  r6_methyl:
    enable_step: true         # -CH3 → -CH2X
    enable_macro: true        # -CF3, -CCl3
    allow_on_methoxy: true

output:
  structure: flat             # Single parquet file (no hierarchical SDF)
  by_rule_csv: true
  parquet: true

qa:
  summary: true
```

### Key Configuration Options

**k_max:**
- Controls maximum substitution depth
- k_max=1: Single substitutions only
- k_max=2: Up to double substitutions (recommended)
- k_max=3: Triple substitutions (much larger output)

**dedup.enable:**
- `true`: Remove structural duplicates (InChI-based)
- `false`: Keep all enumeration paths (raw mode)
- **Recommendation:** Always `true` for production

**sugar_cfg.mode:**
- `heuristic`: Fast pattern-based detection
- `none`: No sugar protection
- **Recommendation:** `heuristic` for natural products with glycosides

**r2.fallback:**
- `true`: Broader prenylation coverage (may include edge cases)
- `false`: Strict prenylation rules only
- **Recommendation:** `true` for discovery, `false` for conservative production

**r6_methyl options:**
- `enable_step`: Stepwise halogenation (-CH3 → -CH2X)
- `enable_macro`: Macro substitution (-CF3, -CCl3)
- Both can be enabled simultaneously

---

## Troubleshooting

### Issue: Excel file not found

**Error:**
```
FileNotFoundError: Excel file not found: data/input/CNPD-ETCM-合并去重.xlsx
```

**Solution:**
1. Verify file exists: `ls data/input/`
2. Check filename encoding (especially Chinese characters)
3. Place Excel file in correct location

---

### Issue: RDKit import error

**Error:**
```
ImportError: No module named 'rdkit'
```

**Solution:**
```bash
# Install RDKit via conda
conda install -c conda-forge rdkit

# Or if using pip
pip install rdkit-pypi
```

---

### Issue: Memory error during enumeration

**Error:**
```
MemoryError: Unable to allocate array
```

**Solutions:**
1. Reduce shard size: `--shard-size 1000`
2. Reduce parallel jobs
3. Reduce RDKit threads: `--rdkit-threads 2`
4. Close other applications
5. Use machine with more RAM

---

### Issue: QC validation failed

**Error:**
```
[FAIL] Field consistency check
```

**Solutions:**
1. Review console output for specific failures
2. Check individual shard logs in `data/output/cnpd_flav_k2_logs/`
3. Identify problematic shards
4. Delete failed shard outputs and re-run: `bash scripts/03_enum_shards.bash 4`
5. If persistent, report issue with example SMILES

---

### Issue: Enumeration extremely slow

**Symptoms:** Single shard taking >2 hours

**Possible Causes:**
1. Very large/complex molecules (e.g., polysaccharides)
2. Too many reactive sites
3. Inefficient RDKit thread configuration

**Solutions:**
1. Check shard log for slow molecules
2. Filter out problematic structures (MW > 1000, >20 rings)
3. Adjust `--rdkit-threads`
4. Consider k_max=1 instead of k_max=2

---

### Issue: Disk space full

**Error:**
```
OSError: [Errno 28] No space left on device
```

**Solutions:**
1. Check disk space: `df -h`
2. Clean temporary files: `rm -rf data/output/cnpd_flav_k2/*_tmp/`
3. Delete pre-dedup file after validation
4. Use external drive for output
5. Compress old outputs: `tar -czf old_run.tar.gz data/output/haloflav_k2/`

---

## Performance Tuning

### CPU and Memory Optimization

**Rule of thumb:**
```
parallel_jobs × rdkit_threads ≈ physical_cpu_cores
```

**Examples:**

| CPU Cores | Parallel Jobs | RDKit Threads | Notes |
|-----------|---------------|---------------|-------|
| 4 | 2 | 2 | Entry-level |
| 8 | 4 | 2 | Standard |
| 16 | 4 | 4 | Recommended |
| 32 | 8 | 4 | High-performance |
| 64 | 8 | 8 | Server |

**Memory requirements:**

| Shard Size | Parallel Jobs | Est. Memory |
|------------|---------------|-------------|
| 1000 | 4 | 8 GB |
| 3000 | 4 | 16 GB |
| 5000 | 4 | 24 GB |
| 3000 | 8 | 32 GB |

### Disk I/O Optimization

1. Use SSD for output directories
2. Use HDD for long-term storage
3. Avoid network drives for enumeration output
4. Compress completed shards: `gzip data/output/cnpd_flav_k2/*/products_k2.parquet`

### Shard Size Selection

| Library Size | Shard Size | Number of Shards | Parallel Jobs |
|--------------|------------|------------------|---------------|
| 1,000 | 500 | 2 | 2 |
| 5,000 | 2,500 | 2 | 2 |
| 10,000 | 3,000 | 3-4 | 4 |
| 20,000 | 3,000 | 6-7 | 4-8 |
| 50,000 | 5,000 | 10 | 8 |

---

## Output Files Reference

### Primary Output

**`data/output/haloflav_k2/products_k2.parquet`**

Main library file containing all halogenated flavonoids (k≤2).

**Schema:**
- `smiles`: Canonical SMILES
- `inchikey`: InChI key (unique structure identifier)
- `k`: Substitution level (0, 1, 2)
- `k_ops`: Number of halogenation operations
- `k_atoms`: Total halogen atoms added
- `rule`: Rule family (R1, R2, R3, R6_methyl)
- `halogen`: Halogen type (F, Cl, Br, I)
- `parent_inchikey`: Parent structure InChI key
- `substitutions`: Array of substitution steps (JSON)
- `macro_label`: Macro type (CF3, CCl3, or null)
- ... (many other fields)

### Statistics Files

**`by_rule.csv`**
Product counts by rule × halogen × k-level

**`parents_coverage.csv`**
Per-parent productivity:
- k=1 product count
- k=2 product count
- Rules triggered
- Halogens used

**`macro_summary.csv`**
Macro substitution (CF3/CCl3) statistics

**`overall_stats.json`**
Global library statistics (JSON format)

**`SUMMARY_REPORT.txt`**
Human-readable summary report

### Metadata

**`RUN_METADATA.json`**
Run metadata including:
- Timestamp
- Configuration used
- Input parameters
- QC validation status
- Summary statistics

---

## Best Practices

### 1. Always Run QC Validation

Never skip QC validation (`--skip-qc`) for production libraries. QC catches:
- Field inconsistencies
- Structural errors
- Parent chain breaks
- Schema violations

### 2. Keep Deduplication Enabled

Global deduplication is essential for clean libraries. Different parents can produce identical structures.

### 3. Monitor Disk Space

Enumeration generates large temporary files. Ensure >100 GB free space.

### 4. Save Intermediate Files

Keep these for debugging:
- `flavonoids_final.parquet`
- Shard manifest
- Enumeration logs
- Pre-dedup parquet file

### 5. Version Control Configurations

Save configuration files with version tags:
```bash
cp configs/flavonoids_k2_prod.yaml configs/flavonoids_k2_prod_v1.0.yaml
```

### 6. Document Changes

Update `RUN_METADATA.json` with notes:
```json
{
  "notes": "Increased R2 fallback coverage for broader prenylation",
  "config_version": "1.1"
}
```

### 7. Backup Final Outputs

```bash
# Create backup
tar -czf haloflav_k2_$(date +%Y%m%d).tar.gz data/output/haloflav_k2/

# Verify backup
tar -tzf haloflav_k2_$(date +%Y%m%d).tar.gz | head
```

---

## Quick Reference Commands

```bash
# Full pipeline execution (Linux/macOS)
python scripts/01_extract_flavonoids.py
python scripts/02_make_shards.py
bash scripts/03_enum_shards.bash 4
python scripts/04_merge_and_qc.py
python scripts/05_summaries.py

# Check progress
tail -f data/output/cnpd_flav_k2_logs/progress.txt

# Check final output
python -c "import pandas as pd; df = pd.read_parquet('data/output/haloflav_k2/products_k2.parquet'); print(f'Products: {len(df):,}')"

# View summary
cat data/output/haloflav_k2/SUMMARY_REPORT.txt
```

---

## Support and Contact

For issues or questions:

1. Check TECHNICAL_GUIDE.md for detailed technical documentation
2. Review troubleshooting section above
3. Check existing test cases in `tests/` directory
4. Consult previous reports: `G1_K2_DISTRIBUTION_ANALYSIS_REPORT.md`, `POST_K2_HARDENING_REPORT.md`

---

## Changelog

**Version 1.0** (Initial Release)
- Complete end-to-end pipeline
- SMARTS-based flavonoid extraction
- Two-stage filtering
- Parallel shard enumeration
- Comprehensive QC validation
- Statistical analysis and reporting
