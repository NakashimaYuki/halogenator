# Halogenator

![CI](https://github.com/NakashimaYuki/halogenator/actions/workflows/ci.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive halogen substitution system for flavonoid natural products with k-dimensional BFS enumeration capabilities.

## Overview

This project implements rule-based halogen substitution (F, Cl, Br, I) for flavonoids with two phases:

- **P0 Phase**: Single substitutions (k=1) with basic enumeration
- **P1 Phase**: K-dimensional BFS enumeration with advanced constraints, performance optimization, and comprehensive QA statistics

## Architecture

The system provides multiple enumeration algorithms:
- **P0**: Simple k=1 enumeration for baseline functionality  
- **P1**: Advanced k-dimensional BFS with constraint validation, ring awareness, and deduplication

## Rules

- **R1**: Aromatic sp2 carbon with H - replace H with halogen
- **R2**: C ring carbon halogenation with subrules:
  - **R2a**: Aromatic sp2 CH in C-ring (flavonoid context) - enabled by default
  - **R2b**: sp3 CH2 at flavanone C3 or oxygenated rings (dual condition) - enabled by default
- **R3**: Replace -OH with halogen (excluding carboxylic acid)
- **R4**: Replace -NHx with halogen (excluding amide)
- **R5**: Replace -C(=O)OH with halogen (whole carboxyl group)

**Note**: R2a/R2b are enabled by default. Use `EnumConfig(..., rules_cfg={'R2': {'sp2_CH_in_C_ring': False, 'sp3_CH2_flavanone': False}})` to opt out.

#### R2b Beta-to-Carbonyl Configuration

R2b site detection uses configurable beta-to-carbonyl semantics for precise control over flavanone C3 position targeting:

```python
# Strict beta=2 bonds (default behavior)
rules_cfg = {'R2': {'allow_alpha_as_beta': False}}

# Alpha-as-beta compatibility mode (includes both alpha=1 and beta=2 positions)
rules_cfg = {'R2': {'allow_alpha_as_beta': True}}
```

**Default Behavior**: Strict beta=2 bond distance ensures only true beta positions are detected.

**Compatibility Mode**: When enabled, also accepts alpha=1 positions for legacy compatibility.

**YAML Configuration**:
```yaml
rules_cfg:
  R2:
    allow_alpha_as_beta: false  # true for compatibility mode
```

#### R2b Fallback Detection Mode

R2b site detection includes a **fallback mechanism** designed for flavanone molecules where strict distance-based detection may miss valid sites. This is particularly relevant for compounds like naringenin that lack aromatic character at the C-ring.

**Default Behavior (Strict Mode):**
- R2b fallback is **disabled by default**
- Only strict aromatic sp2 and beta-carbonyl sites are detected
- Conservative approach minimizes false positives

**Raw Mode (Fallback Enabled):**
- Raw enumeration mode (`--no-constraints --no-sugar-mask --no-sym-fold --no-dedup`) **automatically enables** R2b fallback
- Fallback detection uses relaxed criteria for flavanone-specific patterns
- Detected sites are tagged with `detection='fallback'` in output

**Explicit Control:**
```bash
# Enable R2b fallback explicitly
halogenator enum -c config.yaml --r2-fallback

# Raw mode enables fallback automatically
halogenator enum -c config.yaml --no-constraints --no-sugar-mask --no-sym-fold --no-dedup
```

**YAML Configuration:**
```yaml
rules_cfg:
  R2:
    fallback:
      enable: true  # Explicitly enable R2b fallback detection
```

**Output Fields:**
- `detection`: Detection method used ('strict' or 'fallback')
- `sub_rule`: Sub-rule identifier (R2a, R2b)
- `rule_family`: Grouped family name (R2)

**When to Use Fallback:**
- Testing flavanone-specific enumeration
- Exploratory raw mode runs
- Validating R2b site detection on non-aromatic C-rings

**Note:** Fallback mode increases product count for molecules with flavanone-like structures. Review QA statistics to understand the impact on your specific dataset.

### Rule Execution Order

Rules are executed in a **FIXED** order: **R3, R4, R5** (reaction-type rules) first, then **R1, R2** (site-type rules). 

When configuring rules (e.g., via `--rules` or `EnumConfig(rules=...)`), the parameter only controls which rules are **enabled/disabled**. It does NOT change the execution order, which always follows the pattern above for deterministic and consistent results.

## Quick Start

### Environment Setup

```bash
# Setup environment
conda create -n halogenator python=3.10 -y
conda activate halogenator

# Install RDKit via conda (required - not available via pip on all platforms)
conda install -c conda-forge rdkit -y

# Install other dependencies via pip
pip install pandas pyarrow pyyaml

# Install package (pure Python dependencies only)
pip install -e .
```

### P0 Phase (k=1 enumeration)

```bash
# Run P0 demo
bash scripts/run_p0_demo.sh
```

### P1 Phase (k-dimensional enumeration)

```bash
# Run P1 k=2 enumeration with constraints
halogenator enum -c configs/p1.yaml

# Run P1 performance baseline (10-15 minutes target)
halogenator benchmark --p1-baseline --config configs/p1-baseline.yml

# Compare with DL benchmark results (enhanced with fixed schema)
halogenator compare-dl --ours data/output/p1/products.csv --dl path/to/dl_results.csv --key-strategy inchikey
```

## Usage

### P0 Commands (k=1 enumeration)

```bash
# Ingest flavonoid names and standardize
halogenator ingest -c configs/p0.yaml

# Generate k=1 halogenated products (default: flavonoids only)
halogenator k1 -c configs/p0.yaml

# Generate k=1 products for rule probes (proves R3/R4/R5 functionality)  
halogenator k1 -c configs/p0.yaml --subset probes

# Generate k=1 products for all molecules (flavonoids + probes)
halogenator k1 -c configs/p0.yaml --subset all

# Generate summary report (default: flavonoids)
halogenator report -c configs/p0.yaml

# Use custom output directory (optional)
halogenator k1 -c configs/p0.yaml --outdir tmp/run1
halogenator report -c configs/p0.yaml --outdir tmp/run1
```

### P1 Commands (k-dimensional enumeration)

```bash
# P1 k=2 enumeration with BFS and constraints
halogenator enum -c configs/p1.yaml --k-max 2

# P1 enumeration with custom constraints
halogenator enum -c configs/p1.yaml --constraints '{"per_ring_quota": 2, "min_graph_distance": 2}'

# Performance benchmarking
halogenator benchmark p1-baseline --config configs/p1-baseline.yml --outdir benchmarks/

# Run standard performance suite
halogenator benchmark --standard --output-dir benchmarks/

# Compare results with DL baseline (enhanced mode with low memory option)
halogenator compare-dl --ours data/output/p1/products.csv --dl dl_baseline.csv --key-strategy inchikey --low-mem
```

## Outputs

### Input Files
- `data/input/parents.smi` - Flavonoid molecules (default subset)
- `data/input/rule_probes.smi` - Rule validation probes (phenol, ethylamine, benzoic_acid)
- `examples/input/parents_flavonoids_10.smi` - **Stable reference**: 10 standard flavonoids for reproducible testing
- `configs/p1.yaml` - P1 enumeration configuration
- `configs/p1-baseline.yml` - P1 baseline performance configuration

### Product Files

**P0 Outputs:**
- `data/output/p0/products_k1.sdf` - Flavonoid products (SDF format)
- `data/output/p0/products_k1.parquet` - Flavonoid products (table format)
- `data/output/p0/products_k1_probes.*` - Probe products (with _probes suffix)
- `data/output/p0/products_k1_all.*` - All molecules (with _all suffix)

**P1 Outputs:**
- `data/output/p1/products.parquet` - K-dimensional enumeration results
- `data/output/p1/products.csv` - Products in CSV format for DL comparison
- `data/output/p1/dedup_summary.json` - Deduplication statistics
- `data/output/p1/constraints_violations.csv` - Constraint violation analysis

### Report Files

**P0 Reports:**
- `summary_k1.csv` - Summary statistics (variable schema, backward compatible)
- `summary_k1_pivot.csv` - Rule x halogen matrix

**P1 Reports with Fixed Schema (54 columns for DL workflows):**
- `summary.csv` - **Fixed schema summary** with consistent 54-column structure:
  - 6 metadata columns (RDKit version, platform, timestamp, etc.)
  - 5 overall statistics columns (total parents, products, QC stats)
  - 3 per-parent statistics columns (avg, max, min products)
  - 5 rule breakdown columns (R1-R5 products)
  - 4 halogen breakdown columns (F, Cl, Br, I products)
  - 20 rule x halogen matrix columns (R1-R5 x F,Cl,Br,I)
  - 11 diagnostics columns (QA paths, consistency checks)

**P1 CSV Files:**
- `rule_halogen_k.csv` - **Fixed schema**: parent_key, rule, halogen, k, count
- `constraints_violations.csv` - **Fixed schema**: constraint summary with top 10 violation types
- `pivot_rule_halogen_k.csv` - Rule x halogen x k pivot table

### QA Statistics Files

**Enhanced QA Summary (P1):**
- `qa_summary.json` - Quality assurance statistics with 5-invariant consistency checking
- `qa_summary.metadata.performance` - Performance benchmarking data (when available)

**QA JSON Version Handling:**
- **v2 input (has `version` field)**: Passed through unchanged, preserves pivots structure
- **totals-only input**: 
  - If `metadata.halogens` exists -> generates `version: '1'` with granular slices (`by_rule`, `by_halogen`, `by_rule_halogen`)
  - Otherwise -> generates `version: '1'` with `total` only (legacy compatibility)

**P1 QA Statistics Schema:**
- `attempts`: Total enumeration attempts (rule x halogen x k combinations)  
- `products`: Number of successful attempts
- `no_product_matches`: Failed attempts with supported templates
- `template_unsupported`: Attempts with unsupported templates
- `qa_paths`: Fallback method usage (isotope strategy, AtomMap, heuristics)
- `pivots`: Multi-dimensional breakdowns by rule, halogen, and k value
- `consistency_checks`: 5-invariant validation results

**5-Invariant Consistency Checking:**
1. **H dimension**: Sum_rule,k == halogen_counts[halogen] (existing)
2. **(rule, halogen) dimension**: Sum_k == rule_halogen_counts[rule][halogen] 
3. **k dimension**: Sum_rule,halogen == k_counts[k]
4. **Global**: Sum_rule,halogen,k == attempts_total/pivots.attempts
5. **Parent aggregation**: Sum_parent rule_halogen_k_by_parent == rule_halogen_k_counts

**Core Invariant:** `attempts = products + no_product_matches + template_unsupported`

### DL Comparison Files (Enhanced with Fixed Schema)

**P1 DL Comparison Outputs:**
The `compare-dl` command generates **4 files with fixed schemas** for reliable downstream processing:

1. **`matches.csv`** - Products found in both systems (Fixed Schema):
   - Core columns (fixed order): `key, smiles, inchikey, rule, halogen, k, parent_key_or_smiles, count_ours, count_dl`
   - Extension columns: Additional `ours_*` and `dl_*` columns from input data

2. **`only_ours.csv`** - Products only in our system (Fixed Schema):
   - Core columns (fixed order): `key, smiles, inchikey, rule, halogen, k, parent_key_or_smiles, count`
   - Extension columns: Additional `ours_*` columns from input data

3. **`only_dl.csv`** - Products only in DL system (Fixed Schema):
   - Core columns (fixed order): `key, smiles, inchikey, rule, halogen, k, parent_key_or_smiles, count`
   - Extension columns: Additional `dl_*` columns from input data

4. **`crosstab_rule_halogen_k.csv`** - Dimensional analysis (Fixed Schema):
   - Columns: `rule, halogen, k, ours_count, dl_count, intersect_count, only_ours_count, only_dl_count, match_flag`

5. **`summary.json`** - Comprehensive comparison metadata including:
   - Join key normalization statistics (InChIKey generation, fallback usage)
   - Validation reports (column usage, warnings, errors)
   - Summary statistics (overlap percentages, coverage metrics)
   - Configuration metadata (key strategy, low memory mode)

**Join Key Normalization Strategy:**
- **Default**: InChIKey-based matching (most robust for isomeric differences)
- **Fallback 1**: Canonical SMILES (when InChIKey unavailable)
- **Fallback 2**: Raw SMILES (when RDKit unavailable)
- Automatic deduplication and counting of identical keys

**Advanced Features:**
- `--key-strategy inchikey|canonical_smiles|raw_smiles` - Control join strategy
- `--low-mem` - Enable streaming mode for large datasets (>100K records)
- Strict validation with graceful degradation and detailed error reporting
- Built-in column usage analysis and schema compliance checking

**DL Comparison Commands:**
```bash
# Basic comparison with InChIKey matching (recommended)
halogenator compare-dl --ours our_products.csv --dl dl_products.csv --outdir comparison_results/

# Low memory mode for large datasets
halogenator compare-dl --ours our_products.csv --dl dl_products.csv --outdir comparison_results/ --low-mem

# Custom join key strategy
halogenator compare-dl --ours our_products.csv --dl dl_products.csv --outdir comparison_results/ --key-strategy canonical_smiles

# Legacy compatibility mode
halogenator compare-dl --ours our_products.csv --dl dl_products.csv --outdir comparison_results/ --key smiles
```

**Performance Benchmarking:**
- **Hardware tracking**: Platform, cores, memory (including true peak memory via platform-specific methods)
- **Timing metrics**: P50/P95 per-molecule processing times
- **Target compliance**: Duration <5min, Memory <8GB, Zero failures
- **Baseline molecules**: 10 representative flavonoids for consistent performance testing

**Unified Parent Key Generation:**
The system uses unified parent key generation to ensure consistency between input molecules and output products:

- **Primary**: `IK:<InChIKey>` when RDKit/InChI available (uppercase)
- **Fallback**: `SMI:<normalized_smiles>` when RDKit unavailable  
- **Default**: `SMI:UNKNOWN` for invalid/missing identifiers

**Key Specifications:**
- Case-insensitive "Unknown" handling (`Unknown` = `unknown` = `UNKNOWN`)
- Long SMILES (>120 chars) automatically hashed for storage efficiency
- Whitespace normalization ensures consistent key generation
- InChIKey storage forced to uppercase for standardization
- `normalized_length`: Length after whitespace compression and strip; 0 for unknowns
- `canonical_length`: Present only when RDKit canonicalization succeeds

**RDKit Degradation Matrix:**
- RDKit + InChI working -> `IK:<inchikey>`
- RDKit + InChI failed -> `SMI:<canonical_smiles>`
- RDKit unavailable -> `SMI:<input_smiles>`
- Invalid input -> `SMI:UNKNOWN`

**Legacy Compatibility:**
- Old dedup field names (`statesig_hits`/`inchi_hits`) are read for backward compatibility
- New field names (`dedup_hits_statesig`/`dedup_hits_inchi`) are used in output
- Schema version 2 includes pivot breakdowns for detailed analysis
- Both root-level and nested QA JSON formats supported: `qa_data.get('total', qa_data)`

### Temporary Files
- `parents_all.smi` - Combined parent molecules file (created automatically for 'all' subset)
  - When using `--outdir`, this file is placed in the specified output directory and can be safely deleted after report generation
  - When no `--outdir` is specified, a temporary directory is used and automatically cleaned up

## P1 Advanced Features

### Constraints System
The P1 phase includes comprehensive constraint validation:

- **Per-ring quota**: Maximum substitutions per aromatic ring
- **Min graph distance**: Minimum bond distance between substitution sites  
- **Max per halogen**: Maximum count for each halogen type
- **Ring awareness**: Intelligent ring detection and labeling
- **BFS deduplication**: State signature-based duplicate prevention

### Constraint Configuration
```yaml
constraints:
  per_ring_quota: 2           # Max 2 substitutions per ring
  min_graph_distance: 2       # Min 2 bonds between sites
  max_per_halogen:           # Per-halogen limits
    F: 3
    Cl: 2
    Br: 1
    I: 1
```

### Performance Optimization
- **Lazy RDKit loading**: TYPE_CHECKING pattern prevents import errors
- **Incremental statistics**: Memory-efficient accumulation for large datasets
- **Reaction caching**: Template compilation and reuse
- **Ring labeling cache**: Persistent ring analysis caching
- **Platform-specific memory tracking**: True peak memory measurement

### Sugar Masking & Proximity Guard

For flavonoid glycosides, the system provides advanced sugar masking with proximity guard functionality to prevent halogenation near sugar moieties:

#### Basic Sugar Masking
```bash
# Enable heuristic sugar masking (identifies and masks sugar rings)
halogenator enum quercetin-3-glucoside.smi --sugar.mode=heuristic

# Disable sugar masking (enumerate all sites)
halogenator enum quercetin-3-glucoside.smi --sugar.mode=off
```

#### Proximity Guard System
```bash
# Enable proximity guard with 2-bond radius (prevents halogenation within 2 bonds of sugar)
halogenator enum glycoside.smi --sugar.mode=heuristic --sugar.proximity-guard-radius=2

# Enable proximity guard with 3-bond radius (more conservative filtering)
halogenator enum glycoside.smi --sugar.mode=heuristic --sugar.proximity-guard-radius=3

# Disable proximity guard (only mask sugar ring atoms directly)
halogenator enum glycoside.smi --sugar.mode=heuristic --sugar.proximity-guard-radius=0
```

#### Legacy Compatibility
```bash
# Enable legacy key emission for backward compatibility
halogenator enum glycoside.smi --sugar.mode=heuristic --sugar.emit-legacy-keys

# Batch processing with sugar filtering
halogenator batch config.yml --sugar.mode=heuristic --sugar.proximity-guard-radius=3
```

#### Configuration via YAML
```yaml
sugar_masking:
  mode: heuristic                 # off | heuristic
  proximity_guard_radius: 3       # Bond distance (0=disabled, 2-3=typical)
  emit_legacy_keys: false         # Backward compatibility
  audit: true                     # Enable detailed audit fields
```

#### Advanced Monitoring
```bash
# Enable cache monitoring for debugging performance
export HALOGENATOR_CACHE_MONITOR=1
halogenator enum glycoside.smi --sugar.mode=heuristic

# Run acceptance tests with sugar filtering validation
python scripts/run_pr1_acceptance.py --verbose
```

#### Sugar Events Output
The system tracks three types of sugar-related filtering events:
- **`sugar_mask_filtered`**: Sites filtered by direct sugar mask
- **`sugar_proximity_filtered`**: Sites filtered by proximity guard
- **`post_guard_blocked`**: Sites filtered by post-enumeration guards

These appear in QA statistics and can be monitored for effectiveness:
```json
{
  "qa_paths": {
    "sugar_mask_filtered": 25,
    "sugar_proximity_filtered": 8,
    "post_guard_blocked": 3
  }
}
```

#### Sugar Detection APIs

The system provides two APIs for programmatic sugar detection with backward compatibility:

##### Backward-Compatible API
```python
from halogenator.sugar_mask import get_sugar_mask_with_full_status

# Returns (mask, degraded, audit) - maintains legacy behavior
mask, degraded, audit = get_sugar_mask_with_full_status(mol, sugar_cfg)

# mask: set of atom indices to mask
# degraded: bool indicating if fallback path was used
# audit: dict with detailed sugar detection metadata
```

##### New Recommended API
```python
from halogenator.sugar_mask import get_sugar_mask_status

# Returns (mask, accepted, audit) - clearer semantics
mask, accepted, audit = get_sugar_mask_status(mol, sugar_cfg)

# mask: set of atom indices to mask
# accepted: bool indicating if molecule was accepted as sugar
# audit: dict with detailed sugar detection metadata
```

##### Audit Schema
The audit dictionary contains centralized fields defined in `src.halogenator.schema`:

```python
# Core audit fields (always present)
audit = {
    'path': 'main',                    # 'main' or 'fallback'
    'accepted': True,                  # Sugar acceptance status
    'accepted_via_score': True,        # Accepted through scoring (not fallback)
    'ring_size': 6,                   # Size of detected sugar ring
    'exocyclic_O_count_single_bond': 3, # Exocyclic oxygen count
    'has_cglyco_evidence': False,      # C-glycoside pattern detected
    'ring_score': 12.5                # Evidence-based ring score
}
```

##### Structured Evidence Objects
The system uses `SugarRingEvidence` dataclass for eliminating duplicate scoring calculations:

```python
from halogenator.sugar_mask import SugarRingEvidence

# Evidence contains all scoring data for a sugar ring
evidence = SugarRingEvidence(
    ring=(0, 1, 2, 3, 4, 5),          # Tuple of ring atom indices
    ring_size=6,                       # Ring size
    score=12.5,                        # Calculated evidence score
    exocyclic_o_count_single_bond=3,   # Single-bond exocyclic oxygens
    has_cglyco_evidence=False          # C-glycoside pattern present
)
```

This structured approach ensures consistent scoring across all sugar detection functions and eliminates the need for recalculation in different code paths.

## Testing

### Unit Tests
Run core unit tests (fast, no external dependencies):

```bash
python -m unittest discover -v
```

### Integration Tests  
Run integration tests (slower, requires file I/O):

```bash
HALO_INTEGRATION=1 python -m unittest tests.test_subset_consistency -v
```

### P1 Performance Tests
Run P1-specific performance and validation tests:

```bash
# Run P1 baseline performance test
HALO_BATCH_SIZE=5 python -m halogenator.cli benchmark p1-baseline

# Test consistency checking (5 invariants)  
python -m unittest tests.test_qa_count_semantics tests.test_qa_paths_dimensionality -v

# Test constraint validation
python -m unittest tests.test_p1_constraints -v
```

### Code Quality
Check ASCII compliance of all source files:

```bash
./scripts/check_ascii.sh
```

## Known Issues and Version Compatibility

### RDKit Compatibility
The system uses **TYPE_CHECKING patterns** for robust RDKit integration:
- Module imports work even when RDKit is unavailable
- Runtime RDKit functionality only loaded when needed
- Graceful degradation for unavailable RDKit components
- Some RDKit versions may print non-fatal warnings that do not affect results

### Recommended Versions
The system has been tested with the following dependency versions:
- **RDKit**: 2023.09 or 2024.03+ (recommended via conda-forge)
- **Python**: 3.10 or 3.11
- **Pandas**: 2.0+ 
- **PyArrow**: 12.0+
- **PyYAML**: 6.0+

### Platform Support
- **Windows**: Full support with platform-specific peak memory tracking
- **Linux**: Full support with resource.getrusage memory tracking  
- **macOS**: Full support with macOS-specific memory unit handling
- **Memory tracking**: Automatically selects best method per platform

### Performance Characteristics
- **P0**: Optimized for k=1 enumeration, ~1-5min typical runtime
- **P1**: Target <5min for k=2 enumeration on 4 cores with <8GB RAM
- **Memory efficiency**: Incremental statistics prevent OOM on large datasets
- **Caching**: Reaction templates and ring labeling cached for performance

### Installation Notes
- Install RDKit via `conda install -c conda-forge rdkit` (recommended)
- Avoid installing RDKit via pip as it may lack platform-specific optimizations  
- Different RDKit versions may produce slightly different warning messages but should generate identical molecular outputs
- TYPE_CHECKING imports prevent runtime dependency issues

## Development Standards

### ASCII Only
All source code, documentation, and configuration files must contain only ASCII characters (bytes 0-127). No Unicode characters are allowed.

**Files Checked:**
- Python source files (`.py`)
- Documentation files (`.md`)
- Configuration files (`.yml`, `.yaml`)
- Shell scripts (`.sh`)

**Automated Enforcement:**
- Pre-commit hooks validate ASCII compliance with cross-platform Python implementation
- CI pipeline includes ASCII checking for all platforms
- Unit tests verify repository-wide compliance

**Cross-Platform Implementation:**
The ASCII checks use pure Python implementations for maximum compatibility:
- `scripts/check_ascii_portable.py` - Checks source directories (src/, scripts/, docs/)
- `scripts/check_ascii_root_files.py` - Checks root documentation files (README.md, CHANGELOG.md, etc.)

**Manual Check:**
```bash
# Run ASCII compliance test for all source files
python -m unittest tests.test_repository_ascii_compliance -v

# Check root documentation files specifically
python scripts/check_ascii_root_files.py

# Check source directories
python scripts/check_ascii_portable.py

# Run pre-push checks (includes ASCII validation)
bash scripts/pre-push.sh
```

### Type Safety
- Uses TYPE_CHECKING pattern for optional dependencies
- String literal type annotations for forward references
- Lazy imports for runtime-only dependencies

### Testing Coverage
- Unit tests for all core functionality
- Integration tests for file I/O workflows
- Performance regression tests for P1 baseline
- 5-invariant consistency validation in QA pipeline

## QA Consistency Checks

The report includes QA consistency validation with the following structure:

- Key: `global_product_conservation`
  - Sub-buckets:
    - `global_products`: sum over rule/halogen/k vs `product_count`
    - `product_count_vs_enum`: `product_count` vs enumeration `products`
    - `product_vs_halogen_sum`: `product_count` vs sum of `halogen_counts`
- Key: `three_way_mutex` (separate)
  - Checks the invariant: `attempts == products + no_product_matches + template_unsupported`
  - Runs only when enumeration QA stats have attempts; skipped in degrade mode or when attempts are absent.

The JSON report exposes `qa_source` metadata and boolean `three_way_mutex_checked`.

## V1 Granular QA Statistics Distribution

The P1 phase generates granular QA statistics (version 1) with deterministic distribution algorithms:

### Distribution Algorithm and Remainder Allocation

- **Primary Distribution**: All statistics are first distributed to the 2D structure (`by_rule_halogen`)
- **Remainder Allocation**: Uses lexicographic ordering to distribute remainders deterministically
  - When totals don't divide evenly, remainder values are allocated to the first keys in sorted order
  - Example: 5 items distributed among ['R1', 'R3'] results in R1=3, R3=2 (R1 gets the extra due to sorting)
- **Marginal Consistency**: `by_rule` and `by_halogen` are derived from the 2D structure to ensure:
  - `sum(by_rule_halogen[rule, *]) == by_rule[rule]` for all rules
  - `sum(by_rule_halogen[*, halogen]) == by_halogen[halogen]` for all halogens
  - `sum(by_rule_halogen[*, *]) == totals` for all metrics

### Metadata Documentation

V1 granular outputs include metadata fields documenting the distribution semantics:
- `distribution: "equal_with_lexicographic_remainder"` - Distribution algorithm used
- `marginals_from_2d: true` - Confirms margins are derived from 2D structure

This ensures full transparency and auditability of statistical breakdowns in downstream analysis workflows.

### Non-distributed Overview Counters

V1 granular statistics distinguish between two types of metrics:

- **Distributable Metrics**: QA path counters (isotope_unavailable, atommap_used, etc.) and attempt outcome counters (no_product_matches, template_unsupported) that are distributed across granular dimensions (by_rule, by_halogen, by_rule_halogen)
- **Overview Counters**: High-level statistics (attempts, products, dedup_hits_*) that provide overall enumeration summaries but are NOT distributed to granular dimensions

Overview counters appear only in the `total` section when present in input data, providing context for the distributable metrics while maintaining separation of concerns in the granular breakdown structure.

#### Passthrough Strategy for Overview Counters

- **When `total` is missing**: The system merges top-level overview counters from the input into the constructed `total`, alongside distributable metrics aggregated using multi-path fallback strategy
- **When `total` exists**: The system preserves the existing `total` without modification, maintaining upstream data integrity

#### Multi-Path Fallback for Total Calculation

When `total` is missing in passthrough mode, the system uses a robust fallback strategy to compute distributable metrics:

1. **Priority 1**: Aggregate from `by_rule_halogen` (2D structure - most reliable)
2. **Priority 2**: Fallback to aggregate from `by_rule` (1D fallback)
3. **Priority 3**: Fallback to aggregate from `by_halogen` (1D fallback)
4. **Priority 4**: Default to 0 if no granular data available

This ensures that partial v1 granular inputs (e.g., only `by_rule` provided) can still produce meaningful totals without being forced to zero.

#### Granular Structure Sanitization and Normalization

During passthrough processing, the system applies two cleanup operations to ensure consistent structure:

1. **Sanitization**: Removes non-distributable fields (overview counters) from granular dimensions if they were mistakenly included
2. **Normalization**: Fills missing metrics with 0 in existing granular nodes to ensure consistent shape across all dimensions

These operations maintain data integrity while enforcing the separation between distributable metrics and overview counters.

#### Distribution Order vs Display Order

Note that **distribution order may differ from display order**:

- **Distribution order**: Uses lexicographic (alphabetical) ordering for deterministic remainder allocation during equal distribution
- **Display order**: May follow a different convention (e.g., F, Cl, Br, I for atomic number order)

The exact keys used for distribution are documented in `metadata.distribution_keys` for full algorithmic transparency. Note that distribution semantics primarily apply to construction scenarios (legacy input -> v1 granular) and internal structure completion; passthrough scenarios preserve input granular values and may not undergo distribution.

## Passthrough Completion Modes

The P1 phase supports two completion strategies when processing v1 granular input with missing dimensions:

### Mode Comparison

| Aspect | `zero_fill` (default) | `distribute` |
|--------|----------------------|-------------|
| **Purpose** | Preserve input fidelity | Derive missing dimensions |
| **Missing dimensions** | Filled with zeros | Distributed from marginals |
| **Total source** | Multi-path fallback | Multi-path fallback |
| **marginals_from_2d** | False (input preserved) | True (derived from 2D) |
| **Use case** | Data validation, conservative analysis | Missing data completion |

### 2D as Single Source of Truth (SoT) Execution Order

When completion mode is `distribute` or original input contains 2D structure:

1. **Establish 2D Structure**: Create or validate `by_rule_halogen` as authoritative source
2. **Derive by_rule**: Sum 2D values across halogens for each rule
3. **Derive by_halogen**: Sum 2D values across rules for each halogen
4. **Set metadata.marginals_from_2d = true**: Indicate marginals are derived, not preserved

When completion mode is `zero_fill` and no original 2D structure:

1. **Preserve Input Marginals**: Keep existing `by_rule` and `by_halogen` values
2. **Complete Missing Dimensions**: Fill absent dimensions with all-zero structures
3. **Create Zero 2D**: Generate `by_rule_halogen` with all zeros
4. **Set metadata.marginals_from_2d = false**: Indicate marginals are from input

### Total Calculation Fallback Order

When `total` is missing in passthrough mode, the system uses this priority:

1. **Priority 1**: Sum from 2D structure (if 2D present and valid)
2. **Priority 2**: Use by_rule totals (if by_rule present)
3. **Priority 3**: Use by_halogen totals (if by_halogen present)
4. **Priority 4**: Default to 0 (if no granular data available)

### Base Selection and Cross-Side Fallback Policy

When both `by_rule` and `by_halogen` marginals exist in distribute mode, the system follows strict base priority rules:

**Base Priority**: `by_rule > by_halogen` - when both marginals are present, `by_rule` is always selected as the authoritative base for 2D distribution.

**No Cross-Side Metric Fallback**: When `base=by_rule`, metrics that exist only in `by_halogen` are NOT used to fill missing values in the 2D structure. The 2D distribution uses only metrics present in the selected base marginal.

**Conflict Detection**: While only the base marginal is used for distribution, conflict detection still compares totals between both marginals and generates warnings when differences exceed the tolerance threshold.

#### Example Behavior

```
Input:
  by_rule: {"R1": {"metric_A": 10}}
  by_halogen: {"F": {"metric_A": 8, "metric_B": 5}}

Result (base=by_rule):
  by_rule_halogen: {"R1": {"F": {"metric_A": 5, "metric_B": 0}, "Cl": {"metric_A": 5, "metric_B": 0}}}
```

**Explanation:**
- `by_rule` is missing `metric_B`, `by_halogen` has both metrics
- `by_rule` takes priority as base, so `metric_B` is set to 0 in 2D structure
- `metric_B` is NOT derived from `by_halogen` value (5)
- Conflict warning generated for `metric_A` (10 vs 8) and `metric_B` (0 vs 5)

### Warnings System

The system generates warnings for data quality issues:

#### Empty Sets Warning
```json
{
  "warnings": [{
    "type": "empty_rules_or_halogens"
  }]
}
```
Triggered when rules or halogens lists are empty, preventing division-by-zero errors.

#### Marginal Conflict Warning
```json
{
  "warnings": [{
    "type": "marginal_conflict",
    "base": "by_rule",
    "tolerance": 1,
    "delta": {
      "isotope_unavailable": {
        "lhs": 10,
        "rhs": 13,
        "diff": 3
      }
    }
  }]
}
```
Triggered when by_rule and by_halogen totals differ beyond tolerance threshold. Includes detailed breakdown with left-hand side (by_rule), right-hand side (by_halogen), and difference values. Uses union of metrics to catch inconsistencies in any dimension.

### Overview vs Distributable Metrics

**Distributable Metrics** (distributed across granular dimensions):
- QA path counters: `isotope_unavailable`, `isotope_miss`, `atommap_used`, `heuristic_used`
- Attempt outcome counters: `no_product_matches`, `template_unsupported`

**Overview Counters** (only in total, not distributed):
- Process counters: `attempts`, `products`
- Deduplication counters: `dedup_hits_statesig`, `dedup_hits_inchi`

### CLI Usage

```bash
# Use zero-fill completion (default - preserves input marginals)
halogenator enum -c config.yaml --qa-completion-mode zero_fill

# Use distribute completion (derives missing dimensions from marginals)
halogenator enum -c config.yaml --qa-completion-mode distribute

# Check completion strategy in output metadata
cat data/output/p1/qa_summary.json | jq '.metadata.completion.strategy'
```

The completion strategy and base source are documented in output metadata for full transparency and auditability.

### Minimal Complete Examples

*Note: Examples below demonstrate field structures and completion semantics. Specific domain definitions (rules, halogens) are typically injected by upstream configuration in real usage scenarios.*

#### Zero Fill Mode (Default)
Input with partial granular data:
```json
{
  "version": "1",
  "by_rule": {"R1": {"isotope_unavailable": 5}},
  "total": {"isotope_unavailable": 5, "attempts": 100},
  "metadata": {"rules": ["R1"], "halogens": ["F", "Cl"]}
}
```

Output preserves input marginals and fills missing dimensions with zeros:
```json
{
  "version": "1",
  "total": {"isotope_unavailable": 5, "attempts": 100},
  "by_rule": {"R1": {"isotope_unavailable": 5}},
  "by_halogen": {"F": {"isotope_unavailable": 0}, "Cl": {"isotope_unavailable": 0}},
  "by_rule_halogen": {"R1": {"F": {"isotope_unavailable": 0}, "Cl": {"isotope_unavailable": 0}}},
  "metadata": {
    "completion": {"strategy": "zero_fill"},
    "marginals_from_2d": false,
    "distribution_keys": {"rules": ["R1"], "halogens": ["Cl", "F"]}
  }
}
```

#### Distribute Mode
Input with only one marginal:
```json
{
  "version": "1",
  "by_rule": {"R1": {"isotope_unavailable": 5}},
  "total": {"isotope_unavailable": 5, "attempts": 100},
  "metadata": {"rules": ["R1"], "halogens": ["F", "Cl"]}
}
```

Output derives 2D from by_rule marginal and creates consistent by_halogen:
```json
{
  "version": "1",
  "total": {"isotope_unavailable": 5, "attempts": 100},
  "by_rule": {"R1": {"isotope_unavailable": 5}},
  "by_halogen": {"F": {"isotope_unavailable": 2}, "Cl": {"isotope_unavailable": 3}},
  "by_rule_halogen": {"R1": {"F": {"isotope_unavailable": 2}, "Cl": {"isotope_unavailable": 3}}},
  "metadata": {
    "completion": {"strategy": "distribute", "base": "by_rule"},
    "marginals_from_2d": true,
    "distribution_keys": {"rules": ["R1"], "halogens": ["Cl", "F"]}
  }
}
```

Note: Marginals are derived from 2D distribution (marginals_from_2d: true), and remainder allocation follows lexicographic order.

## Data Compatibility and Behavior Contracts

This section defines explicit contracts for data handling, incompatible input forms, and behavioral guarantees.

### Input Data Validation Contracts

#### QA JSON Version Contracts

**SUPPORTED INPUT FORMS:**
- [OK] Version 1 granular QA JSON with explicit `"version": "1"` field
- [OK] Legacy totals-only JSON without version field (auto-upgraded to v1)
- [OK] Version 2 QA JSON with `"version": "2"` (passthrough mode)

**INCOMPATIBLE INPUT FORMS:**
- [NO] QA JSON with `"version": "0"` or other non-standard version strings
- [NO] QA JSON with version field containing non-string values (integers, null, etc.)
- [NO] Non-JSON input files (XML, CSV, plain text)
- [NO] JSON with missing both `total` and all granular dimensions (`by_rule`, `by_halogen`, `by_rule_halogen`)

**BEHAVIOR ASSERTIONS:**
- Invalid version strings trigger immediate validation error with clear error message
- Missing `total` with empty granular dimensions defaults all metrics to 0
- Non-JSON input files cause parser error before processing begins

#### Granular Structure Validation Contracts

**REQUIRED FIELDS:**
- `metadata.rules` array must be non-empty list of strings
- `metadata.halogens` array must be non-empty list of strings
- All metric values must be non-negative integers (negative values trigger data quality warnings)

**INCOMPATIBLE STRUCTURES:**
- [NO] Granular dimensions with non-string keys (numeric rule IDs, boolean halogen flags)
- [NO] Metric values as strings, floats, or non-numeric types
- [NO] Nested granular structures beyond 2D (3D+ tensor formats not supported)
- [NO] Circular references in metadata (self-referencing objects)

**SANITIZATION CONTRACTS:**
- Non-integer metric values are coerced to integers with data quality warnings
- Coercion method is recorded in `metadata.value_coercion_method` (default: "trunc")
- `value_coercion_method: "trunc"` means truncate toward zero (Python `int()`)
- Warning objects include `coercion_method` field for transparency
- Unknown metrics not in distributable set are removed with structured warnings
- Unknown dimension keys not in metadata arrays are dropped with warnings
- Negative values are preserved but flagged in data quality warnings

**WARNING STRUCTURE EXAMPLES:**

1. **Non-integer value detected** (fields in `where` object):
```json
{
  "type": "non_integer_value_detected",
  "where": {
    "dimension": "by_rule",
    "key": "R1",
    "metric": "products"
  },
  "value": 2.7,
  "coerced_to_int": true
}
```

2. **Unknown metric dropped** (fields at top level, metrics array):
```json
{
  "type": "unknown_metric_dropped",
  "dimension": "by_rule",
  "key": "R1",
  "metrics": ["unknown_metric_a", "unknown_metric_b"]
}
```

3. **Marginal conflict detected** (different structure with delta):
```json
{
  "type": "marginal_conflict",
  "base": "by_rule",
  "tolerance": 2,
  "delta": {"R1": 1, "R2": -1}
}
```

### Processing Mode Contracts

#### Base Selection Deterministic Behavior

**GUARANTEED SELECTION ORDER:**
1. When both `by_rule` and `by_halogen` exist: `by_rule` is ALWAYS selected as base
2. When only `by_rule` exists: `by_rule` is selected as base
3. When only `by_halogen` exists: `by_halogen` is selected as base
4. When neither exists: No base selected (zero-fill or error depending on mode)

**CROSS-SIDE FALLBACK CONTRACT:**
- [NO] Metrics present only in non-base marginal are NOT used for 2D distribution
- [NO] Missing metrics in base marginal do NOT pull values from alternate marginal
- [OK] Missing metrics in base marginal default to 0 in 2D structure
- [OK] Conflict warnings are generated for metrics that differ between marginals

#### 2D Structure Metrics Completeness Guarantees

**METRIC KEY CONSISTENCY:**
- [OK] All 2D cells contain identical metric key sets (complete coverage)
- [OK] Missing metrics in any cell are set to 0 (no sparse representation)
- [OK] Metric keys are union of distributable metrics present across input dimensions after sanitization
- [OK] All cells follow same integer coercion method (metadata.value_coercion_method)

#### Completion Mode Behavioral Guarantees

**ZERO_FILL MODE (Default):**
- Input marginal values are NEVER modified
- Missing dimensions are populated with all-zero structures
- `marginals_from_2d` metadata field is set to `false`
- Total calculation uses multi-path fallback but preserves input fidelity

**DISTRIBUTE MODE:**
- Input marginal values MAY be overwritten by 2D-derived values
- Missing dimensions are computed from base marginal distribution
- `marginals_from_2d` metadata field is set to `true`
- 2D structure becomes single source of truth for marginal consistency

### Error Handling and Degradation Contracts

#### RDKit Unavailability Behavior

**WHEN RDKIT IS MISSING:**
- System continues with limited functionality (no structural analysis)
- Molecular canonicalization falls back to input SMILES preservation
- InChIKey generation is skipped, parent keys use SMILES format
- Warning messages are generated but processing completes

**WHEN RDKIT IS AVAILABLE BUT FAILS:**
- Individual molecule failures are logged but don't halt batch processing
- Failed molecules get fallback parent keys (`SMI:UNKNOWN`)
- Enumeration attempts continue for processable molecules
- Comprehensive error reporting in output logs

#### File I/O Error Contracts

**INPUT FILE ERRORS:**
- Missing configuration files trigger immediate CLI error with usage help
- Corrupted SMILES files skip invalid entries with warning logs
- Permission denied errors halt execution with clear filesystem error messages
- Empty input files generate warnings but produce valid empty outputs

**OUTPUT FILE ERRORS:**
- Existing output files are overwritten without confirmation
- Directory creation failures halt execution with permission error details
- Disk space exhaustion during write triggers cleanup and error reporting
- Partial output files from interrupted runs are left in place for debugging

### ASCII Compliance Contracts

#### Character Set Restrictions

**ENFORCED ASCII-ONLY:**
- All source code files must contain only ASCII characters (0-127)
- Configuration YAML files must use ASCII-only keys and string values
- Documentation and README files must be ASCII-compliant
- Output JSON files use ASCII-safe encoding for all text content

**INCOMPATIBLE NON-ASCII INPUT:**
- [NO] Unicode characters in SMILES input trigger encoding warnings
- [NO] Non-ASCII characters in configuration files cause YAML parse errors
- [NO] UTF-8 encoded source files fail pre-commit ASCII validation checks

### Output Format Guarantees

#### JSON Schema Stability

**GUARANTEED FIELDS (V1 Granular):**
- `version` field is always `"1"` for granular outputs
- `metadata.completion.strategy` documents processing mode used
- `metadata.marginals_from_2d` indicates marginal derivation source
- `metadata.warnings` array contains structured data quality warnings

**SCHEMA EVOLUTION CONTRACTS:**
- New optional fields may be added to metadata without version bump
- Existing field names and types remain stable within major versions
- Removal of fields requires major version increment
- Field semantics remain consistent across point releases

#### CSV Output Determinism

**COLUMN ORDER GUARANTEES:**
- Core columns appear in documented fixed order
- Extension columns (prefixed `ours_`, `dl_`) appear after core columns
- Empty columns are included to maintain schema consistency
- Row order follows lexicographic sorting for deterministic output

**DATA TYPE CONTRACTS:**
- Numeric columns contain only integers or floats (no strings)
- String columns are properly quoted and escaped for CSV compliance
- Boolean columns use literal `true`/`false` values
- Missing values are represented as empty cells, not null strings

