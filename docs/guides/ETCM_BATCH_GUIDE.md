# ETCM 2000 Flavonoids Batch Generation Guide

This guide provides step-by-step instructions for processing the ETCM library of 2000 flavonoids through halogenation enumeration.

## Prerequisites

1. **ETCM Data File**: Ensure `flavonoid_smiles_data_ETCM.xlsx` is in the project root directory
2. **Python Dependencies**: Install required packages (pandas, openpyxl, RDKit, etc.)
3. **Disk Space**: Ensure adequate space (~5-10GB) for products and QA output
4. **Memory**: Recommend 8GB+ RAM for full batch processing

## Quick Start

### Windows Users
```batch
# Run the automated batch script
scripts\run_etcm_batch.bat
```

### Linux/Mac Users
```bash
# Make script executable and run
chmod +x scripts/run_etcm_batch.sh
./scripts/run_etcm_batch.sh
```

## Manual Step-by-Step Process

### Step 1: Data Preparation

Convert ETCM Excel data to SMILES format:

```bash
# Full dataset (2000+ molecules)
python scripts/etcm_to_parents.py -i flavonoid_smiles_data_ETCM.xlsx -o data/input/etcm_flavonoids.smi

# For testing: sample first 50 molecules
python scripts/etcm_to_parents.py -i flavonoid_smiles_data_ETCM.xlsx -o data/input/etcm_flavonoids.smi --sample 50
```

### Step 2: k=2 Enumeration

Run halogenation enumeration with k=2 (di-substituted products):

```bash
python -m halogenator.cli enum \
  -c configs/etcm_flavonoids.yaml \
  --subset flavonoids \
  --outdir data/output/etcm2000 \
  --qa-completion-mode distribute \
  --qa-conflict-tolerance 2 \
  --qa-max-warnings 500
```

### Step 3: Quality Assessment

Review the generated QA summary:

```bash
# View QA statistics
cat data/output/etcm2000/qa_summary.json

# Check product statistics
python -c "
import pandas as pd
df = pd.read_parquet('data/output/etcm2000/products_k2.parquet')
print(f'Total products: {len(df)}')
print(f'Unique parents with products: {df.parent_key.nunique()}')
print(f'Products per parent (avg): {len(df) / df.parent_key.nunique():.1f}')
"
```

### Step 4: Audit Sample Review

Manually review a sample of generated products:

```bash
# Open audit sample in spreadsheet software
# File: data/output/etcm2000/audit/sample_100.csv
#
# Review checklist:
# - Chemical structures look reasonable
# - Substitution patterns follow expected rules
# - No obvious problematic structures (unrealistic bonds, etc.)
# - SMILES strings are valid
```

### Step 5: Extended k=3 Enumeration (Optional)

If k=2 results are satisfactory, run k=3 for tri-substituted products:

```bash
# Create k=3 output directory
mkdir -p data/output/etcm2000/k3

# Run k=3 enumeration (warning: much larger output)
python -m halogenator.cli enum \
  -c configs/etcm_flavonoids_k3.yaml \
  --subset flavonoids \
  --outdir data/output/etcm2000/k3 \
  --qa-completion-mode distribute \
  --qa-conflict-tolerance 3 \
  --qa-max-warnings 1000
```

## Output Structure

After successful completion, expect the following directory structure:

```
data/output/etcm2000/
??? products_k2.parquet          # Main products table
??? products_k2.csv              # CSV version of products
??? summary_k2.csv               # Fixed-schema summary
??? qa_summary.json              # Comprehensive QA analysis
??? pivot_rule_halogen_k.csv     # Dimensional analysis
??? audit/
?   ??? sample_100.csv           # Sample for manual review
?   ??? manual_review_notes.csv  # Your review notes (create)
??? k3/ (if k=3 run)
    ??? products_k3.parquet
    ??? summary_k3.csv
    ??? qa_summary.json
```

## Configuration Details

### `configs/etcm_flavonoids.yaml`
- **Purpose**: k=2 enumeration configuration
- **Target**: di-substituted flavonoid products
- **Constraints**: Balanced speed vs. completeness
- **QA Settings**: Comprehensive analysis with moderate warnings limit

### `configs/etcm_flavonoids_k3.yaml`
- **Purpose**: Extended k=3 enumeration
- **Target**: tri-substituted products
- **Constraints**: More restrictive to manage combinatorial explosion
- **QA Settings**: Higher tolerance for statistical noise

## Performance Expectations

### k=2 Enumeration (2000 molecules)
- **Runtime**: 30-60 minutes (depending on hardware)
- **Products**: ~20,000 - 50,000 (10-25 per parent on average)
- **Memory**: 4-8GB peak usage
- **Disk**: 2-5GB output

### k=3 Enumeration (2000 molecules)
- **Runtime**: 2-6 hours (significantly longer)
- **Products**: ~100,000 - 200,000 (50-100 per parent on average)
- **Memory**: 8-16GB peak usage
- **Disk**: 10-20GB output

## Troubleshooting

### Common Issues

1. **Memory Errors**
   - Reduce batch size by processing subsets
   - Enable state signature pruning
   - Consider running on machine with more RAM

2. **Timeout Errors**
   - Increase timeout settings in configuration
   - Enable fail-fast mode to skip problematic molecules
   - Review molecules causing timeouts

3. **Disk Space Issues**
   - Enable compression in k=3 configuration
   - Process in smaller batches
   - Clean up intermediate files

4. **Quality Issues**
   - Review QA warnings for systematic problems
   - Adjust constraints if too many invalid products
   - Consider stricter PAINS filtering

### Performance Optimization

1. **For Large Batches**
   - Enable checkpointing
   - Use symmetry folding and state signatures
   - Process in chunks if memory limited

2. **For Quality over Quantity**
   - Increase min_graph_distance
   - Enable stricter QC settings
   - Lower per_ring_quota

3. **For Speed**
   - Disable state signature (k=2 only)
   - Increase timeout thresholds
   - Reduce QA analysis depth

## Validation and Quality Control

### Automated Checks
- SMILES validity
- Molecular weight constraints
- PAINS filter compliance
- Structural sanitization

### Manual Review Guidelines
1. Review audit sample (sample_100.csv)
2. Check 5-10 products per parent molecule
3. Verify substitution patterns match expected rules
4. Look for obvious chemical errors
5. Document findings in manual_review_notes.csv

### Success Criteria
- \>95% of parents generate at least one valid product
- QA conflict warnings <10% of total products
- Manual audit shows <5% problematic structures
- Total runtime within expected range

## Next Steps After Completion

1. **Data Analysis**: Use products table for cheminformatics analysis
2. **Database Storage**: Import to chemical database for searching
3. **Comparative Analysis**: Compare with existing halogenated compound libraries
4. **Property Prediction**: Run ADMET predictions on product set
5. **Synthetic Planning**: Prioritize products for synthesis based on properties

## Support and Documentation

- Main documentation: README.md
- Configuration reference: configs/ directory
- Test examples: tests/ directory
- CLI help: `python -m halogenator.cli --help`