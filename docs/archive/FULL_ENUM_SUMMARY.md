# Full-Scale NP Halogenation Library Generation Summary

**Generated:** 2025-12-13 01:30:27

## Overview

- **Classes Processed:** 1
- **Total Enumeration Jobs:** 1
- **Successful Jobs:** 0
- **Failed Jobs:** 1

## Results by Class and K

| Class | k | Parents | Products | Unique | Prod/Parent | Status |
|-------|---|---------|----------|--------|-------------|--------|
| terpenoid | 2 | - | - | - | - | ❌ INFO - halogenator.class_confi |

## Total Statistics

## Output Locations

Generated libraries are stored in: `data/output/nplike/<class>-<k>X/products.parquet`

## Next Steps

1. Perform QC on generated libraries using 09_visualize_library.py
2. Run 05_summaries.py to generate descriptor statistics
3. Export to virtual screening formats using 11_export_for_vs.py
