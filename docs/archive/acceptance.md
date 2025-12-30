# -*- coding: ascii -*-

PR1 Sugar Masking and PR2 Sites: Notes
======================================

Overview
--------
- PR1 provides sugar masking to prevent combinatorial explosion.
- P1-1 adds an evidence-based scoring for sugar ring detection.
- PR2 adds R2a (sp2 CH in C-ring) and R2b (sp3 CH2 in flavanone C-ring) site rules.
- All outputs, logs, and docs are ASCII-only.

Sugar Ring Scoring (P1-1)
-------------------------
- Score factors:
  - Ring size: +2 (6-ring), +1 (5-ring)
  - Exactly one ring oxygen: +3 (required)
  - No inner carbonyl (C in ring double-bonded to O): +2 (required)
  - SP3 ratio on ring carbons: linear bonus up to +3
  - Exocyclic O count (single bonds only): +1 each, capped at +4
  - C-glycoside-like topology: +2
- Threshold: `sugar_ring_score_threshold` (default 8.0). Rings meeting the threshold are accepted.
- Observations written per sample in `*-full.json`:
  - `accepted_via_score` (bool)
  - `accepted_ring_size` (5 or 6)
  - `accepted_ring_score` (float)
  - `exocyclic_O_count_single_bond` (int)
  - `has_cglyco_evidence` (bool)
- Global counters:
  - `accepted_via_score_count`
  - `degraded_due_to_no_score_count`

RDKit Seeding and Logging
-------------------------
- `RDKitSeedManager` centralizes seeding across RDKit subsystems.
- Debounced logging:
  - Always logs one INFO summary: `RDKit seeding summary: seeded=N, failed=M`.
  - DEBUG lists print only once or when the sets change.
  - Set `HALO_RDKitSeedVerbose=1` to force detailed DEBUG lists each call.
- `enumerate_k` logs one INFO summary for Python/NumPy seeds.
- `rdkit_seed_report` is included in acceptance metadata for audit.

PR2 R2b Consistency With Sugar Masking
--------------------------------------
- R2b targets strictly CH2 in C-ring:
  - Carbon, SP3, TotalNumHs==2, in 6-member ring.
- Consistency policy:
  - If masking exists, exclude masked atoms.
  - If masking not available but sugar config is, exclude CH2 that belong to rings accepted by sugar scoring.
  - Otherwise, do not exclude oxygen-containing rings by default.

Usage Notes
-----------
- Run full acceptance:
  - `python scripts/run_pr1_acceptance.py --verbose`
- Verify reproducibility of the two latest runs:
  - `python scripts/verify_reproducibility.py --reuse-latest 2`
- ASCII checks:
  - If available: `python scripts/check_ascii_portable.py`
  - Fallback: `grep -nP "[^\x00-\x7F]" -R src scripts | wc -l`

Configuration Tips
------------------
- Adjust scoring via `sugar_ring_score_threshold` and SP3 thresholds:
  - `sp3_threshold_5_ring` (default 0.4)
  - `sp3_threshold_6_ring` (default 0.5)
- RDKit seeding verbosity:
  - `HALO_RDKitSeedVerbose=1` for detailed DEBUG lists every call.

Limitations
-----------
- Scoring-based sugar detection is heuristic; always validate with domain samples.
- RDKit availability and version may affect internal seed control coverage.
