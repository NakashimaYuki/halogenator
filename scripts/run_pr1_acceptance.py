#!/usr/bin/env python
# -*- coding: ascii -*-
"""
PR1 Sugar Masking Acceptance Test Runner

This script runs acceptance tests for PR1 sugar masking functionality by:
1. Testing 5 glycoside samples + 6 non-sugar controls (11 total)
2. Running k=2 enumeration with sugar.mode=off vs heuristic
3. Collecting statistics on product counts, QA metrics, and symmetry
4. Generating acceptance report for validation

The script validates that:
- Glycosides show significant product reduction with sugar masking
- Non-sugar controls show minimal difference (including lactone/ester filtering)
- QA metrics properly track sugar filtering events
"""

import argparse
import json
import logging
import os
import subprocess
import sys
import hashlib
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import enumerate_products, enumerate_with_stats, EnumConfig
from halogenator.chem_compat import Chem
from halogenator.rdkit_seed_utils import set_rdkit_random_seed

LOG = logging.getLogger(__name__)


def get_git_info() -> Dict[str, str]:
    """Get git repository information."""
    try:
        # Get git commit SHA
        sha_result = subprocess.run(['git', 'rev-parse', 'HEAD'],
                                  capture_output=True, text=True, cwd=Path(__file__).parent.parent)
        git_sha = sha_result.stdout.strip() if sha_result.returncode == 0 else 'unknown'

        # Get short SHA
        short_sha = git_sha[:8] if git_sha != 'unknown' else 'unknown'

        # Get branch name
        branch_result = subprocess.run(['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                                     capture_output=True, text=True, cwd=Path(__file__).parent.parent)
        branch = branch_result.stdout.strip() if branch_result.returncode == 0 else 'unknown'

        # Check if there are uncommitted changes
        status_result = subprocess.run(['git', 'status', '--porcelain'],
                                     capture_output=True, text=True, cwd=Path(__file__).parent.parent)
        has_changes = bool(status_result.stdout.strip()) if status_result.returncode == 0 else True

        return {
            'git_sha': git_sha,
            'git_short_sha': short_sha,
            'git_branch': branch,
            'has_uncommitted_changes': has_changes,
            'git_available': True
        }
    except Exception as e:
        LOG.warning(f"Failed to get git info: {e}")
        return {
            'git_sha': 'unknown',
            'git_short_sha': 'unknown',
            'git_branch': 'unknown',
            'has_uncommitted_changes': True,
            'git_available': False
        }


def calculate_file_hash(file_path: str) -> str:
    """Calculate SHA1 hash of file contents for traceability."""
    try:
        with open(file_path, 'rb') as f:
            return hashlib.sha1(f.read()).hexdigest()[:16]  # First 16 chars for brevity
    except Exception:
        return 'unknown'


def load_samples(samples_file: str) -> List[Dict[str, Any]]:
    """Load test samples from JSON file."""
    with open(samples_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data['samples']


def run_enumeration(smiles: str, sample_id: str, sugar_mode: str, sample_type: str = None, fast_mode: bool = False) -> Dict[str, Any]:
    """
    Run k=2 enumeration for a single sample with specified sugar mode.

    Returns:
        Dictionary with enumeration results and statistics including detailed sugar masking data
    """
    LOG.info(f"Running enumeration for {sample_id} with sugar.mode={sugar_mode}")

    try:
        # Determine proximity guard radius: only enable for glycoside samples with heuristic mode
        is_glycoside = (sample_type == 'glycoside')
        proximity_radius = 3 if (sugar_mode == 'heuristic' and is_glycoside) else 0

        # T28: Unified proximity guard logging for production troubleshooting
        proximity_enabled = proximity_radius > 0
        LOG.info(f"Proximity guard config: {{radius={proximity_radius}, sample_type={sample_type}, sugar_mode={sugar_mode}, enabled={proximity_enabled}}}")

        # H2-A/H2-E: Configure enumeration with optional fast mode
        if fast_mode:
            # Fast mode: k=1, fewer halogens, reduced rules for CI/quick validation
            config = EnumConfig(
                k_max=1,
                halogens=('F', 'Cl'),  # Fewer halogens for speed
                rules=('R1', 'R2'),  # Subset of rules for quick testing
                sugar_cfg={'mode': sugar_mode, 'audit': True, 'proximity_guard_radius': proximity_radius},
                symmetry_cfg={'enabled': True},
                random_seed=12345,
                constraints={'per_ring_quota': 1, 'min_graph_distance': 2, 'max_per_halogen': 2}  # Reduce products
            )
        else:
            # Full mode: k=2 and full rule set for comprehensive testing
            config = EnumConfig(
                k_max=2,
                halogens=('F', 'Cl', 'Br', 'I'),
                rules=('R1', 'R2', 'R3', 'R4', 'R5'),  # Use full rule set for comprehensive testing
                sugar_cfg={'mode': sugar_mode, 'audit': True, 'proximity_guard_radius': proximity_radius},
                symmetry_cfg={'enabled': True},
                random_seed=12345  # Fixed seed for reproducible results
            )

        # Run enumeration with stats collection
        products, qa_stats = enumerate_with_stats(smiles, config)
        pivots = qa_stats.get('pivots', {})

        # Collect detailed sugar masking information
        sugar_mask_info = {}
        if sugar_mode == 'heuristic':
            # Import sugar masking functions to get detailed mask info
            from halogenator.sugar_mask import get_sugar_mask_with_full_status, compute_sugar_audit_fields, _get_default_sugar_cfg
            from halogenator.chem_compat import Chem

            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Use default sugar configuration and capture it for reporting
                sugar_cfg = _get_default_sugar_cfg()
                mask, degraded, status_metadata = get_sugar_mask_with_full_status(mol, mode='heuristic', sugar_cfg=sugar_cfg)

                # Compute audit data directly instead of reading from qa_stats
                audit_data = compute_sugar_audit_fields(mol, mask, sugar_cfg)

                # T2-4: Compute R2a/R2b site counts after masking
                from halogenator.sites import c_ring_sp2_CH_sites, c_ring_sp3_CH2_flavanone_sites
                r2a_sites_after_mask = len(c_ring_sp2_CH_sites(mol, masked_atoms=mask))
                r2b_sites_after_mask = len(c_ring_sp3_CH2_flavanone_sites(mol, masked_atoms=mask, sugar_cfg=sugar_cfg))

                # T2-5: Log R2 summary at sample level to avoid repetition during BFS
                LOG.info("R2 summary: r2a_sites=%d, r2b_sites=%d (after mask)",
                         r2a_sites_after_mask, r2b_sites_after_mask)

                sugar_mask_info = {
                    'mask_size': len(mask),
                    'mask_atoms': sorted(list(mask)) if mask else [],
                    'degraded': degraded,
                    'masked_oh_count': audit_data.get('masked_oh_count', 0),
                    'masked_bridge_o_count': audit_data.get('masked_bridge_o_count', 0),
                    'sugar_rings_count': len(audit_data.get('sugar_rings', [])),
                    'sugar_rings': audit_data.get('sugar_rings', []),
                    'is_c_glycoside_like': audit_data.get('is_c_glycoside_like', False),
                    # P0-H1: Include parameterized configuration values for traceability
                    'sugar_cfg': sugar_cfg,
                    # Observations
                    'accepted_via_score': audit_data.get('accepted_via_score', False),
                    'accepted_ring_size': audit_data.get('accepted_ring_size', 0),
                    'accepted_ring_score': audit_data.get('accepted_ring_score', 0.0),
                    'exocyclic_O_count_single_bond': audit_data.get('exocyclic_O_count_single_bond', 0),
                    'has_cglyco_evidence': audit_data.get('has_cglyco_evidence', False),
                    # T2-4: R2a/R2b site counts after masking
                    'r2a_sites': r2a_sites_after_mask,
                    'r2b_sites': r2b_sites_after_mask
                }

                # Add degraded masking detailed metadata
                if degraded and status_metadata:
                    sugar_mask_info.update({
                        'degraded_reason': status_metadata.get('degraded_reason', []),
                        'fallback_ring_count': status_metadata.get('fallback_ring_count', 0),
                        'cglyco_candidates_count': status_metadata.get('cglyco_candidates_count', 0)
                    })
        else:
            sugar_mask_info = {
                'mask_size': 0,
                'mask_atoms': [],
                'degraded': False,
                'masked_oh_count': 0,
                'masked_bridge_o_count': 0,
                'sugar_rings_count': 0,
                'sugar_rings': [],
                'is_c_glycoside_like': False
            }

        # Build by_rule_halogen_counts from products list (T2-4) - improved robustness
        by_rule_halogen_counts = {}
        if sugar_mode == 'heuristic':
            # Select statistics source (consistent with current sample strategy)
            src_products = products  # products from enumerate_with_stats

            # Derive halogen keys from products, fallback to default if none
            used_halogens = sorted({p.get("halogen") for p in src_products if p.get("halogen")})
            if not used_halogens:
                used_halogens = ["F", "Cl", "Br", "I"]

            by_rule_halogen_counts = {
                "R2a": {h: 0 for h in used_halogens},
                "R2b": {h: 0 for h in used_halogens},
            }

            for p in src_products:
                rule = p.get("rule")
                hal = p.get("halogen")
                if rule in by_rule_halogen_counts and hal in by_rule_halogen_counts[rule]:
                    by_rule_halogen_counts[rule][hal] += 1

        # Extract key statistics
        stats = {
            'sample_id': sample_id,
            'sugar_mode': sugar_mode,
            'smiles': smiles,
            'total_products': len(products),
            'qa_stats': qa_stats,
            'pivots': pivots,
            'products': products,  # Include products for analysis
            'sugar_mask_detail': sugar_mask_info  # Detailed sugar masking data
        }

        # Add by_rule_halogen_counts for heuristic mode
        if sugar_mode == 'heuristic':
            stats['by_rule_halogen_counts'] = by_rule_halogen_counts

        # Extract specific QA metrics of interest
        qa_paths = stats['qa_stats'].get('qa_paths', {})
        stats['sugar_mask_filtered'] = int(qa_paths.get('sugar_mask_filtered', 0))
        stats['sugar_proximity_filtered'] = int(qa_paths.get('sugar_proximity_filtered', 0))
        stats['post_guard_blocked'] = int(qa_paths.get('post_guard_blocked', 0))

        # Backward compatibility: merge legacy keys if present
        legacy_post_guard = int(qa_paths.get('sugar_post_guard_blocked', 0))
        stats['post_guard_blocked'] = max(stats['post_guard_blocked'], legacy_post_guard)

        stats['sugar_mask_degraded'] = qa_paths.get('sugar_mask_degraded', 0)
        stats['sugar_events'] = {
            'total': stats['sugar_mask_filtered'] + stats['sugar_proximity_filtered'] + stats['post_guard_blocked'],
            'sugar_mask_filtered': stats['sugar_mask_filtered'],
            'sugar_proximity_filtered': stats['sugar_proximity_filtered'],
            'post_guard_blocked': stats['post_guard_blocked'],
        }

        LOG.info(
            "Sugar events: mask_filtered=%d, proximity_filtered=%d, post_guard_blocked=%d, total=%d",
            stats['sugar_mask_filtered'], stats['sugar_proximity_filtered'],
            stats['post_guard_blocked'], stats['sugar_events']['total']
        )

        # Extract core enumeration metrics
        stats['attempts'] = stats['qa_stats'].get('attempts', 0)
        stats['template_unsupported'] = stats['qa_stats'].get('template_unsupported', 0)
        stats['no_product_matches'] = stats['qa_stats'].get('no_product_matches', 0)

        # Extract deduplication statistics if available
        dedup_stats = stats['qa_stats'].get('dedup', {})
        stats['dedup_hits_inchi'] = dedup_stats.get('inchi_hits', 0)
        stats['dedup_hits_statesig'] = dedup_stats.get('statesig_hits', 0)
        stats['dedup_total_hits'] = dedup_stats.get('total_hits', 0)

        LOG.info(f"  {sample_id} ({sugar_mode}): {stats['total_products']} products, "
                f"{stats['attempts']} attempts, "
                f"{stats['sugar_mask_filtered']} sugar filtered")

        return stats

    except Exception as e:
        LOG.error(f"Enumeration failed for {sample_id} with {sugar_mode}: {e}")
        return {
            'sample_id': sample_id,
            'sugar_mode': sugar_mode,
            'smiles': smiles,
            'error': str(e),
            'total_products': 0,
            'qa_stats': {},
            'pivots': {},
            'products': []
        }


def validate_p0_h2_assertions(off_stats: Dict, heuristic_stats: Dict, sample_type: str, sample_id: str) -> Dict[str, Any]:
    """
    P0-H2: Validate specific regression test assertions for expanded sample set.

    Positive criteria (glycosides): sugar_events>0 OR attempts_delta>=5 OR reduction>=30%
    Negative criteria (aglycones): reduction<=5% AND mask_size=0 AND degraded=False

    Args:
        off_stats: Statistics from sugar.mode=off
        heuristic_stats: Statistics from sugar.mode=heuristic
        sample_type: 'glycoside' or 'aglycone'
        sample_id: Sample identifier for specific validation rules

    Returns:
        P0-H2 assertion validation results
    """
    off_products = off_stats['total_products']
    heuristic_products = heuristic_stats['total_products']
    off_attempts = off_stats.get('attempts', 0)
    heuristic_attempts = heuristic_stats.get('attempts', 0)

    reduction_pct = (off_products - heuristic_products) / off_products * 100 if off_products > 0 else 0.0
    attempts_delta = off_attempts - heuristic_attempts

    # Extract sugar masking details
    mask_detail = heuristic_stats.get('sugar_mask_detail', {})
    mask_size = mask_detail.get('mask_size', 0)
    degraded = mask_detail.get('degraded', False)
    sugar_events = heuristic_stats.get('sugar_mask_filtered', 0) + heuristic_stats.get('sugar_proximity_filtered', 0) + heuristic_stats.get('sugar_post_guard_blocked', 0)

    p0_h2_result = {
        'sample_id': sample_id,
        'sample_type': sample_type,
        'reduction_pct': reduction_pct,
        'attempts_delta': attempts_delta,
        'mask_size': mask_size,
        'degraded': degraded,
        'sugar_events': sugar_events,
        'assertion_pass': False,
        'assertion_criteria_met': [],
        'assertion_criteria_failed': []
    }

    if sample_type == 'glycoside':
        # Positive criteria: ANY of these conditions must be true
        criteria_sugar_events = sugar_events > 0
        criteria_attempts_delta = attempts_delta >= 5
        criteria_reduction = reduction_pct >= 30.0

        if criteria_sugar_events:
            p0_h2_result['assertion_criteria_met'].append('sugar_events>0')
        else:
            p0_h2_result['assertion_criteria_failed'].append('sugar_events>0')

        if criteria_attempts_delta:
            p0_h2_result['assertion_criteria_met'].append('attempts_delta>=5')
        else:
            p0_h2_result['assertion_criteria_failed'].append('attempts_delta>=5')

        if criteria_reduction:
            p0_h2_result['assertion_criteria_met'].append('reduction>=30%')
        else:
            p0_h2_result['assertion_criteria_failed'].append('reduction>=30%')

        # Pass if ANY criterion is met
        p0_h2_result['assertion_pass'] = criteria_sugar_events or criteria_attempts_delta or criteria_reduction

    else:  # aglycone
        # Negative criteria: ALL of these conditions must be true
        criteria_reduction = reduction_pct <= 5.0
        criteria_mask_size = mask_size == 0
        criteria_degraded = degraded == False

        if criteria_reduction:
            p0_h2_result['assertion_criteria_met'].append('reduction<=5%')
        else:
            p0_h2_result['assertion_criteria_failed'].append('reduction<=5%')

        if criteria_mask_size:
            p0_h2_result['assertion_criteria_met'].append('mask_size=0')
        else:
            p0_h2_result['assertion_criteria_failed'].append('mask_size=0')

        if criteria_degraded:
            p0_h2_result['assertion_criteria_met'].append('degraded=False')
        else:
            p0_h2_result['assertion_criteria_failed'].append('degraded=False')

        # Pass if ALL criteria are met
        p0_h2_result['assertion_pass'] = criteria_reduction and criteria_mask_size and criteria_degraded

    return p0_h2_result


def analyze_comparison(off_stats: Dict, heuristic_stats: Dict, sample_type: str) -> Dict[str, Any]:
    """
    Analyze the comparison between off and heuristic modes with multi-criteria assessment.

    Args:
        off_stats: Statistics from sugar.mode=off
        heuristic_stats: Statistics from sugar.mode=heuristic
        sample_type: 'glycoside' or 'aglycone' for appropriate PASS criteria

    Returns:
        Analysis results with reduction percentages and multi-criteria assessment
    """
    off_products = off_stats['total_products']
    heuristic_products = heuristic_stats['total_products']
    off_attempts = off_stats.get('attempts', 0)
    heuristic_attempts = heuristic_stats.get('attempts', 0)

    if off_products == 0:
        reduction_pct = 0.0
    else:
        reduction_pct = (off_products - heuristic_products) / off_products * 100

    # Extract sugar masking events
    sugar_filtered = heuristic_stats.get('sugar_mask_filtered', 0)
    sugar_proximity_filtered = heuristic_stats.get('sugar_proximity_filtered', 0)
    sugar_blocked = heuristic_stats.get('sugar_post_guard_blocked', 0)
    sugar_degraded = heuristic_stats.get('sugar_mask_degraded', 0)

    # Multi-criteria assessment based on sample type
    pass_criteria = []
    overall_pass = False

    if sample_type == 'glycoside':
        # Glycoside PASS criteria (any one satisfies)

        # Criterion A: Product reduction >= 30%
        product_reduction_pass = reduction_pct >= 30.0
        if product_reduction_pass:
            pass_criteria.append("product_reduction_30pct")

        # Criterion B: QA events indicate sugar masking activity
        qa_events_pass = sugar_filtered > 0 or sugar_blocked > 0
        if qa_events_pass:
            pass_criteria.append("qa_events_detected")

        # Criterion C: Attempts reduction >= 5
        attempts_reduction = off_attempts - heuristic_attempts
        attempts_reduction_pass = attempts_reduction >= 5
        if attempts_reduction_pass:
            pass_criteria.append("attempts_reduction_5plus")

        overall_pass = len(pass_criteria) > 0

    else:  # aglycone/control
        # Control PASS criteria: minimal change expected

        # Criterion: Small product reduction (<= 5%) OR small attempts difference (<= 3)
        small_product_change = reduction_pct <= 5.0
        small_attempts_change = abs(off_attempts - heuristic_attempts) <= 3

        if small_product_change:
            pass_criteria.append("small_product_change_5pct")
        if small_attempts_change:
            pass_criteria.append("small_attempts_change_3")

        overall_pass = small_product_change or small_attempts_change

    analysis = {
        'off_products': off_products,
        'heuristic_products': heuristic_products,
        'off_attempts': off_attempts,
        'heuristic_attempts': heuristic_attempts,
        'reduction_absolute': off_products - heuristic_products,
        'reduction_percentage': reduction_pct,
        'attempts_reduction': off_attempts - heuristic_attempts,
        'significant_reduction': reduction_pct >= 30.0,  # Keep legacy field
        'pass_criteria_met': pass_criteria,
        'overall_pass': overall_pass,
        'assessment_type': sample_type,
        'sugar_events': {
            'filtered_matches': sugar_filtered,
            'proximity_filtered': sugar_proximity_filtered,
            'post_guard_blocked': sugar_blocked,
            'mask_degraded': sugar_degraded,
            'total_sugar_events': sugar_filtered + sugar_proximity_filtered + sugar_blocked
        }
    }

    return analysis


def main():
    """Main acceptance test runner."""
    # H2-E: Parse command line arguments for fast mode support
    parser = argparse.ArgumentParser(description='PR1 Sugar Masking Acceptance Test Runner')
    parser.add_argument('--fast', action='store_true',
                        help='Enable fast mode (k=1, reduced rules) for CI/quick validation')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable verbose logging')
    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')

    # Paths
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    samples_file = project_root / "tests/data/acceptance/pr1_sugar_samples.json"
    output_file = project_root / "benchmarks/pr1_sugar_acceptance.json"

    # Ensure output directory exists
    output_file.parent.mkdir(exist_ok=True)

    LOG.info("Starting PR1 sugar masking acceptance tests")

    # Set RDKit random seed and capture report for metadata
    rdkit_seed_report = set_rdkit_random_seed(12345)

    # Load test samples and calculate file hash for traceability
    samples = load_samples(samples_file)
    samples_hash = calculate_file_hash(samples_file)
    LOG.info(f"Loaded {len(samples)} test samples from {samples_file.resolve()}")
    LOG.info(f"Samples file hash: {samples_hash}")

    # H2-E: Display mode information
    mode_info = "FAST MODE (k=1, F+Cl only, R1+R2 rules)" if args.fast else "FULL MODE (k=2, all halogens, all rules)"
    LOG.info(f"Running in {mode_info}")

    # Run enumeration for each sample in both modes
    results = {}

    for sample in samples:
        sample_id = sample['id']
        smiles = sample['smiles']
        sample_type = sample['type']

        LOG.info(f"Processing {sample_id} ({sample_type}): {sample['name']}")

        # Run with sugar masking off
        off_stats = run_enumeration(smiles, sample_id, 'off', sample_type, args.fast)

        # Run with sugar masking heuristic
        heuristic_stats = run_enumeration(smiles, sample_id, 'heuristic', sample_type, args.fast)

        # Analyze comparison
        analysis = analyze_comparison(off_stats, heuristic_stats, sample_type)

        # P0-H2: Validate regression test assertions for expanded sample set
        p0_h2_validation = validate_p0_h2_assertions(off_stats, heuristic_stats, sample_type, sample_id)

        results[sample_id] = {
            'sample_info': sample,
            'off_mode': off_stats,
            'heuristic_mode': heuristic_stats,
            'analysis': analysis,
            'p0_h2_assertion': p0_h2_validation
        }

    # Get environment and git information
    git_info = get_git_info()
    current_time = datetime.now()

    # Import sugar configuration for reporting
    try:
        from halogenator.sugar_mask import _get_default_sugar_cfg
        default_sugar_cfg = _get_default_sugar_cfg()
    except ImportError:
        default_sugar_cfg = "Unable to import sugar configuration"

    # H2-A: Generate summary report with comprehensive metadata for full traceability
    summary = {
        'test_metadata': {
            'schema_version': 'p1.2',
            'version': '1.0',
            'timestamp': current_time.isoformat(),
            'timestamp_unix': str(Path(__file__).stat().st_mtime),
            'samples_tested': len(samples),
            'k_max': 1 if args.fast else 2,
            'rules': ['R1', 'R2'] if args.fast else ['R1', 'R2', 'R3', 'R4', 'R5'],
            'git_info': git_info,
            # P0-H1: Include parameterized sugar configuration for traceability
            'sugar_masking_config': default_sugar_cfg,
            # H2-A: Add comprehensive traceability information
            'random_seed': 12345,
            'rdkit_seed_report': rdkit_seed_report,
            'samples_file_hash': samples_hash,
            'halogens': ['F', 'Cl'] if args.fast else ['F', 'Cl', 'Br', 'I'],
            'symmetry_enabled': True,
            # H2-E: Add fast mode metadata
            'fast_mode': args.fast,
            'test_mode': 'fast' if args.fast else 'full',
            'environment': {
                'python_version': sys.version,
                'script_path': str(Path(__file__).resolve()),
                'working_directory': str(Path.cwd()),
                'samples_file': str(samples_file.resolve()),
                'samples_file_hash': samples_hash
            }
        },
        'samples': results,
        'overall_assessment': {}
    }

    # Assess overall results
    glycoside_samples = [s for s in samples if s['type'] == 'glycoside']
    control_samples = [s for s in samples if s['type'] == 'aglycone']

    glycoside_reductions = [results[s['id']]['analysis']['reduction_percentage']
                           for s in glycoside_samples]
    control_reductions = [results[s['id']]['analysis']['reduction_percentage']
                         for s in control_samples]

    # Count samples that pass the multi-criteria assessment
    glycoside_passes = sum(1 for s in glycoside_samples
                          if results[s['id']]['analysis']['overall_pass'])
    control_passes = sum(1 for s in control_samples
                        if results[s['id']]['analysis']['overall_pass'])

    # P0-H2: Count samples that pass the regression test assertions
    glycoside_p0_h2_passes = sum(1 for s in glycoside_samples
                                if results[s['id']]['p0_h2_assertion']['assertion_pass'])
    control_p0_h2_passes = sum(1 for s in control_samples
                              if results[s['id']]['p0_h2_assertion']['assertion_pass'])

    # Observation counters across heuristic runs
    accepted_via_score_count = 0
    degraded_due_to_no_score_count = 0
    for sid, item in results.items():
        heur = item.get('heuristic_mode', {})
        detail = heur.get('sugar_mask_detail', {})
        if detail.get('accepted_via_score', False):
            accepted_via_score_count += 1
        if detail.get('degraded', False) and not detail.get('accepted_via_score', False):
            degraded_due_to_no_score_count += 1

    summary['overall_assessment'] = {
        'glycoside_samples': len(glycoside_samples),
        'glycoside_avg_reduction': sum(glycoside_reductions) / len(glycoside_reductions) if glycoside_reductions else 0,
        'glycoside_passes': glycoside_passes,
        'glycoside_significant_reductions': sum(1 for r in glycoside_reductions if r >= 30),  # Keep legacy
        'control_samples': len(control_samples),
        'control_avg_reduction': sum(control_reductions) / len(control_reductions) if control_reductions else 0,
        'control_passes': control_passes,
        'control_significant_reductions': sum(1 for r in control_reductions if r >= 30),  # Keep legacy
        # P0-H2: Regression test assertion results
        'p0_h2_glycoside_passes': glycoside_p0_h2_passes,
        'p0_h2_control_passes': control_p0_h2_passes,
        'p0_h2_total_passes': glycoside_p0_h2_passes + control_p0_h2_passes,
        'p0_h2_total_samples': len(glycoside_samples) + len(control_samples),
        # Observations
        'accepted_via_score_count': accepted_via_score_count,
        'degraded_due_to_no_score_count': degraded_due_to_no_score_count
    }

    # Write results
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    LOG.info(f"Acceptance test results written to {output_file}")

    # B6: Save full mode results to artifacts for traceability (if not in fast mode)
    if not args.fast:
        utc_timestamp = datetime.utcnow().strftime('%Y%m%d_%H%M%S')
        artifacts_dir = project_root / "artifacts/acceptance"
        artifacts_dir.mkdir(parents=True, exist_ok=True)

        full_results_file = artifacts_dir / f"{utc_timestamp}-full.json"
        with open(full_results_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        LOG.info(f"Full mode results archived to {full_results_file}")

    # Print summary
    print("\nPR1 Sugar Masking Acceptance Test Summary")
    print("=" * 50)

    for sample_id, data in results.items():
        sample_info = data['sample_info']
        analysis = data['analysis']
        heuristic_mode = data['heuristic_mode']
        p0_h2_assertion = data['p0_h2_assertion']

        print(f"\n{sample_id} ({sample_info['type']}): {sample_info['name']}")
        print(f"  Products: {analysis['off_products']} (off) -> {analysis['heuristic_products']} (heuristic)")
        print(f"  Attempts: {analysis['off_attempts']} (off) -> {analysis['heuristic_attempts']} (heuristic)")
        print(f"  Reduction: {analysis['reduction_percentage']:.1f}% ({'PASS' if analysis['overall_pass'] else 'FAIL'})")

        if analysis['pass_criteria_met']:
            print(f"  PASS criteria: {', '.join(analysis['pass_criteria_met'])}")
        else:
            print(f"  FAIL: No criteria met for {analysis['assessment_type']} sample")

        # P0-H2: Display regression test assertion results
        p0_h2_status = 'PASS' if p0_h2_assertion['assertion_pass'] else 'FAIL'
        print(f"  P0-H2 Assertion: {p0_h2_status}")
        if p0_h2_assertion['assertion_criteria_met']:
            print(f"    Met: {', '.join(p0_h2_assertion['assertion_criteria_met'])}")
        if p0_h2_assertion['assertion_criteria_failed']:
            print(f"    Failed: {', '.join(p0_h2_assertion['assertion_criteria_failed'])}")

        # Display comprehensive debugging info
        mask_detail = heuristic_mode.get('sugar_mask_detail', {})
        attempts_delta = analysis.get('attempts_reduction', 0)
        sugar_events_total = analysis['sugar_events'].get('total_sugar_events', 0)

        print(f"  Debug info: mask_size={mask_detail.get('mask_size', 0)}, "
              f"degraded={mask_detail.get('degraded', 'N/A')}, "
              f"rings={mask_detail.get('sugar_rings_count', 0)}, "
              f"is_c_glycoside={mask_detail.get('is_c_glycoside_like', False)}")
        print(f"  Events: sugar_total={sugar_events_total}, "
              f"attempts_delta={attempts_delta}, "
              f"dedup={heuristic_mode.get('dedup_total_hits', 0)}")
        print(f"  Details: {analysis['sugar_events']['filtered_matches']} filtered, "
              f"{analysis['sugar_events']['post_guard_blocked']} blocked")

    assessment = summary['overall_assessment']
    print(f"\nOverall Assessment:")
    print(f"  Glycosides: {assessment['glycoside_passes']}/{assessment['glycoside_samples']} "
          f"PASS multi-criteria assessment (avg reduction: {assessment['glycoside_avg_reduction']:.1f}%)")
    print(f"  Controls: {assessment['control_passes']}/{assessment['control_samples']} "
          f"PASS multi-criteria assessment (avg reduction: {assessment['control_avg_reduction']:.1f}%)")
    print(f"  Legacy metric - Significant reductions (>=30%): Glycosides {assessment['glycoside_significant_reductions']}, "
          f"Controls {assessment['control_significant_reductions']}")

    # P0-H2: Display regression test assertion summary
    print(f"\nP0-H2 Regression Test Assertions:")
    print(f"  Glycosides: {assessment['p0_h2_glycoside_passes']}/{assessment['glycoside_samples']} PASS")
    print(f"  Controls: {assessment['p0_h2_control_passes']}/{assessment['control_samples']} PASS")
    print(f"  Total: {assessment['p0_h2_total_passes']}/{assessment['p0_h2_total_samples']} PASS")

    p0_h2_success_rate = (assessment['p0_h2_total_passes'] / assessment['p0_h2_total_samples']) * 100 if assessment['p0_h2_total_samples'] > 0 else 0
    print(f"  Success Rate: {p0_h2_success_rate:.1f}%")

    # P0-H1: Display parameterized configuration values for traceability
    sugar_cfg = summary['test_metadata'].get('sugar_masking_config', {})
    if isinstance(sugar_cfg, dict):
        print(f"\nSugar Masking Configuration (P0-H1 Parameterized):")
        print(f"  SP3 Thresholds: 5-ring >= {sugar_cfg.get('sp3_threshold_5_ring', 'N/A'):.1f}, "
              f"6-ring >= {sugar_cfg.get('sp3_threshold_6_ring', 'N/A'):.1f}")
        print(f"  Few Masked Atoms Threshold: {sugar_cfg.get('few_masked_atoms_threshold', 'N/A')}")
        print(f"  Allow One-Hop C-Glycoside: {sugar_cfg.get('allow_one_hop_cglyco', 'N/A')}")
        print(f"  Include Intermediate in Minimal Mask: {sugar_cfg.get('cglyco_minimal_mask_include_intermediate', 'N/A')}")
        print(f"  Bridge Oxygen Masking: {sugar_cfg.get('mask_glycosidic_bridge_oxygen', 'N/A')}")
        print(f"  Exocyclic Oxygen Masking: {sugar_cfg.get('mask_exocyclic_oxygen', 'N/A')}")

    # H2-A: Display reproducibility information
    test_meta = summary['test_metadata']
    print(f"\nReproducibility Information (H2-A):")
    print(f"  Random Seed: {test_meta.get('random_seed', 'N/A')}")
    print(f"  Samples File Hash: {test_meta.get('samples_file_hash', 'N/A')}")
    print(f"  Rules: {test_meta.get('rules', 'N/A')}")
    print(f"  Halogens: {test_meta.get('halogens', 'N/A')}")
    print(f"  Symmetry Enabled: {test_meta.get('symmetry_enabled', 'N/A')}")

    # H2-B: Determine if P0-H2 assertions pass and implement hard gate
    assertion_failed = assessment['p0_h2_total_passes'] < assessment['p0_h2_total_samples']
    assertion_status = "FAIL" if assertion_failed else "PASS"

    print(f"\n{'='*60}")
    print(f"P0-H2 ASSERTION HARD GATE: {assertion_status}")
    print(f"{'='*60}")

    if assertion_failed:
        print(f"[FAIL] CRITICAL: {assessment['p0_h2_total_samples'] - assessment['p0_h2_total_passes']} assertion(s) failed!")
        print(f"   This indicates regression in sugar masking behavior.")
        print(f"   All assertions must pass before proceeding to P1/PR2 phases.")
        print(f"")
        print(f"Failed samples:")
        for sample_id, data in results.items():
            if not data['p0_h2_assertion']['assertion_pass']:
                sample_info = data['sample_info']
                p0_h2_assertion = data['p0_h2_assertion']
                print(f"  * {sample_id} ({sample_info['type']}): {sample_info['name']}")
                print(f"    Failed: {', '.join(p0_h2_assertion['assertion_criteria_failed'])}")
        print(f"")
        sys.exit(1)  # Exit with non-zero code to fail CI/automation
    else:
        print(f"[PASS] SUCCESS: All {assessment['p0_h2_total_samples']} assertion(s) passed!")
        print(f"   Sugar masking behavior is stable and meets P0-H2 criteria.")
        print(f"   Safe to proceed with P1/PR2 development phases.")

    print(f"{'='*60}")


if __name__ == '__main__':
    main()
