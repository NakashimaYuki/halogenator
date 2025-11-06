# -*- coding: ascii -*-
"""
Combinatorial diagnosis CLI module for halogenation enumeration analysis.

This module provides tools to diagnose combinatorial enumeration issues by:
1. Counting theoretical vs actual vs feasible products for molecules
2. Analyzing site availability and rule applicability
3. Identifying discrepancies in expected vs actual product counts
4. Providing detailed breakdowns by rule and halogen type

Usage:
    python -m halogenator.combo_diagnose [OPTIONS] SMILES

Example:
    python -m halogenator.combo_diagnose "C1=CC(=O)C2=C(C(=C(C=C2C=C1)O)O)O" --config configs/pick_k2_compat.yaml
"""

import argparse
import logging
import sys
from typing import Dict, Any, List, Tuple, Set
from dataclasses import dataclass

from .chem_compat import Chem
from .standardize import std_from_smiles
from .sites import aromatic_CH_indices, c_ring_indices, c_ring_sp2_CH_sites, c_ring_sp3_CH2_flavanone_sites
from .sites_methyl import enumerate_methyl_sites
from .sugar_mask import get_sugar_mask_with_status
from .enumerate_k import enumerate_products, EnumConfig
from .cli import load_config
from .constraints import accept as accept_constraints


LOG = logging.getLogger(__name__)


@dataclass
class SiteAnalysis:
    """Analysis results for sites available for a specific rule."""
    rule: str
    total_sites: int
    unmasked_sites: int
    sites_by_type: Dict[str, int]
    theoretical_combinations: int
    feasible_combinations: int


@dataclass
class ComboAnalysis:
    """Complete combinatorial analysis for a molecule."""
    smiles: str
    inchikey: str
    site_analyses: List[SiteAnalysis]
    actual_products: int
    theoretical_total: int
    feasible_total: int
    discrepancy_ratio: float


def count_theoretical_combinations(sites: int, halogens: int, k_max: int) -> int:
    """
    Count theoretical combinations for k-dimensional enumeration.

    For k_max=2: combinations = sites*halogens + sites*(sites-1)*halogens^2
    For k_max=1: combinations = sites*halogens
    """
    if sites == 0 or halogens == 0:
        return 0

    total = 0

    # k=1 combinations
    total += sites * halogens

    # k=2 combinations (if k_max >= 2)
    if k_max >= 2:
        # Each k=1 product can be further modified at any remaining site
        # This is an approximation - actual enumeration has more complex constraints
        total += sites * (sites - 1) * (halogens ** 2)

    return total


def count_feasible_combinations(mol, sites: List[int], halogens: List[str], k_max: int, cfg: Dict[str, Any]) -> int:
    """
    Count feasible combinations considering constraints like graph distance and ring quotas.

    This provides a more realistic estimate than pure theoretical counting.
    """
    if not sites or not halogens:
        return 0

    # For k=1, all sites are potentially feasible
    k1_feasible = len(sites) * len(halogens)

    if k_max < 2:
        return k1_feasible

    # For k=2, apply constraint filtering
    constraints_cfg = cfg.get('constraints', {})
    min_distance = constraints_cfg.get('min_graph_distance', 0)
    per_ring_quota = constraints_cfg.get('per_ring_quota', None)

    k2_feasible = 0

    if min_distance > 0:
        # Count site pairs that satisfy minimum distance constraint
        try:
            distance_matrix = Chem.GetDistanceMatrix(mol)
            valid_pairs = 0

            for i, site1 in enumerate(sites):
                for j, site2 in enumerate(sites):
                    if i != j and distance_matrix[site1][site2] >= min_distance:
                        valid_pairs += 1

            # Each valid pair can use any halogen combination
            k2_feasible = valid_pairs * (len(halogens) ** 2)

        except Exception:
            # Fallback: assume 50% of theoretical combinations are feasible
            k2_feasible = len(sites) * (len(sites) - 1) * (len(halogens) ** 2) // 2
    else:
        # No distance constraint, all pairs are valid
        k2_feasible = len(sites) * (len(sites) - 1) * (len(halogens) ** 2)

    # Apply ring quota constraints (rough approximation)
    if per_ring_quota and per_ring_quota < len(sites):
        # If ring quota is restrictive, reduce feasible combinations
        reduction_factor = min(1.0, per_ring_quota / len(sites))
        k2_feasible = int(k2_feasible * reduction_factor)

    return k1_feasible + k2_feasible


def analyze_rule_sites(mol, rule: str, cfg: Dict[str, Any], sugar_mask: Set[int]) -> SiteAnalysis:
    """Analyze sites available for a specific rule."""

    if rule in ('R1', 'R2a'):
        # Aromatic CH sites
        all_sites = aromatic_CH_indices(mol)
        rule_halogens = cfg.get('halogens', ['F', 'Cl', 'Br', 'I'])

    elif rule == 'R2b':
        # Flavanone sp3 CH2 sites
        all_sites = c_ring_sp3_CH2_flavanone_sites(mol)
        rule_halogens = cfg.get('rules_cfg', {}).get('R2', {}).get('allowed_halogens', cfg.get('halogens', ['F', 'Cl', 'Br', 'I']))

    elif rule in ('R6', 'R6_methyl'):
        # Methyl sites
        r6_config = cfg.get('rules_cfg', {}).get('R6_methyl', {})
        allow_on_methoxy = r6_config.get('allow_on_methoxy', False)
        allow_allylic_methyl = r6_config.get('allow_allylic_methyl', False)

        site_descriptors = enumerate_methyl_sites(mol, sugar_mask, allow_on_methoxy, allow_allylic_methyl)
        all_sites = [site['idx'] for site in site_descriptors]

        # Count sites by type
        sites_by_type = {}
        for site in site_descriptors:
            site_type = site.get('kind', 'OTHER_CH3')
            sites_by_type[site_type] = sites_by_type.get(site_type, 0) + 1

        rule_halogens = r6_config.get('allowed', ['F', 'Cl'])

    else:
        # R3, R4, R5 - reaction-based rules (harder to analyze)
        all_sites = []
        sites_by_type = {}
        rule_halogens = cfg.get('halogens', ['F', 'Cl', 'Br', 'I'])

    # Filter out masked sites
    unmasked_sites = [site for site in all_sites if site not in sugar_mask]

    # Count combinations
    theoretical = count_theoretical_combinations(len(all_sites), len(rule_halogens), cfg.get('k_max', 1))
    feasible = count_feasible_combinations(mol, unmasked_sites, rule_halogens, cfg.get('k_max', 1), cfg)

    # Default sites_by_type for non-R6 rules
    if rule not in ('R6', 'R6_methyl'):
        sites_by_type = {'total': len(all_sites)}

    return SiteAnalysis(
        rule=rule,
        total_sites=len(all_sites),
        unmasked_sites=len(unmasked_sites),
        sites_by_type=sites_by_type,
        theoretical_combinations=theoretical,
        feasible_combinations=feasible
    )


def run_actual_enumeration(smiles: str, cfg: Dict[str, Any]) -> int:
    """Run actual enumeration and count products."""
    try:
        # Convert dictionary config to EnumConfig object
        enum_cfg = EnumConfig(
            k_max=cfg.get('k_max', 1),
            halogens=tuple(cfg.get('halogens', ['F', 'Cl'])),
            rules=tuple(cfg.get('rules', [])),
            rules_cfg=cfg.get('rules_cfg', {}),
            constraints=cfg.get('constraints', {}),
            engine_cfg=cfg.get('engine', {}),
            std_cfg=cfg.get('standardize', {}),
            qc_cfg=cfg.get('qc', {}),
            pruning_cfg=cfg.get('pruning', {}),
            sugar_cfg=cfg.get('sugar_cfg', {}),
            symmetry_cfg=cfg.get('symmetry', {}),
            random_seed=cfg.get('random_seed', None)
        )

        products = list(enumerate_products(smiles, enum_cfg))
        return len(products)
    except Exception as e:
        LOG.error(f"Enumeration failed: {e}")
        return 0


def analyze_molecule(smiles: str, cfg: Dict[str, Any]) -> ComboAnalysis:
    """Perform complete combinatorial analysis for a molecule."""

    # Standardize molecule
    mol = std_from_smiles(smiles)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")

    # Compute identifiers
    inchikey = Chem.MolToInchiKey(mol) if mol else "UNKNOWN"

    # Get sugar mask
    sugar_cfg = cfg.get('sugar_cfg', {})
    sugar_mask, _ = get_sugar_mask_with_status(mol, mode=sugar_cfg.get('mode', 'heuristic'), sugar_cfg=sugar_cfg)

    # Analyze each rule
    site_analyses = []
    theoretical_total = 0
    feasible_total = 0

    for rule in cfg.get('rules', []):
        analysis = analyze_rule_sites(mol, rule, cfg, sugar_mask)
        site_analyses.append(analysis)
        theoretical_total += analysis.theoretical_combinations
        feasible_total += analysis.feasible_combinations

    # Run actual enumeration
    actual_products = run_actual_enumeration(smiles, cfg)

    # Calculate discrepancy ratio
    discrepancy_ratio = actual_products / feasible_total if feasible_total > 0 else 0.0

    return ComboAnalysis(
        smiles=smiles,
        inchikey=inchikey,
        site_analyses=site_analyses,
        actual_products=actual_products,
        theoretical_total=theoretical_total,
        feasible_total=feasible_total,
        discrepancy_ratio=discrepancy_ratio
    )


def print_analysis_report(analysis: ComboAnalysis) -> None:
    """Print detailed analysis report."""

    print(f"\n=== Combinatorial Analysis Report ===")
    print(f"SMILES: {analysis.smiles}")
    print(f"InChIKey: {analysis.inchikey}")
    print(f"\n--- Site Analysis by Rule ---")

    for site_analysis in analysis.site_analyses:
        print(f"\nRule: {site_analysis.rule}")
        print(f"  Total sites: {site_analysis.total_sites}")
        print(f"  Unmasked sites: {site_analysis.unmasked_sites}")

        if site_analysis.rule in ('R6', 'R6_methyl'):
            print(f"  Sites by type:")
            for site_type, count in site_analysis.sites_by_type.items():
                print(f"    {site_type}: {count}")

        print(f"  Theoretical combinations: {site_analysis.theoretical_combinations}")
        print(f"  Feasible combinations: {site_analysis.feasible_combinations}")

    print(f"\n--- Summary ---")
    print(f"Total theoretical combinations: {analysis.theoretical_total}")
    print(f"Total feasible combinations: {analysis.feasible_total}")
    print(f"Actual products enumerated: {analysis.actual_products}")
    print(f"Discrepancy ratio (actual/feasible): {analysis.discrepancy_ratio:.3f}")

    if analysis.discrepancy_ratio < 0.5:
        print(f"WARNING: LOW EFFICIENCY: Actual products much lower than feasible")
    elif analysis.discrepancy_ratio > 1.1:
        print(f"WARNING: OVER-ENUMERATION: More products than expected")
    else:
        print(f"OK: REASONABLE EFFICIENCY: Actual products match feasible estimates")


def main():
    """CLI entry point for combinatorial diagnosis."""
    parser = argparse.ArgumentParser(description="Diagnose combinatorial enumeration issues")
    parser.add_argument("smiles", help="Input SMILES string")
    parser.add_argument("--config", "-c", required=True, help="Configuration file path")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    parser.add_argument("--quiet", "-q", action="store_true", help="Suppress all output except results")

    args = parser.parse_args()

    # Configure logging
    if args.quiet:
        logging.basicConfig(level=logging.ERROR)
    elif args.verbose:
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    else:
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    try:
        # Load configuration
        cfg = load_config(args.config)

        # Perform analysis
        analysis = analyze_molecule(args.smiles, cfg)

        # Print report
        print_analysis_report(analysis)

        return 0

    except Exception as e:
        LOG.error(f"Analysis failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())