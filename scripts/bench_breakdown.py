# -*- coding: ascii -*-
"""Bench script for rule family breakdown - performance monitoring."""

import sys
import argparse
import pathlib
import time
sys.path.insert(0, 'src')

from halogenator.enumerate_k import enumerate_products, EnumConfig


# Default test molecules for quick benchmarking
DEFAULT_PARENTS = [
    "c1ccccc1",                                           # benzene
    "OC1CCOCC1",                                          # hydroxylated THP
    "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O",        # glucose
    "O=C1CC(c2ccc(O)cc2)Oc2cc(O)cc(O)c12"                # naringenin
]

# Predefined configuration presets
PRESETS = {
    'fast': {
        'k_max': 1,
        'halogens': ('F',),
        'pruning_cfg': {'enable_symmetry_fold': True, 'enable_state_sig': False}
    },
    'standard': {
        'k_max': 2,
        'halogens': ('F', 'Cl'),
        'pruning_cfg': {'enable_symmetry_fold': True, 'enable_state_sig': True}
    },
    'comprehensive': {
        'k_max': 2,
        'halogens': ('F', 'Cl', 'Br', 'I'),
        'pruning_cfg': {'enable_symmetry_fold': True, 'enable_state_sig': True}
    }
}


def load_molecules(file_path):
    """Load molecules from SMILES file."""
    molecules = []
    path = pathlib.Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Molecules file not found: {file_path}")

    with open(path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line and not line.startswith('#'):
                # Handle both bare SMILES and tab-separated format
                parts = line.split('\t')
                smiles = parts[0].strip()
                if smiles:
                    molecules.append(smiles)

    if not molecules:
        raise ValueError(f"No valid SMILES found in {file_path}")

    return molecules


def parse_halogens(halogen_str):
    """Parse comma-separated halogen list."""
    valid_halogens = {'F', 'Cl', 'Br', 'I'}
    halogens = [h.strip() for h in halogen_str.split(',')]

    for h in halogens:
        if h not in valid_halogens:
            raise ValueError(f"Invalid halogen: {h}. Valid options: {', '.join(sorted(valid_halogens))}")

    return tuple(halogens)


def parse_rules(rules_str):
    """Parse comma-separated rules list."""
    valid_rules = {'R1', 'R2', 'R3', 'R4', 'R5'}
    rules = [r.strip() for r in rules_str.split(',')]

    for r in rules:
        if r not in valid_rules:
            raise ValueError(f"Invalid rule: {r}. Valid options: {', '.join(sorted(valid_rules))}")

    return tuple(rules)


def create_config(preset=None, k_max=None, halogens=None, symmetry_fold=None, state_sig=None, rules=None, alpha_as_beta=None):
    """Create EnumConfig from parameters."""
    if preset:
        if preset not in PRESETS:
            raise ValueError(f"Unknown preset: {preset}. Available: {', '.join(PRESETS.keys())}")
        base_cfg = PRESETS[preset].copy()
    else:
        base_cfg = PRESETS['fast'].copy()

    # Override with explicit parameters
    if k_max is not None:
        base_cfg['k_max'] = k_max
    if halogens is not None:
        base_cfg['halogens'] = halogens
    if rules is not None:
        base_cfg['rules'] = rules

    pruning_cfg = base_cfg.get('pruning_cfg', {}).copy()
    if symmetry_fold is not None:
        pruning_cfg['enable_symmetry_fold'] = symmetry_fold
    if state_sig is not None:
        pruning_cfg['enable_state_sig'] = state_sig
    base_cfg['pruning_cfg'] = pruning_cfg

    # Add rules configuration for alpha-as-beta behavior
    if alpha_as_beta is not None:
        rules_cfg = base_cfg.get('rules_cfg', {}).copy()
        r2_cfg = rules_cfg.get('R2', {}).copy()
        r2_cfg['allow_alpha_as_beta'] = bool(alpha_as_beta)
        rules_cfg['R2'] = r2_cfg
        base_cfg['rules_cfg'] = rules_cfg

    return EnumConfig(**base_cfg)


def format_duration(seconds):
    """Format duration in human-readable format."""
    if seconds < 1:
        return f"{seconds*1000:.1f}ms"
    elif seconds < 60:
        return f"{seconds:.1f}s"
    else:
        minutes = int(seconds // 60)
        secs = seconds % 60
        return f"{minutes}m{secs:.1f}s"


def main():
    """Run parametric benchmark to track rule family behavior and product counts."""
    parser = argparse.ArgumentParser(
        description="Rule family breakdown benchmark with parametric options",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick benchmark with default molecules
  python scripts/bench_breakdown.py

  # Use comprehensive preset
  python scripts/bench_breakdown.py --preset comprehensive

  # Custom configuration
  python scripts/bench_breakdown.py --k-max 2 --halogens F,Cl,Br --symmetry-fold

  # Load molecules from file
  python scripts/bench_breakdown.py --molecules data/input/parents.smi --preset standard

  # Performance focus with timing
  python scripts/bench_breakdown.py --timing --max-molecules 10

  # Test only R2 rules with alpha-as-beta compatibility
  python scripts/bench_breakdown.py --rules R2 --alpha-as-beta --timing

  # Compare strict vs compatibility mode for flavanone
  python scripts/bench_breakdown.py --rules R2 --halogens F --quiet  # strict mode
  python scripts/bench_breakdown.py --rules R2 --halogens F --alpha-as-beta --quiet  # compat mode
        """)

    parser.add_argument('--molecules', '-m', type=str,
                       help='Path to SMILES file (default: use built-in test molecules)')
    parser.add_argument('--preset', '-p', choices=list(PRESETS.keys()),
                       help=f'Configuration preset: {", ".join(PRESETS.keys())}')
    parser.add_argument('--k-max', type=int, metavar='N',
                       help='Maximum substitution depth (overrides preset)')
    parser.add_argument('--halogens', type=str, metavar='LIST',
                       help='Comma-separated halogen list: F,Cl,Br,I (overrides preset)')
    parser.add_argument('--symmetry-fold', action='store_true',
                       help='Enable symmetry folding (overrides preset)')
    parser.add_argument('--no-symmetry-fold', dest='symmetry_fold', action='store_false',
                       help='Disable symmetry folding (overrides preset)')
    parser.add_argument('--state-sig', action='store_true',
                       help='Enable state signature deduplication (overrides preset)')
    parser.add_argument('--no-state-sig', dest='state_sig', action='store_false',
                       help='Disable state signature deduplication (overrides preset)')
    parser.add_argument('--max-molecules', type=int, metavar='N',
                       help='Limit number of molecules to process')
    parser.add_argument('--timing', action='store_true',
                       help='Include detailed timing information')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Reduce output verbosity')
    parser.add_argument('--rules', type=str, metavar='LIST',
                       help='Comma-separated rules list: R1,R2,R3,R4,R5 (overrides preset)')
    parser.add_argument('--alpha-as-beta', action='store_true',
                       help='Enable alpha-as-beta compatibility for R2b (allows alpha=1 positions as beta)')

    parser.set_defaults(symmetry_fold=None, state_sig=None)
    args = parser.parse_args()

    try:
        # Load molecules
        if args.molecules:
            parents = load_molecules(args.molecules)
            if not args.quiet:
                print(f"Loaded {len(parents)} molecules from {args.molecules}")
        else:
            parents = DEFAULT_PARENTS
            if not args.quiet:
                print(f"Using {len(parents)} default test molecules")

        # Apply molecule limit
        if args.max_molecules and args.max_molecules < len(parents):
            parents = parents[:args.max_molecules]
            if not args.quiet:
                print(f"Limited to first {len(parents)} molecules")

        # Parse halogens
        halogens = None
        if args.halogens:
            halogens = parse_halogens(args.halogens)

        # Parse rules
        rules = None
        if args.rules:
            rules = parse_rules(args.rules)

        # Create configuration
        cfg = create_config(
            preset=args.preset,
            k_max=args.k_max,
            halogens=halogens,
            symmetry_fold=args.symmetry_fold,
            state_sig=args.state_sig,
            rules=rules,
            alpha_as_beta=args.alpha_as_beta
        )

        if not args.quiet:
            print(f"Configuration: k_max={cfg.k_max}, halogens={cfg.halogens}")
            if hasattr(cfg, 'rules'):
                print(f"Rules: {cfg.rules}")
            print(f"Pruning: symmetry_fold={cfg.pruning_cfg.get('enable_symmetry_fold')}, "
                  f"state_sig={cfg.pruning_cfg.get('enable_state_sig')}")
            if hasattr(cfg, 'rules_cfg') and cfg.rules_cfg.get('R2', {}).get('allow_alpha_as_beta'):
                print(f"R2b alpha-as-beta compatibility: enabled")

    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Initialize counters
    total = 0
    by_family = {}
    by_rule = {}
    by_halogen = {}
    timing_data = []

    if not args.quiet:
        print("Running benchmark...")

    start_time = time.time()

    # Process molecules
    for mol_idx, smi in enumerate(parents, 1):
        if not args.quiet:
            print(f"Processing {mol_idx}/{len(parents)}: {smi[:50]}...")

        mol_start_time = time.time()
        mol_total = 0
        mol_errors = 0

        try:
            for rec in enumerate_products(smi, cfg):
                total += 1
                mol_total += 1

                # Track by rule family (R2a/R2b -> R2)
                family = rec.get('rule_family') or rec.get('rule')
                if family:
                    by_family[family] = by_family.get(family, 0) + 1

                # Track by specific rule (R2a, R2b, etc.)
                rule = rec.get('rule')
                if rule:
                    by_rule[rule] = by_rule.get(rule, 0) + 1

                # Track by halogen
                halogen = rec.get('halogen')
                if halogen:
                    by_halogen[halogen] = by_halogen.get(halogen, 0) + 1

        except Exception as e:
            mol_errors += 1
            if not args.quiet:
                print(f"  Error processing {smi}: {e}")

        mol_duration = time.time() - mol_start_time
        timing_data.append({
            'molecule': smi[:30],
            'products': mol_total,
            'errors': mol_errors,
            'duration': mol_duration
        })

        if not args.quiet:
            duration_str = format_duration(mol_duration)
            print(f"  Products: {mol_total}, Duration: {duration_str}")

    total_duration = time.time() - start_time

    # Output results
    print(f"\n=== BENCHMARK RESULTS ===")
    print(f"Total molecules: {len(parents)}")
    print(f"Total products: {total}")
    print(f"Total duration: {format_duration(total_duration)}")

    if total > 0:
        avg_products_per_mol = total / len(parents)
        avg_time_per_product = total_duration / total
        print(f"Average products/molecule: {avg_products_per_mol:.1f}")
        print(f"Average time/product: {format_duration(avg_time_per_product)}")

    # Rule family breakdown
    if by_family:
        print(f"\nBy rule family:")
        for family in sorted(by_family.keys()):
            count = by_family[family]
            percentage = (count / total * 100) if total > 0 else 0
            print(f"  {family}: {count} ({percentage:.1f}%)")

    # Specific rule breakdown
    if by_rule:
        print(f"\nBy specific rule:")
        for rule in sorted(by_rule.keys()):
            count = by_rule[rule]
            percentage = (count / total * 100) if total > 0 else 0
            print(f"  {rule}: {count} ({percentage:.1f}%)")

    # Halogen breakdown
    if by_halogen:
        print(f"\nBy halogen:")
        for halogen in sorted(by_halogen.keys()):
            count = by_halogen[halogen]
            percentage = (count / total * 100) if total > 0 else 0
            print(f"  {halogen}: {count} ({percentage:.1f}%)")

    # Timing details
    if args.timing and timing_data:
        print(f"\nTiming details:")
        durations = [data['duration'] for data in timing_data]
        durations.sort()

        fastest = min(durations)
        slowest = max(durations)
        median = durations[len(durations) // 2]

        print(f"  Fastest: {format_duration(fastest)}")
        print(f"  Median:  {format_duration(median)}")
        print(f"  Slowest: {format_duration(slowest)}")

        # Show slowest molecules
        timing_data.sort(key=lambda x: x['duration'], reverse=True)
        print(f"  Top 3 slowest molecules:")
        for i, data in enumerate(timing_data[:3], 1):
            print(f"    {i}. {data['molecule']} - {data['products']} products, {format_duration(data['duration'])}")

    # Sanity checks and warnings
    total_errors = sum(data['errors'] for data in timing_data)
    if total_errors > 0:
        print(f"\nWARNING: {total_errors} processing errors encountered")

    if total == 0:
        print("WARNING: No products generated - possible configuration issue")
        return 1

    # R2 family consistency check
    if 'R2' in by_family and ('R2a' in by_rule or 'R2b' in by_rule):
        r2_family_total = by_family['R2']
        r2_rule_total = by_rule.get('R2a', 0) + by_rule.get('R2b', 0)
        if r2_family_total != r2_rule_total:
            print(f"WARNING: R2 family/rule mismatch: family={r2_family_total}, rules={r2_rule_total}")

    if not args.quiet:
        print("\nBenchmark completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())