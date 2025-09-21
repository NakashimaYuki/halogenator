# -*- coding: ascii -*-
"""
Universal k-dimensional BFS halogenation enumerator.

QA Statistics Counting Semantics:
================================
All QA counters follow "per-attempt" semantics unless otherwise noted:

- isotope_unavailable: Rule/halogen/layer attempt where isotope strategy cannot be used
  (no matches, multi-reactant template, or RDKit reaction failure) - count 1 per attempt
- isotope_miss: Isotope tagging succeeded and reaction ran, but tagged site not found 
  in products (rare edge case) - count 1 per failed match  
- atommap_used: Site detection succeeded using AtomMapNum fallback - count 1 per successful attempt
- heuristic_used: Site detection succeeded using heuristic fallback - count 1 per successful attempt
- no_product_matches: Match/attempt that ultimately failed to produce valid products - count 1 per failed match
- template_unsupported: Template has unsupported structure (e.g. multi-reactant) - count 1 per attempt
- dedup_hits_statesig/dedup_hits_inchi: Product deduplicated by state signature/InChIKey - count 1 per duplicate

Note: Counters are aggregated across all rules, halogens, and BFS layers for a given parent molecule.

TODO(v0.X): Remove deprecated _site_proximity_blocked function (replaced by _expand_mask_by_radius).
"""

import logging
import functools
import json
import os
import sys
import types
from typing import Dict, Any, List, Tuple, Iterable, Optional, Set
from collections import deque, defaultdict

from .standardize import std_from_smiles, to_inchikey, to_inchikey_sanitized
from .dedup_util import early_check, commit, compute_stable_key
from .qa_utils import qa_mark, PIVOT_EVENT_KEYS, canonical_event
from .rules import build_reactions
from .sites import aromatic_CH_indices, c_ring_indices, is_c_ring_site_ready, ensure_ready, flavonoid_ring_label, is_carbonyl_carbon, c_ring_sp2_CH_sites, c_ring_sp3_CH2_flavanone_sites, filter_sites_with_mask
from .chem_compat import Chem, canonical_rank_atoms_safe
from .reactions import apply_single_site_halogenation
from .qc import sanitize_ok, basic_descriptors, pains_flags
from .constraints import accept as accept_constraints
from .state_sig import compute_state_signature
from .guard import rdkit_guard
from .schema import empty_qa_paths, QA_PATH_KEYS
from .sugar_mask import get_sugar_mask, get_sugar_mask_with_status, get_sugar_mask_with_full_status, post_guard_blocked, compute_sugar_audit_fields
from .symmetry import canonical_ranks_on_masked_subgraph
from .budget import BudgetState
from .sites_methyl import enumerate_methyl_sites
from .rules_methyl import apply_methyl_step, apply_methyl_macro, validate_methyl_halogen

# Isotope tag used for internal site tracking (no chemical significance)
ISOTOPE_TAG = 109

# Rules version for cache invalidation - increment when reaction rules change
RULES_VERSION = "1.0.0"

# Logger for enumeration process
LOG = logging.getLogger(__name__)

_PYTHON_NUMPY_SEED_INFO_EMITTED = False
_PYTHON_NUMPY_SEED_LAST: Optional[int] = None
_PYTHON_NUMPY_SEED_DEBUG_EMITTED = False

# Global warning deduplication for reaction failures per rule/halogen combination
_reaction_warning_counts = defaultdict(int)
_max_warnings_per_rule_halogen = 1  # Emit at most 1 warning per rule/halogen combination

# Chem alias comes from chem_compat and is patchable via src.halogenator.chem_compat.Chem


def reset_reaction_warning_counts():
    """Reset the global reaction warning counts. For testing and multi-run scenarios."""
    global _reaction_warning_counts
    _reaction_warning_counts.clear()


def print_reaction_warning_summary():
    """Print a single-line summary of deduplicated reaction warnings."""
    if not _reaction_warning_counts:
        return
    
    unique_failures = len(_reaction_warning_counts)
    total_suppressed = sum(count - 1 for count in _reaction_warning_counts.values() if count > 1)
    
    if LOG.isEnabledFor(logging.INFO):
        LOG.info(f"Reaction warnings dedup: {unique_failures} unique rule/halogen failures (suppressed {total_suppressed} duplicates)")
    elif LOG.isEnabledFor(logging.WARNING):
        LOG.warning(f"Reaction warnings dedup: {unique_failures} unique rule/halogen failures (suppressed {total_suppressed} duplicates)")


class QAAggregator:
    """
    Aggregates QA statistics by rule, halogen, and k dimensions.
    Provides real engine-level pivot statistics for M2 reporting.
    """

    def __init__(self, debug_consistency: bool = False):
        self.by_rule = defaultdict(lambda: defaultdict(int))
        self.by_halogen = defaultdict(lambda: defaultdict(int))
        self.by_k = defaultdict(lambda: defaultdict(int))
        self.by_rule_halogen = defaultdict(lambda: defaultdict(int))
        self.by_rule_halogen_k = defaultdict(lambda: defaultdict(int))
        self.by_rule_halogen_ops_atoms = defaultdict(lambda: defaultdict(int))  # Dual-metric tracking
        self.by_rule_ops_atoms = defaultdict(lambda: defaultdict(int))  # Dual-metric tracking
        self.paths = {}  # Dedicated container for qa_paths-like counters
        self.debug_consistency = debug_consistency
        self.qa_paths_pivot_sum = defaultdict(int)  # Track pivot sums for consistency checking
    
    def reset(self):
        """Reset all counters."""
        self.by_rule.clear()
        self.by_halogen.clear()
        self.by_k.clear()
        self.by_rule_halogen.clear()
        self.by_rule_halogen_k.clear()
        self.by_rule_halogen_ops_atoms.clear()
        self.by_rule_ops_atoms.clear()
        self.paths.clear()
        self.qa_paths_pivot_sum.clear()
    
    def record(self, event: str, rule: str = None, halogen: str = None, k: int = None, k_ops: int = None, k_atoms: int = None, amount: int = 1):
        """
        Record a QA event with optional rule/halogen/k/ops/atoms context.

        Event Semantics:
        - attempts: increment ONCE for each rule x halogen x k attempt BEFORE running a reaction or heuristic path
        - products: increment for each attempt that yields at least one product molecule
        - template_unsupported: increment when rxn is None or RunReactants raises
        - isotope_unavailable: increment exactly where qa_paths increments (isotope pattern not available)
        - isotope_miss: increment exactly where qa_paths increments (isotope pattern miss)
        - atommap_used: increment when atom map fallback is engaged
        - heuristic_used: increment when heuristic fallback is engaged
        - no_product_matches: increment when an attempt yields zero valid products

        Args:
            event: Event type (see semantics above)
            rule: Rule identifier (R1, R3, etc.)
            halogen: Halogen symbol (F, Cl, Br, I)
            k: Substitution depth
            k_ops: Number of operations for dual-metric tracking
            k_atoms: Number of atoms for dual-metric tracking
            amount: Number to increment by (default 1)
        """
        # Track qa_paths events in pivot sum for consistency checking
        if event in PIVOT_EVENT_KEYS:
            self.qa_paths_pivot_sum[event] += amount

        if rule:
            self.by_rule[rule][event] += amount
        if halogen:
            self.by_halogen[halogen][event] += amount
        if k is not None:
            self.by_k[k][event] += amount
        if rule and halogen:
            key = f"{rule}_{halogen}"
            self.by_rule_halogen[key][event] += amount
        if rule and halogen and k is not None:
            key = f"{rule}_{halogen}_{k}"
            self.by_rule_halogen_k[key][event] += amount

        # Dual-metric pivot tracking
        if rule and halogen and k_ops is not None and k_atoms is not None:
            key = f"{rule}_{halogen}_{k_ops}_{k_atoms}"
            self.by_rule_halogen_ops_atoms[key][event] += amount
        if rule and k_ops is not None and k_atoms is not None:
            key = f"{rule}_{k_ops}_{k_atoms}"
            self.by_rule_ops_atoms[key][event] += amount
    
    def record_attempt_result(self, rule: str, halogen: str, k: int, produced_count: int, qa_events_dict: dict,
                            *, k_ops: int = None, k_atoms: int = None):
        """
        Record a complete attempt result with proper semantics.

        Args:
            rule: Rule identifier (R1, R3, etc.)
            halogen: Halogen symbol (F, Cl, Br, I)
            k: Substitution depth
            produced_count: Number of products produced by this attempt
            qa_events_dict: Dict of qa event types to counts for this attempt
            k_ops: Operations count for dual-metric pivots (optional)
            k_atoms: Atoms count for dual-metric pivots (optional)
        """
        # Record one attempt
        self.record('attempts', rule=rule, halogen=halogen, k=k, k_ops=k_ops, k_atoms=k_atoms, amount=1)

        # Record products or no_product_matches based on outcome
        if produced_count > 0:
            self.record('products', rule=rule, halogen=halogen, k=k, k_ops=k_ops, k_atoms=k_atoms, amount=1)
        else:
            self.record('no_product_matches', rule=rule, halogen=halogen, k=k, k_ops=k_ops, k_atoms=k_atoms, amount=1)

        # Record QA events for this attempt
        for qa_key, count in (qa_events_dict or {}).items():
            if count > 0:
                self.record(qa_key, rule=rule, halogen=halogen, k=k, k_ops=k_ops, k_atoms=k_atoms, amount=count)
    
    def merge(self, other):
        """Merge another QAAggregator into this one."""
        for rule, events in other.by_rule.items():
            for event, count in events.items():
                self.by_rule[rule][event] += count
        
        for halogen, events in other.by_halogen.items():
            for event, count in events.items():
                self.by_halogen[halogen][event] += count
                
        for k, events in other.by_k.items():
            for event, count in events.items():
                self.by_k[k][event] += count
                
        for key, events in other.by_rule_halogen.items():
            for event, count in events.items():
                self.by_rule_halogen[key][event] += count
                
        for key, events in other.by_rule_halogen_k.items():
            for event, count in events.items():
                self.by_rule_halogen_k[key][event] += count

        # Merge dual-metric pivot dimensions
        for key, events in other.by_rule_halogen_ops_atoms.items():
            for event, count in events.items():
                self.by_rule_halogen_ops_atoms[key][event] += count

        for key, events in other.by_rule_ops_atoms.items():
            for event, count in events.items():
                self.by_rule_ops_atoms[key][event] += count

        # Merge dedicated paths container
        for k, v in other.paths.items():
            self.paths[k] = self.paths.get(k, 0) + int(v)
    
    def record_paths(self, paths: dict):
        """
        Record sample-level QA path events without associating them to rule/halogen/k grids.
        This avoids dimension pollution while ensuring paths appear in final qa_paths output.

        DEFENSIVE: Pivot events are rejected to prevent totals x pivots inconsistency.

        Args:
            paths: Dict of qa_paths event types to counts
        """
        for k, v in (paths or {}).items():
            if v:
                canonical_key = canonical_event(k)
                if canonical_key in PIVOT_EVENT_KEYS:
                    if self.debug_consistency:
                        raise ValueError(f"[qa] STRICT MODE: Pivot event '{canonical_key}' (original: '{k}', value: {v}) "
                                       f"cannot be recorded via record_paths(). Pivot events must come from "
                                       f"rule/halogen/k aggregation only. Use qa_mark() with attempt containers instead.")
                    else:
                        LOG.warning("[qa] record_paths ignored pivot event: %s (original: %s, value: %s)",
                                   canonical_key, k, v)
                    continue
                self.paths[canonical_key] = self.paths.get(canonical_key, 0) + int(v)

    def record_path(self, event: str, amount: int = 1):
        """
        Record a single sample-level QA path event safely.

        This is a safer wrapper around record_paths for single events,
        with built-in pivot validation and canonical naming.

        Args:
            event: QA event name (will be canonicalized)
            amount: Increment amount (default: 1, must be positive)
        """
        if amount <= 0:
            return
        canonical_key = canonical_event(event)
        if canonical_key in PIVOT_EVENT_KEYS:
            if self.debug_consistency:
                raise ValueError(f"[qa] STRICT MODE: Pivot event '{canonical_key}' (original: '{event}') "
                               f"cannot be recorded via record_path(). Pivot events must come from "
                               f"rule/halogen/k aggregation only. Use qa_mark() with attempt containers instead.")
            else:
                LOG.warning("[qa] record_path ignored pivot event: %s (original: %s)",
                           canonical_key, event)
            return
        self.paths[canonical_key] = self.paths.get(canonical_key, 0) + int(amount)

    def pivot_event_sums(self) -> Dict[str, int]:
        """
        Return aggregated pivot event sums across all rule/halogen/k combinations.

        This provides a summary view of pivot events for test validation,
        enabling comparison with qa_paths output from _compute_totals_from_aggregator.

        Returns:
            Dictionary mapping canonical pivot event names to their total counts
        """
        sums = {}
        pivots = self.to_pivots_dict().get('by_rule_halogen_k', {})

        for _grid_key, events in pivots.items():
            for event_key, count in events.items():
                canonical_key = canonical_event(event_key)
                if canonical_key in PIVOT_EVENT_KEYS:
                    sums[canonical_key] = sums.get(canonical_key, 0) + int(count)

        return sums

    def record_guard_failure(self, rule: str, halogen: str, k: int, qa_events_dict: dict):
        """
        Record guard failure as template_unsupported (not no_product_matches).

        Args:
            rule: Rule identifier (usually 'guard')
            halogen: Halogen symbol (usually 'rdkit')
            k: Substitution depth
            qa_events_dict: Dict of qa event types to counts
        """
        # Record one attempt with template_unsupported outcome
        self.record('attempts', rule=rule, halogen=halogen, k=k, amount=1)
        self.record('template_unsupported', rule=rule, halogen=halogen, k=k, amount=1)

        # Record specific qa_paths events (but not template_unsupported again)
        for qa_key, count in qa_events_dict.items():
            if qa_key != 'template_unsupported' and count > 0:
                self.record(qa_key, rule=rule, halogen=halogen, k=k, amount=count)
    
    def to_pivots_dict(self):
        """Convert to dictionary format for JSON serialization with stable schema."""
        def _convert_defaultdict(dd):
            """Recursively convert defaultdict to regular dict."""
            if isinstance(dd, defaultdict):
                return {k: _convert_defaultdict(v) for k, v in dd.items()}
            return dd
        
        # Always include all dimensions for stable schema (allows empty dicts)
        return {
            'by_rule': _convert_defaultdict(self.by_rule),
            'by_halogen': _convert_defaultdict(self.by_halogen),
            'by_k': _convert_defaultdict(self.by_k),
            'by_rule_halogen': _convert_defaultdict(self.by_rule_halogen),
            'by_rule_halogen_k': _convert_defaultdict(self.by_rule_halogen_k),
            'by_rule_halogen_ops_atoms': _convert_defaultdict(self.by_rule_halogen_ops_atoms),
            'by_rule_ops_atoms': _convert_defaultdict(self.by_rule_ops_atoms)
        }


def _check_totals_pivots_consistency(totals_qa_paths: Dict[str, int], pivot_sums: Dict[str, int],
                                   enable_check: bool = False) -> None:
    """
    Check consistency between QA totals and pivot sums.

    Args:
        totals_qa_paths: Dictionary of QA totals from qa_bus
        pivot_sums: Dictionary of pivot sums from aggregator
        enable_check: Whether to perform the check (controlled by config)
    """
    if not enable_check:
        return

    inconsistencies = []
    # Only compare events that belong to pivot-side counting
    keys_to_check = (set(PIVOT_EVENT_KEYS)
                    & (set(totals_qa_paths.keys()) | set(pivot_sums.keys())))

    for key in keys_to_check:
        total_value = totals_qa_paths.get(key, 0)
        pivot_value = pivot_sums.get(key, 0)

        if total_value != pivot_value:
            inconsistencies.append(f"{key}: total={total_value} != pivot_sum={pivot_value}")

    if inconsistencies:
        LOG.warning("QA totals-pivots consistency violations detected: " + "; ".join(inconsistencies))


def _run_reaction_safely(rxn, reactants, rule_id: str, halogen: str, stats_dict: Dict[str, Any] = None, aggregator=None, current_k: int = None):
    """
    Safely execute a reaction with proper error handling and logging.
    
    Args:
        rxn: RDKit reaction object
        reactants: Tuple of reactant molecules
        rule_id: Rule identifier for logging
        halogen: Halogen identifier for logging  
        stats_dict: Dictionary to update with failure counts
        aggregator: DEPRECATED - aggregator recording must happen at attempt boundaries only
        current_k: DEPRECATED - not used without aggregator
    
    Returns:
        List of product sets or empty list on failure
    """
    if rxn is None:
        if stats_dict:
            stats_dict['template_unsupported'] = stats_dict.get('template_unsupported', 0) + 1
        warning_key = f"{rule_id}_{halogen}"
        if _reaction_warning_counts[warning_key] < _max_warnings_per_rule_halogen:
            LOG.warning(f"Rule {rule_id} with {halogen}: reaction is None")
            _reaction_warning_counts[warning_key] += 1
        return []
    
    try:
        reaction_products = rxn.RunReactants(reactants)
        return reaction_products
    except Exception as e:
        if stats_dict:
            stats_dict['template_unsupported'] = stats_dict.get('template_unsupported', 0) + 1
        warning_key = f"{rule_id}_{halogen}"
        if _reaction_warning_counts[warning_key] < _max_warnings_per_rule_halogen:
            LOG.warning(f"Rule {rule_id} with {halogen}: reaction failed - {str(e)}")
            _reaction_warning_counts[warning_key] += 1
        return []


# Cached reaction building for performance
@functools.lru_cache(maxsize=4)
def _build_reactions_cached(rules_version: str):
    """
    Cached version of build_reactions() with version-aware cache key.
    
    Args:
        rules_version: Version string for cache invalidation when rules change
    
    Returns:
        Dictionary of reaction patterns by rule and halogen
    
    Cache size of 4 allows for some version history while preventing unlimited growth.
    """
    LOG.debug(f"Building reaction patterns (cached, version {rules_version})")
    return build_reactions()


def clear_reactions_cache():
    """Clear the reactions cache. Useful when rules change or for testing."""
    _build_reactions_cached.cache_clear()


def get_reactions_cache_info():
    """Get cache statistics for reactions building performance monitoring."""
    return _build_reactions_cached.cache_info()


def _validate_totals_pivots_consistency(qa_stats_dict: Dict[str, Any], aggregator: QAAggregator) -> None:
    """
    Validate that totals in qa_paths match sums across pivot dimensions.
    
    Args:
        qa_stats_dict: Dictionary containing qa_paths totals
        aggregator: QAAggregator containing pivot breakdowns
        
    Raises:
        AssertionError: If totals do not match pivots sums (in debug builds)
        
    Logs warnings for any detected mismatches.
    """
    if not qa_stats_dict.get('qa_paths'):
        return
        
    pivots = aggregator.to_pivots_dict()
    if not pivots.get('by_rule_halogen_k'):
        # No pivot data to validate against
        return
    
    # Calculate sums across the most granular pivot dimension (by_rule_halogen_k)
    pivot_sums = {
        **empty_qa_paths(),
        'attempts': 0,
        'products': 0,
        'no_product_matches': 0,
        'template_unsupported': 0
    }
    
    for key, events in pivots['by_rule_halogen_k'].items():
        for metric, count in events.items():
            if metric in pivot_sums:
                pivot_sums[metric] += count
    
    # Validate qa_paths consistency
    qa_paths_totals = qa_stats_dict['qa_paths']
    mismatches = []
    
    for qa_key in ['isotope_unavailable', 'isotope_miss', 'atommap_used', 'heuristic_used']:
        total_value = qa_paths_totals.get(qa_key, 0)
        pivot_sum = pivot_sums.get(qa_key, 0)
        
        if total_value != pivot_sum:
            mismatches.append(f"{qa_key}: total={total_value} != pivot_sum={pivot_sum}")
    
    # Validate core metrics consistency (full equality constraints)
    for metric in ['attempts', 'products', 'no_product_matches', 'template_unsupported']:
        if metric in qa_stats_dict:
            total_value = qa_stats_dict[metric]
            pivot_sum = pivot_sums.get(metric, 0)
            
            if total_value != pivot_sum:
                mismatches.append(f"{metric}: total={total_value} != pivot_sum={pivot_sum}")
    
    # Validate attempt invariants with three-way exclusive classification
    # Each attempt has exactly one outcome: success, match failure, or template unsupported
    attempts = pivot_sums.get('attempts', 0)
    products = pivot_sums.get('products', 0)
    no_matches = pivot_sums.get('no_product_matches', 0)
    tmpl_unsup = pivot_sums.get('template_unsupported', 0)
    
    if attempts > 0:
        if attempts < products:
            mismatches.append(f"attempts={attempts} < products={products} (invariant violation)")
        
        # Updated invariant: attempts must equal sum of three mutually exclusive outcomes
        if attempts != products + no_matches + tmpl_unsup:
            mismatches.append(
                f"attempts={attempts} != products={products} + "
                f"no_matches={no_matches} + template_unsupported={tmpl_unsup} (invariant violation)"
            )
    
    # Log and assert on mismatches
    if mismatches:
        error_msg = "QA totals-pivots consistency violations detected:\n" + "\n".join(f"  {m}" for m in mismatches)
        LOG.error(error_msg)
        
        # In debug builds, raise assertion error
        debug_env = (os.environ.get('HALO_ASSERT_PIVOT_CONSISTENCY', '') or 
                    os.environ.get('HALOGENATOR_DEBUG', ''))
        if debug_env.lower() in ('1', 'true', 'yes'):
            raise AssertionError(f"Totals-pivots consistency check failed: {mismatches}")
    else:
        LOG.debug("QA totals-pivots consistency validation passed")


def _compute_totals_from_aggregator(aggregator: QAAggregator) -> Dict[str, Any]:
    """
    Compute totals by summing over aggregator pivots.
    This ensures totals are derived from the same source as pivots, eliminating inconsistency.
    
    Args:
        aggregator: QAAggregator containing pivot breakdowns
        
    Returns:
        Dictionary with stable schema totals (always includes core fields, allows zero values)
    """
    pivots = aggregator.to_pivots_dict()
    
    # Initialize stable totals schema (always present, zero if no data)
    totals = {
        'attempts': 0,
        'products': 0, 
        'no_product_matches': 0,
        'template_unsupported': 0,
        'qa_paths': empty_qa_paths()
    }
    
    # Sum across the most granular pivot dimension (by_rule_halogen_k)
    if pivots.get('by_rule_halogen_k'):
        for key, events in pivots['by_rule_halogen_k'].items():
            for metric, count in events.items():
                m = canonical_event(metric)
                if m in PIVOT_EVENT_KEYS:
                    # Add pivot events to qa_paths (they belong here, sourced from pivots)
                    totals['qa_paths'][m] = totals['qa_paths'].get(m, 0) + count
                elif m in totals:
                    # Add other metrics (attempts, products, etc.) to their respective totals
                    totals[m] += count

    # Merge dedicated paths container (only for non-pivot events)
    for k, v in aggregator.paths.items():
        key = canonical_event(k)
        if key in PIVOT_EVENT_KEYS:
            # Defensive: pivot events should not appear in paths (record_paths() should reject them)
            LOG.warning(f"[qa] paths contains pivot key unexpectedly: {key} (v={v})")
            continue
        totals['qa_paths'][key] = totals['qa_paths'].get(key, 0) + int(v)

    return totals


def _snapshot_qa_paths_from_aggregator(agg: QAAggregator) -> Dict[str, int]:
    """
    Build qa_paths snapshot from aggregator (single source of truth).

    Combines:
      - pivot events: sum from by_rule_halogen_k pivots (canonicalized)
      - non-pivot events: from aggregator.paths (already canonicalized & defended)

    This ensures streaming and final snapshots derive from the same source,
    eliminating the T20 streaming vs snapshot inconsistency.

    Args:
        agg: QAAggregator containing both pivot and path events

    Returns:
        Dictionary mapping canonical event names to counts
    """
    qa_paths = dict(agg.pivot_event_sums())  # pivot-only events

    # Add non-pivot events from aggregator.paths
    for k, v in (agg.paths or {}).items():
        ck = canonical_event(k)
        if ck in PIVOT_EVENT_KEYS:
            # Should never happen thanks to record_paths defense, but be defensive
            continue
        qa_paths[ck] = qa_paths.get(ck, 0) + int(v)

    return qa_paths


class QAStats:
    """QA statistics aggregator for enumeration process following per-attempt semantics."""
    
    def __init__(self):
        self.no_product_matches = 0
        self.dedup_hits_statesig = 0  
        self.dedup_hits_inchi = 0
        self.qa_paths = empty_qa_paths()
        self.template_unsupported = 0  # Templates with unsupported structure
        
    def merge(self, other_stats: Dict[str, Any]):
        """Merge stats from a local rule/layer stats dict.

        P0 DEDUP NAMING: Uses unified metric names (dedup_hits_inchi) and provides
        backward compatibility for old field names.
        """
        self.no_product_matches += other_stats.get('no_product_matches', 0)
        self.dedup_hits_statesig += other_stats.get('dedup_hits_statesig', 0)
        self.dedup_hits_inchi += other_stats.get('dedup_hits_inchi', 0)
        self.template_unsupported += other_stats.get('template_unsupported', 0)
        
        # Backward compatibility: merge old dedup field names into new ones
        self.dedup_hits_statesig += other_stats.get('statesig_hits', 0)
        self.dedup_hits_inchi += other_stats.get('inchi_hits', 0)
        
        # Merge QA path stats if present
        if 'qa_paths' in other_stats:
            for key, count in other_stats['qa_paths'].items():
                if key in self.qa_paths:
                    self.qa_paths[key] += count

    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for external consumption."""
        return {
            'no_product_matches': self.no_product_matches,
            'dedup_hits_statesig': self.dedup_hits_statesig,
            'dedup_hits_inchi': self.dedup_hits_inchi,
            'template_unsupported': self.template_unsupported,
            'qa_paths': self.qa_paths.copy()
        }
        
    def log_summary(self):
        """Log QA summary at DEBUG level."""
        LOG.debug(f"QA Summary - No product matches: {self.no_product_matches}, "
                    f"Dedup hits (state sig): {self.dedup_hits_statesig}, "
                    f"Dedup hits (InChI): {self.dedup_hits_inchi}, "
                    f"Template unsupported: {self.template_unsupported}, "
                    f"QA paths: {self.qa_paths}")


class EnumConfig:
    """
    Configuration for k-dimensional enumeration.

    Rule Execution Order:
        Rules are executed in a FIXED order: R3/R4/R5 (reaction-type) first,
        then R1/R2 (site-type). The 'rules' parameter only controls which
        rules are enabled/disabled, NOT their execution order.

    Args:
        k_max: Maximum substitution depth (default: 1)
        halogens: Tuple of halogens to use (default: ('F', 'Cl', 'Br', 'I'))
        rules: Tuple of rules to enable (default: all rules).
               IMPORTANT: Only controls rule enabling, not execution order.
               Execution order is always: R3, R4, R5, then R1, R2.
        constraints: Constraint configuration dict
        std_cfg: Standardization configuration dict
        qc_cfg: Quality control configuration dict
        pruning_cfg: Pruning configuration dict
        sugar_cfg: Sugar masking configuration dict
        symmetry_cfg: Symmetry computation configuration dict
        random_seed: Random seed for reproducible enumeration (default: None for non-deterministic)
    """

    def __init__(self, k_max: int = 1, halogens: Tuple[str, ...] = ('F', 'Cl', 'Br', 'I'),
                 rules: Optional[Tuple[str, ...]] = None,
                 constraints: Optional[Dict[str, Any]] = None,
                 std_cfg: Optional[Dict[str, Any]] = None,
                 qc_cfg: Optional[Dict[str, Any]] = None,
                 pruning_cfg: Optional[Dict[str, Any]] = None,
                 sugar_cfg: Optional[Dict[str, Any]] = None,
                 symmetry_cfg: Optional[Dict[str, Any]] = None,
                 rules_cfg: Optional[Dict[str, Any]] = None,
                 engine_cfg: Optional[Dict[str, Any]] = None,
                 random_seed: Optional[int] = None):
        self.k_max = k_max
        self.halogens = halogens
        self.rules = rules or ('R1', 'R2', 'R3', 'R4', 'R5', 'R6')
        self.constraints = constraints or {'per_ring_quota': 2, 'min_graph_distance': 2, 'max_per_halogen': None}
        self.std_cfg = std_cfg or {'do_tautomer': False}
        self.qc_cfg = qc_cfg or {'sanitize_strict': True, 'pains': True, 'tautomer_canonicalize': False}
        self.pruning_cfg = pruning_cfg or {'enable_symmetry_fold': True, 'enable_state_sig': False}
        self.sugar_cfg = sugar_cfg or {'mode': 'heuristic', 'mask_exocyclic_oxygen': True, 'mask_glycosidic_bridge_oxygen': True, 'audit': False}
        self.symmetry_cfg = symmetry_cfg or {'compute_on_masked_subgraph': True}
        self.rules_cfg = rules_cfg or {
            'R2': {
                'sp2_CH_in_C_ring': True,       # R2a: Enable sp2 CH sites in C-ring
                'sp3_CH2_flavanone': True,      # R2b: Enable sp3 CH2 sites in flavanone C-ring
                'allowed_halogens': ['F', 'Cl', 'Br', 'I']  # Halogens allowed for R2 rules
            },
            'R6_methyl': {
                'enable': False,                # R6: Enable methyl halogenation rules
                'allowed': ['F', 'Cl'],         # Halogens allowed for R6 rules (only F, Cl supported)
                'allow_on_methoxy': False,      # Allow halogenation of O-CH3 groups
                'macro': {
                    'enable': False,            # Enable macro substitution (CF3, CCl3)
                    'labels': ['CF3', 'CCl3']   # Supported macro labels
                }
            }
        }

        # Initialize engine_cfg with backward compatibility migration
        if engine_cfg is None:
            # Check if rules_cfg contains legacy engine settings and migrate
            if 'engine' in self.rules_cfg:
                LOG.debug("Migrating engine settings from rules_cfg to engine_cfg for better organization")
                self.engine_cfg = self.rules_cfg['engine'].copy()
                # Remove from rules_cfg to avoid confusion
                del self.rules_cfg['engine']
            else:
                # Default engine configuration
                self.engine_cfg = {'budget_mode': 'ops', 'emit_legacy_keys': False, 'dedup_stage': 'pre', 'dedup_policy': 'auto'}
        else:
            self.engine_cfg = engine_cfg
            # Ensure dedup_policy has a default value
            if 'dedup_policy' not in self.engine_cfg:
                self.engine_cfg['dedup_policy'] = 'auto'

        self.random_seed = random_seed


def enumerate_products(parent_smi: str, cfg: EnumConfig, 
                      return_qa_stats: bool = False,
                      emit_summary_marker: bool = False,
                      aggregator=None,
                      stream_shape: str = 'legacy') -> Iterable[Dict[str, Any]]:
    """
    Enumerate k-dimensional halogenated products using BFS.
    
    Args:
        parent_smi: Parent molecule SMILES
        cfg: Enumeration configuration
        return_qa_stats: If True, yields (product, qa_stats) tuples instead of just products
        emit_summary_marker: If True and return_qa_stats=True, emits marker record when no products
    
    Yields:
        Product records as dictionaries, or (product, qa_stats) tuples if return_qa_stats=True
    """
    # Standardize parent molecule
    try:
        # Reset ring-label caches to avoid cross-test leakage
        from .sites import reset_ring_cache
        reset_ring_cache()
    except Exception:
        pass
    # Initialize QA stats dict for guard usage
    qa_guard_stats = {
        'qa_paths': {},
        'template_unsupported': 0
    }

    # T28: Ensure aggregator always exists in streaming mode for consistent data flow
    if return_qa_stats and aggregator is None:
        aggregator = QAAggregator(debug_consistency=cfg.engine_cfg.get('qa_debug_consistency', False))

    # H2-A: Set random seed for reproducible enumeration if provided
    if cfg.random_seed is not None:
        import random
        import numpy as np
        random.seed(cfg.random_seed)
        np.random.seed(cfg.random_seed)

        # Set RDKit random seed using centralized seed manager (debounced logging inside)
        from .rdkit_seed_utils import set_rdkit_random_seed
        rdkit_seed_report = set_rdkit_random_seed(cfg.random_seed)

        global _PYTHON_NUMPY_SEED_INFO_EMITTED, _PYTHON_NUMPY_SEED_LAST, _PYTHON_NUMPY_SEED_DEBUG_EMITTED
        seed_changed = (_PYTHON_NUMPY_SEED_LAST != cfg.random_seed)
        if (not _PYTHON_NUMPY_SEED_INFO_EMITTED) or seed_changed:
            LOG.info("Reproducibility seeds set: python=%d, numpy=%d", cfg.random_seed, cfg.random_seed)
            _PYTHON_NUMPY_SEED_INFO_EMITTED = True
            _PYTHON_NUMPY_SEED_LAST = cfg.random_seed
            _PYTHON_NUMPY_SEED_DEBUG_EMITTED = False
        elif not _PYTHON_NUMPY_SEED_DEBUG_EMITTED:
            LOG.debug("Reproducibility seeds unchanged since last report (debounced)")
            _PYTHON_NUMPY_SEED_DEBUG_EMITTED = True

    # Standardize parent molecule with RDKit guard
    parent_mol = None
    with rdkit_guard(qa_guard_stats):
        parent_mol = std_from_smiles(parent_smi, cfg.std_cfg.get('do_tautomer', False))
    
    if not parent_mol:
        # Handle early exit case - need to return guard events if requested
        if return_qa_stats and emit_summary_marker:
            # Create aggregator if needed for guard events
            if aggregator is None:
                aggregator = QAAggregator(debug_consistency=cfg.engine_cfg.get('qa_debug_consistency', False))
            
            # Inject guard events into aggregator as a single attempt
            qa_events = {}
            for metric, count in qa_guard_stats.get('qa_paths', {}).items():
                if count > 0:
                    qa_events[metric] = count
            
            tu = qa_guard_stats.get('template_unsupported', 0)
            if tu > 0:
                qa_events['template_unsupported'] = tu
            
            # Record guard events as template_unsupported failure
            if qa_events:
                aggregator.record_guard_failure('guard', 'rdkit', 1, qa_events)

            # Compute totals from aggregator with guard events
            agg_totals = _compute_totals_from_aggregator(aggregator)

            # Ensure consistent qa_paths key set and dual-write dedup fields
            emit_legacy = bool(cfg.engine_cfg.get('emit_legacy_keys', False))
            qa_paths_snap = agg_totals.get('qa_paths', {})
            from .schema import ensure_full_schema_qa_paths
            qa_paths_snap = ensure_full_schema_qa_paths(qa_paths_snap, emit_legacy_keys=emit_legacy)

            # Emit summary marker with guard events
            # k>1 only supports v2; this produces v2 empty pivots final snapshot
            final_marker = {'is_qa_summary_marker': True}
            snapshot = {
                'version': '2',
                'pivots': aggregator.to_pivots_dict(),
                'attempts': agg_totals.get('attempts', 0),
                'products': agg_totals.get('products', 0),
                'no_product_matches': agg_totals.get('no_product_matches', 0),
                'template_unsupported': agg_totals.get('template_unsupported', 0),
                'qa_paths': qa_paths_snap,
                # Dual-write: top-level dedup = qa_paths values for consistency
                'dedup_hits_statesig': int(qa_paths_snap.get('dedup_hits_statesig', 0)),
                'dedup_hits_inchi': int(qa_paths_snap.get('dedup_hits_inchi', 0))
            }
            yield (final_marker, snapshot)
        return

    # Get parent key and SMILES with RDKit guard
    parent_key = None
    parent_smiles = None
    with rdkit_guard(qa_guard_stats):
        parent_key = to_inchikey(parent_mol)
        parent_smiles = Chem.MolToSmiles(parent_mol, canonical=True)
    
    if not parent_key or not parent_smiles:
        # Handle early exit case - need to return guard events if requested
        if return_qa_stats and emit_summary_marker:
            # Create aggregator if needed for guard events
            if aggregator is None:
                aggregator = QAAggregator(debug_consistency=cfg.engine_cfg.get('qa_debug_consistency', False))
            
            # Inject guard events into aggregator as a single attempt
            qa_events = {}
            for metric, count in qa_guard_stats.get('qa_paths', {}).items():
                if count > 0:
                    qa_events[metric] = count
            
            tu = qa_guard_stats.get('template_unsupported', 0)
            if tu > 0:
                qa_events['template_unsupported'] = tu
            
            # Record guard events as template_unsupported failure
            if qa_events:
                aggregator.record_guard_failure('guard', 'rdkit', 1, qa_events)

            # Compute totals from aggregator with guard events
            agg_totals = _compute_totals_from_aggregator(aggregator)

            # Ensure consistent qa_paths key set and dual-write dedup fields
            emit_legacy = bool(cfg.engine_cfg.get('emit_legacy_keys', False))
            qa_paths_snap = agg_totals.get('qa_paths', {})
            from .schema import ensure_full_schema_qa_paths
            qa_paths_snap = ensure_full_schema_qa_paths(qa_paths_snap, emit_legacy_keys=emit_legacy)

            # Emit summary marker with guard events
            # k>1 only supports v2; this produces v2 empty pivots final snapshot
            final_marker = {'is_qa_summary_marker': True}
            snapshot = {
                'version': '2',
                'pivots': aggregator.to_pivots_dict(),
                'attempts': agg_totals.get('attempts', 0),
                'products': agg_totals.get('products', 0),
                'no_product_matches': agg_totals.get('no_product_matches', 0),
                'template_unsupported': agg_totals.get('template_unsupported', 0),
                'qa_paths': qa_paths_snap,
                # Dual-write: top-level dedup = qa_paths values for consistency
                'dedup_hits_statesig': int(qa_paths_snap.get('dedup_hits_statesig', 0)),
                'dedup_hits_inchi': int(qa_paths_snap.get('dedup_hits_inchi', 0))
            }
            yield (final_marker, snapshot)
        return

    # Global deduplication set for this parent
    seen_keys = {parent_key}

    # Compute sugar mask for this parent molecule
    sugar_mask, mask_degraded, status_metadata = get_sugar_mask_with_full_status(parent_mol, mode=cfg.sugar_cfg.get('mode', 'heuristic'), sugar_cfg=cfg.sugar_cfg)
    LOG.debug(f"Sugar mask computed: {len(sugar_mask)} atoms masked")

    # Sample-level QA event bus (not in product objects, at sample level)
    qa_bus = {
        "sugar_mask_filtered": 0,        # Sites/matches filtered by direct mask
        "sugar_proximity_filtered": 0,   # Sites filtered by proximity guard
        "post_guard_blocked": 0          # Products blocked by post-guard
    }

    # QA statistics tracking
    qa_stats = QAStats()

    # Enable aggregator implicitly for streaming to allow final v2 snapshot, while keeping
    # non-final yields in legacy shape for backward compatibility.
    # Create aggregator when returning QA stats to enable final v2 snapshot.
    # Per-record legacy shape does not depend on aggregator.
    if return_qa_stats and aggregator is None:
        aggregator = QAAggregator(debug_consistency=cfg.engine_cfg.get('qa_debug_consistency', False))

    # Record sugar mask degradation using unified counting
    if mask_degraded:
        qa_stats.qa_paths['sugar_mask_degraded'] = qa_stats.qa_paths.get('sugar_mask_degraded', 0) + 1
        if aggregator:
            aggregator.record_paths({'sugar_mask_degraded': 1})

    # BFS frontier: (mol, depth, history, sugar_mask, budget_payload)
    # Initialize with budget payload from engine configuration
    # P0 CONFIG FIX: Uses cfg.engine_cfg.budget_mode instead of hardcoded 'ops'
    initial_budget_payload = {
        'mode': cfg.engine_cfg.get('budget_mode', 'ops'),
        'k_ops': 0, 'k_atoms': 0,
        'site_tokens': {}, 'site_state': {}
    }
    frontier = deque([(parent_mol, 0, [], sugar_mask, initial_budget_payload)])
    
    # State signature tracking for pruning (optional)
    state_signatures = set() if cfg.pruning_cfg.get('enable_state_sig', False) else None
    
    while frontier:
        if len(frontier[0]) == 4:
            # Legacy format: (mol, depth, history, current_mask)
            mol, depth, history, current_mask = frontier.popleft()
            budget_payload = initial_budget_payload.copy()
        else:
            # New format: (mol, depth, history, current_mask, budget_payload)
            mol, depth, history, current_mask, budget_payload = frontier.popleft()

        # Stop if reached maximum depth
        if depth >= cfg.k_max:
            continue

        # Expand current molecule
        next_frontier, product_records, layer_qa_stats = _layer_expand(
            mol, depth, history, cfg, seen_keys, state_signatures, aggregator, current_mask, qa_bus, budget_payload
        )
        
        # Merge QA statistics from this layer
        qa_stats.merge(layer_qa_stats)
        
        # Yield all products from this expansion
        for record in product_records:
            if return_qa_stats:
                if stream_shape == 'v2':
                    agg_totals = _compute_totals_from_aggregator(aggregator) if aggregator else {
                        'attempts': 0, 'products': 0, 'no_product_matches': qa_stats.no_product_matches,
                        'template_unsupported': qa_stats.template_unsupported,
                        'qa_paths': qa_stats.qa_paths.copy()
                    }

                    # Ensure consistent qa_paths key set and dual-write dedup fields
                    emit_legacy = bool(cfg.engine_cfg.get('emit_legacy_keys', False))
                    if aggregator:
                        qa_paths_snap = _snapshot_qa_paths_from_aggregator(aggregator)
                    else:
                        qa_paths_snap = qa_stats.qa_paths.copy()

                    from .schema import ensure_full_schema_qa_paths
                    qa_paths_snap = ensure_full_schema_qa_paths(qa_paths_snap, emit_legacy_keys=emit_legacy)

                    snapshot = {
                        'version': '2',
                        'pivots': aggregator.to_pivots_dict() if aggregator else {'by_rule': {}, 'by_halogen': {}, 'by_k': {}, 'by_rule_halogen': {}, 'by_rule_halogen_k': {}},
                        'attempts': agg_totals.get('attempts', 0),
                        'products': agg_totals.get('products', 0),
                        'no_product_matches': agg_totals.get('no_product_matches', 0),
                        'template_unsupported': agg_totals.get('template_unsupported', 0),
                        'qa_paths': qa_paths_snap,
                        # Dual-write: top-level dedup = qa_paths values for consistency
                        'dedup_hits_statesig': int(qa_paths_snap.get('dedup_hits_statesig', 0)),
                        'dedup_hits_inchi': int(qa_paths_snap.get('dedup_hits_inchi', 0))
                    }
                    yield (record, snapshot)
                else:
                    # For non-v2 streaming, also use aggregator if available
                    emit_legacy = bool(cfg.engine_cfg.get('emit_legacy_keys', False))
                    qa_snapshot = qa_stats.to_dict()

                    if aggregator is not None:
                        qa_paths_snap = _snapshot_qa_paths_from_aggregator(aggregator)
                    else:
                        qa_paths_snap = qa_snapshot.get('qa_paths', {}).copy()

                    from .schema import ensure_full_schema_qa_paths
                    qa_paths_snap = ensure_full_schema_qa_paths(qa_paths_snap, emit_legacy_keys=emit_legacy)
                    qa_snapshot['qa_paths'] = qa_paths_snap

                    # Dual-write: top-level dedup = qa_paths values for consistency
                    qa_snapshot['dedup_hits_statesig'] = int(qa_paths_snap.get('dedup_hits_statesig', 0))
                    qa_snapshot['dedup_hits_inchi'] = int(qa_paths_snap.get('dedup_hits_inchi', 0))

                    yield (record, qa_snapshot)
            else:
                yield record
            
        # Add new molecules to frontier for next iteration
        frontier.extend(next_frontier)
    
    # If return_qa_stats was requested but no products were yielded,
    # yield a dummy record with final QA stats to ensure they're available
    if return_qa_stats and qa_stats.to_dict():
        # Check if we have any non-zero stats (indicating enumeration was attempted)
        stats_dict = qa_stats.to_dict()
        has_activity = any(v > 0 for v in stats_dict.values() if isinstance(v, int))
        has_qa_activity = any(v > 0 for v in stats_dict.get('qa_paths', {}).values())
        
        if emit_summary_marker and (has_activity or has_qa_activity):
            # Create a special marker record to carry final QA stats (v2 snapshot)
            final_marker = {
                'is_qa_summary_marker': True,
                'parent_smiles': parent_smiles,
                'parent_inchikey': parent_key
            }
            agg_totals = _compute_totals_from_aggregator(aggregator) if aggregator else {
                'attempts': 0, 'products': 0, 'no_product_matches': qa_stats.no_product_matches,
                'template_unsupported': qa_stats.template_unsupported,
                'qa_paths': qa_stats.qa_paths.copy()
            }

            # Ensure consistent qa_paths key set and dual-write dedup fields
            emit_legacy = bool(cfg.engine_cfg.get('emit_legacy_keys', False))
            if aggregator:
                qa_paths_snap = _snapshot_qa_paths_from_aggregator(aggregator)
            else:
                qa_paths_snap = qa_stats.qa_paths.copy()

            from .schema import ensure_full_schema_qa_paths
            qa_paths_snap = ensure_full_schema_qa_paths(qa_paths_snap, emit_legacy_keys=emit_legacy)

            snapshot = {
                'version': '2',
                'pivots': aggregator.to_pivots_dict() if aggregator else {'by_rule': {}, 'by_halogen': {}, 'by_k': {}, 'by_rule_halogen': {}, 'by_rule_halogen_k': {}},
                'attempts': agg_totals.get('attempts', 0),
                'products': agg_totals.get('products', 0),
                'no_product_matches': agg_totals.get('no_product_matches', 0),
                'template_unsupported': agg_totals.get('template_unsupported', 0),
                'qa_paths': qa_paths_snap,
                # Dual-write: top-level dedup = qa_paths values for consistency
                'dedup_hits_statesig': int(qa_paths_snap.get('dedup_hits_statesig', 0)),
                'dedup_hits_inchi': int(qa_paths_snap.get('dedup_hits_inchi', 0))
            }
            yield (final_marker, snapshot)
    
    # Add guard events to aggregator before final logging
    if aggregator:
        # Combine all guard events into a single attempt
        qa_events = {}
        for metric, count in qa_guard_stats.get('qa_paths', {}).items():
            if count > 0:
                qa_events[metric] = count
        
        tu = qa_guard_stats.get('template_unsupported', 0)
        if tu > 0:
            qa_events['template_unsupported'] = tu
        
        # Record guard events as template_unsupported failure
        if qa_events:
            aggregator.record_guard_failure('guard', 'rdkit', 1, qa_events)
    
    # Merge qa_bus events into qa_stats before final logging/return
    if qa_stats:
        qa_paths = qa_stats.qa_paths
        for k, v in qa_bus.items():
            qa_paths[k] = qa_paths.get(k, 0) + int(v)

        # Backward compatibility: maintain sugar_post_guard_blocked for legacy tools
        # Controlled by engine_cfg.emit_legacy_keys (default False for cleaner new installations)
        emit_legacy = bool(cfg.engine_cfg.get('emit_legacy_keys', False))
        if emit_legacy and "post_guard_blocked" in qa_paths:
            qa_paths["sugar_post_guard_blocked"] = (
                qa_paths.get("sugar_post_guard_blocked", 0) + qa_paths["post_guard_blocked"]
            )


    # Also record qa_bus events in the aggregator for enumerate_with_stats
    # Only record non-pivot events - pivot events should come from pivot aggregation
    if aggregator and qa_bus:
        qa_bus_events = {
            k: v for k, v in qa_bus.items()
            if v > 0 and canonical_event(k) not in PIVOT_EVENT_KEYS
        }
        if qa_bus_events:
            aggregator.record_paths(qa_bus_events)

    # Log final QA summary
    qa_stats.log_summary()


def enumerate_with_stats(parent_smi: str, cfg: EnumConfig) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Enumerate products and return final QA statistics separately.
    
    This provides a stable interface where products and QA stats are cleanly separated,
    avoiding the streaming cumulative snapshots of enumerate_products(..., return_qa_stats=True).
    
    Args:
        parent_smi: Parent molecule SMILES
        cfg: Enumeration configuration
        
    Returns:
        (products_list, final_qa_stats_dict) - stats dict includes version and pivots for M2
    """
    # Reset ring-label caches to avoid cross-test leakage
    try:
        from .sites import reset_ring_cache
        reset_ring_cache()
    except Exception:
        pass
    products = []
    final_qa_stats = {
        'no_product_matches': 0,
        'dedup_hits_statesig': 0,
        'dedup_hits_inchi': 0,
        'template_unsupported': 0,
        'qa_paths': empty_qa_paths()
    }  # Initialize with stable schema
    
    # Initialize pivot aggregator for M2 statistics
    aggregator = QAAggregator(debug_consistency=cfg.engine_cfg.get('qa_debug_consistency', False))
    
    for item in enumerate_products(parent_smi, cfg, return_qa_stats=True, emit_summary_marker=True, aggregator=aggregator):
        if isinstance(item, tuple) and len(item) == 2:
            product, qa_stats_snapshot = item
            # Check if this is a QA summary marker (not a real product)
            if product.get('is_qa_summary_marker', False):
                # This is just a marker to carry final QA stats, don't add to products
                final_qa_stats = qa_stats_snapshot
            else:
                # Regular product - aggregator events already recorded in enumerate_products
                products.append(product)
                final_qa_stats = qa_stats_snapshot  # Keep overwriting to get final state
        else:
            # Fallback for non-QA mode (shouldn't happen with return_qa_stats=True)
            products.append(item)
    
    # Compute final stats with totals derived from aggregator (eliminates totals-pivots inconsistency)
    aggregator_totals = _compute_totals_from_aggregator(aggregator)
    
    # Build final stats dict with aggregator-derived totals and version 2 format
    # Ensure consistent qa_paths key set and backward compatibility before final return
    from .schema import ensure_full_schema_qa_paths
    compatible_qa_paths = ensure_full_schema_qa_paths(
        aggregator_totals.get('qa_paths', empty_qa_paths()),
        emit_legacy_keys=cfg.engine_cfg.get('emit_legacy_keys', False)
    )

    final_qa_stats_v2 = {
        'version': '2',
        'pivots': aggregator.to_pivots_dict(),
        # Core totals from aggregator (stable schema - always present)
        'attempts': aggregator_totals.get('attempts', 0),
        'products': aggregator_totals.get('products', 0),
        'no_product_matches': aggregator_totals.get('no_product_matches', 0),
        'template_unsupported': aggregator_totals.get('template_unsupported', 0),
        'qa_paths': compatible_qa_paths,
        # Dedup stats: prefer aggregator values for consistency, fall back to streaming stats
        'dedup_hits_statesig': compatible_qa_paths.get('dedup_hits_statesig', final_qa_stats.get('dedup_hits_statesig', 0)),
        'dedup_hits_inchi': compatible_qa_paths.get('dedup_hits_inchi', final_qa_stats.get('dedup_hits_inchi', 0))
    }
    
    # Check totals-pivots consistency with new unified approach
    _check_totals_pivots_consistency(
        totals_qa_paths=final_qa_stats_v2.get('qa_paths', {}),
        pivot_sums=dict(aggregator.qa_paths_pivot_sum),
        enable_check=cfg.engine_cfg.get('qa_debug_consistency', False)  # Default disabled for production
    )

    # T29: Optional cache hit monitoring for debugging
    try:
        from .sites import log_sugar_ring_cache_stats
        log_sugar_ring_cache_stats()
    except Exception:
        pass

    return products, final_qa_stats_v2


def _layer_expand(mol, depth: int, history: List[Dict[str, Any]],
                  cfg: EnumConfig, seen_global: set,
                  state_sigs: Optional[set] = None, aggregator=None, sugar_mask: Optional[set] = None, qa_bus: Optional[Dict] = None, budget_payload: Optional[dict] = None) -> Tuple[List[Tuple], List[Dict[str, Any]], Dict[str, Any]]:
    """
    Expand one layer of the BFS tree.
    
    Returns:
        (next_frontier, product_records, layer_qa_stats)
    """
    next_frontier = []
    product_records = []
    
    # Aggregate QA statistics for this layer
    layer_qa_stats = {
        'no_product_matches': 0,
        'dedup_hits_statesig': 0,
        'dedup_hits_inchi': 0,
        'template_unsupported': 0,
        'qa_paths': empty_qa_paths()
    }
    
    # Build reactions once (cached for performance)
    reactions = _build_reactions_cached(RULES_VERSION)
    
    # Default to empty mask if not provided
    if sugar_mask is None:
        sugar_mask = set()

    # Default to empty qa_bus if not provided
    if qa_bus is None:
        qa_bus = {"sugar_mask_filtered": 0, "post_guard_blocked": 0}

    # 1) Apply reaction-type rules (R3, R4, R5)
    reaction_rules = [r for r in cfg.rules if r in ('R3', 'R4', 'R5')]
    for rule_id in reaction_rules:
        reaction_products, rule_qa_stats = _apply_reaction_rule(
            mol, rule_id, reactions, cfg, history, depth, seen_global, aggregator, sugar_mask, qa_bus, budget_payload
        )
        for prod_mol, new_history, record, payload in reaction_products:
            # Compute new sugar mask for the product molecule
            new_mask, mask_degraded = get_sugar_mask_with_status(prod_mol, mode=cfg.sugar_cfg.get('mode', 'heuristic'), sugar_cfg=cfg.sugar_cfg)
            if mask_degraded:
                layer_qa_stats['qa_paths']['sugar_mask_degraded'] = layer_qa_stats['qa_paths'].get('sugar_mask_degraded', 0) + 1

            # Apply post-guard check for sugar masking violations using product mask
            if post_guard_blocked(prod_mol, new_mask):
                LOG.debug(f"Product blocked by sugar mask post-guard: {rule_id}")
                qa_bus["post_guard_blocked"] += 1
                continue

            product_records.append(record)
            # For reaction rules (R3,R4,R5), pass budget payload unchanged
            next_frontier.append((prod_mol, depth + 1, new_history, new_mask, payload))

        # Merge rule QA stats into layer stats (for backward compatibility with totals)
        for key in ['no_product_matches', 'dedup_hits_statesig', 'dedup_hits_inchi', 'template_unsupported']:
            count = rule_qa_stats.get(key, 0)
            layer_qa_stats[key] += count

        if 'qa_paths' in rule_qa_stats:
            for qa_key, qa_count in rule_qa_stats['qa_paths'].items():
                if qa_key in layer_qa_stats['qa_paths']:
                    layer_qa_stats['qa_paths'][qa_key] += qa_count

    # 2) Apply site-type rules with symmetry folding (R1, R2, R6)
    site_products = _apply_site_rules(mol, cfg, history, depth, seen_global, aggregator, sugar_mask, qa_bus, budget_payload)
    for product_tuple in site_products:
        if len(product_tuple) == 4:
            # R6 format: (prod_mol, new_history, record, new_budget_payload)
            prod_mol, new_history, record, new_budget_payload = product_tuple
        else:
            # R1/R2 format: (prod_mol, new_history, record)
            prod_mol, new_history, record = product_tuple
            new_budget_payload = budget_payload  # Pass through unchanged

        # Compute new sugar mask for the product molecule
        new_mask, mask_degraded = get_sugar_mask_with_status(prod_mol, mode=cfg.sugar_cfg.get('mode', 'heuristic'), sugar_cfg=cfg.sugar_cfg)
        if mask_degraded:
            layer_qa_stats['qa_paths']['sugar_mask_degraded'] = layer_qa_stats['qa_paths'].get('sugar_mask_degraded', 0) + 1

        # Apply post-guard check for sugar masking violations using product mask
        if post_guard_blocked(prod_mol, new_mask):
            LOG.debug(f"Site product blocked by sugar mask post-guard")
            qa_bus["post_guard_blocked"] += 1
            continue

        product_records.append(record)
        next_frontier.append((prod_mol, depth + 1, new_history, new_mask, new_budget_payload))
    
    return next_frontier, product_records, layer_qa_stats


def _find_reaction_matches(mol, rule_id: str, reactions: Dict, halogen: str, sugar_mask: Optional[set] = None) -> Tuple[List[Tuple[int, Tuple]], int, int]:
    """
    Find reaction matches with site index and full match tuple.

    Assumption: All R3/R4/R5 templates have exactly one reactant template.
    For multi-reactant templates, this function will issue a warning and return empty matches.

    Returns:
        (matches_with_sites, template_unsupported_count, sugar_filtered_count)
    """
    matches_with_sites = []
    template_unsupported = 0
    sugar_filtered_count = 0

    if rule_id not in reactions or halogen not in reactions[rule_id]:
        return matches_with_sites, template_unsupported, sugar_filtered_count
    
    rxn = reactions[rule_id][halogen]
    try:
        # Get the reactant template (first template)
        reactant_templates = rxn.GetReactants()
        if not reactant_templates:
            return matches_with_sites, template_unsupported
        
        # Robustness check: Warn if template has multiple reactants
        if len(reactant_templates) != 1:
            LOG.warning(f"Rule {rule_id} with {halogen} has {len(reactant_templates)} reactant templates, "
                          f"expected exactly 1. Skipping isotope-based matching, will use fallback.")
            template_unsupported = 1  # Count this as unsupported template attempt
            return matches_with_sites, template_unsupported
            
        reactant_template = reactant_templates[0]
        
        # Find all matches for the reactant pattern
        matches = mol.GetSubstructMatches(reactant_template, uniquify=True)
        
        for match in matches:
            # Find the atom with mapping number 1 (the substitution site)
            site_atom_idx = None
            for template_idx, mol_idx in enumerate(match):
                template_atom = reactant_template.GetAtomWithIdx(template_idx)
                if template_atom.GetAtomMapNum() == 1:
                    site_atom_idx = mol_idx
                    break
            
            if site_atom_idx is not None:
                matches_with_sites.append((site_atom_idx, match))

    except Exception:
        # If SMIRKS matching fails, return empty list to trigger fallback
        pass

    # Filter out matches targeting masked atoms
    if sugar_mask:
        original_count = len(matches_with_sites)
        matches_with_sites = [(site_idx, match) for site_idx, match in matches_with_sites
                             if not _match_hits_mask(match, sugar_mask)]
        filtered_count = len(matches_with_sites)
        if original_count > filtered_count:
            sugar_filtered_count = original_count - filtered_count
            LOG.debug(f"Reaction {rule_id}_{halogen}: filtered {sugar_filtered_count} "
                        f"matches intersecting with masked atoms")

    return matches_with_sites, template_unsupported, sugar_filtered_count


def _match_hits_mask(match_tuple, mask_set: Set[int]) -> bool:
    """
    Check if any atom in a reaction match intersects with the mask.

    Args:
        match_tuple: Tuple of atom indices from reaction pattern matching
        mask_set: Set of atom indices that are masked

    Returns:
        True if any atom in the match is in the mask
    """
    try:
        # Convert match tuple to atom indices if needed
        atom_indices = []
        for item in match_tuple:
            if isinstance(item, int):
                atom_indices.append(item)
            else:
                # Handle RDKit atom objects
                try:
                    atom_indices.append(item.GetIdx())
                except AttributeError:
                    # If it's not an atom object, try to convert directly
                    atom_indices.append(int(item))

        # Check if any atom in match is in mask
        return any(idx in mask_set for idx in atom_indices)

    except Exception as e:
        LOG.debug(f"Match-mask intersection check failed: {e}")
        # Conservative: block on failure
        return True


def _find_isotope_tagged_site(product_mol, halogen: str, isotope_tag: int) -> Optional[int]:
    """Find carbon with isotope tag that is adjacent to the new halogen."""
    try:
        # Look for halogen atoms in the product
        for atom in product_mol.GetAtoms():
            if atom.GetSymbol() == halogen:
                # Check neighbors of this halogen
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C' and neighbor.GetIsotope() == isotope_tag:
                        return neighbor.GetIdx()
    except Exception:
        pass
    return None


def _clear_isotope_tags(mol) -> None:
    """Clear only our internal isotope tags from molecule, preserving real isotopes."""
    try:
        for atom in mol.GetAtoms():
            if atom.GetIsotope() == ISOTOPE_TAG:
                atom.SetIsotope(0)
    except Exception:
        pass


def _detect_site_from_product_with_mapping_fallback(mol, product_mol, 
                                                    rule_id: str, halogen: str, 
                                                    fallback_stats: Optional[Dict[str, int]] = None) -> Tuple[Optional[int], Optional[str]]:
    """Enhanced site detection: try AtomMapNum first, then heuristic fallback."""
    site_atom_idx = None
    ring_tag = None
    
    try:
        # Primary: Look for atoms with mapping number 1 (before clearing)
        for atom in product_mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num == 1:
                # This is the carbon atom where substitution occurred
                site_atom_idx = atom.GetIdx()
                # Get ring tag from original molecule using same atom index
                if site_atom_idx is not None and site_atom_idx < mol.GetNumAtoms():
                    ring_tag = flavonoid_ring_label(mol, site_atom_idx)
                # Clear mapping number after reading
                atom.SetAtomMapNum(0)
                break
        
        # Clear all remaining mapping numbers
        for atom in product_mol.GetAtoms():
            atom.SetAtomMapNum(0)
            
    except Exception:
        pass
    
    # Fallback: Use heuristic detection if mapping failed
    if site_atom_idx is None:
        if fallback_stats is not None:
            fallback_stats['heuristic_used'] = fallback_stats.get('heuristic_used', 0) + 1
        site_atom_idx, ring_tag = _detect_site_from_product(mol, product_mol, rule_id, halogen)
    else:
        # AtomMapNum was successful
        if fallback_stats is not None:
            fallback_stats['atommap_used'] = fallback_stats.get('atommap_used', 0) + 1
    
    return site_atom_idx, ring_tag



def _detect_site_from_product(mol, product_mol, rule_id: str, halogen: str) -> Tuple[int, str]:
    """Fallback: detect substitution site from product molecule (heuristic approach)."""
    site_atom_idx = None
    ring_tag = None
    
    # Find carbon atoms connected to the new halogen atom
    for atom in product_mol.GetAtoms():
        if atom.GetSymbol() == halogen:
            # Find the carbon this halogen is bonded to
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    candidate_idx = neighbor.GetIdx()
                    # Verify this carbon exists in original molecule
                    if candidate_idx < mol.GetNumAtoms():
                        if rule_id == 'R3':
                            # For R3, find the carbon that was bonded to -OH in original
                            oh_carbons = []
                            for orig_idx in range(mol.GetNumAtoms()):
                                orig_at = mol.GetAtomWithIdx(orig_idx)
                                if orig_at.GetSymbol() == 'C':
                                    has_oh = any(n.GetSymbol() == 'O' and n.GetTotalNumHs() > 0 
                                                for n in orig_at.GetNeighbors())
                                    if has_oh:
                                        oh_carbons.append(orig_idx)
                            if oh_carbons:
                                site_atom_idx = oh_carbons[0]
                                ring_tag = flavonoid_ring_label(mol, site_atom_idx)
                                break
                        elif rule_id == 'R4':
                            # For R4, find the carbon that was bonded to -NH in original
                            nh_carbons = []
                            for orig_idx in range(mol.GetNumAtoms()):
                                orig_at = mol.GetAtomWithIdx(orig_idx)
                                if orig_at.GetSymbol() == 'C':
                                    has_nh = any(n.GetSymbol() == 'N' for n in orig_at.GetNeighbors())
                                    if has_nh:
                                        nh_carbons.append(orig_idx)
                            if nh_carbons:
                                site_atom_idx = nh_carbons[0]
                                ring_tag = flavonoid_ring_label(mol, site_atom_idx)
                                break
                        elif rule_id == 'R5':
                            # For R5, find the carbon that was bonded to carboxyl group
                            cooh_carbons = []
                            for orig_idx in range(mol.GetNumAtoms()):
                                orig_at = mol.GetAtomWithIdx(orig_idx)
                                if orig_at.GetSymbol() == 'C':
                                    # Check if connected to carboxyl carbon
                                    for neighbor in orig_at.GetNeighbors():
                                        if neighbor.GetSymbol() == 'C':
                                            # Check if this neighbor carbon has both =O and -OH
                                            has_carbonyl = False
                                            has_oh = False
                                            for n2 in neighbor.GetNeighbors():
                                                if n2.GetSymbol() == 'O':
                                                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                                                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                                                        has_carbonyl = True
                                                    elif n2.GetTotalNumHs() > 0:
                                                        has_oh = True
                                            if has_carbonyl and has_oh:
                                                cooh_carbons.append(orig_idx)
                                                break
                            if cooh_carbons:
                                site_atom_idx = cooh_carbons[0]
                                ring_tag = flavonoid_ring_label(mol, site_atom_idx)
                                break
                    if site_atom_idx is not None:
                        break
            if site_atom_idx is not None:
                break
                
    return site_atom_idx, ring_tag


def _expand_mask_by_radius(mol, mask: set, radius: int) -> set:
    """
    Precompute all atoms within 'radius' bonds of any atom in 'mask'.
    This enables O(1) proximity checks instead of per-site BFS.
    """
    if radius <= 0 or not mask:
        return set(mask or ())
    visited = set(mask)
    frontier = list(mask)
    for _ in range(radius):
        nxt = []
        for aidx in frontier:
            atom = mol.GetAtomWithIdx(aidx)
            for nbr in atom.GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx in visited:
                    continue
                visited.add(nidx)
                nxt.append(nidx)
        frontier = nxt
        if not frontier:
            break
    return visited


def _site_proximity_blocked(mol, site_idx: int, sugar_mask: set, radius: int) -> bool:
    """
    DEPRECATED: Return True if 'site_idx' is within 'radius' bonds of any sugar-masked atom.

    This function is deprecated and replaced by precomputed expanded_mask membership checks
    in _apply_proximity_guard for better performance. Will be removed in a future release.

    Use _expand_mask_by_radius() + membership check instead.
    """
    import warnings
    warnings.warn(
        "_site_proximity_blocked is deprecated. Use _expand_mask_by_radius + membership check instead.",
        DeprecationWarning,
        stacklevel=2
    )

    if radius <= 0 or not sugar_mask:
        return False
    if site_idx in sugar_mask:
        return True
    visited = {site_idx}
    frontier = [site_idx]
    for _ in range(radius):
        nxt = []
        for aidx in frontier:
            atom = mol.GetAtomWithIdx(aidx)
            for nbr in atom.GetNeighbors():
                nidx = nbr.GetIdx()
                if nidx in visited:
                    continue
                if nidx in sugar_mask:
                    return True
                visited.add(nidx)
                nxt.append(nidx)
        frontier = nxt
        if not frontier:
            break
    return False


def _apply_proximity_guard(sites: List[int], expanded_mask: set, qa_bus: dict, rule_id: str, key: str = "sugar_proximity_filtered") -> List[int]:
    """Apply proximity guard to filter sites near expanded mask atoms."""
    if not expanded_mask:
        return sites

    kept, dropped = [], 0
    dropped_sample = []
    for idx in sites:
        if idx in expanded_mask:
            qa_bus[key] += 1
            dropped += 1
            if len(dropped_sample) < 5:  # Sample up to 5 dropped sites for debugging
                dropped_sample.append(idx)
        else:
            kept.append(idx)

    LOG.debug(f"{rule_id} proximity-guard: dropped={dropped}, kept={len(kept)}")
    if dropped and LOG.isEnabledFor(logging.DEBUG) and dropped_sample:
        LOG.debug(f"{rule_id} proximity-guard dropped sample sites: {dropped_sample}")
    return kept


def _sites_with_mask_drop(site_fn, *, mol, sugar_mask, qa_bus, **kwargs):
    """Call site function twice: raw(no masking) vs masked(with masking), record difference to sugar_mask_filtered."""
    masked_sites = site_fn(mol, masked_atoms=sugar_mask, **kwargs)
    if sugar_mask:
        raw_sites = site_fn(mol, masked_atoms=set(), **kwargs)
        drop = max(0, len(raw_sites) - len(masked_sites))
        qa_bus["sugar_mask_filtered"] += drop
        LOG.debug(f"Site masking: {site_fn.__name__} raw={len(raw_sites)}, masked={len(masked_sites)}, drop={drop}")
    return masked_sites


def _apply_reaction_rule(mol, rule_id: str, reactions: Dict, cfg: EnumConfig,
                        history: List[Dict[str, Any]], depth: int, seen_global: set, aggregator=None, sugar_mask: Optional[set] = None, qa_bus: Optional[Dict] = None, budget_payload: Optional[dict] = None) -> Tuple[List[Tuple], Dict[str, Any]]:
    """
    Apply reaction-based rule (R3, R4, R5) with proper attempt boundary semantics.
    
    ATTEMPT BOUNDARY: Each rule x halogen x k combination is exactly ONE attempt.
    - attempts=1 per rule x halogen x k (regardless of number of matches)
    - products=1 if attempt yields >=1 valid products, 0 otherwise  
    - no_product_matches=1 if attempt yields 0 valid products, 0 otherwise
    
    Returns (results, qa_stats).
    """
    results = []

    # Default to empty mask if not provided
    if sugar_mask is None:
        sugar_mask = set()

    # Default to empty qa_bus if not provided
    if qa_bus is None:
        qa_bus = {"sugar_mask_filtered": 0, "post_guard_blocked": 0}

    if rule_id not in reactions:
        return results, {
            'dedup_hits_statesig': 0, 
            'dedup_hits_inchi': 0, 
            'no_product_matches': 0,
            'template_unsupported': 0,
            'qa_paths': empty_qa_paths()
        }
    
    # Track aggregate attempt outcomes for return structure
    no_product_matches = 0
    template_unsupported = 0
    
    for halogen in cfg.halogens:
        if halogen not in reactions[rule_id]:
            # Record attempt with template_unsupported outcome if aggregator available
            if aggregator:
                qa_events = {'template_unsupported': 1}
                aggregator.record_attempt_result(rule_id, halogen, depth + 1, produced_count=0, qa_events_dict=qa_events,
                                               k_ops=None, k_atoms=None)
            template_unsupported += 1
            continue
            
        # This is a real attempt for this rule/halogen combination
        rxn = reactions[rule_id][halogen]
        attempt_qa_events = {}
        attempt_produced_any = False
        isotope_unavailable_occurred = False
        atommap_used_in_attempt = False
        heuristic_used_in_attempt = False
        
        try:
            # Get all matches with site information
            matches_with_sites, template_unsupported_count, sugar_filtered_count = _find_reaction_matches(mol, rule_id, reactions, halogen, sugar_mask)
            if template_unsupported_count > 0:
                attempt_qa_events['template_unsupported'] = template_unsupported_count
                template_unsupported += template_unsupported_count

            if sugar_filtered_count > 0:
                qa_mark(qa_bus, attempt_qa_events, 'sugar_mask_filtered', sugar_filtered_count)
            
            if not matches_with_sites:
                # Isotope strategy unavailable for this rule/halogen attempt
                isotope_unavailable_occurred = True
                
                # Fallback: Use traditional site detection if pre-matching failed
                reaction_products = _run_reaction_safely(rxn, (mol,), rule_id, halogen, qa_bus, aggregator=None, current_k=depth+1)
                for product_mol in _iter_reaction_mols(reaction_products):
                    if product_mol is None:
                        continue
                        
                    # Use fallback site detection (isotope tagging failed)
                    # Create temporary dict to capture qa events for this product
                    product_qa_events = {}
                    site_atom_idx, ring_tag = _detect_site_from_product_with_mapping_fallback(
                        mol, product_mol, rule_id, halogen, product_qa_events
                    )
                    
                    # Merge product QA events into attempt tracking (no sentinel key pollution)
                    for qa_key, count in product_qa_events.items():
                        if count > 0:
                            if qa_key == 'atommap_used':
                                atommap_used_in_attempt = True
                            elif qa_key == 'heuristic_used':
                                heuristic_used_in_attempt = True
                            else:
                                # Accumulative semantics for other QA events
                                qa_mark(qa_bus, attempt_qa_events, qa_key, count)
                    
                    if site_atom_idx is None:
                        continue
                    
                    # Process this product
                    final_prod = _process_reaction_product(
                        product_mol, mol, site_atom_idx, ring_tag, rule_id, halogen,
                        history, depth, cfg, seen_global, qa_bus, budget_payload,
                        attempt=attempt_qa_events
                    )
                    if final_prod:
                        results.append(final_prod)
                        attempt_produced_any = True
            else:
                # Process each match independently with isotope tagging
                for site_atom_idx, match in matches_with_sites:
                    ring_tag = flavonoid_ring_label(mol, site_atom_idx)
                    
                    # Create a copy and add isotope tag to the site carbon
                    mol_copy = Chem.RWMol(mol)
                    site_carbon = mol_copy.GetAtomWithIdx(site_atom_idx)
                    site_carbon.SetIsotope(ISOTOPE_TAG)
                    
                    # Run reaction on tagged molecule
                    try:
                        reaction_products = rxn.RunReactants((mol_copy,))
                        product_found = False
                        
                        for product_mol in _iter_reaction_mols(reaction_products):
                            if product_mol is None:
                                continue
                                
                            # Find product corresponding to this tagged site
                            tagged_site_in_product = _find_isotope_tagged_site(
                                product_mol, halogen, ISOTOPE_TAG
                            )
                            
                            if tagged_site_in_product is not None:
                                # Clear isotope tags from product
                                _clear_isotope_tags(product_mol)
                                
                                # Process this product with the known site
                                final_prod = _process_reaction_product(
                                    product_mol, mol, site_atom_idx, ring_tag, rule_id, halogen,
                                    history, depth, cfg, seen_global, qa_bus, budget_payload,
                                    attempt=attempt_qa_events
                                )
                                if final_prod:
                                    results.append(final_prod)
                                    attempt_produced_any = True
                                    product_found = True
                                    break
                        
                        if not product_found:
                            # Isotope match failed to produce valid product
                            qa_mark(qa_bus, attempt_qa_events, 'isotope_miss')
                            no_product_matches += 1
                            
                    except Exception:
                        # If isotope tagging fails, treat as unavailable strategy
                        # Note: This is NOT a product match failure, just strategy unavailable
                        isotope_unavailable_occurred = True
                        continue
                        
        except Exception as e:
            # Log high-level rule failures (already deduplicated in _run_reaction_safely for reaction-specific failures)
            warning_key = f"{rule_id}_{halogen}_toplevel"
            if _reaction_warning_counts[warning_key] < _max_warnings_per_rule_halogen:
                # Use module-level logger
                LOG.warning(f"Rule {rule_id} with {halogen} failed at top level: {str(e)}")
                _reaction_warning_counts[warning_key] += 1
            attempt_qa_events['template_unsupported'] = attempt_qa_events.get('template_unsupported', 0) + 1
            template_unsupported += 1
        
        # Apply binary semantics for all binary QA events (at most 1 per attempt)
        if isotope_unavailable_occurred:
            qa_mark(qa_bus, attempt_qa_events, 'isotope_unavailable')
        
        if atommap_used_in_attempt:
            qa_mark(qa_bus, attempt_qa_events, 'atommap_used')

        if heuristic_used_in_attempt:
            qa_mark(qa_bus, attempt_qa_events, 'heuristic_used')
        
        # Record this attempt result with aggregator if available
        if aggregator:
            aggregator.record_attempt_result(rule_id, halogen, depth + 1,
                                           produced_count=(1 if attempt_produced_any else 0),
                                           qa_events_dict=attempt_qa_events,
                                           k_ops=None, k_atoms=None)
    
    # Build final qa_stats from qa_bus only (unified counting)
    final_qa_stats = {
        'dedup_hits_statesig': qa_bus.get('dedup_hits_statesig', 0),
        'dedup_hits_inchi': qa_bus.get('dedup_hits_inchi', 0),
        'no_product_matches': no_product_matches,
        'template_unsupported': template_unsupported,
        'qa_paths': {
            'isotope_unavailable': qa_bus.get('isotope_unavailable', 0),
            'isotope_miss': qa_bus.get('isotope_miss', 0),
            'atommap_used': qa_bus.get('atommap_used', 0),
            'heuristic_used': qa_bus.get('heuristic_used', 0),
            'sugar_mask_filtered': qa_bus.get('sugar_mask_filtered', 0),
            'sugar_mask_degraded': qa_bus.get('sugar_mask_degraded', 0),
            'post_guard_blocked': qa_bus.get('post_guard_blocked', 0),
            'rdkit_error': qa_bus.get('rdkit_error', 0),
        }
    }
    return results, final_qa_stats


def _process_reaction_product(product_mol, original_mol, site_atom_idx: int,
                            ring_tag: str, rule_id: str, halogen: str, history: List[Dict[str, Any]],
                            depth: int, cfg: EnumConfig, seen_global: set,
                            qa_bus: Optional[Dict[str, int]] = None, budget_payload: Optional[dict] = None,
                            attempt: Optional[Dict[str, int]] = None) -> Optional[Tuple]:
    """Process a single reaction product with correct deduplication strategy.

    P0 RETURN FORMAT: Returns 4-tuple (mol, history, record, budget_payload)
    for consistency with site rules.
    """
    try:
        # Sanitize the product  
        Chem.SanitizeMol(product_mol)
        final_prod = product_mol
    except (Chem.AtomValenceException, Chem.KekulizeException, ValueError) as e:
        # Common RDKit sanitization failures
        LOG.debug(f"Molecule sanitization failed: {type(e).__name__}: {e}")
        return None
    except Exception as e:
        # Catch any other unexpected errors during sanitization
        LOG.debug(f"Unexpected sanitization error: {type(e).__name__}: {e}")
        return None
    
    # Build history item with site and ring_tag
    history_item = {
        'rule': rule_id,
        'site': site_atom_idx,
        'sym': None,
        'ring_tag': ring_tag if ring_tag is not None else '',
        'halogen': halogen,
        'depth': depth + 1
    }
    new_history = history + [history_item]
    
    # Check constraints
    ok, violations = accept_constraints(final_prod, new_history, cfg.constraints)
    if not ok:
        return None

    # Configurable deduplication policy
    policy = cfg.engine_cfg.get('dedup_policy', 'auto')
    if policy == 'auto':
        policy = 'stable_key' if cfg.pruning_cfg.get('enable_symmetry_fold', True) else 'state_sig'

    # Always compute sanitized key for the record, independent of dedup policy
    product_inchikey = None

    # Prepare state signature only if needed
    state_sig = None
    if policy == 'state_sig':
        state_sig = compute_state_signature(new_history)

    # Check for duplicates using the configured policy
    dedup_key, is_dup = early_check(
        final_prod, seen_global, qa_bus,
        metric='dedup_hits_inchi',
        policy=policy,
        state_sig=state_sig,
        attempt=attempt
    )
    if is_dup:
        # When using state_sig, also bump statesig counter for clarity
        if policy == 'state_sig':
            qa_mark(qa_bus, attempt, 'dedup_hits_statesig')
        return None

    # Commit deduplication key after all validations pass
    commit(dedup_key, seen_global)

    # Compute record inchikey (stable) if still None
    if product_inchikey is None:
        product_inchikey = compute_stable_key(final_prod) or "UNKNOWN"

    # Create product record with sugar audit (if enabled)
    product_mask, mask_degraded = get_sugar_mask_with_status(final_prod, mode=cfg.sugar_cfg.get('mode', 'heuristic'), sugar_cfg=cfg.sugar_cfg)
    if mask_degraded and qa_bus is not None:
        qa_bus['sugar_mask_degraded'] = qa_bus.get('sugar_mask_degraded', 0) + 1
    if cfg.sugar_cfg.get('audit', False):
        sugar_audit = compute_sugar_audit_fields(final_prod, product_mask, cfg.sugar_cfg)
    else:
        sugar_audit = None
    record = _make_record(final_prod, original_mol, depth + 1, rule_id, halogen, new_history, True, sugar_audit, inchikey=product_inchikey)
    return (final_prod, new_history, record, budget_payload)


def _apply_site_rules(mol, cfg: EnumConfig, history: List[Dict[str, Any]],
                     depth: int, seen_global: set, aggregator=None, sugar_mask: Optional[set] = None, qa_bus: Optional[Dict] = None, budget_payload: Optional[dict] = None) -> List[Tuple]:
    """
    Apply site-based rules (R1, R2) with proper attempt boundary semantics.

    ATTEMPT BOUNDARY: Each site x halogen x k combination is exactly ONE attempt.
    - attempts=1 per site x halogen x k
    - products=1 if attempt yields >=1 valid products, 0 otherwise
    - no_product_matches=1 if attempt yields 0 valid products, 0 otherwise

    Uses symmetry folding to reduce representative sites while maintaining correct attempt counting.
    """
    results = []

    # Reset ring label cache at entry to avoid cross-test state leakage when
    # site rules are invoked directly by tests.
    try:
        from .sites import reset_ring_cache
        reset_ring_cache()
    except Exception:
        pass

    # Default to empty mask if not provided
    if sugar_mask is None:
        sugar_mask = set()

    # Default to empty qa_bus if not provided
    if qa_bus is None:
        qa_bus = {"sugar_mask_filtered": 0, "sugar_proximity_filtered": 0, "post_guard_blocked": 0}

    # Precompute expanded mask once for proximity guard (reused across R1/R2a/R2b)
    radius = int(cfg.sugar_cfg.get('proximity_guard_radius', 0) or 0)
    expanded_mask = _expand_mask_by_radius(mol, sugar_mask, radius) if (radius > 0 and sugar_mask) else set()

    # Prepare symmetry classes and ring labels
    ensure_ready(mol)
    # Use masked subgraph symmetry if enabled, otherwise standard symmetry
    if cfg.symmetry_cfg.get('compute_on_masked_subgraph', True) and sugar_mask:
        ranks = canonical_ranks_on_masked_subgraph(mol, sugar_mask)
    else:
        ranks = canonical_rank_atoms_safe(mol, breakTies=False)
    ring_label_map = _build_ring_label_map(mol)
    
    # Process rules with proper attempt boundaries (site x halogen x k)
    site_rules = [r for r in cfg.rules if r in ('R1', 'R2')]
    for rule_id in site_rules:
        # Initialize carbonyl_unknown_count for each rule
        carbonyl_unknown_count = 0
        
        # Get sites for this rule
        if rule_id == 'R1':
            sites = _sites_with_mask_drop(
                aromatic_CH_indices, mol=mol, sugar_mask=sugar_mask, qa_bus=qa_bus
            )
            LOG.debug(f"R1 sites after masking: {len(sites)}")

            # Apply proximity guard to filter sites near sugar mask
            sites = _apply_proximity_guard(sites, expanded_mask, qa_bus, 'R1')

            if cfg.pruning_cfg.get('enable_symmetry_fold', True):
                representatives = _pick_one_per_sym_group(sites, ranks, 'R1')
            else:
                representatives = sites
        else:  # R2 (legacy - only when R2a/R2b are not enabled)
            # Check if PR2 extensions (R2a/R2b) are configured
            r2_config = cfg.rules_cfg.get('R2', {})
            has_pr2_config = ('sp2_CH_in_C_ring' in r2_config or 'sp3_CH2_flavanone' in r2_config)

            # Skip legacy R2 if PR2 extensions are configured (even if disabled)
            if has_pr2_config:
                LOG.debug("PR2 extensions (R2a/R2b) configured, skipping legacy R2 rule")
                continue

            # Process all C ring indices and collect carbonyl_unknown events
            all_c_ring_sites = c_ring_indices(mol)
            # Filter out masked atoms
            filtered_c_ring_sites = filter_sites_with_mask(all_c_ring_sites, sugar_mask)
            LOG.debug(f"R2 sites: {len(all_c_ring_sites)} total, {len(filtered_c_ring_sites)} after masking")
            
            # Count carbonyl_unknown events for sites that couldn't be determined
            for site_idx in filtered_c_ring_sites:
                try:
                    atom = mol.GetAtomWithIdx(site_idx)
                    carbonyl_result = is_carbonyl_carbon(atom)

                    if carbonyl_result is None:
                        # Count carbonyl_unknown events
                        carbonyl_unknown_count += 1
                except Exception:
                    # If we can't check carbonyl status, count as unknown
                    carbonyl_unknown_count += 1

            # Filter for ready sites (already filtered for masking)
            sites = [i for i in filtered_c_ring_sites if is_c_ring_site_ready(mol, i)]

            if cfg.pruning_cfg.get('enable_symmetry_fold', True):
                representatives = _pick_one_per_group_with_ring(sites, ranks, ring_label_map, 'R2')
            else:
                representatives = sites
        
        # Each site per halogen is one attempt (correct attempt boundaries)
        first_r2_attempt = True
        for site in representatives:
            for halogen in cfg.halogens:
                produced_any = False
                qa_events = {}
                
                # For R2 rule, add carbonyl_unknown events to the first attempt only
                if rule_id == 'R2' and first_r2_attempt and carbonyl_unknown_count > 0:
                    qa_events['carbonyl_unknown'] = carbonyl_unknown_count
                    first_r2_attempt = False
                
                site_result = _apply_single_site(
                    mol, site, halogen, rule_id, ranks, ring_label_map,
                    history, depth, cfg, seen_global, budget_payload, qa_bus,
                    attempt=qa_events
                )
                if site_result:
                    results.append(site_result)
                    produced_any = True
                
                # Record this single site x halogen attempt
                if aggregator:
                    aggregator.record_attempt_result(rule_id, halogen, depth + 1,
                                                   produced_count=(1 if produced_any else 0),
                                                   qa_events_dict=qa_events,
                                                   k_ops=None, k_atoms=None)

    # PR2 Extensions: R2a and R2b (only when explicitly enabled)
    # ========================================================

    # T2-5: Initialize R2 halogen counters
    # Get R2 allowed halogens once (used by both R2a and R2b)
    # Intersect with main config halogens to respect user choice
    r2_allowed_halogens_raw = cfg.rules_cfg.get('R2', {}).get('allowed_halogens', cfg.halogens)
    r2_allowed_halogens = [h for h in r2_allowed_halogens_raw if h in cfg.halogens]

    by_rule_halogen_counts = {
        "R2a": {halogen: 0 for halogen in r2_allowed_halogens},
        "R2b": {halogen: 0 for halogen in r2_allowed_halogens}
    }

    # R2a: sp2 CH sites in C-ring (controlled by rules_cfg.R2.sp2_CH_in_C_ring)
    r2a_sites = []
    if cfg.rules_cfg.get('R2', {}).get('sp2_CH_in_C_ring', True):
        r2a_sites = _sites_with_mask_drop(
            c_ring_sp2_CH_sites, mol=mol, sugar_mask=sugar_mask, qa_bus=qa_bus
        )

        # Apply proximity guard to filter sites near sugar mask
        r2a_sites = _apply_proximity_guard(r2a_sites, expanded_mask, qa_bus, 'R2a')

        if cfg.pruning_cfg.get('enable_symmetry_fold', True):
            r2a_representatives = _pick_one_per_group_with_ring(r2a_sites, ranks, ring_label_map, 'R2a')
        else:
            r2a_representatives = r2a_sites

        for site in r2a_representatives:
            for halogen in r2_allowed_halogens:
                attempt_qa = {}
                site_result = _apply_single_site(
                    mol, site, halogen, 'R2a', ranks, ring_label_map,
                    history, depth, cfg, seen_global, budget_payload, qa_bus,
                    attempt=attempt_qa
                )
                if site_result:
                    results.append(site_result)
                    # T2-5: Increment halogen counter for successful R2a products
                    by_rule_halogen_counts["R2a"][halogen] += 1

                # Record attempt for R2a
                if aggregator:
                    aggregator.record_attempt_result('R2a', halogen, depth + 1,
                                                   produced_count=(1 if site_result else 0),
                                                   qa_events_dict=attempt_qa,
                                                   k_ops=None, k_atoms=None)

    # R2b: sp3 CH2 sites in flavanone C-ring (controlled by rules_cfg.R2.sp3_CH2_flavanone)
    r2b_sites = []
    if cfg.rules_cfg.get('R2', {}).get('sp3_CH2_flavanone', True):
        r2b_sites = _sites_with_mask_drop(
            c_ring_sp3_CH2_flavanone_sites, mol=mol, sugar_mask=sugar_mask, qa_bus=qa_bus, sugar_cfg=cfg.sugar_cfg
        )

        # Apply proximity guard to filter sites near sugar mask
        r2b_sites = _apply_proximity_guard(r2b_sites, expanded_mask, qa_bus, 'R2b')

        if cfg.pruning_cfg.get('enable_symmetry_fold', True):
            # R2b symmetry with 'CH2' tag for proper grouping: ('R2b', sym_class, ring_tag, 'CH2')
            r2b_representatives = _pick_one_per_group_with_ring(r2b_sites, ranks, ring_label_map, 'R2b', extra_tag='CH2')
        else:
            r2b_representatives = r2b_sites

        for site in r2b_representatives:
            for halogen in r2_allowed_halogens:
                attempt_qa = {}
                site_result = _apply_single_site(
                    mol, site, halogen, 'R2b', ranks, ring_label_map,
                    history, depth, cfg, seen_global, budget_payload, qa_bus,
                    attempt=attempt_qa
                )
                if site_result:
                    results.append(site_result)
                    # T2-5: Increment halogen counter for successful R2b products
                    by_rule_halogen_counts["R2b"][halogen] += 1

                # Record attempt for R2b
                if aggregator:
                    aggregator.record_attempt_result('R2b', halogen, depth + 1,
                                                   produced_count=(1 if site_result else 0),
                                                   qa_events_dict=attempt_qa,
                                                   k_ops=None, k_atoms=None)

    # R6: Methyl halogenation rules (controlled by rules_cfg.R6_methyl.enable)
    # =====================================================================
    #
    # P0 ARCHITECTURAL NOTE: R6 implementation uses temporary budget copies
    # to prevent cross-branch contamination during enumeration. Each branch
    # gets its own budget_tmp = budget.copy() to ensure failed paths don't
    # corrupt the budget state for successful paths.
    #
    r6_config = cfg.rules_cfg.get('R6_methyl', {})
    if r6_config.get('enable', False):
        LOG.debug("Processing R6 methyl halogenation rules")

        # Get methyl sites
        allow_on_methoxy = r6_config.get('allow_on_methoxy', False)
        r6_sites = enumerate_methyl_sites(mol, sugar_mask or set(), allow_on_methoxy)
        LOG.debug(f"R6 sites found: {len(r6_sites)}")

        # Apply proximity guard to filter sites near sugar mask
        r6_sites = _apply_proximity_guard(r6_sites, expanded_mask, qa_bus, 'R6')

        # Initialize budget state from payload or create new one
        if budget_payload is None:
            # Create initial budget for depth 0
            budget_mode = cfg.engine_cfg.get('budget_mode', 'ops')
            budget = BudgetState(budget_mode, cfg.k_max)
        else:
            # Restore budget state from payload
            from .budget import BudgetState
            budget = BudgetState.from_payload(budget_payload, cfg.k_max)

        # Track local state dedup per site
        seen_local_states = {}  # cidx -> set of state_keys

        # Get allowed halogens for R6
        r6_allowed = r6_config.get('allowed', ['F', 'Cl'])

        # Process each methyl site
        for cidx in r6_sites:
            seen_local_states[cidx] = set()

            # Try step halogenation for each allowed halogen
            for X in r6_allowed:
                if not validate_methyl_halogen(X):
                    continue

                # Create per-attempt QA events container
                attempt_qa = {}

                # Calculate costs
                op_cost = budget.token_cost(cidx)
                atom_cost = 1

                # Create per-branch budget copy
                budget_tmp = budget.copy()
                if not budget_tmp.can_afford(op_cost, atom_cost):
                    continue

                # Check local state dedup
                current_state = budget.get_state_key(cidx)
                # Simulate next state after adding X
                next_state_dict = dict(current_state)
                next_state_dict[X] = next_state_dict.get(X, 0) + 1
                next_state_key = tuple(sorted(next_state_dict.items()))

                if next_state_key in seen_local_states[cidx]:
                    continue  # Skip duplicate local state

                seen_local_states[cidx].add(next_state_key)

                # Apply step halogenation
                new_mol = apply_methyl_step(mol, cidx, X)
                if new_mol is None:
                    # Record failed attempt
                    if aggregator:
                        aggregator.record_attempt_result('R6', X, depth + 1,
                                                       produced_count=0,
                                                       qa_events_dict=attempt_qa,
                                                       k_ops=getattr(budget_tmp, 'k_ops', None),
                                                       k_atoms=getattr(budget_tmp, 'k_atoms', None))
                    continue

                # Configurable deduplication stage and policy
                dedup_stage = cfg.engine_cfg.get('dedup_stage', 'pre')
                policy = cfg.engine_cfg.get('dedup_policy', 'auto')
                if policy == 'auto':
                    policy = 'stable_key' if cfg.pruning_cfg.get('enable_symmetry_fold', True) else 'state_sig'

                # Always compute sanitized key for the record, independent of dedup policy
                product_inchikey = None

                # Prepare state signature only if needed
                state_sig = None
                if policy == 'state_sig':
                    # For R6, create minimal history for state signature computation
                    temp_history = history + [{
                        'rule': 'R6',
                        'site': cidx,
                        'halogen': X,
                        'type': 'step'
                    }]
                    state_sig = compute_state_signature(temp_history)

                dedup_key = None

                # Pre-constraint deduplication check
                if dedup_stage == 'pre':
                    dedup_key, is_dup = early_check(
                        new_mol, seen_global, qa_bus,
                        metric='dedup_hits_inchi',
                        policy=policy,
                        state_sig=state_sig,
                        attempt=attempt_qa
                    )
                    if is_dup:
                        # When using state_sig, also bump statesig counter for clarity
                        if policy == 'state_sig':
                            qa_mark(qa_bus, attempt_qa, 'dedup_hits_statesig')
                        # Record failed attempt due to dedup
                        if aggregator:
                            aggregator.record_attempt_result('R6', X, depth + 1,
                                                           produced_count=0,
                                                           qa_events_dict=attempt_qa,
                                                           k_ops=getattr(budget_tmp, 'k_ops', None),
                                                           k_atoms=getattr(budget_tmp, 'k_atoms', None))
                        continue

                # Commit on tmp budget only
                budget_tmp.charge(op_cost, atom_cost)
                budget_tmp.mark_first_touch(cidx)
                budget_tmp.bump_state(cidx, X)

                # Create product history record
                new_history = history + [{
                    'rule': 'R6',
                    'site': cidx,
                    'halogen': X,
                    'type': 'step',
                    'k_ops': budget_tmp.k_ops,
                    'k_atoms': budget_tmp.k_atoms,
                    'budget_mode': budget_tmp.budget_mode
                }]

                # Apply constraints and create product record
                if not accept_constraints(new_mol, new_history, cfg.constraints):
                    # Record failed attempt
                    if aggregator:
                        aggregator.record_attempt_result('R6', X, depth + 1,
                                                       produced_count=0,
                                                       qa_events_dict=attempt_qa,
                                                       k_ops=getattr(budget_tmp, 'k_ops', None),
                                                       k_atoms=getattr(budget_tmp, 'k_atoms', None))
                    continue

                # Post-constraint deduplication check
                if dedup_stage == 'post':
                    if dedup_key is None:
                        dedup_key, is_dup = early_check(
                            new_mol, seen_global, qa_bus,
                            metric='dedup_hits_inchi',
                            policy=policy,
                            state_sig=state_sig,
                            attempt=attempt_qa
                        )
                        if is_dup:
                            if policy == 'state_sig':
                                qa_mark(qa_bus, attempt_qa, 'dedup_hits_statesig')
                            # Record failed attempt
                            if aggregator:
                                aggregator.record_attempt_result('R6', X, depth + 1,
                                                               produced_count=0,
                                                               qa_events_dict=attempt_qa,
                                                               k_ops=getattr(budget_tmp, 'k_ops', None),
                                                               k_atoms=getattr(budget_tmp, 'k_atoms', None))
                            continue

                # Apply sugar mask post-guard
                if post_guard_blocked(new_mol, cfg.sugar_cfg):
                    qa_mark(qa_bus, attempt_qa, 'sugar_post_guard_blocked')
                    # Record failed attempt
                    if aggregator:
                        aggregator.record_attempt_result('R6', X, depth + 1,
                                                       produced_count=0,
                                                       qa_events_dict=attempt_qa,
                                                       k_ops=getattr(budget_tmp, 'k_ops', None),
                                                       k_atoms=getattr(budget_tmp, 'k_atoms', None))
                    continue

                # Commit deduplication key after all validations pass
                commit(dedup_key, seen_global)

                # Compute record inchikey (stable) if still None
                if product_inchikey is None:
                    product_inchikey = compute_stable_key(new_mol) or "UNKNOWN"

                # Create product record
                record = {
                    'smiles': Chem.MolToSmiles(new_mol),
                    'inchikey': product_inchikey,
                    'rule': 'R6',
                    'rule_family': 'R6',  # Rule family for grouping
                    'halogen': X,
                    'k': budget_tmp.k_atoms,
                    'k_ops': budget_tmp.k_ops,
                    'k_atoms': budget_tmp.k_atoms,
                    'budget_mode': budget_tmp.budget_mode,
                    'parent_smi': Chem.MolToSmiles(mol),
                    'site_tokens_json': json.dumps(getattr(budget_tmp, 'site_tokens', {}) or {}, separators=(',', ':'))
                }

                # Include updated budget payload for BFS continuation
                new_budget_payload = budget_tmp.to_payload()
                results.append((new_mol, new_history, record, new_budget_payload))

                # Record successful attempt
                if aggregator:
                    aggregator.record_attempt_result('R6', X, depth + 1,
                                                   produced_count=1,
                                                   qa_events_dict=attempt_qa,
                                                   k_ops=getattr(budget_tmp, 'k_ops', None),
                                                   k_atoms=getattr(budget_tmp, 'k_atoms', None))

        # Process macro halogenation if enabled
        macro_cfg = r6_config.get('macro', {})
        if macro_cfg.get('enable', False):
            macro_labels = macro_cfg.get('labels', ['CF3', 'CCl3'])

            for cidx in r6_sites:
                for label in macro_labels:
                    X = 'F' if label == 'CF3' else 'Cl'

                    if not validate_methyl_halogen(X):
                        continue

                    # Create per-attempt QA events container
                    attempt_qa = {}

                    # Calculate costs for macro
                    op_cost = budget.token_cost(cidx)
                    atom_cost = 3  # CF3/CCl3 adds 3 atoms

                    # Create per-branch budget copy
                    budget_tmp = budget.copy()
                    if not budget_tmp.can_afford(op_cost, atom_cost):
                        continue

                    # Apply macro halogenation
                    new_mol = apply_methyl_macro(mol, cidx, label)
                    if new_mol is None:
                        # Record failed attempt
                        if aggregator:
                            aggregator.record_attempt_result('R6', X, depth + 1,
                                                           produced_count=0,
                                                           qa_events_dict=attempt_qa,
                                                           k_ops=getattr(budget_tmp, 'k_ops', None),
                                                           k_atoms=getattr(budget_tmp, 'k_atoms', None))
                        continue

                    # Configurable deduplication stage and policy
                    dedup_stage = cfg.engine_cfg.get('dedup_stage', 'pre')
                    policy = cfg.engine_cfg.get('dedup_policy', 'auto')
                    if policy == 'auto':
                        policy = 'stable_key' if cfg.pruning_cfg.get('enable_symmetry_fold', True) else 'state_sig'

                    # Always compute sanitized key for the record, independent of dedup policy
                    product_inchikey = None

                    # Prepare state signature only if needed
                    state_sig = None
                    if policy == 'state_sig':
                        # For R6 macro, create minimal history for state signature computation
                        temp_history = history + [{
                            'rule': 'R6',
                            'site': cidx,
                            'halogen': X,
                            'type': 'macro',
                            'label': label
                        }]
                        state_sig = compute_state_signature(temp_history)

                    dedup_key = None

                    # Pre-constraint deduplication check
                    if dedup_stage == 'pre':
                        dedup_key, is_dup = early_check(
                            new_mol, seen_global, qa_bus,
                            metric='dedup_hits_inchi',
                            policy=policy,
                            state_sig=state_sig,
                            attempt=attempt_qa
                        )
                        if is_dup:
                            # When using state_sig, also bump statesig counter for clarity
                            if policy == 'state_sig':
                                qa_mark(qa_bus, attempt_qa, 'dedup_hits_statesig')
                            # Record failed attempt due to dedup
                            if aggregator:
                                aggregator.record_attempt_result('R6', X, depth + 1,
                                                               produced_count=0,
                                                               qa_events_dict=attempt_qa,
                                                               k_ops=getattr(budget_tmp, 'k_ops', None),
                                                               k_atoms=getattr(budget_tmp, 'k_atoms', None))
                            continue

                    # Commit on tmp budget only
                    budget_tmp.charge(op_cost, atom_cost)
                    budget_tmp.mark_first_touch(cidx)
                    for _ in range(3):
                        budget_tmp.bump_state(cidx, X)

                    # Create product history record
                    new_history = history + [{
                        'rule': 'R6',
                        'site': cidx,
                        'halogen': X,
                        'type': 'macro',
                        'label': label,
                        'k_ops': budget_tmp.k_ops,
                        'k_atoms': budget_tmp.k_atoms,
                        'budget_mode': budget_tmp.budget_mode
                    }]

                    # Apply constraints and create product record
                    if not accept_constraints(new_mol, new_history, cfg.constraints):
                        # Record failed attempt
                        if aggregator:
                            aggregator.record_attempt_result('R6', X, depth + 1,
                                                           produced_count=0,
                                                           qa_events_dict=attempt_qa,
                                                           k_ops=getattr(budget_tmp, 'k_ops', None),
                                                           k_atoms=getattr(budget_tmp, 'k_atoms', None))
                        continue

                    # Post-constraint deduplication check
                    if dedup_stage == 'post':
                        if dedup_key is None:
                            dedup_key, is_dup = early_check(
                                new_mol, seen_global, qa_bus,
                                metric='dedup_hits_inchi',
                                policy=policy,
                                state_sig=state_sig,
                                attempt=attempt_qa
                            )
                            if is_dup:
                                if policy == 'state_sig':
                                    qa_mark(qa_bus, attempt_qa, 'dedup_hits_statesig')
                                # Record failed attempt
                                if aggregator:
                                    aggregator.record_attempt_result('R6', X, depth + 1,
                                                                   produced_count=0,
                                                                   qa_events_dict=attempt_qa,
                                                                   k_ops=getattr(budget_tmp, 'k_ops', None),
                                                                   k_atoms=getattr(budget_tmp, 'k_atoms', None))
                                continue

                    # Apply sugar mask post-guard
                    if post_guard_blocked(new_mol, cfg.sugar_cfg):
                        qa_mark(qa_bus, attempt_qa, 'sugar_post_guard_blocked')
                        # Record failed attempt
                        if aggregator:
                            aggregator.record_attempt_result('R6', X, depth + 1,
                                                           produced_count=0,
                                                           qa_events_dict=attempt_qa,
                                                           k_ops=getattr(budget_tmp, 'k_ops', None),
                                                           k_atoms=getattr(budget_tmp, 'k_atoms', None))
                        continue

                    # Commit deduplication key after all validations pass
                    commit(dedup_key, seen_global)

                    # Compute record inchikey (stable) if still None
                    if product_inchikey is None:
                        product_inchikey = compute_stable_key(new_mol) or "UNKNOWN"

                    # Create product record
                    record = {
                        'smiles': Chem.MolToSmiles(new_mol),
                        'inchikey': product_inchikey,
                        'rule': 'R6',
                        'rule_family': 'R6',  # Rule family for grouping
                        'halogen': X,
                        'macro_label': label,
                        'k': budget_tmp.k_atoms,
                        'k_ops': budget_tmp.k_ops,
                        'k_atoms': budget_tmp.k_atoms,
                        'budget_mode': budget_tmp.budget_mode,
                        'parent_smi': Chem.MolToSmiles(mol),
                        'site_tokens_json': json.dumps(getattr(budget_tmp, 'site_tokens', {}) or {}, separators=(',', ':'))
                    }

                    # Include updated budget payload for BFS continuation
                    new_budget_payload = budget_tmp.to_payload()
                    results.append((new_mol, new_history, record, new_budget_payload))

                    # Record successful attempt
                    if aggregator:
                        aggregator.record_attempt_result('R6', X, depth + 1,
                                                       produced_count=1,
                                                       qa_events_dict=attempt_qa,
                                                       k_ops=getattr(budget_tmp, 'k_ops', None),
                                                       k_atoms=getattr(budget_tmp, 'k_atoms', None))

    return results


def _apply_single_site(mol, site: int, halogen: str, rule: str,
                      ranks: List[int], ring_label_map: Dict[int, Any],
                      history: List[Dict[str, Any]], depth: int,
                      cfg: EnumConfig, seen_global: set, budget_payload: Optional[dict] = None,
                      qa_bus: Optional[Dict[str, int]] = None,
                      attempt: Optional[Dict[str, int]] = None) -> Optional[Tuple]:
    """Apply halogenation to a single site."""

    # Apply halogenation
    product_mol = apply_single_site_halogenation(mol, site, halogen)
    if product_mol is None:
        return None
    
    # Build history item
    history_item = {
        'rule': rule,
        'site': int(site),
        'sym': int(ranks[site]),
        'ring_tag': ring_label_map.get(site),
        'halogen': halogen,
        'depth': depth + 1
    }
    new_history = history + [history_item]

    # Configurable deduplication stage and policy
    dedup_stage = cfg.engine_cfg.get('dedup_stage', 'pre')
    policy = cfg.engine_cfg.get('dedup_policy', 'auto')
    if policy == 'auto':
        policy = 'stable_key' if cfg.pruning_cfg.get('enable_symmetry_fold', True) else 'state_sig'

    # Always compute sanitized key for the record, independent of dedup policy
    product_inchikey = None

    # Prepare state signature only if needed
    state_sig = None
    if policy == 'state_sig':
        state_sig = compute_state_signature(new_history)

    dedup_key = None

    # Pre-constraint deduplication check
    if dedup_stage == 'pre':
        dedup_key, is_dup = early_check(
            product_mol, seen_global, qa_bus,
            metric='dedup_hits_inchi',
            policy=policy,
            state_sig=state_sig,
            attempt=attempt
        )
        if is_dup:
            # When using state_sig, also bump statesig counter for clarity
            if policy == 'state_sig':
                qa_mark(qa_bus, attempt, 'dedup_hits_statesig')
            return None

    # Check constraints
    ok, violations = accept_constraints(product_mol, new_history, cfg.constraints)
    if not ok:
        return None

    # Post-constraint deduplication check
    if dedup_stage == 'post':
        if dedup_key is None:
            dedup_key, is_dup = early_check(
                product_mol, seen_global, qa_bus,
                metric='dedup_hits_inchi',
                policy=policy,
                state_sig=state_sig,
                attempt=attempt
            )
            if is_dup:
                if policy == 'state_sig':
                    qa_mark(qa_bus, attempt, 'dedup_hits_statesig')
                return None

    # Commit deduplication key after all validations pass
    commit(dedup_key, seen_global)

    # Compute record inchikey (stable) if still None
    if product_inchikey is None:
        product_inchikey = compute_stable_key(product_mol) or "UNKNOWN"

    # Create product record with sugar audit (if enabled)
    product_mask, _ = get_sugar_mask_with_status(product_mol, mode=cfg.sugar_cfg.get('mode', 'heuristic'), sugar_cfg=cfg.sugar_cfg)
    if cfg.sugar_cfg.get('audit', False):
        sugar_audit = compute_sugar_audit_fields(product_mol, product_mask, cfg.sugar_cfg)
    else:
        sugar_audit = None
    record = _make_record(product_mol, mol, depth + 1, rule, halogen, new_history, True, sugar_audit, inchikey=product_inchikey)
    return (product_mol, new_history, record, budget_payload)


def _replace_star_and_sanitize(mol, halogen: str):
    """Replace * atom with halogen and sanitize."""
    try:
        from .chem_compat import Chem
        editable = Chem.EditableMol(mol)
        halogen_atomic_num = {'F': 9, 'Cl': 17, 'Br': 35, 'I': 53}[halogen]
        
        # Find and replace star atoms
        for atom_idx in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 0:  # Star atom
                editable.ReplaceAtom(atom_idx, Chem.Atom(halogen_atomic_num))
        
        result_mol = editable.GetMol()
        
        # Sanitize
        try:
            Chem.SanitizeMol(result_mol)
        except Exception:
            pass  # Continue even if sanitization fails
        ensure_ready(result_mol)
        return result_mol
        
    except Exception:
        return None


def _build_ring_label_map(mol) -> Dict[int, Any]:
    """Build atom index to flavonoid ring label mapping (A/B/C)."""
    ring_label_map = {}
    try:
        # Use flavonoid-specific ring labeling for stable per-ring quotas
        for atom_idx in range(mol.GetNumAtoms()):
            label = flavonoid_ring_label(mol, atom_idx)
            if label:  # Only store non-empty labels
                ring_label_map[atom_idx] = label
                    
    except Exception:
        pass
    
    return ring_label_map


def _pick_one_per_sym_group(sites: List[int], ranks: List[int], rule: str) -> List[int]:
    """Pick one representative per symmetry group."""
    if not sites:
        return []
    
    # Group by symmetry class
    groups = {}
    for site in sites:
        sym_class = ranks[site]
        if sym_class not in groups:
            groups[sym_class] = []
        groups[sym_class].append(site)
    
    # Pick first representative from each group
    return [sites[0] for sites in groups.values()]


def _pick_one_per_group_with_ring(sites: List[int], ranks: List[int],
                                 ring_label_map: Dict[int, Any], rule: str, extra_tag: str = None) -> List[int]:
    """Pick one representative per (symmetry, ring_tag[, extra_tag]) group."""
    if not sites:
        return []

    # Group by (sym_class, ring_tag) or (rule, sym_class, ring_tag, extra_tag) for enhanced grouping
    groups = {}
    for site in sites:
        sym_class = ranks[site]
        ring_tag = ring_label_map.get(site)

        if extra_tag:
            # Enhanced grouping key for rules like R2b: (rule, sym_class, ring_tag, extra_tag)
            key = (rule, sym_class, ring_tag, extra_tag)
        else:
            # Standard grouping key: (sym_class, ring_tag)
            key = (sym_class, ring_tag)

        if key not in groups:
            groups[key] = []
        groups[key].append(site)

    # Pick first representative from each group
    return [sites[0] for sites in groups.values()]


def _make_record(product_mol, parent_mol, depth: int, rule: str,
                halogen: str, history: List[Dict[str, Any]], constraints_ok: bool,
                sugar_audit: Optional[Dict[str, Any]] = None,
                inchikey: Optional[str] = None) -> Dict[str, Any]:
    """Create product record dictionary."""
    # Basic identifiers (use chem_compat safe helpers)
    try:
        smiles_prod = Chem.MolToSmiles(product_mol, canonical=True)
    except Exception:
        smiles_prod = ''
    try:
        smiles_parent = Chem.MolToSmiles(parent_mol, canonical=True)
    except Exception:
        smiles_parent = ''

    # Use provided inchikey if available, otherwise compute sanitized version as fallback
    if inchikey is None:
        inchikey = to_inchikey_sanitized(product_mol) or to_inchikey(product_mol)

    # Use sanitized parent inchikey for consistency
    parent_key = to_inchikey_sanitized(parent_mol) or to_inchikey(parent_mol)

    # Determine rule family for grouping (R2a/R2b both belong to R2 family)
    rule_family = 'R2' if rule in ('R2a', 'R2b') else rule

    record = {
        'inchikey': inchikey or 'UNKNOWN',
        'smiles': smiles_prod,
        'parent_inchikey': parent_key or 'UNKNOWN',
        'parent_smiles': smiles_parent,
        'k': depth,
        'rule': rule,  # Last rule applied (for compatibility)
        'rule_family': rule_family,  # Rule family for grouping
        'halogen': halogen,  # Last halogen applied (for compatibility)
        'substitutions': history[:],  # Complete history
        'constraints_ok': constraints_ok,
        'constraints_violations': {} if constraints_ok else {}
    }
    
    # Add descriptors
    try:
        descriptors = basic_descriptors(product_mol)
        record.update(descriptors)
    except Exception:
        pass
    
    # Add QC flags
    try:
        record['sanitize_ok'] = sanitize_ok(product_mol)
        record['pains_flags'] = pains_flags(product_mol)
    except Exception:
        record['sanitize_ok'] = False
        record['pains_flags'] = []

    # Add sugar masking audit fields
    if sugar_audit:
        record.update(sugar_audit)
    else:
        # Default empty sugar audit fields
        record.update({
            'sugar_mask_atoms': [],
            'sugar_rings': [],
            'masked_oh_count': 0,
            'masked_bridge_o_count': 0
        })

    return record


def _iter_reaction_mols(reaction_products):
    """Yield product molecules from possibly nested reaction outputs (lists/tuples)."""
    if not reaction_products:
        return
    stack = [reaction_products]
    while stack:
        x = stack.pop()
        if x is None:
            continue
        if isinstance(x, (list, tuple)):
            for item in reversed(list(x)):
                stack.append(item)
        else:
            yield x
