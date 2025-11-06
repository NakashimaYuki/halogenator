# -*- coding: ascii -*-
"""k=1 enumeration core logic."""

import logging
from typing import List, Dict, Any, Optional, Tuple

from .rules import build_reactions
from .sites import (
    aromatic_CH_indices, c_ring_indices, symmetry_groups,
    get_ring_tag_for_atom, ensure_ready, is_c_ring_site_ready,
    flavonoid_ring_label, is_carbonyl_carbon, c_ring_sp2_CH_sites, c_ring_sp3_CH2_flavanone_sites
)
from .reactions import apply_single_site_halogenation
from .qc import sanitize_ok, basic_descriptors, pains_flags
from .dedupe import dedupe_with_props
from .standardize import to_inchikey
from .enumerate_k import _run_reaction_safely, _validate_totals_pivots_consistency, QAAggregator, _compute_totals_from_aggregator, emit_product
from .sites_methyl import enumerate_methyl_sites
from .rules_methyl import apply_methyl_step, validate_methyl_halogen, apply_methyl_macro
from .guard import rdkit_guard
from .sugar_mask import get_sugar_mask_with_full_status
from .schema import empty_qa_paths

# Logger for k=1 enumeration
LOG = logging.getLogger(__name__)


def enumerate_k1_with_stats(parent_smi: str, cfg, stream_shape: str = 'legacy') -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Enumerate k=1 halogenated products with QA statistics.
    Returns products and QA stats in the same schema as enumerate_with_stats.
    Interface matches enumerate_with_stats for CLI compatibility.
    
    Args:
        parent_smi: Parent molecule SMILES string
        cfg: EnumConfig object with halogens, constraints, etc.
    
    Returns:
        (products_list, qa_stats_dict) where qa_stats_dict matches enumerate_with_stats schema
    """
    try:
        from .chem_compat import Chem
        # Test if RDKit is actually available
        if not hasattr(Chem, 'MolFromSmiles') or not callable(Chem.MolFromSmiles):
            raise ImportError("RDKit not fully available")
    except (ImportError, Exception):
        LOG.debug("RDKit not available")
        return [], {
            'no_product_matches': 1,
            'template_unsupported': 0,
            'dedup_hits_statesig': 0,
            'dedup_hits_inchi': 0,
            'qa_paths': empty_qa_paths()
        }
    
    # Initialize QA stats early for guard usage
    qa_stats = {
        'no_product_matches': 0,
        'template_unsupported': 0,
        'dedup_hits_statesig': 0,
        'dedup_hits_inchi': 0,
        'qa_paths': empty_qa_paths()
    }
    
    # Convert SMILES to molecule with RDKit guard
    parent_mol = None
    with rdkit_guard(qa_stats):
        parent_mol = Chem.MolFromSmiles(parent_smi)
    
    if parent_mol is None:
        # Return empty results for invalid SMILES with proper version based on stream_shape
        empty_qa_stats = {
            'no_product_matches': 1,
            'template_unsupported': qa_stats.get('template_unsupported', 0),
            'dedup_hits_statesig': 0,
            'dedup_hits_inchi': 0,
            'qa_paths': qa_stats.get('qa_paths', empty_qa_paths())
        }
        
        # Set version based on stream_shape
        if stream_shape == 'v2':
            empty_qa_stats['version'] = '2'
            empty_qa_stats['pivots'] = {}  # Empty pivots for v2
        else:
            empty_qa_stats['version'] = '1'
        
        return [], empty_qa_stats
    
    # Extract parameters from EnumConfig
    halogens = list(cfg.halogens)
    rules = list(cfg.rules)  # Use configured rules, don't force R2
    config = {
        'constraints': cfg.constraints,
        'qc': cfg.qc_cfg,
        'standardize': cfg.std_cfg,
        'rules_cfg': cfg.rules_cfg
    }
    
    # Initialize pivot aggregator for M2 statistics (version 2 format)
    aggregator = QAAggregator()
    
    # Get products using modified function with stats tracking
    product_tuples = _enumerate_k1_halogenation_with_stats_tracking(parent_mol, halogens, rules, config, qa_stats, aggregator)
    
    # Convert to records format (list of dicts) and track stats
    products = []
    for product_mol, props in product_tuples:
        try:
            # Build history for k=1 product (single substitution step)
            history_entry = {
                'rule': props['rule'],
                'site': props.get('site_index'),
                'sym': props.get('sym_class'),
                'ring_tag': props.get('ring_tag', ''),
                'halogen': props['halogen'],
                'depth': 1
            }

            # Add macro substitution metadata to history if present
            if props.get('substitution_type') == 'macro':
                history_entry['type'] = 'macro'
                history_entry['label'] = props.get('macro_label')

            # Add detection metadata to history if present (for R2b fallback tracking)
            if props.get('detection'):
                history_entry['detection'] = props.get('detection')

            history = [history_entry]

            # Prepare extra fields (include macro_label if present)
            extra_fields = {
                'site_index': props.get('site_index'),
                'sym_class': props.get('sym_class'),
                'ring_tag': props.get('ring_tag')
            }

            # Add macro label to record if this is a macro substitution
            if props.get('macro_label'):
                extra_fields['macro_label'] = props['macro_label']

            # Add detection to extra_fields if present (for R2b fallback tracking)
            if props.get('detection'):
                extra_fields['detection'] = props.get('detection')

            # Add sub_rule for R2 series (R2a, R2b) to enable detailed rule tracking
            if props['rule'] in ('R2a', 'R2b'):
                extra_fields['sub_rule'] = props['rule']

            # For macro substitutions (CF3/CCl3), override k metrics
            # Macro: k_ops=1 (single operation), k_atoms=3 (three atoms replaced)
            metrics_override = None
            if props.get('substitution_type') == 'macro':
                metrics_override = {'k_ops': 1, 'k_atoms': 3}

            # Use emit_product for unified record generation
            # For k=1: parent_mol is the root parent, so root_parent = parent
            product_record = emit_product(
                product_mol=product_mol,
                parent_mol=parent_mol,  # immediate parent
                rule=props['rule'],
                halogen=props['halogen'],
                history=history,
                depth=1,
                budget_state=None,  # k=1 has no budget state
                extra_fields=extra_fields,
                root_parent_mol=parent_mol,  # For k=1, root = immediate parent
                root_parent_inchikey=props['parent_inchikey'],
                root_parent_smiles=props['parent_smiles'],
                metrics_override=metrics_override  # Override for macro: k_ops=1, k_atoms=3
            )

            # Check if product generation succeeded (emit_product returns error indicators)
            if product_record.get('inchikey') == 'ERROR' or not product_record.get('smiles'):
                qa_stats['no_product_matches'] += 1
                continue

            products.append(product_record)

        except Exception as e:
            LOG.warning(f"Failed to convert product to record: {e}")
            qa_stats['no_product_matches'] += 1
            continue
    
    # Compute final stats with totals derived from aggregator (eliminates totals-pivots inconsistency)
    aggregator_totals = _compute_totals_from_aggregator(aggregator)
    
    # Merge guard events from qa_stats into aggregator_totals before final output
    merged_qa_paths = dict(aggregator_totals.get('qa_paths', {}))
    for key, value in qa_stats.get('qa_paths', {}).items():
        merged_qa_paths[key] = merged_qa_paths.get(key, 0) + value
    
    # Merge guard top-level counts
    template_unsupported = aggregator_totals.get('template_unsupported', 0) + qa_stats.get('template_unsupported', 0)
    no_product_matches = aggregator_totals.get('no_product_matches', 0) + qa_stats.get('no_product_matches', 0)
    
    # Build final stats dict based on stream_shape parameter
    if stream_shape == 'v2':
        # Version 2 format with pivots
        final_qa_stats = {
            'version': '2',
            'pivots': aggregator.to_pivots_dict()
        }
        
        # Include aggregator totals (only if they exist - no sentinel pollution)
        for key in ['attempts', 'products']:
            if key in aggregator_totals:
                final_qa_stats[key] = aggregator_totals[key]
        
        # Use merged counts for guard-sensitive metrics
        final_qa_stats['no_product_matches'] = no_product_matches
        final_qa_stats['template_unsupported'] = template_unsupported
        
        # Use merged qa_paths for v2 format
        final_qa_stats['qa_paths'] = merged_qa_paths
                
        # Legacy dedup stats (not tracked in aggregator) - only if they exist
        if qa_stats.get('dedup_hits_statesig', 0) > 0:
            final_qa_stats['dedup_hits_statesig'] = qa_stats['dedup_hits_statesig']
        if qa_stats.get('dedup_hits_inchi', 0) > 0:
            final_qa_stats['dedup_hits_inchi'] = qa_stats['dedup_hits_inchi']
            
        # Validate totals-pivots consistency
        _validate_totals_pivots_consistency(final_qa_stats, aggregator)
    else:
        # Legacy format (stream_shape == 'legacy')
        final_qa_stats = {
            'version': '1',
            'no_product_matches': no_product_matches,
            'template_unsupported': template_unsupported,
            'dedup_hits_statesig': qa_stats.get('dedup_hits_statesig', 0),
            'dedup_hits_inchi': qa_stats.get('dedup_hits_inchi', 0),
            'qa_paths': merged_qa_paths
        }
        
        # Basic consistency validation for legacy format
        attempts = aggregator_totals.get('attempts', 0)
        products_count = aggregator_totals.get('products', 0)
        if attempts < 0 or products_count < 0:
            LOG.warning(f"Legacy format consistency check: negative counts detected (attempts={attempts}, products={products_count})")
        if products_count > attempts:
            LOG.warning(f"Legacy format consistency check: products ({products_count}) > attempts ({attempts})")
    
    return products, final_qa_stats


def _enumerate_k1_halogenation_with_stats_tracking(
    parent_mol, 
    halogens: List[str], 
    rules: List[str], 
    config: Dict[str, Any], 
    stats_dict: Dict[str, Any], 
    aggregator=None
):
    """Internal function for k=1 enumeration with stats tracking aligned with k>1 semantics."""
    if parent_mol is None:
        return []
    
    try:
        from .chem_compat import Chem
        # SanitizeFlags are optional - try to import but don't fail if unavailable
        try:
            from rdkit.Chem import SanitizeFlags
        except ImportError:
            SanitizeFlags = None
    except Exception:
        LOG.debug("RDKit not available")
        return []
    
    # Ensure molecule is properly sanitized before any operations
    try:
        Chem.SanitizeMol(parent_mol)
    except Exception:
        parent_mol.UpdatePropertyCache(strict=False)
    
    parent_smiles = Chem.MolToSmiles(parent_mol, canonical=True)
    parent_inchikey = to_inchikey(parent_mol)

    # Compute sugar mask for R6 enumeration
    sugar_cfg = config.get('sugar_cfg', {})
    sugar_mode = sugar_cfg.get('mode', 'heuristic')
    sugar_mask, mask_degraded, status_metadata = get_sugar_mask_with_full_status(
        parent_mol,
        mode=sugar_mode,
        sugar_cfg=sugar_cfg
    )
    LOG.debug(f"k=1 sugar mask: {len(sugar_mask)} atoms masked, degraded={mask_degraded}, mode={sugar_mode}")

    # Emit QA event if sugar masking is disabled (raw mode)
    if aggregator is not None and sugar_mode == 'off':
        aggregator.record('sugar_mask_mode_off', amount=1)

    products = []
    reactions = build_reactions()
    
    # Step 1: Apply reaction-based rules (R1, R3, R4, R5) with proper attempt boundaries
    for rule in ['R1', 'R3', 'R4', 'R5']:
        if rule in rules and rule in reactions:
            for halogen in halogens:
                # Each (rule, halogen) combination is one attempt
                attempt_products_count = 0
                attempt_qa_events = {}
                
                if halogen not in reactions[rule]:
                    # Record failed attempt due to missing reaction template
                    if aggregator:
                        attempt_qa_events['template_unsupported'] = 1
                        aggregator.record_attempt_result(rule, halogen, 1, 0, attempt_qa_events,
                                                       k_ops=None, k_atoms=None)
                    stats_dict['template_unsupported'] += 1
                    continue
                
                rxn = reactions[rule][halogen]
                local_seen = set()  # Local deduplication per (rule, halogen)
                
                # Run reaction safely (without passing aggregator to avoid double counting)
                reaction_products = _run_reaction_safely(rxn, (parent_mol,), rule, halogen, stats_dict, None, current_k=1)
                
                for product_set in reaction_products:
                    for product_mol in product_set:
                        if product_mol is not None:
                            # Sanitize product molecule (lightweight first, full if needed)
                            try:
                                Chem.SanitizeMol(product_mol, sanitizeOps=SanitizeFlags.SANITIZE_PROPERTIES)
                            except Exception:
                                try:
                                    Chem.SanitizeMol(product_mol)
                                except Exception:
                                    product_mol.UpdatePropertyCache(strict=False)
                            
                            # Ensure ring info and properties are ready for downstream operations
                            try:
                                ensure_ready(product_mol)
                            except Exception:
                                pass  # Continue with potentially incomplete initialization
                            
                            try:
                                product_smiles = Chem.MolToSmiles(product_mol, canonical=True)
                            except Exception:
                                # fallback: skip local dedup on this product, but keep pipeline alive
                                product_smiles = None
                            
                            if product_smiles is None or product_smiles not in local_seen:
                                if product_smiles is not None:
                                    local_seen.add(product_smiles)
                                props = {
                                    'parent_smiles': parent_smiles,
                                    'parent_inchikey': parent_inchikey,
                                    'rule': rule,
                                    'site_index': None,
                                    'sym_class': None,
                                    'ring_tag': None,
                                    'halogen': halogen,
                                    'depth': 1
                                }
                                products.append((product_mol, props))
                                attempt_products_count += 1
                
                # Record this attempt result with proper semantics
                if aggregator:
                    aggregator.record_attempt_result(rule, halogen, 1, attempt_products_count, attempt_qa_events,
                                                   k_ops=None, k_atoms=None)
    
    # Step 2: Apply R2a/R2b rules (separate processing to match k>=2 path)
    if 'R2' in rules:
        # Get R2 configuration
        r2_config = config.get('rules_cfg', {}).get('R2', {})

        # [DIAGNOSTIC] Entry gating for R2
        r2_enabled = 'R2' in rules
        r2a_gate = r2_config.get('sp2_CH_in_C_ring', True)
        r2b_gate = r2_config.get('sp3_CH2_flavanone', True)
        LOG.info("[R2] Entry gates: R2=%s, R2a=%s, R2b=%s; rules=%s",
                 r2_enabled, r2a_gate, r2b_gate, rules)
        LOG.debug("Effective R2 configuration (k=1): %s", r2_config)

        # Get R2 allowed halogens (intersect with global halogens)
        r2_allowed_halogens_raw = r2_config.get('allowed_halogens', halogens)
        r2_allowed_halogens = [h for h in r2_allowed_halogens_raw if h in halogens]
        LOG.debug("[R2] Allowed halogens: %s (from config: %s, global: %s)",
                 r2_allowed_halogens, r2_allowed_halogens_raw, halogens)

        # R2a: sp2 CH sites in C-ring (controlled by rules_cfg.R2.sp2_CH_in_C_ring)
        if r2_config.get('sp2_CH_in_C_ring', True):
            try:
                r2a_sites = c_ring_sp2_CH_sites(parent_mol, set())
                LOG.info("[R2a] Site enumeration: found %d candidates", len(r2a_sites))
                if r2a_sites:
                    LOG.debug("[R2a] Sites: %s", r2a_sites[:10])  # Show first 10

                for site in r2a_sites:
                    for halogen in r2_allowed_halogens:
                        produced_any = False
                        qa_events = {}

                        try:
                            product_mol = apply_single_site_halogenation(parent_mol, site, halogen)
                            if product_mol is not None:
                                ring_tag = flavonoid_ring_label(parent_mol, site)
                                props = {
                                    'parent_smiles': parent_smiles,
                                    'parent_inchikey': parent_inchikey,
                                    'rule': 'R2a',
                                    'site_index': site,
                                    'sym_class': None,
                                    'ring_tag': ring_tag,
                                    'halogen': halogen,
                                    'depth': 1
                                }
                                products.append((product_mol, props))
                                produced_any = True
                        except Exception as e:
                            LOG.debug(f"R2a site halogenation failed for site {site} + {halogen}: {e}")
                            qa_events['no_product_matches'] = 1

                        # Record R2a attempt
                        if aggregator:
                            aggregator.record_attempt_result('R2a', halogen, 1,
                                                           produced_count=(1 if produced_any else 0),
                                                           qa_events_dict=qa_events,
                                                           k_ops=None, k_atoms=None)
            except Exception as e:
                LOG.warning(f"R2a processing failed: {e}")
                for halogen in r2_allowed_halogens:
                    if aggregator:
                        aggregator.record_attempt_result('R2a', halogen, 1,
                                                       produced_count=0,
                                                       qa_events_dict={'template_unsupported': 1},
                                                       k_ops=None, k_atoms=None)

        # R2b: sp3 CH2 sites in flavanone C-ring (controlled by rules_cfg.R2.sp3_CH2_flavanone)
        if r2_config.get('sp3_CH2_flavanone', True):
            try:
                sugar_cfg = config.get('sugar_cfg', {})
                # Request detection information to track strict vs fallback enumeration
                r2b_sites, r2b_used_fallback = c_ring_sp3_CH2_flavanone_sites(
                    parent_mol, set(), sugar_cfg, config.get('rules_cfg'),
                    return_detection=True
                )
                LOG.info("[R2b] Site enumeration: found %d candidates (sugar_cfg.mode=%s, detection=%s)",
                         len(r2b_sites), sugar_cfg.get('mode', 'heuristic'),
                         'fallback' if r2b_used_fallback else 'strict')
                if r2b_sites:
                    LOG.debug("[R2b] Sites: %s", r2b_sites[:10])  # Show first 10

                for site in r2b_sites:
                    for halogen in r2_allowed_halogens:
                        produced_any = False
                        qa_events = {}

                        try:
                            product_mol = apply_single_site_halogenation(parent_mol, site, halogen)
                            if product_mol is not None:
                                ring_tag = flavonoid_ring_label(parent_mol, site)
                                props = {
                                    'parent_smiles': parent_smiles,
                                    'parent_inchikey': parent_inchikey,
                                    'rule': 'R2b',
                                    'site_index': site,
                                    'sym_class': None,
                                    'ring_tag': ring_tag,
                                    'halogen': halogen,
                                    'depth': 1,
                                    'detection': 'fallback' if r2b_used_fallback else 'strict'  # Track detection method
                                }
                                products.append((product_mol, props))
                                produced_any = True
                        except Exception as e:
                            LOG.debug(f"R2b site halogenation failed for site {site} + {halogen}: {e}")
                            qa_events['no_product_matches'] = 1

                        # Record R2b attempt
                        if aggregator:
                            aggregator.record_attempt_result('R2b', halogen, 1,
                                                           produced_count=(1 if produced_any else 0),
                                                           qa_events_dict=qa_events,
                                                           k_ops=None, k_atoms=None)
            except Exception as e:
                LOG.warning(f"R2b processing failed: {e}")
                for halogen in r2_allowed_halogens:
                    if aggregator:
                        aggregator.record_attempt_result('R2b', halogen, 1,
                                                       produced_count=0,
                                                       qa_events_dict={'template_unsupported': 1},
                                                       k_ops=None, k_atoms=None)

    # Step 3: Apply R6 methyl halogenation rules (if enabled)
    if 'R6_methyl' in rules:
        # Get R6 configuration from config
        r6_config = config.get('rules_cfg', {}).get('R6_methyl', {})

        if r6_config.get('enable', False):
            try:
                # Get allowed halogens for R6 (intersect with global halogens)
                r6_allowed_raw = r6_config.get('allowed', ['F', 'Cl'])
                r6_allowed = [h for h in r6_allowed_raw if h in halogens]

                # Get R6 configuration parameters
                allow_on_methoxy = r6_config.get('allow_on_methoxy', False)
                allow_allylic_methyl = r6_config.get('allow_allylic_methyl', False)

                # Use computed sugar_mask for consistency with k>=2 path
                r6_sites = enumerate_methyl_sites(parent_mol, sugar_mask, allow_on_methoxy, allow_allylic_methyl)
                LOG.debug(f"ENUM_K1 DEBUG: Found {len(r6_sites)} methyl sites with allow_on_methoxy={allow_on_methoxy}, allow_allylic_methyl={allow_allylic_methyl}")
                for i, site in enumerate(r6_sites):
                    if isinstance(site, dict):
                        LOG.debug(f"ENUM_K1 DEBUG: Site {i}: idx={site['idx']}, kind={site.get('kind', 'UNKNOWN')}")
                    else:
                        LOG.debug(f"ENUM_K1 DEBUG: Site {i}: {site} (legacy format)")

                # Each site x halogen is one attempt
                for site in r6_sites:
                    # Extract site index from structured descriptor (backward compatibility)
                    if isinstance(site, dict):
                        site_idx = site['idx']
                    else:
                        site_idx = site  # Legacy integer format

                    for halogen in r6_allowed:
                        produced_any = False
                        qa_events = {}

                        try:
                            # Validate halogen is supported for R6
                            if not validate_methyl_halogen(halogen):
                                qa_events['template_unsupported'] = 1
                            else:
                                # Apply methyl halogenation (pass full site descriptor)
                                product_mol = apply_methyl_step(parent_mol, site, halogen)

                                if product_mol is not None:
                                    # Get ring tag for the site (use integer index)
                                    ring_tag = flavonoid_ring_label(parent_mol, site_idx)

                                    props = {
                                        'parent_smiles': parent_smiles,
                                        'parent_inchikey': parent_inchikey,
                                        'rule': 'R6_methyl',  # Use R6_methyl for consistency with config
                                        'site_index': site_idx,  # Store integer index
                                        'sym_class': None,
                                        'ring_tag': ring_tag,
                                        'halogen': halogen,
                                        'depth': 1
                                    }
                                    products.append((product_mol, props))
                                    produced_any = True
                                else:
                                    qa_events['no_product_matches'] = 1

                        except Exception as e:
                            LOG.debug(f"R6 methyl halogenation failed for site {site} + {halogen}: {e}")
                            qa_events['no_product_matches'] = 1

                        # Record this single site x halogen attempt
                        if aggregator:
                            aggregator.record_attempt_result('R6_methyl', halogen, 1,
                                                           produced_count=(1 if produced_any else 0),
                                                           qa_events_dict=qa_events,
                                                           k_ops=None, k_atoms=None)

            except Exception as e:
                # R6 processing failed entirely - record template_unsupported attempts for all allowed halogens
                LOG.warning(f"R6 methyl processing failed entirely: {e}")
                r6_allowed_fallback = r6_config.get('allowed', ['F', 'Cl'])
                r6_allowed_fallback = [h for h in r6_allowed_fallback if h in halogens]

                for halogen in r6_allowed_fallback:
                    if aggregator:
                        aggregator.record_attempt_result('R6_methyl', halogen, 1,
                                                       produced_count=0,
                                                       qa_events_dict={'template_unsupported': 1},
                                                       k_ops=None, k_atoms=None)

            # Process R6_methyl macro substitution if enabled
            macro_cfg = r6_config.get('macro', {})
            if macro_cfg.get('enable', False):
                try:
                    macro_labels = macro_cfg.get('labels', ['CF3', 'CCl3'])
                    LOG.debug(f"ENUM_K1 MACRO DEBUG: Processing macro substitution with labels={macro_labels}")

                    # Re-enumerate sites for macro substitution (same as step substitution)
                    r6_macro_sites = enumerate_methyl_sites(parent_mol, sugar_mask, allow_on_methoxy, allow_allylic_methyl)

                    # Each site x label is one attempt
                    for site in r6_macro_sites:
                        # Extract site index from structured descriptor
                        if isinstance(site, dict):
                            site_idx = site['idx']
                        else:
                            site_idx = site  # Legacy integer format

                        for label in macro_labels:
                            # Determine halogen from label (CF3 -> F, CCl3 -> Cl)
                            halogen = 'F' if label == 'CF3' else 'Cl'

                            # Check if halogen is in allowed list
                            if halogen not in r6_allowed:
                                continue

                            produced_any = False
                            qa_events = {}

                            try:
                                # Validate halogen is supported
                                if not validate_methyl_halogen(halogen):
                                    qa_events['template_unsupported'] = 1
                                else:
                                    # Apply macro halogenation
                                    product_mol = apply_methyl_macro(parent_mol, site_idx, label)

                                    if product_mol is not None:
                                        # Get ring tag for the site
                                        ring_tag = flavonoid_ring_label(parent_mol, site_idx)

                                        # Create props dict with macro metadata
                                        props = {
                                            'parent_smiles': parent_smiles,
                                            'parent_inchikey': parent_inchikey,
                                            'rule': 'R6_methyl',
                                            'site_index': site_idx,
                                            'sym_class': None,
                                            'ring_tag': ring_tag,
                                            'halogen': halogen,
                                            'depth': 1,
                                            'substitution_type': 'macro',  # Mark as macro substitution
                                            'macro_label': label  # Store macro label
                                        }
                                        products.append((product_mol, props))
                                        produced_any = True
                                        LOG.debug(f"ENUM_K1 MACRO DEBUG: Generated macro product for site {site_idx} + {label}")
                                    else:
                                        qa_events['no_product_matches'] = 1

                            except Exception as e:
                                LOG.debug(f"R6 macro substitution failed for site {site} + {label}: {e}")
                                qa_events['no_product_matches'] = 1

                            # Record this single site x label attempt
                            if aggregator:
                                aggregator.record_attempt_result('R6_methyl', halogen, 1,
                                                               produced_count=(1 if produced_any else 0),
                                                               qa_events_dict=qa_events,
                                                               k_ops=1,  # Macro is 1 operation
                                                               k_atoms=3)  # Macro is 3 atoms (CF3 or CCl3)

                except Exception as e:
                    # Macro processing failed entirely
                    LOG.warning(f"R6 macro processing failed entirely: {e}")
                    # No need to record attempts here as we couldn't enumerate sites

    return products


def enumerate_k1_halogenation(
    parent_mol,
    halogens: List[str],
    rules: List[str],
    config: Dict[str, Any]
):
    """Enumerate k=1 halogenated products with properties."""
    if parent_mol is None:
        return []

    try:
        from .chem_compat import Chem
        # SanitizeFlags are optional - try to import but don't fail if unavailable
        try:
            from rdkit.Chem import SanitizeFlags
        except ImportError:
            SanitizeFlags = None
    except Exception:
        LOG.debug("RDKit not available")
        return []

    # Ensure molecule is properly sanitized before any operations
    try:
        Chem.SanitizeMol(parent_mol)
    except Exception:
        parent_mol.UpdatePropertyCache(strict=False)

    parent_smiles = Chem.MolToSmiles(parent_mol, canonical=True)
    parent_inchikey = to_inchikey(parent_mol)

    # CRITICAL FIX: Compute sugar mask for consistency with k>=2 path (moved before R6 section)
    sugar_cfg = config.get('sugar_cfg', {})
    sugar_mask, mask_degraded, status_metadata = get_sugar_mask_with_full_status(
        parent_mol,
        mode=sugar_cfg.get('mode', 'heuristic'),
        sugar_cfg=sugar_cfg
    )
    LOG.debug(f"k=1 sugar mask computed: {len(sugar_mask)} atoms masked, degraded={mask_degraded}")

    products = []
    reactions = build_reactions()
    
    # Step 1: Apply reaction-based rules (R1, R3, R4, R5)
    for rule in ['R1', 'R3', 'R4', 'R5']:
        if rule in rules and rule in reactions:
            for halogen in halogens:
                if halogen in reactions[rule]:
                    rxn = reactions[rule][halogen]
                    local_seen = set()  # Local deduplication per (rule, halogen)
                    reaction_products = _run_reaction_safely(rxn, (parent_mol,), rule, halogen, None, aggregator=None, current_k=1)
                    for product_set in reaction_products:
                        for product_mol in product_set:
                            if product_mol is not None:
                                # Sanitize product molecule (lightweight first, full if needed)
                                try:
                                    Chem.SanitizeMol(product_mol, sanitizeOps=SanitizeFlags.SANITIZE_PROPERTIES)
                                except Exception:
                                    try:
                                        Chem.SanitizeMol(product_mol)
                                    except Exception:
                                        product_mol.UpdatePropertyCache(strict=False)
                                
                                # Ensure ring info and properties are ready for downstream operations
                                try:
                                    ensure_ready(product_mol)
                                except Exception:
                                    pass  # Continue with potentially incomplete initialization
                                
                                try:
                                    product_smiles = Chem.MolToSmiles(product_mol, canonical=True)
                                except Exception:
                                    # fallback: skip local dedup on this product, but keep pipeline alive
                                    product_smiles = None
                                
                                if product_smiles is None or product_smiles not in local_seen:
                                    if product_smiles is not None:
                                        local_seen.add(product_smiles)
                                    props = {
                                        'parent_smiles': parent_smiles,
                                        'parent_inchikey': parent_inchikey,
                                        'rule': rule,
                                        'site_index': None,
                                        'sym_class': None,
                                        'ring_tag': None,
                                        'halogen': halogen,
                                        'depth': 1
                                    }
                                    products.append((product_mol, props))
    
    # Step 2: Apply site-based rules with symmetry folding (R2)  
    # Note: R1 now uses reaction engine above, may add back symmetry folding in future
    
    if 'R2' in rules:
        r2_products = _apply_r2_rule(
            parent_mol, parent_smiles, parent_inchikey, halogens
        )
        products.extend(r2_products)
    
    # Step 3: QC filtering
    if config.get('qc', {}).get('sanitize_strict', True):
        products = [(mol, props) for mol, props in products if sanitize_ok(mol)]
    
    # Step 4: Deduplication
    if config.get('dedupe', {}).get('method', 'inchikey') == 'inchikey':
        products = dedupe_with_props(products)
    
    # Step 5: Add descriptors and final properties
    final_products = []
    for mol, props in products:
        if mol is None:
            continue
        
        # Add descriptors
        descriptors = basic_descriptors(mol)
        props.update(descriptors)
        
        # Add QC flags
        props['sanitize_ok'] = sanitize_ok(mol)
        props['pains_flags'] = pains_flags(mol)
        
        final_products.append((mol, props))
    
    return final_products






def _apply_r2_rule(parent_mol, parent_smiles: str, parent_inchikey: str,
                   halogens: List[str]):
    """Apply R2 rule with symmetry folding."""
    try:
        from .chem_compat import Chem
    except Exception:
        LOG.debug("RDKit not available")
        return []
    
    products = []
    
    # Ensure parent molecule is ready for ring operations
    ensure_ready(parent_mol)
    
    # Get C ring indices
    cring_indices = c_ring_indices(parent_mol)
    if not cring_indices:
        return products
    
    # Group by symmetry with ring tags
    indexed_sites = []
    for idx in cring_indices:
        ring_tag = get_ring_tag_for_atom(parent_mol, idx)
        indexed_sites.append((idx, ring_tag))
    
    # Group by (rank, ring_tag)
    ranks = Chem.CanonicalRankAtoms(parent_mol, breakTies=False)
    groups = {}
    
    for idx, ring_tag in indexed_sites:
        rank = ranks[idx]
        key = (rank, ring_tag)
        if key not in groups:
            groups[key] = []
        groups[key].append(idx)
    
    # For each group, take representative
    for (rank, ring_tag), indices in groups.items():
        representative_idx = indices[0]
        
        # Apply each halogen
        for halogen in halogens:
            product_mol = apply_single_site_halogenation(parent_mol, representative_idx, halogen)
            if product_mol is not None:
                props = {
                    'parent_smiles': parent_smiles,
                    'parent_inchikey': parent_inchikey,
                    'rule': 'R2',
                    'site_index': representative_idx,
                    'sym_class': str(rank),
                    'ring_tag': str(ring_tag) if ring_tag is not None else None,
                    'halogen': halogen,
                    'depth': 1
                }
                products.append((product_mol, props))
    
    return products