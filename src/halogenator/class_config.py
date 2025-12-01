# -*- coding: ascii -*-
"""
Per-class halogenation configuration loader.

Loads and applies NP class-specific halogenation rules from
configs/halogen_rules_by_class.yaml.
"""

import logging
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List

LOG = logging.getLogger(__name__)

# Mapping from semantic rule IDs to legacy rule IDs
# The enumerate functions currently only recognize legacy IDs
SEMANTIC_TO_LEGACY_RULES = {
    'RING_SP2__CH__TO__X': 'R1',
    'COOH__TO__CX': 'R5',
    'ALPHA_CARBONYL__CH2__TO__X': 'ALPHA_CARBONYL__CH2__TO__X',  # Keep as-is
    'RING_SP3__CH__TO__X': 'RING_SP3__CH__TO__X',  # Keep as-is
    'PRIMARY_OH__CH2OH__TO__X': 'PRIMARY_OH__CH2OH__TO__X',  # Keep as-is
    'R1': 'R1',
    'R3': 'R3',
    'R4': 'R4',
    'R5': 'R5',
}


def load_class_halogen_config(
    np_class: str,
    k: int,
    config_path: Optional[Path] = None
) -> Dict[str, Any]:
    """
    Load halogenation configuration for a specific NP class and k value.

    Args:
        np_class: NP class name (e.g., 'terpenoid', 'alkaloid', 'polyphenol')
        k: Halogenation depth (1 or 2)
        config_path: Optional path to config YAML (defaults to configs/halogen_rules_by_class.yaml)

    Returns:
        Configuration dict with keys:
        - enabled: bool
        - rules: List[str] - rule IDs to use
        - halogens: List[str] - halogens to use
        - max_sites_per_parent: int - max sites per molecule (-1 = unlimited)
        - sugar_mask: bool - whether to apply sugar masking
        - per_rule_overrides: Dict[str, Any] - per-rule specific overrides
        - batch_size: int
        - rdkit_threads: int
    """
    if config_path is None:
        # Default to configs/halogen_rules_by_class.yaml in repo root
        config_path = Path(__file__).parent.parent.parent / "configs" / "halogen_rules_by_class.yaml"

    if not config_path.exists():
        LOG.warning(f"Class config not found: {config_path}, using defaults")
        return _default_config()

    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)
    except Exception as e:
        LOG.error(f"Failed to load config from {config_path}: {e}")
        return _default_config()

    classes = data.get('classes', {})
    class_cfg = classes.get(np_class, {})

    if not class_cfg:
        LOG.warning(f"Class '{np_class}' not found in config, using defaults")
        return _default_config()

    # Check if class is enabled
    if class_cfg.get('enabled') == False:
        LOG.warning(f"Class '{np_class}' is disabled in config")
        return {'enabled': False}

    # Get k-specific config
    k_key = f'k{k}'
    k_cfg = class_cfg.get(k_key, {})

    if not k_cfg:
        LOG.warning(f"No k{k} config for class '{np_class}', using defaults")
        return _default_config()

    # Build result
    result = {
        'enabled': True,
        'rules': k_cfg.get('include_rules', []),
        'halogens': data.get('defaults', {}).get('halogens', ['F', 'Cl', 'Br', 'I']),
        'max_sites_per_parent': k_cfg.get('max_sites_per_parent', -1),
        'sugar_mask': k_cfg.get('sugar_mask', False),
        'per_rule_overrides': class_cfg.get('per_rule_overrides', {}),
        'batch_size': data.get('defaults', {}).get('batch_size', 5000),
        'rdkit_threads': data.get('defaults', {}).get('rdkit_threads', 8),
    }

    LOG.info(
        f"Loaded config for class='{np_class}' k={k}: "
        f"{len(result['rules'])} rules, max_sites={result['max_sites_per_parent']}, "
        f"sugar_mask={result['sugar_mask']}"
    )

    return result


def _default_config() -> Dict[str, Any]:
    """Return default configuration when class config is not available."""
    return {
        'enabled': True,
        'rules': ['HALO_ARYL_CH__H__TO__X', 'HALO_CARBOXYL__COOH__TO__COX'],
        'halogens': ['F', 'Cl', 'Br', 'I'],
        'max_sites_per_parent': -1,
        'sugar_mask': False,
        'per_rule_overrides': {},
        'batch_size': 5000,
        'rdkit_threads': 8,
    }


def apply_rule_overrides(
    rules: List[str],
    overrides: Dict[str, Any]
) -> List[str]:
    """
    Apply per-rule overrides to filter and modify rule list.

    Args:
        rules: Original list of rule IDs
        overrides: Dict mapping rule_id -> {enabled: bool, ...}

    Returns:
        Filtered list of rules with disabled rules removed
    """
    filtered = []
    for rule in rules:
        override = overrides.get(rule, {})
        if override.get('enabled', True):  # Default to enabled
            filtered.append(rule)
        else:
            LOG.info(f"Rule '{rule}' disabled by override")

    return filtered


def map_semantic_to_legacy_rules(rules: List[str]) -> List[str]:
    """
    Map semantic rule IDs to legacy rule IDs for enumerate compatibility.

    Args:
        rules: List of semantic or legacy rule IDs

    Returns:
        List of legacy rule IDs
    """
    mapped = []
    for rule in rules:
        legacy = SEMANTIC_TO_LEGACY_RULES.get(rule, rule)
        mapped.append(legacy)
        if legacy != rule:
            LOG.debug(f"Mapped rule '{rule}' -> '{legacy}'")

    return mapped
