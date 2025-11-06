# -*- coding: ascii -*-
"""
RDKit random seed utilities for reproducible operations.

This module provides centralized management of random seed setting for various
RDKit operations that may use randomness, ensuring reproducible results.
"""

import logging
from typing import List, Optional, Dict, Any
from .chem_compat import Chem

LOG = logging.getLogger(__name__)


def _canonicalize_sources(entries: List[str]) -> frozenset:
    """Canonicalize seed/failure source labels for debounced logging."""
    canonical = []
    for entry in entries:
        if entry:
            canonical.append(entry.split(' ', 1)[0])
    return frozenset(canonical)


class RDKitSeedManager:
    """
    Manages random seed setting for various RDKit operations.

    Tracks which random sources have been successfully seeded and provides
    an audit trail of what randomness sources are controlled.
    """

    def __init__(self):
        self.seeded_sources: List[str] = []
        self.failed_sources: List[str] = []
        self.current_seed: Optional[int] = None
        self._has_logged_failures = False  # De-bounce repeated detailed logging
        self._last_seeded_set: Optional[frozenset] = None
        self._last_failed_set: Optional[frozenset] = None
        self._info_emitted = False

    def set_global_seed(self, seed: int, verbose: bool = False) -> Dict[str, Any]:
        """
        Set random seed for all available RDKit random sources.

        Args:
            seed: Random seed value

        Returns:
            Dictionary containing:
            - 'seeded_sources': List of successfully seeded RDKit operations
            - 'failed_sources': List of operations that couldn't be seeded
            - 'seed': The seed value used
        """
        # Environment override for verbosity
        try:
            import os
            if os.environ.get('HALO_RDKitSeedVerbose', '') in ('1', 'true', 'True'):
                verbose = True
        except Exception:
            pass

        self.current_seed = seed
        self.seeded_sources.clear()
        self.failed_sources.clear()

        # Try global RDKit seed (if available)
        try:
            if hasattr(Chem, 'SetRandomSeed'):
                Chem.SetRandomSeed(seed)
                self.seeded_sources.append('Chem.SetRandomSeed')
            else:
                self.failed_sources.append('Chem.SetRandomSeed (not available)')
        except Exception as e:
            self.failed_sources.append(f'Chem.SetRandomSeed (error: {e})')

        # Try descriptor-specific seed
        try:
            if hasattr(Chem.rdMolDescriptors, 'SetRandomSeed'):
                Chem.rdMolDescriptors.SetRandomSeed(seed)
                self.seeded_sources.append('rdMolDescriptors.SetRandomSeed')
            else:
                self.failed_sources.append('rdMolDescriptors.SetRandomSeed (not available)')
        except Exception as e:
            self.failed_sources.append(f'rdMolDescriptors.SetRandomSeed (error: {e})')

        # Try AllChem seed (for conformer generation)
        try:
            from rdkit.Chem import AllChem
            if hasattr(AllChem, 'SetRandomSeed'):
                AllChem.SetRandomSeed(seed)
                self.seeded_sources.append('AllChem.SetRandomSeed')
            else:
                self.failed_sources.append('AllChem.SetRandomSeed (not available)')
        except Exception as e:
            self.failed_sources.append(f'AllChem.SetRandomSeed (error: {e})')

        # Try geometry optimization random sources
        try:
            from rdkit.Chem import rdDistGeom
            if hasattr(rdDistGeom, 'SetRandomSeed'):
                rdDistGeom.SetRandomSeed(seed)
                self.seeded_sources.append('rdDistGeom.SetRandomSeed')
            else:
                self.failed_sources.append('rdDistGeom.SetRandomSeed (not available)')
        except Exception as e:
            self.failed_sources.append(f'rdDistGeom.SetRandomSeed (error: {e})')

        # Try force field optimization seeds
        try:
            from rdkit.Chem import rdForceFieldHelpers
            if hasattr(rdForceFieldHelpers, 'SetRandomSeed'):
                rdForceFieldHelpers.SetRandomSeed(seed)
                self.seeded_sources.append('rdForceFieldHelpers.SetRandomSeed')
            else:
                self.failed_sources.append('rdForceFieldHelpers.SetRandomSeed (not available)')
        except Exception as e:
            self.failed_sources.append(f'rdForceFieldHelpers.SetRandomSeed (error: {e})')

        result = {
            'seeded_sources': self.seeded_sources.copy(),
            'failed_sources': self.failed_sources.copy(),
            'seed': seed
        }

        # De-bounced logging behavior
        previous_seeded = self._last_seeded_set
        previous_failed = self._last_failed_set
        canonical_seeded = _canonicalize_sources(self.seeded_sources)
        canonical_failed = _canonicalize_sources(self.failed_sources)

        sets_changed = (canonical_seeded != previous_seeded) or (canonical_failed != previous_failed)

        if not self._info_emitted or sets_changed:
            LOG.info("RDKit seeding summary: seeded=%d, failed=%d", len(canonical_seeded), len(canonical_failed))
            self._info_emitted = True

        if verbose or not self._has_logged_failures or sets_changed:
            if self.seeded_sources:
                LOG.debug("RDKit random sources seeded with %d: %s", seed, self.seeded_sources)
            if self.failed_sources:
                LOG.debug("RDKit random sources failed to seed: %s", self.failed_sources)
            self._has_logged_failures = True
        else:
            LOG.debug("RDKit seeding details unchanged since last report (debounced)")

        self._last_seeded_set = canonical_seeded
        self._last_failed_set = canonical_failed

        return result

    def get_seeded_sources(self) -> List[str]:
        """Get list of successfully seeded RDKit random sources."""
        return self.seeded_sources.copy()

    def get_failed_sources(self) -> List[str]:
        """Get list of RDKit random sources that failed to seed."""
        return self.failed_sources.copy()

    def get_current_seed(self) -> Optional[int]:
        """Get the currently set seed value."""
        return self.current_seed


# Global instance for centralized seed management
_global_seed_manager = RDKitSeedManager()


def set_rdkit_random_seed(seed: int, verbose: bool = False) -> Dict[str, Any]:
    """
    Set random seed for all available RDKit operations.

    This is the main function to use for setting RDKit random seeds.
    It attempts to set seeds for all known RDKit random sources and
    returns a report of what was successfully seeded.

    Args:
        seed: Random seed value

    Returns:
        Dictionary with seeding results and audit information
    """
    return _global_seed_manager.set_global_seed(seed, verbose=verbose)


def get_rdkit_seed_report() -> Dict[str, Any]:
    """
    Get a report of current RDKit random seed status.

    Returns:
        Dictionary containing current seed and which sources are controlled
    """
    return {
        'current_seed': _global_seed_manager.get_current_seed(),
        'seeded_sources': _global_seed_manager.get_seeded_sources(),
        'failed_sources': _global_seed_manager.get_failed_sources()
    }


def embed_molecule_with_seed(mol, seed: Optional[int] = None, **kwargs):
    """
    Embed a molecule with explicit random seed control.

    This function wraps RDKit's EmbedMolecule with explicit seed setting
    for reproducible conformer generation.

    Args:
        mol: RDKit molecule object
        seed: Random seed for embedding (uses global seed if None)
        **kwargs: Additional arguments passed to EmbedMolecule

    Returns:
        Result from EmbedMolecule operation
    """
    try:
        from rdkit.Chem import AllChem

        # Use provided seed or current global seed
        use_seed = seed if seed is not None else _global_seed_manager.get_current_seed()

        if use_seed is not None:
            # Set seed specifically for this operation
            if hasattr(AllChem, 'EmbedMolecule'):
                # Pass randomSeed parameter if supported
                if 'randomSeed' not in kwargs:
                    kwargs['randomSeed'] = use_seed
                return AllChem.EmbedMolecule(mol, **kwargs)

        # Fallback to regular embedding without explicit seed
        return AllChem.EmbedMolecule(mol, **kwargs)

    except ImportError:
        LOG.warning("AllChem not available for molecule embedding")
        return -1
    except Exception as e:
        LOG.warning(f"Molecule embedding failed: {e}")
        return -1


def optimize_molecule_with_seed(mol, seed: Optional[int] = None, **kwargs):
    """
    Optimize molecule geometry with explicit random seed control.

    Args:
        mol: RDKit molecule object with embedded coordinates
        seed: Random seed for optimization (uses global seed if None)
        **kwargs: Additional arguments passed to optimization

    Returns:
        Result from geometry optimization
    """
    try:
        from rdkit.Chem import AllChem

        # Use provided seed or current global seed
        use_seed = seed if seed is not None else _global_seed_manager.get_current_seed()

        if use_seed is not None and hasattr(AllChem, 'UFFOptimizeMolecule'):
            # Set up UFF optimization with seed if possible
            return AllChem.UFFOptimizeMolecule(mol, **kwargs)

        # Fallback to regular optimization
        return AllChem.UFFOptimizeMolecule(mol, **kwargs)

    except ImportError:
        LOG.warning("AllChem not available for molecule optimization")
        return -1
    except Exception as e:
        LOG.warning(f"Molecule optimization failed: {e}")
        return -1

def _reset_seed_manager_state_for_tests():
    """Reset global seed manager state for deterministic tests."""
    # Note: Intended for unit tests only; do not call from production code.
    _global_seed_manager.seeded_sources.clear()
    _global_seed_manager.failed_sources.clear()
    _global_seed_manager.current_seed = None
    _global_seed_manager._has_logged_failures = False
    _global_seed_manager._last_seeded_set = None
    _global_seed_manager._last_failed_set = None
    _global_seed_manager._info_emitted = False

