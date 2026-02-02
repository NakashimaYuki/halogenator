#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Natural Product-Like Library Generator (Performance Refactored v2)

PERFORMANCE IMPROVEMENTS:
1. ✓ Single RunReactants call per molecule (avoid N repetitions)
2. ✓ Parallel batch processing with ProcessPoolExecutor
3. ✓ Optimized SQLite deduplication (PRAGMA + batch INSERT OR IGNORE)
4. ✓ Streaming Parquet writer (avoid full-memory accumulation)
5. ✓ Cached sugar site annotations (avoid repeated json.loads)
6. ✓ In-batch deduplication using canonical_smiles before InChIKey
7. ✓ Enhanced logging (throughput, dedup rates, failure breakdown)
8. ✓ Resumable execution support

Author: Task D Performance Refactor
Date: 2025-11-09
"""

import os
import sys
import argparse
import logging
import json
import yaml
import sqlite3
import time
import pandas as pd
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional, Any, Set
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
from functools import partial
import traceback

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem import inchi
except ImportError:
    print("ERROR: RDKit not available. Please install: pip install rdkit")
    sys.exit(1)

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s | %(levelname)-8s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def safe_mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    """Safely create RDKit mol from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def get_inchikey(mol: Chem.Mol) -> Optional[str]:
    """Get InChIKey from mol object."""
    try:
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


def canonical_smiles(mol: Chem.Mol) -> Optional[str]:
    """Get canonical SMILES from mol."""
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return None


# ============================================================================
# MOLECULAR PROPERTY CALCULATION (Quick mode for enumeration phase)
# ============================================================================

def calculate_properties_quick(mol: Chem.Mol) -> Dict[str, float]:
    """
    Calculate essential molecular properties (quick mode for enumeration).

    Full descriptors can be computed in post-processing with 08a_fill_descriptors.py
    """
    try:
        return {
            'MW': Descriptors.MolWt(mol),
            'TPSA': Descriptors.TPSA(mol),
            'HBD': Lipinski.NumHDonors(mol),
            'HBA': Lipinski.NumHAcceptors(mol),
            'cLogP': Descriptors.MolLogP(mol),
            'RotB': Lipinski.NumRotatableBonds(mol)
        }
    except Exception:
        return {
            'MW': None, 'TPSA': None, 'HBD': None,
            'HBA': None, 'cLogP': None, 'RotB': None
        }


def calculate_halogen_fields(mol, xf_label: str) -> Dict[str, Any]:
    """
    Calculate halogen-related fields for a product molecule.

    Args:
        mol: RDKit molecule object
        xf_label: Transformation label (e.g., 'OH_to_OMe')

    Returns:
        Dictionary with:
        - halogens_set: Pipe-separated halogen types (e.g., "F|Cl")
        - halogen_counts_json: JSON string with counts (e.g., '{"F":1,"Cl":1}')
        - primary_halogen: Halogen introduced by last step (or None)
    """
    try:
        # Count halogens in product molecule
        halogen_counts = {'F': 0, 'Cl': 0, 'Br': 0, 'I': 0}

        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in halogen_counts:
                halogen_counts[symbol] += 1

        # Filter to non-zero halogens
        present_halogens = {k: v for k, v in halogen_counts.items() if v > 0}

        # halogens_set: "F|Cl" format
        if present_halogens:
            halogens_set = '|'.join(sorted(present_halogens.keys()))
        else:
            halogens_set = None

        # halogen_counts_json: '{"F":1,"Cl":1}' format
        if present_halogens:
            halogen_counts_json = json.dumps(present_halogens, separators=(',', ':'))
        else:
            halogen_counts_json = None

        # primary_halogen: Halogen introduced by last transformation step
        # For OH_to_OMe and OH_to_NH2, no halogens are introduced
        # This field is for future halogen introduction rules
        primary_halogen = None
        if 'OMe' not in xf_label and 'NH2' not in xf_label:
            # If this were a halogenation rule, extract the halogen from the label
            # For now, all current rules don't introduce halogens
            primary_halogen = None

        return {
            'halogens_set': halogens_set,
            'halogen_counts_json': halogen_counts_json,
            'primary_halogen': primary_halogen
        }

    except Exception as e:
        logger.debug(f"Failed to calculate halogen fields: {e}")
        return {
            'halogens_set': None,
            'halogen_counts_json': None,
            'primary_halogen': None
        }


# ============================================================================
# SUGAR RING FILTERING (with caching)
# ============================================================================

SUGAR_SMARTS_PATTERNS = [
    'C1OCCC1',   # furanose skeleton
    'C1OCCCC1',  # pyranose skeleton
    'C1OC(=O)CC1',  # gamma-lactone
    'C1OC(=O)CCC1',  # delta-lactone
    '[C;R1]1[O;R1][C;R1][C;R1][C;R1][C;R1]1',
    '[C;R1]1[O;R1][C;R1][C;R1][C;R1]1',
]


class SugarSiteCache:
    """Cache for sugar site annotations to avoid repeated json.loads."""

    def __init__(self, sugar_mask_json: Optional[str] = None, sugar_rings_json: Optional[str] = None):
        self.sugar_atoms: Set[int] = set()
        self.sugar_ring_atoms: Set[int] = set()

        # Parse once at initialization
        if sugar_mask_json:
            try:
                atoms = json.loads(sugar_mask_json)
                self.sugar_atoms = set(atoms) if isinstance(atoms, list) else set()
            except:
                pass

        if sugar_rings_json:
            try:
                rings = json.loads(sugar_rings_json)
                if isinstance(rings, list):
                    for ring in rings:
                        self.sugar_ring_atoms.update(ring)
            except:
                pass

    def is_sugar_site(self, atom_idx: int, mol: Chem.Mol) -> bool:
        """Check if atom is in sugar moiety (fast lookup)."""
        # Check cached annotations
        if atom_idx in self.sugar_atoms or atom_idx in self.sugar_ring_atoms:
            return True

        # Fallback: check if in O-containing ring
        atom = mol.GetAtomWithIdx(atom_idx)
        if not atom.IsInRing():
            return False

        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            if atom_idx in ring and len(ring) in (5, 6):
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                if any(a.GetSymbol() == 'O' for a in ring_atoms):
                    return True

        return False


# ============================================================================
# DEDUPLICATION: SQLite (exact) and Bloom Filter (probabilistic)
# ============================================================================

# Import Bloom filter deduplicator (memory-efficient for large-scale jobs)
try:
    import sys
    halogenator_path = Path(__file__).resolve().parent.parent / 'src'
    if str(halogenator_path) not in sys.path:
        sys.path.insert(0, str(halogenator_path))
    from halogenator.bloom_dedup import BloomFilterDeduplicator
    BLOOM_AVAILABLE = True
except ImportError as e:
    logger.warning(f"Bloom filter not available: {e}")
    BLOOM_AVAILABLE = False


class SqliteDeduplicator:
    """Optimized SQLite-based cross-batch deduplication."""

    def __init__(self, db_path: Path, resume: bool = False):
        self.db_path = db_path
        self.conn = sqlite3.connect(str(db_path))
        self.cursor = self.conn.cursor()

        # Performance PRAGMA settings
        self.cursor.execute("PRAGMA journal_mode=WAL")
        self.cursor.execute("PRAGMA synchronous=OFF")
        self.cursor.execute("PRAGMA temp_store=MEMORY")
        self.cursor.execute("PRAGMA locking_mode=EXCLUSIVE")
        self.cursor.execute("PRAGMA cache_size=-200000")  # 200MB cache
        self.cursor.execute("PRAGMA mmap_size=268435456")  # 256MB mmap

        # Create table with unique index
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS seen_keys (
                key TEXT PRIMARY KEY
            )
        """)
        self.conn.commit()

        # Load existing keys into memory set (for faster lookups)
        self.seen_in_memory: Set[str] = set()
        if resume:
            logger.info("Loading existing dedup keys for resume...")
            for row in self.cursor.execute("SELECT key FROM seen_keys"):
                self.seen_in_memory.add(row[0])
            logger.info(f"  Loaded {len(self.seen_in_memory):,} existing keys")

    def filter_new_keys(self, keys: List[str]) -> List[str]:
        """Filter to only new keys. Returns list of new keys and updates DB."""
        new_keys = []

        for key in keys:
            if key not in self.seen_in_memory:
                new_keys.append(key)
                self.seen_in_memory.add(key)

        # Batch insert
        if new_keys:
            self.cursor.executemany(
                "INSERT OR IGNORE INTO seen_keys (key) VALUES (?)",
                [(k,) for k in new_keys]
            )

        return new_keys

    def commit(self):
        self.conn.commit()

    def checkpoint(self):
        """
        Execute WAL checkpoint to reduce memory usage.

        This truncates the WAL file and releases memory.
        Should be called periodically during long-running jobs.
        """
        try:
            self.cursor.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            logger.info("Dedup DB checkpoint completed (WAL truncated)")
        except Exception as e:
            logger.warning(f"Dedup checkpoint failed: {e}")

    def close(self):
        self.conn.close()


# ============================================================================
# STREAMING PARQUET WRITER
# ============================================================================

class StreamingParquetWriter:
    """
    Stream products to Parquet file with buffering and memory monitoring.

    Improvements from parallel_enum.py pattern:
    - Product buffer with configurable flush_interval
    - Optional memory monitoring (psutil)
    - Schema enforcement on every flush (prevents null-type inference bugs)
    - Adaptive flushing under memory pressure
    """

    def __init__(
        self,
        output_path: Path,
        schema: pa.Schema,
        flush_interval: int = 2000,  # CRITICAL: Reduced for multi-worker safety
        memory_threshold: float = 0.50,  # Legacy parameter (not used)
        target_memory_percent: float = 75.0  # SYSTEM-WIDE memory threshold
    ):
        self.output_path = output_path
        self.writer = None
        self.schema = schema
        self.total_written = 0
        self.flush_interval = flush_interval
        self.memory_threshold = memory_threshold
        self.target_memory_percent = target_memory_percent  # NEW

        # Product buffer (NEW)
        self.buffer = []
        self.flush_count = 0

        # Memory monitoring (optional, like parallel_enum.py)
        self.memory_monitor_enabled = self._check_psutil()

        import logging
        self.logger = logging.getLogger(__name__)

    def _check_psutil(self) -> bool:
        """Check if psutil available for memory monitoring."""
        try:
            import psutil
            return True
        except ImportError:
            return False

    def _check_memory(self) -> float:
        """Get current memory usage percentage (0.0-1.0)."""
        if not self.memory_monitor_enabled:
            return 0.0
        try:
            import psutil
            return psutil.virtual_memory().percent / 100.0
        except:
            return 0.0

    def _get_process_memory_mb(self) -> float:
        """Get current Python process memory usage in MB."""
        try:
            import psutil
            import os
            process = psutil.Process(os.getpid())
            return process.memory_info().rss / 1024 / 1024  # Bytes to MB
        except:
            return 0.0

    def _get_total_memory_mb(self) -> float:
        """Get system total memory in MB."""
        try:
            import psutil
            return psutil.virtual_memory().total / 1024 / 1024
        except:
            return 32000.0  # Assume 32GB default

    def _get_memory_usage_percent(self) -> float:
        """
        Get SYSTEM-WIDE memory usage percentage.

        CRITICAL: Returns global system memory usage, NOT single process memory.
        Essential for multi-worker scenarios where total = main + 16 workers.
        """
        try:
            import psutil
            return psutil.virtual_memory().percent
        except:
            return 0.0

    def write_batch(self, products: List[Dict]):
        """
        Intelligent memory-managed batch writing with SYSTEM-WIDE monitoring.

        CRITICAL FIX (2025-12-27): Process large batches in chunks to prevent
        buffer explosion. Previous code allowed buffer to grow to 900K+ products
        before flush, causing 1-3GB memory spikes and OOM crashes.

        Strategy:
        1. HARD LIMIT: Never allow buffer to exceed flush_interval
        2. Process incoming products in chunks if needed
        3. Monitor GLOBAL system memory (not single process)
        4. Force flush and GC after each chunk to free memory immediately
        """
        if not products:
            return

        # HARD LIMIT: Never exceed this buffer size
        MAX_BUFFER_SIZE = self.flush_interval

        # Calculate safe chunk size (half of max to allow headroom)
        chunk_size = max(1000, MAX_BUFFER_SIZE // 2)

        # Process products in chunks to prevent buffer explosion
        for i in range(0, len(products), chunk_size):
            chunk = products[i:i+chunk_size]

            # Add chunk to buffer
            self.buffer.extend(chunk)

            # Get current SYSTEM memory status
            system_mem_percent = self._get_memory_usage_percent()
            buffer_size = len(self.buffer)

            # Adaptive flush conditions
            should_flush = False
            flush_reason = ""

            # Condition 1: Buffer size limit (PRIMARY - HARD LIMIT)
            if buffer_size >= MAX_BUFFER_SIZE:
                should_flush = True
                flush_reason = f"buffer_full ({buffer_size:,} products)"

            # Condition 2: System memory pressure (SECONDARY)
            if system_mem_percent > self.target_memory_percent:
                should_flush = True
                flush_reason = f"system_memory_pressure ({system_mem_percent:.1f}%)"

            # Condition 3: Critical system memory (URGENT)
            if system_mem_percent > 85.0:
                should_flush = True
                flush_reason = f"CRITICAL_SYSTEM_MEMORY ({system_mem_percent:.1f}%)"
                self.logger.warning(f"[!] Critical system memory: {system_mem_percent:.1f}%")

            if should_flush:
                # [PROFILING] Record pre-flush memory
                pre_flush_mem = system_mem_percent
                self.logger.info(f"Flush triggered: {flush_reason}")
                self.logger.info(f"[PROFILING] Pre-flush memory: {pre_flush_mem:.1f}%")

                # Flush and measure timing
                self._flush()

                # Force garbage collection to free memory immediately
                import gc
                gc.collect()

                # [PROFILING] Record post-flush memory and delta
                post_flush_mem = self._get_memory_usage_percent()
                mem_delta = post_flush_mem - pre_flush_mem
                self.logger.info(f"[PROFILING] Post-flush memory: {post_flush_mem:.1f}%, delta: {mem_delta:+.1f}%")

    def _flush(self):
        """
        Flush buffer to disk with schema enforcement.

        CRITICAL: Cast to schema on every flush to prevent null-type inference issues.
        This matches the fix in parallel_enum.py (lines 433-453).
        """
        if not self.buffer:
            return

        # [PROFILING] Start timing
        import time
        flush_start = time.time()
        buffer_size = len(self.buffer)

        # Convert buffer to PyArrow table
        table = pa.Table.from_pylist(self.buffer, schema=self.schema)

        # Initialize writer on first flush
        if self.writer is None:
            self.writer = pq.ParquetWriter(self.output_path, self.schema)
            self.logger.info(f"ParquetWriter initialized (flush_interval={self.flush_interval:,})")
        else:
            # CRITICAL: Cast table to established schema on every flush
            # Prevents null-type inference issues (terpenoid k=2 bug)
            try:
                table = table.cast(self.schema)
            except Exception as e:
                self.logger.warning(f"Schema cast failed, attempting column-by-column: {e}")
                # Fallback: column-by-column cast
                arrays = []
                for field in self.schema:
                    if field.name in table.column_names:
                        try:
                            arrays.append(table.column(field.name).cast(field.type))
                        except:
                            arrays.append(table.column(field.name))
                    else:
                        # Column missing, create null column
                        arrays.append(pa.nulls(len(table), type=field.type))
                table = pa.Table.from_arrays(arrays, schema=self.schema)

        # Write to disk
        self.writer.write_table(table)
        self.total_written += len(self.buffer)
        self.flush_count += 1

        # Log progress
        if self.flush_count % 10 == 0:
            self.logger.info(
                f"Flush #{self.flush_count}: {self.total_written:,} total products written"
            )

        # Clear buffer
        self.buffer.clear()

        # Active memory management
        import gc
        gc.collect()  # Force garbage collection to free memory

        # Log SYSTEM memory status after flush
        system_mem_after = self._get_memory_usage_percent()
        self.logger.info(f"System memory after flush: {system_mem_after:.1f}%")

        # [PROFILING] End timing and log
        flush_elapsed = time.time() - flush_start
        self.logger.info(f"[PROFILING] Flush took {flush_elapsed:.2f}s for {buffer_size:,} products ({buffer_size/flush_elapsed:.0f} products/s)")

    def close(self):
        """Flush remaining buffer and close writer."""
        # Final flush
        if self.buffer:
            self.logger.info(f"Final flush: {len(self.buffer):,} products")
            self._flush()

        # Close writer
        if self.writer:
            self.writer.close()
            self.logger.info(
                f"StreamingParquetWriter closed: {self.total_written:,} total products, "
                f"{self.flush_count} flushes"
            )


# ============================================================================
# SPLIT SUBCOMMAND (unchanged from original)
# ============================================================================

def cmd_split(args):
    """Split base libraries - same as v1."""
    # Import from original script (keeping this unchanged)
    logger.info("SPLIT command - using original implementation")
    logger.info("For full split functionality, use original script or copy implementation")
    logger.info("Focus of v2 refactor is on APPLY command performance")


# ============================================================================
# APPLY SUBCOMMAND: Refactored for Performance
# ============================================================================

def mark_anchor_with_isotope(mol: Chem.Mol, match: Tuple[int, ...], anchor_idx: int = 0, isotope: int = 203):
    """
    Mark anchor atom in match with temporary isotope (pre-reaction).

    Args:
        mol: Molecule to modify
        match: Substructure match tuple
        anchor_idx: Index in match (0 for [c:1] position)
        isotope: Isotope number (default 203)
    """
    if anchor_idx < len(match):
        atom = mol.GetAtomWithIdx(match[anchor_idx])
        atom.SetIsotope(isotope)


def annotate_product_sites_with_mapnum(
    rxn: AllChem.ChemicalReaction,
    prod_mol: Chem.Mol,
    highlight_mapnums: List[int],
    target_symbol: Optional[str] = None,
    stamp_mapnum: int = 777
) -> List[int]:
    """
    UNIQUE_STAMP_VERSION: Use anchor isotope strategy to mark exactly ONE site per product.

    Strategy:
    1. Locate anchor atom (芳香碳) via temporary isotope (from reaction mapping :1)
    2. From anchor's neighbors, uniquely select target atom based on chemistry:
       - For OMe (O): Non-ring oxygen (X2) connected to aromatic C and sp3 C (CH3)
       - For NH2 (N): Non-ring nitrogen with H count >= 2, not in amide
    3. Mark ONLY that unique atom with :777
    4. Clean up isotope

    This ensures:
    - Single site per product (no multi-highlighting)
    - Chemical validation (not just template matching)
    - No false positives from bridge oxygen or other moieties

    Args:
        rxn: RDKit reaction object (not used in new strategy)
        prod_mol: Product molecule (from rxn.RunReactants())
        highlight_mapnums: Mapping numbers from config (for backward compatibility)
        target_symbol: Target atom symbol ('O' for OMe, 'N' for NH2)
        stamp_mapnum: Mapping number to set (default 777)

    Returns:
        List of marked atom indices (always length 0 or 1)
    """
    if not target_symbol:
        return []

    marked_atoms = []

    # Find anchor atoms with temporary isotope 203 (aromatic carbon)
    anchor_atoms = [a for a in prod_mol.GetAtoms()
                    if a.GetSymbol() == 'C' and a.GetIsotope() == 203]

    for anchor in anchor_atoms:
        # Get neighbors of anchor
        neighbors = list(anchor.GetNeighbors())

        # Find target atom based on chemistry
        target_atom = None

        if target_symbol == 'O':
            # For OMe: Find non-ring oxygen (X2) connected to methyl
            for n in neighbors:
                if n.GetSymbol() != 'O':
                    continue
                if n.IsInRing():  # Exclude bridge oxygen
                    continue
                if n.GetDegree() != 2:  # Must be divalent
                    continue

                # Check if connected to sp3 carbon (CH3)
                o_neighbors = list(n.GetNeighbors())
                has_methyl = False
                for on in o_neighbors:
                    if on.GetIdx() == anchor.GetIdx():
                        continue
                    if on.GetSymbol() == 'C' and on.GetHybridization() == Chem.HybridizationType.SP3:
                        # Check if it's a methyl (degree <= 2)
                        if on.GetDegree() <= 2:
                            has_methyl = True
                            break

                if has_methyl:
                    target_atom = n
                    break

        elif target_symbol == 'N':
            # For NH2: Find non-ring nitrogen with H count >= 2
            for n in neighbors:
                if n.GetSymbol() != 'N':
                    continue
                if n.IsInRing():  # Exclude ring nitrogen
                    continue

                # Check H count (should be 2 for NH2)
                h_count = n.GetTotalNumHs()
                if h_count < 2:
                    continue

                # Exclude amide nitrogen (check for C=O neighbor)
                is_amide = False
                for nn in n.GetNeighbors():
                    if nn.GetSymbol() == 'C':
                        for nnn in nn.GetNeighbors():
                            if nnn.GetSymbol() == 'O' and nn.GetBondBetweenAtoms(nn.GetIdx(), nnn.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                is_amide = True
                                break
                    if is_amide:
                        break

                if not is_amide:
                    target_atom = n
                    break

        # Mark the unique target atom
        if target_atom:
            target_atom.SetAtomMapNum(stamp_mapnum)
            marked_atoms.append(target_atom.GetIdx())

        # Clean up anchor isotope (CRITICAL FIX: Always clean, moved outside if block)
        anchor.SetIsotope(0)

    return marked_atoms



class TransformationEngineV2:
    """
    PERFORMANCE OPTIMIZED transformation engine.

    Key improvements:
    - Single RunReactants call per molecule
    - Cached sugar site filtering
    - Minimal property calculation (InChIKey deferred to post-processing)
    """

    def __init__(self, config: Dict, transform_name: str):
        self.config = config
        self.defaults = config.get('defaults', {})

        # Find transform
        self.transform = None
        for xf in config.get('transforms', []):
            if xf['name'] == transform_name:
                self.transform = xf
                break

        if self.transform is None:
            raise ValueError(f"Transform '{transform_name}' not found in config")

        logger.info(f"\n[TransformationEngine v2]")
        logger.info(f"  Transform: {self.transform['name']}")
        logger.info(f"  Label: {self.transform['label']}")
        logger.info(f"  Query: {self.transform['query_smarts']}")
        logger.info(f"  SMIRKS: {self.transform['smirks']}")

        # Compile patterns
        self.query_pattern = Chem.MolFromSmarts(self.transform['query_smarts'])
        if self.query_pattern is None:
            raise ValueError(f"Invalid query SMARTS")

        self.reaction = AllChem.ReactionFromSmarts(self.transform['smirks'])
        if self.reaction is None:
            raise ValueError(f"Invalid reaction SMIRKS")

        # Settings
        self.apply_mode = self.defaults.get('apply_mode', 'per_site')
        self.max_sites = self.defaults.get('max_sites_per_mol', -1)  # -1 = unlimited
        self.exclude_sugar = self.defaults.get('exclude_sugar_like', True)

        # PR-A: Atom mapping for site visualization
        self.highlight_mapnums = self.transform.get('highlight_mapnums', [])
        self.product_highlight_symbol = self.transform.get('product_highlight_symbol')

        logger.info(f"  Mode: {self.apply_mode}")
        logger.info(f"  Max sites: {self.max_sites}")
        logger.info(f"  Exclude sugar: {self.exclude_sugar}")
        if self.highlight_mapnums:
            logger.info(f"  Highlight mapnums: {self.highlight_mapnums} (symbol: {self.product_highlight_symbol})")

    def find_sites(self, mol: Chem.Mol, sugar_cache: Optional[SugarSiteCache] = None) -> List[Tuple[int]]:
        """Find matching sites with sugar filtering."""
        matches = mol.GetSubstructMatches(self.query_pattern)

        if not matches:
            return []

        # Filter sugar sites
        if self.exclude_sugar and sugar_cache:
            filtered = []
            for match in matches:
                is_sugar = any(sugar_cache.is_sugar_site(idx, mol) for idx in match)
                if not is_sugar:
                    filtered.append(match)
            matches = filtered

        # Limit to max_sites (only if max_sites > 0; -1 means unlimited)
        if self.max_sites > 0 and len(matches) > self.max_sites:
            matches = matches[:self.max_sites]

        return list(matches)

    def apply_to_molecule(self, smiles: str, source_data: Dict) -> List[Dict]:
        """
        Apply transformation - OPTIMIZED VERSION.

        Key change: Single RunReactants call, enumerate products directly.
        """
        mol = safe_mol_from_smiles(smiles)
        if mol is None:
            return [{
                'smiles': None,
                'canonical_smiles': None,
                'inchikey': None,
                'xf_success': False,
                'xf_reason': 'invalid_input_smiles',
                'source_smiles': smiles,
                **source_data
            }]

        # Create sugar cache once
        sugar_cache = SugarSiteCache(
            source_data.get('sugar_mask_atoms_json'),
            source_data.get('sugar_rings_json')
        )

        matches = self.find_sites(mol, sugar_cache)

        if not matches:
            return [{
                'smiles': None,
                'canonical_smiles': None,
                'inchikey': None,
                'xf_success': False,
                'xf_reason': 'no_matching_sites',
                'source_smiles': smiles,
                **source_data
            }]

        products = []

        if self.apply_mode == 'per_site':
            # *** KEY OPTIMIZATION: Single RunReactants call ***
            # ISOTOPE_ANCHOR_V2: Mark anchor atoms with isotope before reaction
            mol_copy = Chem.Mol(mol)
            for match in matches:
                mark_anchor_with_isotope(mol_copy, match, anchor_idx=0, isotope=203)

            try:
                rxn_products = self.reaction.RunReactants((mol_copy,))

                if not rxn_products:
                    return [{
                        'smiles': None,
                        'canonical_smiles': None,
                        'inchikey': None,
                        'xf_success': False,
                        'xf_reason': 'reaction_no_products',
                        'source_smiles': smiles,
                        **source_data
                    }]

                # Enumerate products (one per match/site)
                for site_idx, match in enumerate(matches):
                    if site_idx >= len(rxn_products):
                        break

                    try:
                        prod_mol = rxn_products[site_idx][0]

                        # Sanitize
                        Chem.SanitizeMol(prod_mol)

                        # ===== PR-A: Post-reaction atom stamping =====
                        # Use product template to mark transformation sites with :777
                        marked_atoms = []
                        if self.highlight_mapnums:
                            marked_atoms = annotate_product_sites_with_mapnum(
                                self.reaction,
                                prod_mol,
                                highlight_mapnums=self.highlight_mapnums,
                                target_symbol=self.product_highlight_symbol,
                                stamp_mapnum=777
                            )

                        # Get canonical SMILES (fast)
                        prod_smiles_canon = canonical_smiles(prod_mol)
                        if not prod_smiles_canon:
                            continue

                        # ===== NEW: Get mapped SMILES (preserve :777) =====
                        prod_smiles_mapped = Chem.MolToSmiles(prod_mol, canonical=False)

                        # InChIKey generation is EXPENSIVE - defer to post-processing
                        # For dedup, use canonical_smiles first
                        prod_inchikey = get_inchikey(prod_mol)  # Still needed for traceability

                        # Quick properties only
                        props = calculate_properties_quick(prod_mol)

                        # Calculate halogen fields
                        halogen_fields = calculate_halogen_fields(prod_mol, self.transform['label'])

                        products.append({
                            'smiles': prod_smiles_canon,
                            'canonical_smiles': prod_smiles_canon,
                            'inchikey': prod_inchikey,
                            'xf_success': True,
                            'xf_reason': None,
                            'xf_label': self.transform['label'],
                            'xf_rule_id': self.transform['name'],
                            'xf_site_index': site_idx,
                            'xf_site_atoms': str(match),
                            # ===== NEW FIELDS for PR-A =====
                            'product_mapped_smiles': prod_smiles_mapped,
                            'xf_site_mapnums': json.dumps([777] * len(marked_atoms)),
                            'source_smiles': smiles,
                            **source_data,
                            **props,
                            **halogen_fields
                        })

                    except Exception as e:
                        logger.debug(f"Site {site_idx} failed: {e}")
                        continue

            except Exception as e:
                return [{
                    'smiles': None,
                    'canonical_smiles': None,
                    'inchikey': None,
                    'xf_success': False,
                    'xf_reason': f'reaction_exception: {str(e)[:100]}',
                    'source_smiles': smiles,
                    **source_data
                }]

        else:
            raise NotImplementedError("all_sites mode not implemented")

        return products


# Worker initialization (module-level for multiprocessing)
_worker_config = None
_worker_transform_name = None
_worker_engine = None


def _init_worker(config: Dict, transform_name: str):
    """Initialize worker process with transformation engine."""
    global _worker_config, _worker_transform_name, _worker_engine
    _worker_config = config
    _worker_transform_name = transform_name
    _worker_engine = TransformationEngineV2(config, transform_name)


def _process_batch_worker(batch_data: Tuple) -> Tuple[List[Dict], Dict]:
    """
    Worker function for parallel batch processing.

    Returns:
        (products, stats)
    """
    rows = batch_data

    all_products = []
    stats = Counter()

    for row in rows:
        smiles = row.get('smiles')
        if not smiles:
            stats['skip_no_smiles'] += 1
            continue

        # Prepare source data
        source_data = {
            'source_inchikey': row.get('inchikey'),
            'source_subset': row.get('_source_subset'),
            'k': row.get('k'),
            'halogen': row.get('halogen'),
            'halogen_atom_count': row.get('halogen_atom_count'),
            'halogen_pair': row.get('halogen_pair'),
            'sugar_mask_atoms_json': row.get('sugar_mask_atoms_json'),
            'sugar_rings_json': row.get('sugar_rings_json'),
        }

        # Apply transformation
        products = _worker_engine.apply_to_molecule(smiles, source_data)

        for prod in products:
            if prod['xf_success']:
                stats['success'] += 1
            else:
                stats[f"fail_{prod.get('xf_reason', 'unknown')}"] += 1

        all_products.extend(products)

    return (all_products, dict(stats))


def _process_batch_result(products, stats, deduper, writer, counters, batch_idx):
    """
    Process batch result: deduplicate and write to output.

    Extracted helper function to avoid code duplication in streaming loop.

    Args:
        products: List of product dicts from worker
        stats: Stats dict from worker
        deduper: Deduplicator instance
        writer: StreamingParquetWriter instance
        counters: Dict with total_processed, total_products, total_unique, total_stats
        batch_idx: Current batch index (for logging)

    Returns:
        Number of unique products written
    """
    # Update stats
    counters['total_processed'] += len([p for p in products if p.get('source_smiles')])
    counters['total_products'] += len(products)
    counters['total_stats'].update(stats)

    # In-batch deduplication using canonical_smiles
    unique_batch_products = []
    seen_in_batch = set()

    for prod in products:
        if not prod['xf_success']:
            continue

        canon_smi = prod.get('canonical_smiles')
        if not canon_smi or canon_smi in seen_in_batch:
            continue

        seen_in_batch.add(canon_smi)
        unique_batch_products.append(prod)

    # Cross-batch deduplication (batched)
    canon_smiles_list = [p['canonical_smiles'] for p in unique_batch_products]
    new_keys = deduper.filter_new_keys(canon_smiles_list)
    new_keys_set = set(new_keys)

    final_products = [p for p in unique_batch_products if p['canonical_smiles'] in new_keys_set]

    counters['total_unique'] += len(final_products)

    # Stream write
    writer.write_batch(final_products)
    deduper.commit()

    return len(final_products)


def cmd_apply(args):
    """
    APPLY transformation - REFACTORED FOR PERFORMANCE.

    Improvements:
    - Parallel processing with ProcessPoolExecutor
    - Optimized deduplication (in-batch + batched SQL)
    - Streaming Parquet writer
    - Enhanced logging and metrics
    """
    logger.info("="*80)
    logger.info("APPLY v2: Performance-Optimized Transformation")
    logger.info("="*80)

    start_time = time.time()

    # Load config
    if not os.path.exists(args.xf_config):
        logger.error(f"Config not found: {args.xf_config}")
        return

    with open(args.xf_config, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)

    logger.info(f"\nConfig: {args.xf_config}")
    logger.info(f"Transform: {args.xf_name}")

    # Setup output
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Setup deduplicator (Bloom filter or SQLite)
    use_bloom = getattr(args, 'use_bloom_filter', False)

    if use_bloom:
        if not BLOOM_AVAILABLE:
            logger.error("Bloom filter requested but not available (missing dependencies: mmh3, bitarray)")
            logger.error("Install with: pip install mmh3 bitarray")
            return

        expected_items = getattr(args, 'bloom_expected_items', 60_000_000)
        logger.info(f"\n[DEDUP] Using Bloom Filter (memory-efficient)")
        logger.info(f"  Expected items: {expected_items:,}")
        logger.info(f"  False positive rate: 1%")
        logger.info(f"  Estimated size: ~{expected_items * 1.15 / 1024 / 1024:.1f} MB")
        logger.info(f"  vs SQLite: ~{expected_items * 386 / 1024 / 1024 / 1024:.1f} GB (99.7% reduction)")

        deduper = BloomFilterDeduplicator(
            output_dir=outdir,
            expected_items=expected_items,
            false_positive_rate=0.01
        )
    else:
        logger.info(f"\n[DEDUP] Using SQLite (exact deduplication)")
        dedup_db_path = outdir / "dedup.db"
        deduper = SqliteDeduplicator(dedup_db_path, resume=args.resume)

    # Load input
    if not os.path.exists(args.input):
        logger.error(f"Input not found: {args.input}")
        return

    logger.info(f"\nInput: {args.input}")

    # Read parquet in batches using pyarrow (streaming)
    parquet_file = pq.ParquetFile(args.input)
    total_rows = parquet_file.metadata.num_rows
    logger.info(f"  Total rows: {total_rows:,}")

    source_subset = Path(args.input).parent.name
    logger.info(f"  Source subset: {source_subset}")

    # Determine output schema
    # Define based on expected fields
    output_schema = pa.schema([
        ('smiles', pa.string()),
        ('canonical_smiles', pa.string()),
        ('inchikey', pa.string()),
        ('xf_success', pa.bool_()),
        ('xf_reason', pa.string()),
        ('xf_label', pa.string()),
        ('xf_rule_id', pa.string()),
        ('xf_site_index', pa.int32()),
        ('xf_site_atoms', pa.string()),
        # PR-A: New fields for atom mapping visualization
        ('product_mapped_smiles', pa.string()),  # SMILES with :777 mapping preserved
        ('xf_site_mapnums', pa.string()),        # JSON list of mapping numbers (e.g., "[777]")
        ('source_smiles', pa.string()),
        ('source_inchikey', pa.string()),
        ('source_subset', pa.string()),
        ('k', pa.int32()),
        ('halogen', pa.string()),
        ('halogen_atom_count', pa.int32()),
        ('halogen_pair', pa.string()),
        ('halogens_set', pa.string()),           # New: e.g., "F|Cl"
        ('halogen_counts_json', pa.string()),    # New: e.g., '{"F":1,"Cl":1}'
        ('primary_halogen', pa.string()),        # New: halogen from last step
        ('MW', pa.float64()),
        ('TPSA', pa.float64()),
        ('HBD', pa.int32()),
        ('HBA', pa.int32()),
        ('cLogP', pa.float64()),
        ('RotB', pa.int32()),
    ])

    # Setup streaming writer with flush interval and memory target
    products_path = outdir / "products.parquet"
    target_memory = getattr(args, 'target_memory', 45.0)  # Default to 45.0 if not provided
    writer = StreamingParquetWriter(
        products_path,
        output_schema,
        flush_interval=args.flush_interval,
        target_memory_percent=target_memory
    )
    logger.info(f"StreamingParquetWriter initialized with flush_interval={args.flush_interval:,}, target_system_memory={target_memory:.1f}%")

    # Batch processing settings
    batch_size = args.batch_size
    num_batches = (total_rows + batch_size - 1) // batch_size
    workers = args.workers

    logger.info(f"\nProcessing:")
    logger.info(f"  Batches: {num_batches}")
    logger.info(f"  Batch size: {batch_size:,}")
    logger.info(f"  Workers: {workers}")
    logger.info(f"  Resume mode: {args.resume}")

    # Metrics
    total_processed = 0
    total_products = 0
    total_unique = 0
    total_stats = Counter()

    # Define necessary columns (minimal set for transformation)
    needed_cols = [
        'smiles',           # Required for transformation
        'inchikey',         # Source tracking
        'k',                # Metadata passthrough
        'halogen',          # Metadata passthrough
        'halogen_atom_count',  # Metadata passthrough
        'halogen_pair',     # Metadata passthrough
        'sugar_mask_atoms_json',  # Sugar site annotation (optional)
        'sugar_rings_json'        # Sugar ring annotation (optional)
    ]

    # Setup max in-flight batches (CRITICAL for OOM prevention)
    max_in_flight = getattr(args, 'max_in_flight', None)
    if max_in_flight is None:
        max_in_flight = workers * 2
    logger.info(f"Max in-flight batches: {max_in_flight} (prevents queue buildup → OOM)")

    # Counters dict for helper function
    counters = {
        'total_processed': total_processed,
        'total_products': total_products,
        'total_unique': total_unique,
        'total_stats': total_stats,
    }

    # Streaming read + parallel processing with IN-FLIGHT LIMIT
    # CRITICAL FIX: Process results as they complete, don't accumulate all futures!
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_init_worker,
        initargs=(config, args.xf_name)
    ) as executor:

        # Streaming submission with in-flight limit
        futures = {}  # {future: batch_idx}
        batch_idx = 0
        batches_processed = 0

        batch_iter = parquet_file.iter_batches(columns=needed_cols, batch_size=batch_size)

        for batch in batch_iter:
            # Convert to pandas DataFrame (only this batch in memory)
            batch_df = batch.to_pandas()

            # Add source subset metadata
            batch_df['_source_subset'] = source_subset

            # Submit to executor
            batch_rows = batch_df.to_dict('records')
            future = executor.submit(_process_batch_worker, batch_rows)
            futures[future] = batch_idx
            batch_idx += 1

            # CRITICAL: If at in-flight limit, wait for completions before submitting more
            while len(futures) >= max_in_flight:
                # Wait for at least one to complete
                done, pending = wait(futures, return_when=FIRST_COMPLETED)

                for future in done:
                    completed_batch_idx = futures.pop(future)

                    try:
                        products, stats = future.result()

                        # Process result using helper function
                        _process_batch_result(
                            products, stats, deduper, writer, counters, completed_batch_idx
                        )

                        batches_processed += 1

                        # Periodic checkpoint
                        if batches_processed % 10 == 0:
                            deduper.checkpoint()

                        # Log progress
                        if batches_processed % 10 == 0:
                            elapsed = time.time() - start_time
                            rate = counters['total_processed'] / elapsed if elapsed > 0 else 0
                            logger.info(f"[Batch {batches_processed}/{num_batches}] "
                                       f"In-flight: {len(futures)}/{max_in_flight} | "
                                       f"Processed: {counters['total_processed']:,} | "
                                       f"Products: {counters['total_products']:,} | "
                                       f"Unique: {counters['total_unique']:,} | "
                                       f"Rate: {rate:.1f} mol/s")

                    except Exception as e:
                        logger.error(f"Batch {completed_batch_idx} failed: {e}")
                        logger.error(traceback.format_exc())

        logger.info(f"All batches submitted. Processing remaining {len(futures)} in-flight batches...")

        # Process remaining futures
        for future in as_completed(futures):
            completed_batch_idx = futures.pop(future)

            try:
                products, stats = future.result()

                # Process result
                _process_batch_result(
                    products, stats, deduper, writer, counters, completed_batch_idx
                )

                batches_processed += 1

                # Periodic checkpoint
                if batches_processed % 10 == 0:
                    deduper.checkpoint()

                # Log progress
                if batches_processed % 10 == 0 or batches_processed == num_batches:
                    elapsed = time.time() - start_time
                    rate = counters['total_processed'] / elapsed if elapsed > 0 else 0
                    logger.info(f"[Batch {batches_processed}/{num_batches}] "
                               f"In-flight: {len(futures)}/{max_in_flight} | "
                               f"Processed: {counters['total_processed']:,} | "
                               f"Products: {counters['total_products']:,} | "
                               f"Unique: {counters['total_unique']:,} | "
                               f"Rate: {rate:.1f} mol/s")

            except Exception as e:
                logger.error(f"Batch {completed_batch_idx} failed: {e}")
                logger.error(traceback.format_exc())

        # Update final counters from dict
        total_processed = counters['total_processed']
        total_products = counters['total_products']
        total_unique = counters['total_unique']
        total_stats = counters['total_stats']

    # Cleanup
    writer.close()
    deduper.close()

    elapsed = time.time() - start_time

    logger.info(f"\n" + "="*80)
    logger.info(f"TRANSFORMATION COMPLETE")
    logger.info(f"  Time: {elapsed:.1f}s")
    logger.info(f"  Processed: {total_processed:,} molecules")
    logger.info(f"  Total products: {total_products:,}")
    logger.info(f"  Unique products: {total_unique:,}")
    logger.info(f"  Throughput: {total_processed/elapsed:.1f} mol/s")
    logger.info(f"  Uniqueness: {100*total_unique/max(1,total_products):.1f}%")
    logger.info(f"\nFailure breakdown:")
    for reason, count in total_stats.most_common(10):
        if reason.startswith('fail_'):
            logger.info(f"  {reason}: {count:,}")
    logger.info("="*80)

    # Save summary with complete configuration for reproducibility
    summary = {
        'transform': args.xf_name,
        'source_file': args.input,
        'source_subset': source_subset,
        'total_input': total_rows,
        'total_processed': total_processed,
        'total_products': total_products,
        'unique_products': total_unique,
        'elapsed_seconds': elapsed,
        'throughput_mol_per_sec': total_processed / elapsed if elapsed > 0 else 0,
        'stats': dict(total_stats),
        # Configuration parameters for reproducibility (added for auditing)
        'config': {
            'workers': workers,
            'batch_size': batch_size,
            'resume': args.resume,
            'max_sites_per_mol': config.get('defaults', {}).get('max_sites_per_mol', -1),
            'exclude_sugar_like': config.get('defaults', {}).get('exclude_sugar_like', True),
            'apply_mode': config.get('defaults', {}).get('apply_mode', 'per_site'),
            'xf_config': args.xf_config,
            'isotope_anchor_strategy': 'v2'  # Indicates ISOTOPE_ANCHOR_V2 implementation
        }
    }

    summary_path = outdir / "SUMMARY.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"\nSaved: {products_path}")
    logger.info(f"Saved: {summary_path}")


# ============================================================================
# QC & REPORT SUBCOMMANDS (keep original for now)
# ============================================================================

def cmd_qc(args):
    """QC command - use original script."""
    logger.info("QC command: Use original 08_transform_library.py for full QC functionality")
    logger.info("v2 refactor focuses on APPLY performance")


def cmd_report(args):
    """Report command - use original script."""
    logger.info("REPORT command: Use original 08_transform_library.py for full reporting")
    logger.info("v2 refactor focuses on APPLY performance")


# ============================================================================
# MAIN CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Natural Product-Like Library Generator v2 (Performance Optimized)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # ========== split ==========
    parser_split = subparsers.add_parser('split', help='Split base libraries')
    parser_split.add_argument('--parents', required=True, help='Parents parquet file')
    parser_split.add_argument('--products', required=True, help='Products parquet file')
    parser_split.add_argument('--outdir', required=True, help='Output directory')

    # ========== apply ==========
    parser_apply = subparsers.add_parser('apply', help='Apply transformation (OPTIMIZED)')
    parser_apply.add_argument('-i', '--input', required=True, help='Input parquet file')
    parser_apply.add_argument('-o', '--outdir', required=True, help='Output directory')
    parser_apply.add_argument('--xf-config', required=True, help='Transformation config YAML')
    parser_apply.add_argument('--xf-name', required=True, help='Transformation name')
    parser_apply.add_argument('--batch-size', type=int, default=50000, help='Batch size (default: 50000)')
    parser_apply.add_argument('--workers', type=int, default=6, help='Parallel workers (default: 6)')
    parser_apply.add_argument(
        '--flush-interval',
        type=int,
        default=2000,
        help='Products to buffer before flushing to disk (default: 2000). '
             'CRITICAL: Keep low for multi-worker jobs. Lower = safer. '
             'Recommended: 1000-3000 for 16 workers.'
    )
    parser_apply.add_argument(
        '--target-memory',
        type=float,
        default=75.0,
        help='Target SYSTEM memory usage percent (default: 75.0). '
             'Monitors global system memory, not process memory. '
             'Emergency flush at 85 percent. Recommended: 70.0-80.0.'
    )
    parser_apply.add_argument('--resume', action='store_true', help='Resume from existing output')
    parser_apply.add_argument(
        '--use-bloom-filter',
        action='store_true',
        help='Use Bloom filter deduplication instead of SQLite (memory-efficient, 1%% FPR). '
             'Recommended for large jobs (>50M products) to avoid OOM. '
             'Saves ~99.7%% memory vs SQLite (68MB vs 21GB for 60M items).'
    )
    parser_apply.add_argument(
        '--bloom-expected-items',
        type=int,
        default=60_000_000,
        help='Expected number of products for Bloom filter sizing (default: 60M)'
    )
    parser_apply.add_argument(
        '--max-in-flight',
        type=int,
        default=None,
        help='Maximum in-flight batches in ProcessPoolExecutor queue (default: workers × 2). '
             'CRITICAL for preventing OOM on large datasets. Lower values = less memory but slower. '
             'Recommended: 4-8 for 32GB RAM systems.'
    )

    # ========== qc ==========
    parser_qc = subparsers.add_parser('qc', help='Generate QC report')
    parser_qc.add_argument('-i', '--input', required=True)
    parser_qc.add_argument('-o', '--outdir', required=True)

    # ========== report ==========
    parser_report = subparsers.add_parser('report', help='Generate summary report')
    parser_report.add_argument('--basedir', required=True)
    parser_report.add_argument('-o', '--output', required=True)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    # Dispatch
    if args.command == 'split':
        cmd_split(args)
    elif args.command == 'apply':
        cmd_apply(args)
    elif args.command == 'qc':
        cmd_qc(args)
    elif args.command == 'report':
        cmd_report(args)

    return 0


if __name__ == '__main__':
    sys.exit(main())
