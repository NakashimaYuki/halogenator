# -*- coding: ascii -*-
"""Multi-process parallelization for k-dimensional enumeration."""

import os
import sys
from typing import List, Dict, Any, Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import repeat
import tempfile
import sqlite3

from .enumerate_k import enumerate_products, EnumConfig
from .io_utils import write_table, read_table


def run_parallel(parents: List[str], cfg: EnumConfig, workers: int = None,
                output_path: str = None, batch_size: int = 1000) -> None:
    """
    Run k-dimensional enumeration in parallel across parent molecules.
    
    Args:
        parents: List of parent SMILES strings
        cfg: Enumeration configuration
        workers: Number of worker processes (default: CPU count)
        output_path: Output file path for products
        batch_size: Batch size for streaming writes
    """
    if workers is None:
        workers = min(len(parents), os.cpu_count() or 4)
    
    print(f"Starting parallel enumeration with {workers} workers")
    print(f"Processing {len(parents)} parent molecules")
    
    # Create output directory if needed
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Use SQLite for global deduplication if specified
    use_sqlite = cfg.pruning_cfg.get('use_sqlite_dedup', False)
    dedup_db = None
    
    if use_sqlite and output_path:
        dedup_db = _create_dedup_database(output_path)
    
    try:
        # Process parents in parallel
        all_records = []
        processed_count = 0
        
        with ProcessPoolExecutor(max_workers=workers) as executor:
            # Submit all parent enumeration tasks
            futures = {
                executor.submit(_enumerate_single_parent, parent_smi, cfg): parent_smi
                for parent_smi in parents
            }
            
            # Collect results as they complete
            for future in as_completed(futures):
                parent_smi = futures[future]
                try:
                    parent_records = future.result()
                    
                    # Apply global deduplication
                    if use_sqlite and dedup_db:
                        unique_records = _dedupe_with_sqlite(parent_records, dedup_db)
                    else:
                        unique_records = _dedupe_in_memory(parent_records, all_records)
                    
                    all_records.extend(unique_records)
                    processed_count += 1
                    
                    print(f"Completed {processed_count}/{len(parents)}: {parent_smi} "
                          f"({len(unique_records)} unique products)")
                    
                    # Batch write to avoid memory issues
                    if len(all_records) >= batch_size and output_path:
                        _append_records_to_file(all_records, output_path)
                        all_records.clear()
                        
                except Exception as e:
                    print(f"Error processing {parent_smi}: {e}")
                    continue
        
        # Write remaining records
        if all_records and output_path:
            _append_records_to_file(all_records, output_path)
            
    finally:
        if dedup_db:
            dedup_db.close()
    
    print(f"Parallel enumeration completed")


def _enumerate_single_parent(parent_smi: str, cfg: EnumConfig) -> List[Dict[str, Any]]:
    """
    Enumerate products for a single parent molecule.
    
    This function runs in a worker process, so it must be pickleable
    and cannot share state with the main process.
    """
    try:
        records = list(enumerate_products(parent_smi, cfg))
        return records
    except Exception as e:
        print(f"Error in worker process for {parent_smi}: {e}")
        return []


def _create_dedup_database(output_path: str) -> sqlite3.Connection:
    """Create SQLite database for global deduplication."""
    db_path = output_path.replace('.parquet', '_dedup.db')
    conn = sqlite3.connect(db_path)
    
    # Create deduplication table
    conn.execute('''
        CREATE TABLE IF NOT EXISTS seen_products (
            inchikey TEXT PRIMARY KEY,
            first_seen_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    conn.commit()
    
    return conn


def _dedupe_with_sqlite(records: List[Dict[str, Any]], db: sqlite3.Connection) -> List[Dict[str, Any]]:
    """Deduplicate records using SQLite database."""
    unique_records = []
    
    for record in records:
        inchikey = record.get('inchikey')
        if not inchikey:
            continue
        
        try:
            # Try to insert the key
            db.execute('INSERT INTO seen_products (inchikey) VALUES (?)', (inchikey,))
            db.commit()
            
            # If insert succeeded, this is a new product
            unique_records.append(record)
            
        except sqlite3.IntegrityError:
            # Key already exists, skip this record
            continue
    
    return unique_records


def _dedupe_in_memory(new_records: List[Dict[str, Any]], 
                     existing_records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Deduplicate records using in-memory set."""
    # Build set of existing keys
    existing_keys = {r.get('inchikey') for r in existing_records if r.get('inchikey')}
    
    # Filter new records
    unique_records = []
    for record in new_records:
        inchikey = record.get('inchikey')
        if inchikey and inchikey not in existing_keys:
            unique_records.append(record)
            existing_keys.add(inchikey)  # Add to set to avoid duplicates within new_records
    
    return unique_records


def _append_records_to_file(records: List[Dict[str, Any]], output_path: str) -> None:
    """Append records to output file."""
    if not records:
        return
    
    try:
        # Check if file exists to determine write mode
        file_exists = os.path.exists(output_path)
        
        if file_exists:
            # Read existing records and merge
            existing_records = read_table(output_path)
            all_records = existing_records + records
        else:
            all_records = records
        
        # Write combined records
        write_table(all_records, output_path)
        
    except Exception as e:
        print(f"Error writing to {output_path}: {e}")


class ParallelConfig:
    """Configuration for parallel enumeration."""
    
    def __init__(self, workers: int = None, batch_size: int = 1000,
                 use_sqlite_dedup: bool = False, temp_dir: str = None):
        self.workers = workers or min(8, os.cpu_count() or 4)
        self.batch_size = batch_size
        self.use_sqlite_dedup = use_sqlite_dedup
        self.temp_dir = temp_dir or tempfile.gettempdir()


def create_parallel_config(**kwargs) -> ParallelConfig:
    """Create parallel configuration with defaults."""
    return ParallelConfig(**kwargs)