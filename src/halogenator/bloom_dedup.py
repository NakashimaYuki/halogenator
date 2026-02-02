"""
Memory-efficient Bloom filter-based deduplication.

Replaces SQLite dedup (24GB) with Bloom filter (~115MB) for large-scale transforms.
Trade-off: 1% false positive rate (some duplicates may pass through).
"""

import mmh3  # MurmurHash3 for fast hashing
from bitarray import bitarray
import math
from typing import Set, Optional
import pickle
from pathlib import Path


class BloomFilter:
    """
    Probabilistic data structure for membership testing.

    Memory-efficient alternative to exact deduplication.
    False positive rate: ~1% at optimal parameters.
    False negative rate: 0% (guaranteed).

    Size calculation:
        For n items, p% false positive rate:
        m = -(n * ln(p)) / (ln(2)^2)
        k = (m / n) * ln(2)

    Example:
        60M items, 1% FPR:
        m = 575M bits = 71.9 MB
        k = 7 hash functions
    """

    def __init__(self, expected_items: int, false_positive_rate: float = 0.01):
        """
        Initialize Bloom filter.

        Args:
            expected_items: Expected number of unique items
            false_positive_rate: Target false positive rate (default 1%)
        """

        self.expected_items = expected_items
        self.false_positive_rate = false_positive_rate

        # Calculate optimal parameters
        self.size_bits = self._optimal_num_bits(expected_items, false_positive_rate)
        self.num_hashes = self._optimal_num_hashes(self.size_bits, expected_items)

        # Initialize bit array
        self.bit_array = bitarray(self.size_bits)
        self.bit_array.setall(0)

        # Track statistics
        self.items_added = 0
        self.estimated_duplicates = 0

    @staticmethod
    def _optimal_num_bits(n: int, p: float) -> int:
        """
        Calculate optimal number of bits.

        Formula: m = -(n * ln(p)) / (ln(2)^2)
        """
        m = -(n * math.log(p)) / (math.log(2) ** 2)
        return int(m)

    @staticmethod
    def _optimal_num_hashes(m: int, n: int) -> int:
        """
        Calculate optimal number of hash functions.

        Formula: k = (m / n) * ln(2)
        """
        k = (m / n) * math.log(2)
        return int(k)

    def _hash(self, item: str) -> list:
        """
        Generate k hash values for item using MurmurHash3.

        Uses double hashing technique:
        h_i(x) = (h1(x) + i * h2(x)) % m
        """
        # Two independent hash values
        h1 = mmh3.hash(item, 0) % self.size_bits
        h2 = mmh3.hash(item, h1) % self.size_bits

        # Generate k hashes using double hashing
        hashes = []
        for i in range(self.num_hashes):
            hash_val = (h1 + i * h2) % self.size_bits
            hashes.append(hash_val)

        return hashes

    def add(self, item: str):
        """
        Add item to Bloom filter.

        Args:
            item: String to add (e.g., SMILES, InChIKey)
        """
        for hash_val in self._hash(item):
            self.bit_array[hash_val] = 1

        self.items_added += 1

    def contains(self, item: str) -> bool:
        """
        Test if item is in the filter.

        Args:
            item: String to test

        Returns:
            True if item is POSSIBLY in set (may be false positive)
            False if item is DEFINITELY NOT in set (guaranteed)
        """
        for hash_val in self._hash(item):
            if self.bit_array[hash_val] == 0:
                return False
        return True

    def add_if_new(self, item: str) -> bool:
        """
        Add item if not already present.

        Args:
            item: String to add

        Returns:
            True if item was new (added)
            False if item was duplicate (not added, possibly false positive)
        """
        if self.contains(item):
            self.estimated_duplicates += 1
            return False
        else:
            self.add(item)
            return True

    def save(self, path: Path):
        """Save Bloom filter to disk."""
        with open(path, 'wb') as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(path: Path) -> 'BloomFilter':
        """Load Bloom filter from disk."""
        with open(path, 'rb') as f:
            return pickle.load(f)

    def get_stats(self) -> dict:
        """Get filter statistics."""
        size_bytes = self.size_bits // 8
        size_mb = size_bytes / 1024 / 1024

        # Estimate actual false positive rate
        # p ≈ (1 - e^(-kn/m))^k
        if self.items_added > 0:
            actual_fpr = (1 - math.exp(-self.num_hashes * self.items_added / self.size_bits)) ** self.num_hashes
        else:
            actual_fpr = 0

        return {
            'expected_items': self.expected_items,
            'target_fpr': self.false_positive_rate,
            'actual_fpr': actual_fpr,
            'size_bits': self.size_bits,
            'size_mb': size_mb,
            'num_hashes': self.num_hashes,
            'items_added': self.items_added,
            'estimated_duplicates': self.estimated_duplicates,
            'fill_ratio': sum(self.bit_array) / self.size_bits
        }

    def __repr__(self):
        stats = self.get_stats()
        return (f"BloomFilter(size={stats['size_mb']:.1f}MB, "
                f"items={stats['items_added']:,}, "
                f"fpr={stats['actual_fpr']:.2%}, "
                f"fill={stats['fill_ratio']:.2%})")


class BloomFilterDeduplicator:
    """
    Drop-in replacement for SqliteDeduplicator using Bloom filter.

    Compatible interface with transform pipeline but uses memory-efficient
    Bloom filter instead of SQLite database.
    """

    def __init__(self,
                 output_dir: Path,
                 expected_items: int = 60_000_000,
                 false_positive_rate: float = 0.01):
        """
        Initialize Bloom filter deduplicator.

        Args:
            output_dir: Directory to save Bloom filter
            expected_items: Expected number of products (default 60M)
            false_positive_rate: Target FPR (default 1%)
        """

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.bloom_path = self.output_dir / "dedup_bloom.pkl"

        # Initialize or load Bloom filter
        if self.bloom_path.exists():
            self.bloom = BloomFilter.load(self.bloom_path)
        else:
            self.bloom = BloomFilter(expected_items, false_positive_rate)

        # Statistics
        self.total_checked = 0
        self.duplicates_found = 0

    def filter_new_keys(self, keys: list) -> list:
        """
        Filter to only new keys (compatible with SqliteDeduplicator).

        Args:
            keys: List of deduplication keys (SMILES, InChIKey, etc.)

        Returns:
            List of keys that are new (not seen before)
        """
        new_keys = []

        for key in keys:
            self.total_checked += 1

            if self.bloom.contains(key):
                # Already seen (or false positive)
                self.duplicates_found += 1
            else:
                # New item
                self.bloom.add(key)
                new_keys.append(key)

        return new_keys

    def is_duplicate(self, key: str) -> bool:
        """
        Check if key has been seen before.

        Args:
            key: Deduplication key (SMILES, InChIKey, etc.)

        Returns:
            True if duplicate (or false positive)
            False if new item
        """
        self.total_checked += 1

        if self.bloom.contains(key):
            self.duplicates_found += 1
            return True
        else:
            self.bloom.add(key)
            return False

    def commit(self):
        """Commit (no-op for Bloom filter, kept for compatibility)."""
        pass

    def checkpoint(self):
        """
        Save Bloom filter to disk.

        Called periodically to persist state.
        Much faster than SQLite checkpoint (no WAL truncation).
        """
        self.bloom.save(self.bloom_path)

    def get_stats(self) -> dict:
        """Get deduplication statistics."""
        bloom_stats = self.bloom.get_stats()

        return {
            'total_checked': self.total_checked,
            'duplicates_found': self.duplicates_found,
            'unique_items': self.total_checked - self.duplicates_found,
            'duplication_rate': self.duplicates_found / self.total_checked if self.total_checked > 0 else 0,
            'bloom_size_mb': bloom_stats['size_mb'],
            'bloom_fill': bloom_stats['fill_ratio'],
            'estimated_fpr': bloom_stats['actual_fpr']
        }

    def close(self):
        """Final checkpoint before closing."""
        self.checkpoint()

    def __repr__(self):
        stats = self.get_stats()
        return (f"BloomFilterDeduplicator(checked={stats['total_checked']:,}, "
                f"unique={stats['unique_items']:,}, "
                f"dupes={stats['duplicates_found']:,}, "
                f"size={stats['bloom_size_mb']:.1f}MB)")


def estimate_bloom_size(num_items: int, fpr: float = 0.01) -> dict:
    """
    Estimate Bloom filter size for given parameters.

    Args:
        num_items: Expected number of items
        fpr: False positive rate

    Returns:
        Dict with size estimates
    """

    # Calculate optimal parameters
    size_bits = BloomFilter._optimal_num_bits(num_items, fpr)
    num_hashes = BloomFilter._optimal_num_hashes(size_bits, num_items)

    size_bytes = size_bits // 8
    size_mb = size_bytes / 1024 / 1024
    size_gb = size_mb / 1024

    return {
        'num_items': num_items,
        'false_positive_rate': fpr,
        'size_bits': size_bits,
        'size_bytes': size_bytes,
        'size_mb': size_mb,
        'size_gb': size_gb,
        'num_hashes': num_hashes,
        'comparison_sqlite_gb': num_items * 386 / 1024 / 1024 / 1024  # Based on observed overhead
    }


if __name__ == "__main__":
    # Test and demonstrate Bloom filter

    print("="*80)
    print("BLOOM FILTER MEMORY ESTIMATES")
    print("="*80)

    test_cases = [
        (1_000_000, 0.01),     # 1M items, 1% FPR
        (10_000_000, 0.01),    # 10M items, 1% FPR
        (60_000_000, 0.01),    # 60M items, 1% FPR (polyphenol-2X)
        (60_000_000, 0.001),   # 60M items, 0.1% FPR
        (100_000_000, 0.01),   # 100M items, 1% FPR
    ]

    for num_items, fpr in test_cases:
        est = estimate_bloom_size(num_items, fpr)

        print(f"\n{num_items:,} items @ {fpr:.1%} FPR:")
        print(f"  Bloom filter: {est['size_mb']:,.1f} MB ({est['size_gb']:.2f} GB)")
        print(f"  Hash functions: {est['num_hashes']}")
        print(f"  vs SQLite: {est['comparison_sqlite_gb']:.2f} GB")
        print(f"  Reduction: {(1 - est['size_gb'] / est['comparison_sqlite_gb']) * 100:.1f}%")

    print("\n" + "="*80)
    print("QUICK TEST")
    print("="*80)

    # Quick test
    bf = BloomFilter(expected_items=10000, false_positive_rate=0.01)

    # Add items
    test_smiles = [
        "CCO",
        "CCCO",
        "C1CCCCC1",
        "c1ccccc1",
        "CC(C)C",
    ]

    print(f"\nAdding {len(test_smiles)} unique SMILES...")
    for smi in test_smiles:
        bf.add(smi)

    print(f"Testing membership...")
    for smi in test_smiles:
        assert bf.contains(smi), f"Should contain {smi}"
    print(f"  ✓ All {len(test_smiles)} items found correctly")

    print(f"\nTesting non-membership...")
    new_items = ["CCCCO", "CC(C)(C)C", "c1ccc(O)cc1"]
    for smi in new_items:
        is_present = bf.contains(smi)
        print(f"  {smi}: {'FOUND (false positive)' if is_present else 'not found (correct)'}")

    print(f"\nFilter stats:")
    for key, value in bf.get_stats().items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value:,}")

    print("\n" + "="*80)
