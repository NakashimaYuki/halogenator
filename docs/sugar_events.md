# Sugar Events Pipeline Documentation

## Overview

The Sugar Events Pipeline provides comprehensive tracking and filtering of halogenation sites near sugar moieties in polyphenol molecules, particularly flavonoid glycosides. This system implements a three-component observability architecture that tracks different types of sugar-related filtering events.

## Architecture

### Three-Component Event Classification

1. **`sugar_mask_filtered`**: Sites filtered during initial site identification due to being within the computed sugar mask
2. **`sugar_proximity_filtered`**: Sites filtered by the proximity guard that lie within the configured radius of sugar mask atoms
3. **`post_guard_blocked`**: Sites filtered by post-enumeration guards due to sugar-related constraints

### Core Components

#### 1. Sugar Mask Computation
- **Location**: `src/halogenator/sugar_mask.py`
- **Purpose**: Identifies sugar ring atoms and associated functional groups
- **Key Function**: `get_sugar_mask_with_full_status()`
- **Caching**: Hybrid WeakKeyDictionary + SMILES-LRU strategy for performance

#### 2. Proximity Guard System
- **Location**: `src/halogenator/enumerate_k.py:_apply_proximity_guard()`
- **Purpose**: Filters sites within N bonds of sugar mask atoms
- **Performance**: O(1) membership checks using precomputed expanded masks
- **Configuration**: `sugar_cfg.proximity_guard_radius` (default: 0=disabled)

#### 3. QA Event Aggregation
- **Location**: `src/halogenator/enumerate_k.py:QAAggregator`
- **Purpose**: Collects and aggregates sugar events across enumeration
- **API**: Dedicated `record_paths()` method prevents semantic pollution

## Configuration

### CLI Parameters

```bash
# Enable proximity guard with 2-bond radius
python -m halogenator.cli enum --sugar.mode=heuristic --sugar.proximity-guard-radius=2

# Enable legacy key emission for backward compatibility
python -m halogenator.cli enum --sugar.emit-legacy-keys

# Enable cache monitoring for debugging
export HALOGENATOR_CACHE_MONITOR=1
python -m halogenator.cli enum --sugar.mode=heuristic
```

### EnumConfig Options

```python
from halogenator.enumerate_k import EnumConfig

# Standard configuration
config = EnumConfig(
    k_max=2,
    halogens=('F', 'Cl', 'Br'),
    rules=('R1', 'R2'),
    sugar_cfg={
        'mode': 'heuristic',
        'proximity_guard_radius': 3,  # Enable 3-bond proximity guard
        'emit_legacy_keys': True      # Emit legacy key aliases
    }
)

# Performance-optimized configuration
config = EnumConfig(
    k_max=1,
    sugar_cfg={
        'mode': 'heuristic',
        'proximity_guard_radius': 2,  # Lighter radius for performance
        'cache_mode': 'weakkey'       # Force WeakKeyDictionary if available
    }
)
```

## Usage Examples

### Basic Usage

```python
from halogenator.enumerate_k import enumerate_with_stats, EnumConfig

# Enumerate quercetin 3-O-glucoside with sugar filtering
smiles = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12.OC1OC(CO)C(O)C(O)C1O"
config = EnumConfig(
    k_max=2,
    halogens=('F', 'Cl'),
    rules=('R1', 'R2'),
    sugar_cfg={'mode': 'heuristic', 'proximity_guard_radius': 3}
)

products, qa_stats = enumerate_with_stats(smiles, config)

# Extract sugar events
qa_paths = qa_stats.get('qa_paths', {})
sugar_events = {
    'mask_filtered': qa_paths.get('sugar_mask_filtered', 0),
    'proximity_filtered': qa_paths.get('sugar_proximity_filtered', 0),
    'post_guard_blocked': qa_paths.get('post_guard_blocked', 0),
    'total': (qa_paths.get('sugar_mask_filtered', 0) +
             qa_paths.get('sugar_proximity_filtered', 0) +
             qa_paths.get('post_guard_blocked', 0))
}

print(f"Sugar events: {sugar_events}")
```

### Advanced Monitoring

```python
# Enable comprehensive monitoring
import os
os.environ['HALOGENATOR_CACHE_MONITOR'] = '1'

from halogenator.sites import get_sugar_ring_cache_stats, reset_sugar_ring_cache_stats

# Reset cache stats before test
reset_sugar_ring_cache_stats()

# Run enumeration
products, qa_stats = enumerate_with_stats(smiles, config)

# Check cache performance
cache_stats = get_sugar_ring_cache_stats()
print(f"Cache performance: {cache_stats}")
```

## Implementation Details

### Proximity Guard Algorithm

The proximity guard implements an optimized BFS-based algorithm:

1. **Mask Expansion**: Compute all atoms within radius R of sugar mask atoms once per enumeration
2. **Site Filtering**: Apply O(1) set membership checks during site identification
3. **Event Recording**: Track filtered sites in dedicated `sugar_proximity_filtered` counter

```python
def _apply_proximity_guard(sites: List[int], expanded_mask: set,
                          qa_bus: dict, rule_id: str,
                          key: str = "sugar_proximity_filtered") -> List[int]:
    """Apply proximity guard to filter sites near expanded mask atoms."""
    if not expanded_mask:
        return sites

    filtered_sites = []
    blocked_count = 0

    for site in sites:
        if site in expanded_mask:
            blocked_count += 1
        else:
            filtered_sites.append(site)

    # Record events in QA bus
    if blocked_count > 0:
        qa_bus[key] = qa_bus.get(key, 0) + blocked_count

    return filtered_sites
```

### Caching Strategy

The system implements a hybrid caching approach for sugar ring identification:

```python
# Primary: WeakKeyDictionary for automatic cleanup
try:
    _sugar_rings_cache = weakref.WeakKeyDictionary()
    _CACHE_MODE = "weakkey"
except TypeError:
    # Fallback: SMILES-based LRU cache
    _sugar_rings_cache = None
    _CACHE_MODE = "smiles-lru"
```

### Legacy Compatibility

For backward compatibility, the system optionally emits legacy key aliases:

- `sugar_post_guard_blocked` ? `post_guard_blocked`
- Controlled via `sugar_cfg.emit_legacy_keys`

## Performance Considerations

### Optimization Strategies

1. **Mask Expansion Caching**: Compute expanded masks once per enumeration, reuse across rules
2. **Sugar Ring Caching**: Cache ring identification results with hybrid WeakKey/SMILES strategy
3. **O(1) Proximity Checks**: Use set membership instead of iterative distance calculations
4. **Selective Activation**: Only enable proximity guard for glycoside samples in heuristic mode

### Benchmarking

Typical performance improvements with optimized implementation:

- **Cache Hit Rate**: 85-95% for repeated molecules
- **Proximity Guard Overhead**: <5% when properly configured
- **Memory Usage**: WeakKeyDictionary automatically manages cleanup

## Troubleshooting

### Common Issues

#### Zero Sugar Events
```
Symptoms: sugar_events = 0 for known glycosides
Diagnosis: Check proximity_guard_radius > 0 and sugar.mode = heuristic
Solution: Set --sugar.proximity-guard-radius=2 or higher
```

#### Poor Cache Performance
```
Symptoms: High miss rate in cache monitoring logs
Diagnosis: WeakKeyDictionary may not be working
Solution: Enable monitoring and check fallback to SMILES-LRU
```

#### Legacy Key Conflicts
```
Symptoms: Double counting in unit tests
Diagnosis: Both new and legacy keys being summed
Solution: Use unified calculation matching acceptance script
```

### Debug Logging

#### Proximity Guard Configuration
```
LOG: Proximity guard config: {radius=3, sample_type=glycoside, sugar_mode=heuristic, enabled=True}
```

#### Cache Performance (when HALOGENATOR_CACHE_MONITOR=1)
```
LOG: Sugar ring cache: mode=weakkey, hits=15, misses=3, hit_rate=83.3%, cache_size=8
```

#### Sugar Events Summary
```
LOG: Sugar events: mask_filtered=42, proximity_filtered=15, post_guard_blocked=3, total=60
```

## Testing

### Unit Tests

Key test files:
- `tests/test_sugar_events.py`: Core functionality tests
- `tests/test_aggregator_paths.py`: QA aggregation tests
- `tests/test_sugar_mask_basic.py`: Basic sugar masking tests

### Acceptance Tests

```bash
# Run full acceptance test suite
python scripts/run_pr1_acceptance.py

# Fast mode for CI
python scripts/run_pr1_acceptance.py --fast

# With verbose logging
python scripts/run_pr1_acceptance.py --verbose
```

### Sample Test Data

The system is tested against 11 standardized samples:
- 5 glycoside samples (should show sugar events > 0)
- 6 aglycone controls (should show minimal events)

## Integration

### With CLI

The sugar events system integrates seamlessly with the CLI:

```bash
# Standard glycoside enumeration
halogenator enum "glycoside.smi" --sugar.mode=heuristic --sugar.proximity-guard-radius=3

# Batch processing with sugar filtering
halogenator batch config.yml --sugar.mode=heuristic
```

### With External Scripts

```python
# Integration example for external tools
from halogenator.enumerate_k import enumerate_with_stats, EnumConfig

def analyze_glycoside(smiles):
    config = EnumConfig(
        k_max=2,
        sugar_cfg={'mode': 'heuristic', 'proximity_guard_radius': 3}
    )

    products, qa_stats = enumerate_with_stats(smiles, config)

    return {
        'products': len(products),
        'sugar_events': qa_stats.get('qa_paths', {})
    }
```

## API Reference

### Key Functions

#### enumerate_with_stats()
```python
def enumerate_with_stats(parent_smi: str, cfg: EnumConfig) -> Tuple[List[Dict], Dict]:
    """
    Enumerate products with comprehensive QA statistics.

    Returns:
        (products_list, qa_stats_dict) with sugar events in qa_stats['qa_paths']
    """
```

#### get_sugar_ring_cache_stats()
```python
def get_sugar_ring_cache_stats() -> Dict[str, Any]:
    """
    Get current cache statistics for debugging.

    Returns:
        Dict with keys: mode, hit, miss, hit_rate, cache_size
    """
```

#### log_sugar_ring_cache_stats()
```python
def log_sugar_ring_cache_stats(logger=None, prefix="Sugar ring cache"):
    """
    Log cache statistics if HALOGENATOR_CACHE_MONITOR=1.
    """
```

### Configuration Schema

```python
sugar_cfg = {
    'mode': str,                      # 'off', 'heuristic', 'manual'
    'proximity_guard_radius': int,    # Bond distance radius (0=disabled)
    'emit_legacy_keys': bool,         # Enable legacy key aliases
    'audit': bool,                    # Enable detailed audit fields
    # ... other sugar masking options
}
```

## Version History

- **v1.0**: Initial implementation with basic sugar masking
- **v1.1**: Added proximity guard system with configurable radius
- **v1.2**: Implemented three-component observability architecture
- **v1.3**: Added hybrid caching strategy and performance optimizations
- **v1.4**: Added comprehensive monitoring and debugging capabilities

## Related Documentation

- [Sugar Masking Algorithm](sugar_mask.md)
- [CLI Reference](cli.md)
- [Performance Tuning](performance.md)
- [Testing Guide](testing.md)