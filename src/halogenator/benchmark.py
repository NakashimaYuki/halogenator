# -*- coding: ascii -*-
"""Performance benchmark harness for halogenator operations."""

import time
import functools
import logging
from typing import Dict, Any, List, Callable, Optional, Tuple
from collections import defaultdict, namedtuple
import json
import os
from datetime import datetime

# Benchmark result structure
BenchmarkResult = namedtuple('BenchmarkResult', [
    'operation', 'duration_ms', 'success', 'result_count', 
    'memory_usage_mb', 'metadata'
])

# Module-level logger
LOG = logging.getLogger(__name__)

# Global benchmark data collection
_benchmark_data = defaultdict(list)
_current_session = None


class BenchmarkSession:
    """Context manager for benchmark sessions."""
    
    def __init__(self, session_name: str):
        self.session_name = session_name
        self.start_time = None
        self.end_time = None
        self.results = []
    
    def __enter__(self):
        global _current_session
        _current_session = self
        self.start_time = time.perf_counter()  # Use perf_counter for session timing
        LOG.info(f"Starting benchmark session: {self.session_name}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        global _current_session
        self.end_time = time.perf_counter()
        total_duration = (self.end_time - self.start_time) * 1000
        LOG.info(f"Finished benchmark session: {self.session_name} ({total_duration:.2f}ms)")
        _current_session = None
    
    def add_result(self, result: BenchmarkResult):
        """Add a benchmark result to this session."""
        self.results.append(result)
        _benchmark_data[self.session_name].append(result)


def benchmark_operation(operation_name: str, include_memory: bool = False):
    """
    Decorator to benchmark function execution time and optionally memory usage.
    
    Args:
        operation_name: Name of the operation being benchmarked
        include_memory: Whether to measure memory usage (requires psutil)
    """
    def decorator(func: Callable):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start_time = time.perf_counter()  # Use perf_counter for accurate benchmarking
            memory_before = 0
            memory_after = 0
            
            # Measure memory usage if requested
            if include_memory:
                try:
                    import psutil
                    process = psutil.Process()
                    memory_before = process.memory_info().rss / 1024 / 1024  # MB
                except ImportError:
                    LOG.warning("psutil not available for memory monitoring")
            
            try:
                result = func(*args, **kwargs)
                success = True
                
                # Try to count results if it's a list/generator
                result_count = 0
                if hasattr(result, '__len__'):
                    result_count = len(result)
                elif hasattr(result, '__iter__') and not isinstance(result, (str, bytes)):
                    # Convert generator to list to count (be careful with large datasets)
                    if hasattr(result, '__next__'):
                        result = list(result)
                        result_count = len(result)
                
            except Exception as e:
                LOG.error(f"Benchmark operation {operation_name} failed: {e}")
                result = None
                success = False
                result_count = 0
            
            end_time = time.perf_counter()
            duration_ms = (end_time - start_time) * 1000
            
            # Measure memory usage after
            if include_memory:
                try:
                    memory_after = process.memory_info().rss / 1024 / 1024  # MB
                except (ImportError, NameError):
                    pass
            
            memory_usage_mb = memory_after - memory_before if include_memory else 0
            
            # Create benchmark result
            benchmark_result = BenchmarkResult(
                operation=operation_name,
                duration_ms=duration_ms,
                success=success,
                result_count=result_count,
                memory_usage_mb=memory_usage_mb,
                metadata={
                    'function': func.__name__,
                    'args_count': len(args),
                    'kwargs_count': len(kwargs)
                }
            )
            
            # Add to current session if available
            if _current_session:
                _current_session.add_result(benchmark_result)
            else:
                _benchmark_data['default'].append(benchmark_result)
            
            LOG.debug(f"Benchmark {operation_name}: {duration_ms:.2f}ms, "
                        f"success={success}, results={result_count}")
            
            return result
        
        return wrapper
    return decorator


def get_cache_statistics() -> Dict[str, Any]:
    """Collect cache statistics from various components."""
    stats = {}
    
    try:
        from .enumerate_k import get_reaction_cache_stats
        stats['reaction_cache'] = get_reaction_cache_stats()
    except ImportError:
        pass
    
    try:
        from .sites import get_ring_labeling_cache_stats
        stats['ring_labeling_cache'] = get_ring_labeling_cache_stats()
    except ImportError:
        pass
    
    return stats


def benchmark_enumeration(smiles: str, config: Any, operation_name: str = "enumeration") -> List[Dict[str, Any]]:
    """
    Benchmark enumeration operation with detailed metrics.
    
    Args:
        smiles: Parent SMILES string
        config: Enumeration configuration
        operation_name: Name for this benchmark operation
        
    Returns:
        List of enumeration products
    """
    from .enumerate_k import enumerate_products
    
    @benchmark_operation(operation_name, include_memory=True)
    def _benchmark_enumerate():
        return list(enumerate_products(smiles, config))
    
    return _benchmark_enumerate()


def generate_benchmark_report(session_name: str = None, output_file: str = None) -> Dict[str, Any]:
    """
    Generate a comprehensive benchmark report.
    
    Args:
        session_name: Specific session to report on (default: all sessions)
        output_file: Optional file to write JSON report to
        
    Returns:
        Dictionary containing benchmark report
    """
    report = {
        'timestamp': datetime.now().isoformat(),
        'sessions': {},
        'cache_statistics': get_cache_statistics(),
        'summary': {}
    }
    
    # Collect data from specified session or all sessions
    if session_name:
        sessions_to_report = {session_name: _benchmark_data[session_name]}
    else:
        sessions_to_report = dict(_benchmark_data)
    
    total_operations = 0
    total_duration = 0
    success_count = 0
    
    for session, results in sessions_to_report.items():
        session_stats = {
            'operation_count': len(results),
            'total_duration_ms': sum(r.duration_ms for r in results),
            'average_duration_ms': 0,
            'success_rate': 0,
            'operations': {}
        }
        
        if results:
            session_stats['average_duration_ms'] = session_stats['total_duration_ms'] / len(results)
            session_stats['success_rate'] = sum(1 for r in results if r.success) / len(results)
        
        # Group by operation type
        by_operation = defaultdict(list)
        for result in results:
            by_operation[result.operation].append(result)
        
        for op_name, op_results in by_operation.items():
            op_stats = {
                'count': len(op_results),
                'total_duration_ms': sum(r.duration_ms for r in op_results),
                'average_duration_ms': sum(r.duration_ms for r in op_results) / len(op_results),
                'min_duration_ms': min(r.duration_ms for r in op_results),
                'max_duration_ms': max(r.duration_ms for r in op_results),
                'success_rate': sum(1 for r in op_results if r.success) / len(op_results),
                'total_results': sum(r.result_count for r in op_results),
                'average_results': sum(r.result_count for r in op_results) / len(op_results)
            }
            
            if any(r.memory_usage_mb for r in op_results):
                memory_values = [r.memory_usage_mb for r in op_results if r.memory_usage_mb]
                if memory_values:
                    op_stats['memory'] = {
                        'total_mb': sum(memory_values),
                        'average_mb': sum(memory_values) / len(memory_values),
                        'max_mb': max(memory_values)
                    }
            
            session_stats['operations'][op_name] = op_stats
        
        report['sessions'][session] = session_stats
        
        # Update global totals
        total_operations += session_stats['operation_count']
        total_duration += session_stats['total_duration_ms']
        success_count += sum(1 for r in results if r.success)
    
    # Generate summary
    report['summary'] = {
        'total_operations': total_operations,
        'total_duration_ms': total_duration,
        'average_duration_ms': total_duration / total_operations if total_operations else 0,
        'overall_success_rate': success_count / total_operations if total_operations else 0
    }
    
    # Write to file if requested
    if output_file:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, indent=2)
        LOG.info(f"Benchmark report written to {output_file}")
    
    return report


def clear_benchmark_data(session_name: str = None):
    """Clear benchmark data for specified session or all sessions."""
    if session_name:
        _benchmark_data.pop(session_name, None)
    else:
        _benchmark_data.clear()
    LOG.info(f"Cleared benchmark data for session: {session_name or 'all'}")


def run_standard_benchmarks(output_dir: str = "benchmarks") -> Dict[str, Any]:
    """
    Run a standard set of benchmarks for common operations.
    
    Args:
        output_dir: Directory to write benchmark reports to
        
    Returns:
        Combined benchmark report
    """
    from .enumerate_k import EnumConfig
    
    # Standard test molecules
    test_molecules = [
        ("benzene", "c1ccccc1"),
        ("toluene", "Cc1ccccc1"),
        ("phenol", "Oc1ccccc1"),
        ("quercetin", "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12")
    ]
    
    # Standard configurations
    configs = [
        ("k1_basic", EnumConfig(k_max=1, halogens=('F', 'Cl'))),
        ("k2_constrained", EnumConfig(
            k_max=2, 
            halogens=('F', 'Cl'), 
            constraints={'per_ring_quota': 2, 'min_graph_distance': 2}
        )),
        ("k3_complex", EnumConfig(
            k_max=3, 
            halogens=('F', 'Cl', 'Br'), 
            constraints={'per_ring_quota': 3, 'min_graph_distance': 1}
        ))
    ]
    
    all_reports = {}
    
    for config_name, config in configs:
        with BenchmarkSession(f"standard_benchmark_{config_name}") as session:
            for mol_name, smiles in test_molecules:
                try:
                    products = benchmark_enumeration(
                        smiles, config, f"{mol_name}_{config_name}"
                    )
                    LOG.info(f"Benchmarked {mol_name} with {config_name}: {len(products)} products")
                except Exception as e:
                    LOG.error(f"Failed to benchmark {mol_name} with {config_name}: {e}")
        
        # Generate report for this configuration
        report_file = os.path.join(output_dir, f"{config_name}_benchmark.json")
        report = generate_benchmark_report(f"standard_benchmark_{config_name}", report_file)
        all_reports[config_name] = report
    
    # Generate combined report
    combined_report_file = os.path.join(output_dir, "combined_benchmark_report.json")
    combined_report = generate_benchmark_report(output_file=combined_report_file)
    
    LOG.info(f"Standard benchmarks completed. Reports written to {output_dir}/")
    return combined_report


def run_p1_baseline_k2(output_dir: str = "benchmarks") -> Dict[str, Any]:
    """
    Run P1 baseline k=2 performance benchmark with 10 representative flavonoids.
    
    Target: Complete in 10-15 minutes on 4 cores with <=8GB RAM.
    
    Args:
        output_dir: Directory to write benchmark reports to
        
    Returns:
        P1 baseline benchmark report
    """
    import yaml
    import os
    from .enumerate_k import EnumConfig
    
    # Load baseline configuration
    config_path = os.path.join("configs", "p1-baseline.yml")
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"P1 baseline config not found: {config_path}")
        
    with open(config_path, 'r') as f:
        baseline_config = yaml.safe_load(f)
    
    # Extract molecules and configuration
    molecules = [(mol['name'], mol['smiles']) for mol in baseline_config['baseline_molecules']]
    
    # Extract warmup/repeat configuration
    warmup_runs = baseline_config.get('performance_targets', {}).get('warmup_runs', 0)
    repeat_runs = baseline_config.get('performance_targets', {}).get('repeat_runs', 1)
    
    # Create enumeration config
    enum_config = EnumConfig(
        k_max=baseline_config['k_max'],
        halogens=tuple(baseline_config['halogens']),
        constraints=baseline_config.get('constraints', {}),
        qc_cfg=baseline_config.get('qc', {}),
        std_cfg=baseline_config.get('standardization', {})
    )
    
    # Record hardware info
    hardware_info = _get_hardware_info()
    
    # Run benchmark session
    session_name = "p1_baseline_k2"
    with BenchmarkSession(session_name) as session:
        session_start_time = time.perf_counter()  # Use perf_counter for baseline timing
        
        # Track per-molecule performance
        molecule_times = []  # Only repeat runs, no warmup
        total_products = 0
        failure_count = 0
        total_runs_executed = 0  # Track for reporting
        warmup_runs_executed = 0
        
        for mol_name, smiles in molecules:
            # Warmup runs (not counted in statistics)
            for warmup_i in range(warmup_runs):
                try:
                    LOG.debug(f"P1 baseline warmup {warmup_i+1}/{warmup_runs}: {mol_name}")
                    benchmark_enumeration(
                        smiles, enum_config, f"{mol_name}_warmup_{warmup_i}"
                    )
                    warmup_runs_executed += 1
                except Exception as e:
                    LOG.warning(f"P1 baseline warmup failed for {mol_name} (run {warmup_i+1}): {e}")
            
            # Repeat runs (counted in statistics)
            for repeat_i in range(repeat_runs):
                try:
                    mol_start = time.perf_counter()  # Use perf_counter for accurate timing
                    products = benchmark_enumeration(
                        smiles, enum_config, f"{mol_name}_repeat_{repeat_i}"
                    )
                    mol_duration = time.perf_counter() - mol_start
                    
                    molecule_times.append(mol_duration * 1000)  # Convert to ms
                    total_products += len(products)
                    total_runs_executed += 1
                    
                    LOG.info(f"P1 baseline repeat {repeat_i+1}/{repeat_runs}: {mol_name} -> {len(products)} products in {mol_duration:.2f}s")
                    
                except Exception as e:
                    failure_count += 1
                    LOG.error(f"P1 baseline repeat failed for {mol_name} (run {repeat_i+1}): {e}")
        
        total_duration = time.perf_counter() - session_start_time
    
    # Calculate statistics with correct percentile computation
    if molecule_times:
        molecule_times.sort()
        n = len(molecule_times)
        
        # Correct percentile calculation
        # p50: median
        p50_time = molecule_times[n // 2]
        
        # p95: Use ceil(0.95 * n) - 1 for proper p95 calculation
        import math
        p95_index = min(math.ceil(0.95 * n) - 1, n - 1)  # Ensure within bounds
        p95_time = molecule_times[p95_index]
        
        # total_samples should be repeat_runs * num_molecules (as per spec)
        expected_samples = repeat_runs * len(molecules)
        total_samples = n  # Actual samples collected (may be less due to failures)
    else:
        p50_time = p95_time = 0
        total_samples = 0
        expected_samples = repeat_runs * len(molecules)
    
    # Generate report
    report_data = {
        'benchmark_name': 'p1-baseline-k2',
        'timestamp': time.time(),
        'config': baseline_config,
        'hardware_info': hardware_info,
        'performance': {
            'total_duration_seconds': total_duration,
            'total_duration_minutes': total_duration / 60,
            'molecules_processed': len(molecules),
            'total_products': total_products,
            'failure_count': failure_count,
            'sampling_strategy': {
                'warmup_runs': warmup_runs,
                'repeat_runs': repeat_runs,
                'warmup_runs_executed': warmup_runs_executed,
                'repeat_runs_executed': total_runs_executed,
                'expected_samples': expected_samples,
                'warmup_excluded_from_stats': True
            },
            'per_molecule_times_ms': {
                'p50': p50_time,
                'p95': p95_time,
                'min': min(molecule_times) if molecule_times else 0,
                'max': max(molecule_times) if molecule_times else 0,
                'total_samples': total_samples,
                'all_times': molecule_times  # Keep all samples for analysis (repeat runs only)
            }
        },
        'targets': {
            'max_duration_minutes': baseline_config.get('performance_targets', {}).get('max_duration_minutes', 15),
            'max_memory_gb': baseline_config.get('performance_targets', {}).get('max_memory_gb', 8),
            'target_cores': baseline_config.get('performance_targets', {}).get('target_cores', 4)
        },
        'target_compliance': {
            'duration_ok': total_duration / 60 <= baseline_config.get('performance_targets', {}).get('max_duration_minutes', 15),
            'memory_ok': _compute_memory_ok(hardware_info.get('peak_memory_gb'), baseline_config.get('performance_targets', {}).get('max_memory_gb', 8)),
            'no_failures': failure_count == 0
        }
    }
    
    # Write report to file
    os.makedirs(output_dir, exist_ok=True)
    report_file = os.path.join(output_dir, "p1-baseline-k2.json")
    
    with open(report_file, 'w') as f:
        json.dump(report_data, f, indent=2)
    
    # Try to locate and update product output directory qa_summary
    product_output_dir = _try_locate_product_output_dir(baseline_config)
    qa_summary_updated = False
    
    if product_output_dir and os.path.exists(product_output_dir):
        product_qa_summary_path = os.path.join(product_output_dir, 'qa_summary.json')
        qa_summary_updated = update_qa_summary_with_performance(product_qa_summary_path, report_data)
        if qa_summary_updated:
            LOG.info(f"Updated product qa_summary with performance: {product_qa_summary_path}")
            report_data['product_output_dir'] = product_output_dir
        else:
            LOG.warning(f"Failed to update product qa_summary: {product_qa_summary_path}")
    
    # Also update/create benchmarks qa_summary as backup/central location
    benchmarks_qa_summary_path = os.path.join(output_dir, 'qa_summary.json')
    benchmarks_updated = update_qa_summary_with_performance(benchmarks_qa_summary_path, report_data)
    
    if not qa_summary_updated:
        # Record reason for not finding product directory
        report_data['product_output_dir'] = product_output_dir or 'not_located'
        report_data['qa_summary_location'] = 'benchmarks_only' if benchmarks_updated else 'none'
        LOG.info(f"Performance data stored in benchmarks qa_summary: {benchmarks_qa_summary_path}")
    else:
        report_data['qa_summary_location'] = 'both' if benchmarks_updated else 'product_only'
    
    # Re-write report with updated metadata
    with open(report_file, 'w') as f:
        json.dump(report_data, f, indent=2)
    
    LOG.info(f"P1 baseline k=2 completed in {total_duration/60:.2f} min. Report: {report_file}")
    
    return report_data


def _get_hardware_info() -> Dict[str, Any]:
    """Get hardware information for benchmark reporting."""
    import platform
    hardware = {
        'platform': platform.platform(),
        'processor': platform.processor(),
        'python_version': platform.python_version(),
        'logical_cores': None,
        'total_memory_gb': None,
        'peak_memory_gb': None
    }
    
    try:
        import os
        hardware['logical_cores'] = os.cpu_count()
    except Exception:
        pass
    
    try:
        import psutil
        memory = psutil.virtual_memory()
        hardware['total_memory_gb'] = memory.total / (1024**3)
        
        # Track true peak memory using platform-specific methods
        peak_memory_gb, peak_memory_method = _get_peak_memory_gb()
        hardware['peak_memory_gb'] = peak_memory_gb
        hardware['peak_memory_method'] = peak_memory_method
        
    except ImportError:
        LOG.debug("psutil not available - memory tracking disabled")
    except Exception as e:
        LOG.debug(f"Failed to get memory info: {e}")
    
    return hardware


def _get_peak_memory_gb() -> Tuple[Optional[float], str]:
    """
    Get true peak memory usage in GB using platform-specific methods.
    
    Returns:
        Tuple of (peak_memory_gb, method_name) where:
        - peak_memory_gb: Peak memory in GB, or None if unavailable
        - method_name: String describing the method used
    """
    import platform
    import sys
    
    try:
        system = platform.system().lower()
        
        if system == 'windows':
            # Try psutil memory tracking on Windows - prioritize memory_full_info()
            try:
                import psutil
                process = psutil.Process()
                
                # Priority 1: memory_full_info().peak_wset (most accurate on Windows)
                try:
                    full_memory_info = process.memory_full_info()
                    if hasattr(full_memory_info, 'peak_wset'):
                        return full_memory_info.peak_wset / (1024**3), 'full_info_peak_wset'
                except Exception:
                    LOG.debug("memory_full_info().peak_wset not available")
                
                # Priority 2: memory_info().peak_wset (fallback)
                try:
                    memory_info = process.memory_info()
                    if hasattr(memory_info, 'peak_wset'):
                        return memory_info.peak_wset / (1024**3), 'info_peak_wset'
                except Exception:
                    LOG.debug("memory_info().peak_wset not available")
                
                # Priority 3: current RSS (last resort)
                try:
                    memory_info = process.memory_info()
                    return memory_info.rss / (1024**3), 'rss_fallback'
                except Exception:
                    LOG.debug("memory_info().rss not available")
                    
            except ImportError:
                LOG.debug("psutil not available on Windows")
            except Exception as e:
                LOG.debug(f"Windows memory detection failed: {e}")
                
        elif system in ['linux', 'darwin', 'unix']:
            # Try resource.getrusage on Unix-like systems
            try:
                import resource
                usage = resource.getrusage(resource.RUSAGE_SELF)
                # Note: ru_maxrss is in KB on Linux, bytes on macOS
                if system == 'darwin':  # macOS returns bytes
                    return usage.ru_maxrss / (1024**3), 'ru_maxrss'
                else:  # Linux returns KB
                    return usage.ru_maxrss / (1024**2), 'ru_maxrss'  # KB to GB
            except Exception:
                pass
        
        # Fallback: Try psutil current RSS if available
        try:
            import psutil
            process = psutil.Process()
            memory_info = process.memory_info()
            return memory_info.rss / (1024**3), 'rss_fallback'
        except Exception:
            pass
            
    except Exception as e:
        LOG.debug(f"Failed to get peak memory: {e}")
    
    # Return None if all methods fail
    return None, 'none'


def _compute_memory_ok(peak_memory_gb: Optional[float], target_memory_gb: float) -> Optional[bool]:
    """
    Compute memory compliance using three-state logic.
    
    Args:
        peak_memory_gb: Peak memory usage in GB, or None if unavailable
        target_memory_gb: Target memory limit in GB
        
    Returns:
        True if within limit, False if exceeds limit, None if undetermined
    """
    if peak_memory_gb is None:
        return None  # Cannot determine - neither pass nor fail
    
    return peak_memory_gb <= target_memory_gb


def _try_locate_product_output_dir(baseline_config: Dict[str, Any]) -> Optional[str]:
    """
    Try to locate the product output directory from baseline configuration.
    
    Args:
        baseline_config: Baseline configuration dictionary
        
    Returns:
        Path to product output directory, or None if not found
    """
    try:
        # Look for common output directory patterns in config
        io_config = baseline_config.get('io', {})
        
        # Try various output directory keys
        potential_keys = [
            'output_dir', 'output_directory', 'products_dir',
            'enumeration_output_dir', 'results_dir'
        ]
        
        for key in potential_keys:
            output_dir = io_config.get(key)
            if output_dir and os.path.exists(output_dir):
                LOG.debug(f"Located product output directory via {key}: {output_dir}")
                return output_dir
        
        # Try to infer from products_table path
        products_table = io_config.get('products_table')
        if products_table:
            products_dir = os.path.dirname(products_table)
            if os.path.exists(products_dir):
                LOG.debug(f"Inferred product output directory from products_table: {products_dir}")
                return products_dir
        
        # Try to infer from summary_csv path
        summary_csv = io_config.get('summary_csv')
        if summary_csv:
            summary_dir = os.path.dirname(summary_csv)
            if os.path.exists(summary_dir):
                LOG.debug(f"Inferred product output directory from summary_csv: {summary_dir}")
                return summary_dir
        
        # Try common default patterns
        default_patterns = [
            'data/output/p1',
            'data/output/p0', 
            'output',
            'results'
        ]
        
        for pattern in default_patterns:
            if os.path.exists(pattern):
                LOG.debug(f"Located product output directory via default pattern: {pattern}")
                return pattern
        
        LOG.debug("Could not locate product output directory")
        return None
        
    except Exception as e:
        LOG.debug(f"Error locating product output directory: {e}")
        return None


def update_qa_summary_with_performance(qa_summary_path: str, performance_data: Dict[str, Any]) -> bool:
    """
    Update qa_summary.json with performance benchmark data.
    
    Args:
        qa_summary_path: Path to qa_summary.json file
        performance_data: Performance data from benchmark run
        
    Returns:
        True if successfully updated, False otherwise
    """
    try:
        import json
        import os
        
        # Load existing qa_summary if it exists
        if os.path.exists(qa_summary_path):
            with open(qa_summary_path, 'r') as f:
                qa_data = json.load(f)
        else:
            qa_data = {}
        
        # Ensure metadata structure exists
        if 'metadata' not in qa_data:
            qa_data['metadata'] = {}
        
        # Create performance summary (without large arrays)
        performance_summary = {
            'benchmark_name': performance_data.get('benchmark_name', 'unknown'),
            'duration_seconds': performance_data.get('performance', {}).get('total_duration_seconds', 0),
            'duration_minutes': performance_data.get('performance', {}).get('total_duration_minutes', 0),
            'molecules_processed': performance_data.get('performance', {}).get('molecules_processed', 0),
            'total_products': performance_data.get('performance', {}).get('total_products', 0),
            'failure_count': performance_data.get('performance', {}).get('failure_count', 0),
            'per_molecule_times_ms': {
                'p50': performance_data.get('performance', {}).get('per_molecule_times_ms', {}).get('p50', 0),
                'p95': performance_data.get('performance', {}).get('per_molecule_times_ms', {}).get('p95', 0),
                'min': performance_data.get('performance', {}).get('per_molecule_times_ms', {}).get('min', 0),
                'max': performance_data.get('performance', {}).get('per_molecule_times_ms', {}).get('max', 0)
            },
            'hardware_info': {
                'logical_cores': performance_data.get('hardware_info', {}).get('logical_cores'),
                'total_memory_gb': performance_data.get('hardware_info', {}).get('total_memory_gb'),
                'peak_memory_gb': performance_data.get('hardware_info', {}).get('peak_memory_gb'),
                'peak_memory_method': performance_data.get('hardware_info', {}).get('peak_memory_method', 'unknown'),
                'platform': performance_data.get('hardware_info', {}).get('platform')
            },
            'targets': performance_data.get('targets', {}),
            'target_compliance': performance_data.get('target_compliance', {}),
            'timestamp': performance_data.get('timestamp')
        }
        
        # Write performance summary to metadata.performance
        qa_data['metadata']['performance'] = performance_summary
        
        # Write back to file
        os.makedirs(os.path.dirname(qa_summary_path), exist_ok=True)
        with open(qa_summary_path, 'w') as f:
            json.dump(qa_data, f, indent=2)
        
        LOG.info(f"Updated qa_summary with performance data: {qa_summary_path}")
        return True
        
    except Exception as e:
        LOG.error(f"Failed to update qa_summary with performance data: {e}")
        return False