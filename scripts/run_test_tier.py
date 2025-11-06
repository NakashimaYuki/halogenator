#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
Test tier runner for halogenator CI.

Organizes tests into tiers for better CI feedback and resource management:
- Tier 0 (Smoke): Fast basic functionality tests (< 2 minutes)
- Tier 1 (Unit): Comprehensive unit tests (< 5 minutes) 
- Tier 2 (Integration): End-to-end workflow tests (< 15 minutes)
- Tier 3 (Full): All tests including long-running ones (< 30 minutes)
"""

import sys
import subprocess
import time
import argparse
from pathlib import Path


# Test tier definitions
TEST_TIERS = {
    'smoke': {
        'description': 'Fast smoke tests for basic functionality',
        'timeout_minutes': 2,
        'patterns': [
            'tests.test_schema_records',
            'tests.test_parent_key_metadata',
            'tests.test_qa_json_compatibility'
        ]
    },
    'unit': {
        'description': 'Comprehensive unit tests',
        'timeout_minutes': 5,
        'patterns': [
            'tests.test_*',
            '!tests.test_subset_consistency',  # Exclude integration tests
            '!tests.test_cli_*',               # Exclude CLI integration tests
            '!tests.test_fallback_*',          # Exclude fallback integration tests
            '!tests.test_qa_report_integration'
        ]
    },
    'integration': {
        'description': 'Integration and end-to-end tests',
        'timeout_minutes': 15,
        'patterns': [
            'tests.test_subset_consistency',
            'tests.test_cli_end_to_end_consistency',
            'tests.test_cli_pivots_merge',
            'tests.test_qa_report_integration',
            'tests.test_fallback_integration'
        ]
    },
    'full': {
        'description': 'All tests including long-running ones',
        'timeout_minutes': 30,
        'patterns': ['tests']  # Run everything
    }
}


def run_command(cmd, timeout_seconds=None):
    """Run a command with optional timeout."""
    print(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            timeout=timeout_seconds,
            cwd=Path(__file__).parent.parent
        )
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        print(f"Command timed out after {timeout_seconds} seconds")
        return 124, "", "Timeout"


def expand_test_patterns(patterns):
    """Expand test patterns into specific test modules."""
    # For now, return the patterns as-is
    # In a more sophisticated implementation, this could:
    # - Discover all test files
    # - Apply include/exclude patterns
    # - Return specific test modules to run
    
    expanded = []
    excludes = set()
    
    for pattern in patterns:
        if pattern.startswith('!'):
            excludes.add(pattern[1:])
        else:
            expanded.append(pattern)
    
    # Filter out excludes
    final_patterns = []
    for pattern in expanded:
        if pattern not in excludes:
            final_patterns.append(pattern)
    
    return final_patterns


def run_test_tier(tier_name, verbose=False, dry_run=False):
    """Run tests for a specific tier."""
    if tier_name not in TEST_TIERS:
        print(f"Unknown test tier: {tier_name}")
        print(f"Available tiers: {', '.join(TEST_TIERS.keys())}")
        return 1
    
    tier = TEST_TIERS[tier_name]
    print(f"Running {tier_name} tier tests: {tier['description']}")
    print(f"Timeout: {tier['timeout_minutes']} minutes")
    
    # Expand patterns
    test_patterns = expand_test_patterns(tier['patterns'])
    
    if dry_run:
        print("Dry run mode - would run:")
        for pattern in test_patterns:
            print(f"  {pattern}")
        return 0
    
    # Run tests
    start_time = time.time()
    total_failures = 0
    
    for pattern in test_patterns:
        print(f"\n{'='*60}")
        print(f"Running test pattern: {pattern}")
        print('='*60)
        
        cmd = ['python', '-m', 'unittest']
        if verbose:
            cmd.append('-v')
        
        if pattern == 'tests':
            # Special case for full test discovery
            cmd.extend(['discover', '-s', 'tests', '-p', 'test_*.py'])
        else:
            cmd.append(pattern)
        
        timeout_seconds = tier['timeout_minutes'] * 60
        returncode, stdout, stderr = run_command(cmd, timeout_seconds)
        
        print(stdout)
        if stderr:
            print("STDERR:", stderr)
        
        if returncode != 0:
            total_failures += 1
            print(f"FAILED: {pattern} (exit code: {returncode})")
        else:
            print(f"PASSED: {pattern}")
    
    # Summary
    elapsed = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"Test tier '{tier_name}' completed in {elapsed:.1f} seconds")
    print(f"Patterns run: {len(test_patterns)}")
    print(f"Failures: {total_failures}")
    
    if total_failures > 0:
        print(f"TIER FAILED: {total_failures} test pattern(s) failed")
        return 1
    else:
        print("TIER PASSED: All test patterns passed")
        return 0


def main():
    parser = argparse.ArgumentParser(description='Run test tiers for halogenator CI')
    parser.add_argument('tier', choices=list(TEST_TIERS.keys()) + ['list'], 
                       help='Test tier to run, or "list" to show available tiers')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Run tests in verbose mode')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be run without executing')
    
    args = parser.parse_args()
    
    if args.tier == 'list':
        print("Available test tiers:")
        for tier_name, tier_info in TEST_TIERS.items():
            print(f"  {tier_name:12} - {tier_info['description']} ({tier_info['timeout_minutes']}min)")
        return 0
    
    return run_test_tier(args.tier, args.verbose, args.dry_run)


if __name__ == '__main__':
    sys.exit(main())