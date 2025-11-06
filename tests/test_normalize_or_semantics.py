#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Specific tests for OR semantics in normalize_rules_cfg_keys.

Tests edge cases where boolean OR logic is crucial for correct behavior.
"""

import sys
from pathlib import Path

# Add the src directory to Python path for imports
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

from halogenator.cli import normalize_rules_cfg_keys


def test_boolean_or_semantics_explicit():
    """
    Test explicit boolean OR semantics when R2a and R2b both specify boolean values.

    This test checks if overlapping boolean configuration uses OR logic:
    - R2a: sp2_CH_in_C_ring=False, sp3_CH2_flavanone=True
    - R2b: sp2_CH_in_C_ring=True, sp3_CH2_flavanone=False

    Expected result with OR semantics:
    - sp2_CH_in_C_ring=True (False OR True = True)
    - sp3_CH2_flavanone=True (True OR False = True)
    """
    config_with_boolean_conflict = {
        'R2a': {
            'sp2_CH_in_C_ring': False,      # R2a says False
            'sp3_CH2_flavanone': True,      # R2a says True
            'allowed_halogens': ['F']
        },
        'R2b': {
            'sp2_CH_in_C_ring': True,       # R2b says True
            'sp3_CH2_flavanone': False,     # R2b says False
            'allowed_halogens': ['Cl']
        }
    }

    normalized = normalize_rules_cfg_keys(config_with_boolean_conflict)

    print("Input config:")
    print(f"  R2a: sp2_CH_in_C_ring={config_with_boolean_conflict['R2a']['sp2_CH_in_C_ring']}")
    print(f"  R2a: sp3_CH2_flavanone={config_with_boolean_conflict['R2a']['sp3_CH2_flavanone']}")
    print(f"  R2b: sp2_CH_in_C_ring={config_with_boolean_conflict['R2b']['sp2_CH_in_C_ring']}")
    print(f"  R2b: sp3_CH2_flavanone={config_with_boolean_conflict['R2b']['sp3_CH2_flavanone']}")

    print("\\nNormalized result:")
    print(f"  R2: sp2_CH_in_C_ring={normalized['R2']['sp2_CH_in_C_ring']}")
    print(f"  R2: sp3_CH2_flavanone={normalized['R2']['sp3_CH2_flavanone']}")

    # With OR semantics, both should be True (False OR True = True, True OR False = True)
    expected_sp2 = True   # False OR True = True
    expected_sp3 = True   # True OR False = True

    print(f"\\nExpected with OR semantics:")
    print(f"  R2: sp2_CH_in_C_ring={expected_sp2} (False OR True)")
    print(f"  R2: sp3_CH2_flavanone={expected_sp3} (True OR False)")

    # Check if current implementation uses OR semantics
    actual_sp2 = normalized['R2']['sp2_CH_in_C_ring']
    actual_sp3 = normalized['R2']['sp3_CH2_flavanone']

    print(f"\\nActual result:")
    print(f"  R2: sp2_CH_in_C_ring={actual_sp2}")
    print(f"  R2: sp3_CH2_flavanone={actual_sp3}")

    # Analysis
    or_semantics_working = (actual_sp2 == expected_sp2) and (actual_sp3 == expected_sp3)

    if or_semantics_working:
        print("\\n[PASS] OR semantics working correctly")
        return True
    else:
        print("\\n[ISSUE] OR semantics NOT working - current implementation uses 'last writer wins'")
        print("Current behavior: .update() overwrites previous values instead of OR-ing booleans")
        return False


def test_auto_enable_vs_explicit_false():
    """
    Test interaction between auto-enable logic and explicit False values.

    R2a config without explicit sp2_CH_in_C_ring should auto-enable it to True.
    But if another config explicitly sets it to False, what happens?
    """
    config = {
        'R2a': {
            # sp2_CH_in_C_ring not specified, should auto-enable to True
            'allowed_halogens': ['F']
        },
        'R2': {
            'sp2_CH_in_C_ring': False,  # Explicitly set to False
            'sp3_CH2_flavanone': False
        }
    }

    normalized = normalize_rules_cfg_keys(config)

    print("\\nTest: Auto-enable vs Explicit False")
    print("Input:")
    print(f"  R2a: sp2_CH_in_C_ring not specified (should auto-enable to True)")
    print(f"  R2: sp2_CH_in_C_ring=False (explicitly False)")

    result_sp2 = normalized['R2']['sp2_CH_in_C_ring']
    print(f"Result: sp2_CH_in_C_ring={result_sp2}")

    # The question is: should auto-enable override explicit False, or should explicit False win?
    # OR semantics would suggest: False OR True = True
    if result_sp2:
        print("[PASS] Auto-enable logic working (OR semantics: False OR True = True)")
        return True
    else:
        print("[INFO] Auto-enable did not override explicit False (explicit wins)")
        return False


def run_or_semantics_verification():
    """Run all OR semantics verification tests."""
    print("=== Verifying OR Semantics in normalize_rules_cfg_keys ===")

    result1 = test_boolean_or_semantics_explicit()
    result2 = test_auto_enable_vs_explicit_false()

    print(f"\\n=== Summary ===")
    print(f"Boolean OR semantics: {'WORKING' if result1 else 'NEEDS IMPROVEMENT'}")
    print(f"Auto-enable logic: {'WORKING' if result2 else 'EXPLICIT VALUES WIN'}")

    if result1 and result2:
        print("\\n[CONCLUSION] OR semantics are working correctly")
    elif not result1:
        print("\\n[CONCLUSION] OR semantics need to be implemented for explicit boolean conflicts")
    else:
        print("\\n[CONCLUSION] Mixed results - review implementation")

    return result1, result2


if __name__ == "__main__":
    run_or_semantics_verification()