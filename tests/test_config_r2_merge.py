#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Unit tests for R2 configuration merging to prevent regression issues.

These tests ensure that:
1. YAML-only R2 configuration values are preserved
2. YAML True + CLI unspecified preserves True values
3. YAML True + CLI explicit disable affects only target key
4. Historical R2a/R2b key normalization uses OR semantics for boolean values
5. Configuration merge order follows: defaults <- YAML <- CLI

Tests the fix for the configuration pipeline bug where:
- CLI parameters with action='store_true' defaulted to False, overriding YAML True values
- R2a/R2b normalization didn't use OR semantics for boolean merging
"""

import tempfile
import yaml
import argparse
import sys
from pathlib import Path
from unittest.mock import Mock

# Add the src directory to Python path for imports
src_path = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(src_path))

from halogenator.cli import normalize_rules_cfg_keys, deep_merge


class TestConfigurationMerge:
    """Test configuration merge pipeline components."""

    def test_yaml_only_r2_configuration_preservation(self):
        """
        Test case 1: YAML-only R2 configuration should preserve True values.

        User config with R2.sp2=true, R2.sp3ch2=true, allow_alpha_as_beta=true
        should result in effective configuration keeping all True.
        """
        yaml_config = {
            'R2': {
                'sp2_CH_in_C_ring': True,
                'sp3_CH2_flavanone': True,
                'allow_alpha_as_beta': True,
                'allowed_halogens': ['F', 'Cl']
            }
        }

        # Test normalization (should be pass-through for R2)
        normalized = normalize_rules_cfg_keys(yaml_config)

        assert 'R2' in normalized
        assert normalized['R2']['sp2_CH_in_C_ring'] == True
        assert normalized['R2']['sp3_CH2_flavanone'] == True
        assert normalized['R2']['allow_alpha_as_beta'] == True
        assert normalized['R2']['allowed_halogens'] == ['F', 'Cl']

    def test_yaml_true_cli_unspecified_preserves_true(self):
        """
        Test case 2: YAML True + CLI unspecified should preserve True values.

        When YAML has True and CLI arguments are None (unspecified),
        the merge should preserve the YAML True values.
        """
        defaults = {
            'R2': {
                'sp2_CH_in_C_ring': False,
                'sp3_CH2_flavanone': False,
                'allow_alpha_as_beta': False
            }
        }

        yaml_config = {
            'R2': {
                'sp2_CH_in_C_ring': True,
                'sp3_CH2_flavanone': True,
                'allow_alpha_as_beta': True
            }
        }

        # Simulate CLI with None values (unspecified)
        cli_overrides = {}  # Empty means no CLI overrides

        # Merge: defaults <- YAML <- CLI
        result = deep_merge(defaults, yaml_config)
        result = deep_merge(result, cli_overrides)

        assert result['R2']['sp2_CH_in_C_ring'] == True
        assert result['R2']['sp3_CH2_flavanone'] == True
        assert result['R2']['allow_alpha_as_beta'] == True

    def test_yaml_true_cli_explicit_disable_selective(self):
        """
        Test case 3: YAML True + CLI explicit disable should affect only target key.

        When YAML has True values and CLI explicitly disables one key,
        only that key should become False, others remain True.
        """
        defaults = {
            'R2': {
                'sp2_CH_in_C_ring': False,
                'sp3_CH2_flavanone': False,
                'allow_alpha_as_beta': False
            }
        }

        yaml_config = {
            'R2': {
                'sp2_CH_in_C_ring': True,
                'sp3_CH2_flavanone': True,
                'allow_alpha_as_beta': True
            }
        }

        # Simulate CLI explicit disable of only sp3_CH2_flavanone
        cli_overrides = {
            'R2': {
                'sp3_CH2_flavanone': False  # CLI explicitly set to False
            }
        }

        # Merge: defaults <- YAML <- CLI
        result = deep_merge(defaults, yaml_config)
        result = deep_merge(result, cli_overrides)

        # Only the CLI-overridden key should be False
        assert result['R2']['sp2_CH_in_C_ring'] == True      # YAML True preserved
        assert result['R2']['sp3_CH2_flavanone'] == False    # CLI override applied
        assert result['R2']['allow_alpha_as_beta'] == True   # YAML True preserved

    def test_historical_r2a_r2b_key_normalization_individual(self):
        """
        Test case 4a: Historical R2a/R2b individual key normalization.

        Legacy keys R2a and R2b should normalize to R2 with appropriate sub-settings.
        Each should enable their respective boolean flags when not explicitly provided.
        """
        # Test R2a normalization
        r2a_config = {
            'R2a': {
                'allowed_halogens': ['F', 'Cl']
            }
        }

        normalized_r2a = normalize_rules_cfg_keys(r2a_config)

        assert 'R2' in normalized_r2a
        assert 'R2a' not in normalized_r2a
        assert normalized_r2a['R2']['sp2_CH_in_C_ring'] == True  # Auto-enabled
        assert normalized_r2a['R2']['allowed_halogens'] == ['F', 'Cl']

        # Test R2b normalization
        r2b_config = {
            'R2b': {
                'allowed_halogens': ['Br', 'I'],
                'allow_alpha_as_beta': True
            }
        }

        normalized_r2b = normalize_rules_cfg_keys(r2b_config)

        assert 'R2' in normalized_r2b
        assert 'R2b' not in normalized_r2b
        assert normalized_r2b['R2']['sp3_CH2_flavanone'] == True  # Auto-enabled
        assert normalized_r2b['R2']['allow_alpha_as_beta'] == True
        assert normalized_r2b['R2']['allowed_halogens'] == ['Br', 'I']

    def test_historical_r2a_r2b_coexistence_needs_or_semantics(self):
        """
        Test case 4b: Historical R2a/R2b coexistence should use OR semantics.

        When both R2a and R2b configurations are present, boolean values
        should be OR-ed together, not overwritten.

        NOTE: This test demonstrates the current limitation and the expected behavior.
        """
        # Configuration with both R2a and R2b present
        mixed_config = {
            'R2a': {
                'sp2_CH_in_C_ring': True,
                'allowed_halogens': ['F']
            },
            'R2b': {
                'sp3_CH2_flavanone': True,
                'allow_alpha_as_beta': True,
                'allowed_halogens': ['Cl']  # This will overwrite R2a's setting currently
            }
        }

        normalized = normalize_rules_cfg_keys(mixed_config)

        # Both should be enabled (this tests the auto-enable logic)
        assert normalized['R2']['sp2_CH_in_C_ring'] == True
        assert normalized['R2']['sp3_CH2_flavanone'] == True
        assert normalized['R2']['allow_alpha_as_beta'] == True

        # LIMITATION: Current implementation overwrites allowed_halogens
        # Expected behavior would be to merge/union the lists
        # For now, we document the current behavior:
        assert normalized['R2']['allowed_halogens'] == ['Cl']  # Last writer wins

    def test_configuration_merge_order_precedence(self):
        """
        Test case 5: Configuration merge order should follow defaults <- YAML <- CLI.

        Verify that the deep_merge function preserves the correct precedence.
        """
        defaults = {
            'R2': {
                'sp2_CH_in_C_ring': False,
                'sp3_CH2_flavanone': False,
                'allow_alpha_as_beta': False,
                'allowed_halogens': ['F']
            },
            'other_setting': 'default_value'
        }

        yaml_config = {
            'R2': {
                'sp2_CH_in_C_ring': True,  # Override default
                'sp3_CH2_flavanone': True,  # Override default
                # allow_alpha_as_beta not specified, should keep default
                'allowed_halogens': ['F', 'Cl']  # Override default
            },
            'yaml_only_setting': 'yaml_value'
        }

        cli_config = {
            'R2': {
                'allow_alpha_as_beta': True,  # CLI override
                # Other R2 values not specified, should keep YAML values
            },
            'cli_only_setting': 'cli_value'
        }

        # Apply merge order: defaults <- YAML <- CLI
        step1 = deep_merge(defaults, yaml_config)
        final_result = deep_merge(step1, cli_config)

        # Check precedence
        assert final_result['R2']['sp2_CH_in_C_ring'] == True      # YAML overrides default
        assert final_result['R2']['sp3_CH2_flavanone'] == True     # YAML overrides default
        assert final_result['R2']['allow_alpha_as_beta'] == True   # CLI overrides default
        assert final_result['R2']['allowed_halogens'] == ['F', 'Cl']  # YAML overrides default

        assert final_result['other_setting'] == 'default_value'
        assert final_result['yaml_only_setting'] == 'yaml_value'
        assert final_result['cli_only_setting'] == 'cli_value'

    def test_deep_merge_preserves_user_values(self):
        """
        Test that deep_merge correctly preserves user values over defaults.

        This tests the core merge function used in the configuration pipeline.
        """
        defaults = {
            'level1': {
                'bool_false_default': False,
                'bool_true_default': True,
                'string_default': 'default',
                'nested': {
                    'deep_bool': False,
                    'deep_string': 'deep_default'
                }
            }
        }

        user_config = {
            'level1': {
                'bool_false_default': True,     # User override
                # bool_true_default not specified, keep default
                'string_default': 'user_value',  # User override
                'nested': {
                    'deep_bool': True,          # User override
                    # deep_string not specified, keep default
                },
                'user_only': 'user_added'       # User addition
            }
        }

        result = deep_merge(defaults, user_config)

        assert result['level1']['bool_false_default'] == True      # User override
        assert result['level1']['bool_true_default'] == True       # Default preserved
        assert result['level1']['string_default'] == 'user_value' # User override
        assert result['level1']['nested']['deep_bool'] == True     # User override
        assert result['level1']['nested']['deep_string'] == 'deep_default'  # Default preserved
        assert result['level1']['user_only'] == 'user_added'      # User addition


class TestNormalizeRulesCfgKeysORSemantics:
    """Test OR semantics implementation for normalize_rules_cfg_keys."""

    def test_current_implementation_limitation(self):
        """
        Document current limitation: boolean values get overwritten instead of OR-ed.

        This test demonstrates the issue that needs to be fixed for complete OR semantics.
        """
        # Scenario: R2a enables sp2, R2b enables sp3, both should result in both enabled
        config = {
            'R2': {
                'sp2_CH_in_C_ring': True,
                'sp3_CH2_flavanone': False,  # Will be updated by R2a/R2b logic
            },
            'R2a': {
                'some_r2a_setting': True
            },
            'R2b': {
                'some_r2b_setting': True
            }
        }

        normalized = normalize_rules_cfg_keys(config)

        # The auto-enable logic should make both True
        assert normalized['R2']['sp2_CH_in_C_ring'] == True   # From base + R2a auto-enable
        assert normalized['R2']['sp3_CH2_flavanone'] == True  # From R2b auto-enable


def run_all_tests():
    """Run all test cases and report results."""
    test_obj = TestConfigurationMerge()
    test_obj2 = TestNormalizeRulesCfgKeysORSemantics()

    tests = [
        ("YAML-only R2 configuration preservation", test_obj.test_yaml_only_r2_configuration_preservation),
        ("YAML True + CLI unspecified preserves True", test_obj.test_yaml_true_cli_unspecified_preserves_true),
        ("YAML True + CLI explicit disable selective", test_obj.test_yaml_true_cli_explicit_disable_selective),
        ("Historical R2a/R2b individual normalization", test_obj.test_historical_r2a_r2b_key_normalization_individual),
        ("Historical R2a/R2b coexistence OR semantics", test_obj.test_historical_r2a_r2b_coexistence_needs_or_semantics),
        ("Configuration merge order precedence", test_obj.test_configuration_merge_order_precedence),
        ("Deep merge preserves user values", test_obj.test_deep_merge_preserves_user_values),
        ("Current implementation OR semantics", test_obj2.test_current_implementation_limitation),
    ]

    passed = 0
    failed = 0

    for test_name, test_func in tests:
        try:
            test_func()
            print(f'[PASS] {test_name}')
            passed += 1
        except Exception as e:
            print(f'[FAIL] {test_name} - {str(e)}')
            failed += 1

    print(f'\n--- Test Results ---')
    print(f'Passed: {passed}')
    print(f'Failed: {failed}')
    print(f'Total: {passed + failed}')

    return failed == 0


if __name__ == "__main__":
    run_all_tests()