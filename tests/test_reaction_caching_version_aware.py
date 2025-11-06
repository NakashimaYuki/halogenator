# -*- coding: ascii -*-
"""Tests for version-aware reaction caching functionality."""

import unittest
from unittest.mock import patch
from src.halogenator.enumerate_k import (
    _build_reactions_cached, 
    clear_reactions_cache, 
    get_reactions_cache_info,
    RULES_VERSION
)


class TestReactionCachingVersionAware(unittest.TestCase):
    """Test that reaction caching respects version changes and configuration."""
    
    def setUp(self):
        """Clear cache before each test."""
        clear_reactions_cache()
    
    def tearDown(self):
        """Clear cache after each test."""
        clear_reactions_cache()
    
    def test_same_version_uses_cache(self):
        """Test that repeated calls with same version use cache."""
        # Clear cache and get initial state
        clear_reactions_cache()
        initial_info = get_reactions_cache_info()
        self.assertEqual(initial_info.hits, 0)
        self.assertEqual(initial_info.misses, 0)
        
        # First call should be a cache miss
        reactions1 = _build_reactions_cached(RULES_VERSION)
        info_after_first = get_reactions_cache_info()
        self.assertEqual(info_after_first.misses, 1)
        self.assertEqual(info_after_first.hits, 0)
        
        # Second call with same version should be a cache hit
        reactions2 = _build_reactions_cached(RULES_VERSION)
        info_after_second = get_reactions_cache_info()
        self.assertEqual(info_after_second.misses, 1)
        self.assertEqual(info_after_second.hits, 1)
        
        # Results should be identical (same object due to caching)
        self.assertIs(reactions1, reactions2)
    
    def test_different_version_triggers_rebuild(self):
        """Test that different version strings trigger cache miss."""
        # First call with version 1
        reactions1 = _build_reactions_cached("1.0.0")
        info_after_first = get_reactions_cache_info()
        self.assertEqual(info_after_first.misses, 1)
        
        # Second call with version 2 should be another miss
        reactions2 = _build_reactions_cached("1.1.0")
        info_after_second = get_reactions_cache_info()
        self.assertEqual(info_after_second.misses, 2)
        self.assertEqual(info_after_second.hits, 0)
        
        # Results should be different objects
        self.assertIsNot(reactions1, reactions2)
        # But should have same structure (both are valid reaction sets)
        self.assertEqual(set(reactions1.keys()), set(reactions2.keys()))
    
    def test_cache_clear_functionality(self):
        """Test that cache clearing works correctly."""
        # Populate cache
        _build_reactions_cached(RULES_VERSION)
        info_populated = get_reactions_cache_info()
        self.assertGreater(info_populated.currsize, 0)
        
        # Clear cache
        clear_reactions_cache()
        info_cleared = get_reactions_cache_info()
        self.assertEqual(info_cleared.currsize, 0)
        self.assertEqual(info_cleared.hits, 0)
        self.assertEqual(info_cleared.misses, 0)
        
        # Next call should be miss again
        _build_reactions_cached(RULES_VERSION)
        info_after_clear = get_reactions_cache_info()
        self.assertEqual(info_after_clear.misses, 1)
    
    def test_multiple_versions_cached_simultaneously(self):
        """Test that multiple versions can be cached simultaneously."""
        # Build reactions for different versions
        reactions_v1 = _build_reactions_cached("1.0.0")
        reactions_v2 = _build_reactions_cached("1.1.0") 
        reactions_v3 = _build_reactions_cached("1.2.0")
        
        # All should be cache misses initially
        info_after_builds = get_reactions_cache_info()
        self.assertEqual(info_after_builds.misses, 3)
        self.assertEqual(info_after_builds.hits, 0)
        
        # Repeat calls should hit cache
        reactions_v1_repeat = _build_reactions_cached("1.0.0")
        reactions_v2_repeat = _build_reactions_cached("1.1.0")
        
        info_after_repeats = get_reactions_cache_info()
        self.assertEqual(info_after_repeats.misses, 3)
        self.assertEqual(info_after_repeats.hits, 2)
        
        # Same objects should be returned
        self.assertIs(reactions_v1, reactions_v1_repeat)
        self.assertIs(reactions_v2, reactions_v2_repeat)
    
    def test_cache_size_limit_respected(self):
        """Test that cache respects maxsize=4 limit."""
        # Build reactions for 5 different versions (more than cache size of 4)
        versions = ["1.0.0", "1.1.0", "1.2.0", "1.3.0", "1.4.0"]
        reactions = {}
        
        for version in versions:
            reactions[version] = _build_reactions_cached(version)
        
        # Should have 5 misses initially
        info_after_all = get_reactions_cache_info()
        self.assertEqual(info_after_all.misses, 5)
        
        # Current cache size should be at most 4 due to LRU eviction
        self.assertLessEqual(info_after_all.currsize, 4)
        
        # The first version might have been evicted, so calling it again 
        # might be a cache miss
        _build_reactions_cached(versions[0])  # Might be miss if evicted
        info_final = get_reactions_cache_info()
        
        # Total cache size should never exceed 4
        self.assertLessEqual(info_final.currsize, 4)
    
    def test_current_rules_version_works(self):
        """Test that the current RULES_VERSION constant works correctly."""
        # Should be able to build reactions with current version
        reactions = _build_reactions_cached(RULES_VERSION)
        self.assertIsNotNone(reactions)
        self.assertIsInstance(reactions, dict)
        
        # Should have expected rule types
        expected_rules = {"R1", "R3", "R4", "R5"}
        self.assertEqual(set(reactions.keys()), expected_rules)
        
        # Each rule should have halogen reactions
        for rule in expected_rules:
            self.assertIsInstance(reactions[rule], dict)
            self.assertGreater(len(reactions[rule]), 0)  # Should have some halogens
    
    def test_reactions_structure_consistent_across_versions(self):
        """Test that different versions produce consistent reaction structure."""
        v1_reactions = _build_reactions_cached("1.0.0")
        v2_reactions = _build_reactions_cached("2.0.0")
        
        # Both should have same top-level structure
        self.assertEqual(set(v1_reactions.keys()), set(v2_reactions.keys()))
        
        # Each rule should have similar structure
        for rule in v1_reactions:
            self.assertIsInstance(v1_reactions[rule], dict)
            self.assertIsInstance(v2_reactions[rule], dict)
            # Should have same halogen keys (assuming same HALOGENS constant)
            self.assertEqual(set(v1_reactions[rule].keys()), set(v2_reactions[rule].keys()))
    
    def test_cache_info_provides_useful_statistics(self):
        """Test that cache info provides useful statistics for monitoring."""
        # Start fresh
        clear_reactions_cache()
        
        # Build some reactions
        _build_reactions_cached("1.0.0")
        _build_reactions_cached("1.0.0")  # Hit
        _build_reactions_cached("1.1.0")  # Miss
        
        info = get_reactions_cache_info()
        
        # Check that all expected attributes exist
        self.assertTrue(hasattr(info, 'hits'))
        self.assertTrue(hasattr(info, 'misses'))
        self.assertTrue(hasattr(info, 'maxsize'))
        self.assertTrue(hasattr(info, 'currsize'))
        
        # Check expected values
        self.assertEqual(info.hits, 1)
        self.assertEqual(info.misses, 2)
        self.assertEqual(info.maxsize, 4)
        self.assertEqual(info.currsize, 2)
    
    @patch('src.halogenator.enumerate_k.build_reactions')
    def test_cache_reduces_build_reactions_calls(self, mock_build_reactions):
        """Test that caching actually reduces calls to build_reactions."""
        # Mock the underlying build_reactions function
        mock_build_reactions.return_value = {"R1": {}, "R3": {}, "R4": {}, "R5": {}}
        
        # Make multiple calls with same version
        _build_reactions_cached("test_version")
        _build_reactions_cached("test_version")
        _build_reactions_cached("test_version")
        
        # build_reactions should only be called once due to caching
        self.assertEqual(mock_build_reactions.call_count, 1)
        
        # Different version should trigger another call
        _build_reactions_cached("test_version_2")
        self.assertEqual(mock_build_reactions.call_count, 2)


if __name__ == '__main__':
    unittest.main()