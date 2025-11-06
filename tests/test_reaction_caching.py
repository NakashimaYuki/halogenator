# -*- coding: ascii -*-
"""Tests for reaction building caching."""

import unittest
from unittest.mock import patch, call
from src.halogenator.enumerate_k import (
    _build_reactions_cached, enumerate_products, EnumConfig, 
    RULES_VERSION, clear_reactions_cache, get_reactions_cache_info
)


class TestReactionCaching(unittest.TestCase):
    """Test reaction building is properly cached."""
    
    def test_reactions_built_once(self):
        """Test that build_reactions is called only once when cached."""
        # Clear any existing cache
        clear_reactions_cache()
        
        with patch('src.halogenator.enumerate_k.build_reactions') as mock_build:
            # Mock return value
            mock_reactions = {
                'R3': {'F': None, 'Cl': None},
                'R4': {'F': None, 'Cl': None},
                'R5': {'F': None, 'Cl': None}
            }
            mock_build.return_value = mock_reactions
            
            # Call cached function multiple times with same version
            result1 = _build_reactions_cached(RULES_VERSION)
            result2 = _build_reactions_cached(RULES_VERSION)
            result3 = _build_reactions_cached(RULES_VERSION)
            
            # Should all return the same cached result
            self.assertEqual(result1, mock_reactions)
            self.assertEqual(result2, mock_reactions)
            self.assertEqual(result3, mock_reactions)
            
            # build_reactions should only be called once
            mock_build.assert_called_once()
    
    def test_caching_in_enumeration(self):
        """Test that caching works in actual enumeration process."""
        # Clear cache
        clear_reactions_cache()
        
        with patch('src.halogenator.enumerate_k.build_reactions') as mock_build:
            # Mock return value with minimal reactions
            mock_reactions = {
                'R3': {'F': None, 'Cl': None}, 
                'R4': {'F': None, 'Cl': None},
                'R5': {'F': None, 'Cl': None}
            }
            mock_build.return_value = mock_reactions
            
            cfg = EnumConfig(k_max=1, halogens=('F', 'Cl'))
            
            # Run enumeration (which internally calls _layer_expand -> _build_reactions_cached)
            list(enumerate_products('c1ccccc1O', cfg))
            
            # build_reactions should be called only once despite multiple layer expansions
            self.assertEqual(mock_build.call_count, 1)
    
    def test_cache_info_tracking(self):
        """Test that cache statistics are tracked correctly."""
        # Clear cache and check initial state
        clear_reactions_cache()
        cache_info = get_reactions_cache_info()
        self.assertEqual(cache_info.hits, 0)
        self.assertEqual(cache_info.misses, 0)
        
        with patch('src.halogenator.enumerate_k.build_reactions') as mock_build:
            mock_reactions = {'R3': {'F': None}, 'R4': {'F': None}, 'R5': {'F': None}}
            mock_build.return_value = mock_reactions
            
            # First call should be a miss
            _build_reactions_cached(RULES_VERSION)
            cache_info = get_reactions_cache_info()
            self.assertEqual(cache_info.hits, 0)
            self.assertEqual(cache_info.misses, 1)
            
            # Second call should be a hit
            _build_reactions_cached(RULES_VERSION)
            cache_info = get_reactions_cache_info()
            self.assertEqual(cache_info.hits, 1)
            self.assertEqual(cache_info.misses, 1)
            
            # Third call should be another hit
            _build_reactions_cached(RULES_VERSION)
            cache_info = get_reactions_cache_info()
            self.assertEqual(cache_info.hits, 2)
            self.assertEqual(cache_info.misses, 1)


if __name__ == '__main__':
    unittest.main()