# -*- coding: ascii -*-
"""Test stream_shape parameter detection efficiency."""

import unittest
from unittest.mock import patch, MagicMock
from src.halogenator.cli import cmd_enum


class TestStreamShapeDetection(unittest.TestCase):
    """Test stream_shape parameter detection efficiency."""
    
    def test_stream_shape_detection_happens_once(self):
        """Stream shape detection should happen only once, not per molecule."""
        
        # Mock config and args
        config = {
            'io': {'smiles_file': 'test.smi'},
            'halogens': ['Cl'],
            'k_max': 1
        }
        args = MagicMock()
        args.k = 1
        args.subset = 'flavonoids'
        args.stream_shape = 'v2'
        hasattr(args, 'k')  # Make sure hasattr returns True
        
        # Mock enumerate functions that support stream_shape
        mock_k1_fn = MagicMock()
        mock_k1_fn.return_value = ([], {'no_product_matches': 0, 'qa_paths': {}})
        
        mock_k_fn = MagicMock()
        mock_k_fn.return_value = ([], {'no_product_matches': 0, 'qa_paths': {}})
        
        # Mock read_smi to return test data
        test_records = [('CCO', 'ethanol'), ('CCC', 'propane')]
        
        with patch('src.halogenator.cli.read_smi', return_value=test_records):
            with patch('inspect.signature') as mock_signature:
                # Mock signature to show both functions support stream_shape
                mock_sig = MagicMock()
                mock_sig.parameters = {'parent_smi': None, 'cfg': None, 'stream_shape': None}
                mock_signature.return_value = mock_sig
                
                # Run the command
                cmd_enum(config, args, enumerate_k_fn=mock_k_fn, enumerate_k1_fn=mock_k1_fn)
                
                # Verify that inspect.signature was called exactly twice (once for each function)
                # Not 2 * number_of_molecules times
                self.assertEqual(mock_signature.call_count, 2)
                
                # Verify that k1 function was called with stream_shape parameter
                self.assertEqual(mock_k1_fn.call_count, len(test_records))
                for call in mock_k1_fn.call_args_list:
                    self.assertIn('stream_shape', call.kwargs)
                    self.assertEqual(call.kwargs['stream_shape'], 'v2')
    
    def test_stream_shape_detection_covers_k1_and_k_paths(self):
        """Both k=1 and k>1 paths should support stream_shape when available."""
        
        config = {
            'io': {'smiles_file': 'test.smi'},
            'halogens': ['Cl'],
            'k_max': 2
        }
        args = MagicMock()
        args.k = 2
        args.subset = 'flavonoids'
        args.stream_shape = 'legacy'
        
        # Mock enumerate functions
        mock_k1_fn = MagicMock()
        mock_k1_fn.return_value = ([], {'no_product_matches': 0, 'qa_paths': {}})
        
        mock_k_fn = MagicMock()
        mock_k_fn.return_value = ([], {'no_product_matches': 0, 'qa_paths': {}})
        
        test_records = [('CCO', 'ethanol')]
        
        with patch('src.halogenator.cli.read_smi', return_value=test_records):
            with patch('inspect.signature') as mock_signature:
                # Mock signature to show both functions support stream_shape
                mock_sig = MagicMock()
                mock_sig.parameters = {'parent_smi': None, 'cfg': None, 'stream_shape': None}
                mock_signature.return_value = mock_sig
                
                # Run with k>1
                cmd_enum(config, args, enumerate_k_fn=mock_k_fn, enumerate_k1_fn=mock_k1_fn)
                
                # Should have called k_fn (not k1_fn) with stream_shape
                self.assertEqual(mock_k_fn.call_count, len(test_records))
                self.assertEqual(mock_k1_fn.call_count, 0)
                
                for call in mock_k_fn.call_args_list:
                    self.assertIn('stream_shape', call.kwargs)
                    self.assertEqual(call.kwargs['stream_shape'], 'legacy')
    
    def test_stream_shape_graceful_fallback_when_not_supported(self):
        """Functions without stream_shape support should be called without the parameter."""
        
        config = {
            'io': {'smiles_file': 'test.smi'},
            'halogens': ['Cl'],
            'k_max': 1
        }
        args = MagicMock()
        args.k = 1
        args.subset = 'flavonoids'
        args.stream_shape = 'v2'
        
        # Mock enumerate function that does NOT support stream_shape
        mock_k1_fn = MagicMock()
        mock_k1_fn.return_value = ([], {'no_product_matches': 0, 'qa_paths': {}})
        
        test_records = [('CCO', 'ethanol')]
        
        with patch('src.halogenator.cli.read_smi', return_value=test_records):
            with patch('inspect.signature') as mock_signature:
                # Mock signature to show function does NOT support stream_shape
                mock_sig = MagicMock()
                mock_sig.parameters = {'parent_smi': None, 'cfg': None}  # No stream_shape
                mock_signature.return_value = mock_sig
                
                # Run the command
                cmd_enum(config, args, enumerate_k_fn=None, enumerate_k1_fn=mock_k1_fn)
                
                # Should have called k1_fn without stream_shape parameter
                self.assertEqual(mock_k1_fn.call_count, len(test_records))
                for call in mock_k1_fn.call_args_list:
                    self.assertNotIn('stream_shape', call.kwargs)


if __name__ == '__main__':
    unittest.main()