# -*- coding: ascii -*-
"""
Test sugar audit fields presence and type in writer output.

This test ensures that the 4 sugar audit fields are properly written to
output files with correct types: JSON for lists, numeric for counts.
"""

import json
import tempfile
import unittest
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from halogenator.enumerate_k import enumerate_with_stats, EnumConfig
from halogenator.io_utils import write_table, read_table
from halogenator.chem_compat import Chem


class TestSugarAuditFieldsWriter(unittest.TestCase):
    """Test sugar audit fields are properly written to output files."""

    def test_sugar_audit_fields_presence_and_type_csv(self):
        """
        Test that sugar audit fields are present in CSV output with correct types.
        """
        # Use quercetin (flavonoid) which should have some maskable sites but still produce products
        quercetin = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"

        config = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic', 'audit': True}  # Enable audit
        )

        products, qa_stats = enumerate_with_stats(quercetin, config)
        self.assertGreater(len(products), 0, "Should generate products for testing")

        # Write to temporary CSV file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            csv_path = f.name

        try:
            # Write products to CSV
            write_table(products, csv_path)

            # Read back and check fields
            records = read_table(csv_path)
            self.assertGreater(len(records), 0, "Should read back products from CSV")

            # Check first record for required fields
            record = records[0]

            # Check presence of all 4 audit fields
            expected_fields = ['sugar_mask_atoms', 'sugar_rings', 'masked_oh_count', 'masked_bridge_o_count']
            for field in expected_fields:
                self.assertIn(field, record, f"Field '{field}' should be present in CSV record")

            # Check types
            self.assertIsInstance(record['sugar_mask_atoms'], list,
                                "sugar_mask_atoms should be a list")
            self.assertIsInstance(record['sugar_rings'], list,
                                "sugar_rings should be a list")
            self.assertIsInstance(record['masked_oh_count'], int,
                                "masked_oh_count should be an integer")
            self.assertIsInstance(record['masked_bridge_o_count'], int,
                                "masked_bridge_o_count should be an integer")

            # Check that values are reasonable
            self.assertGreaterEqual(record['masked_oh_count'], 0,
                                  "masked_oh_count should be non-negative")
            self.assertGreaterEqual(record['masked_bridge_o_count'], 0,
                                  "masked_bridge_o_count should be non-negative")

            # For quercetin (non-sugar) with sugar masking, masked atoms may be 0 or few
            self.assertGreaterEqual(len(record['sugar_mask_atoms']), 0,
                                  "quercetin should have non-negative masked atoms count")

        finally:
            Path(csv_path).unlink(missing_ok=True)

    def test_sugar_audit_fields_presence_and_type_parquet(self):
        """
        Test that sugar audit fields are present in Parquet output with correct types.
        """
        # Use a simple aromatic molecule that should produce products
        test_molecule = "c1ccc(O)cc1"  # Phenol

        config = EnumConfig(
            k_max=1,
            halogens=('Cl',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic', 'audit': True}
        )

        products, qa_stats = enumerate_with_stats(test_molecule, config)
        self.assertGreater(len(products), 0, "Should generate products for testing")

        # Write to temporary Parquet file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.parquet', delete=False) as f:
            parquet_path = f.name

        try:
            # Write products to Parquet
            write_table(products, parquet_path)

            # Read back and check fields
            records = read_table(parquet_path)
            self.assertGreater(len(records), 0, "Should read back products from Parquet")

            # Check first record for required fields
            record = records[0]

            # Check presence of all 4 audit fields
            audit_fields = {
                'sugar_mask_atoms': list,
                'sugar_rings': list,
                'masked_oh_count': int,
                'masked_bridge_o_count': int
            }

            for field, expected_type in audit_fields.items():
                self.assertIn(field, record, f"Field '{field}' should be present in Parquet record")
                self.assertIsInstance(record[field], expected_type,
                                    f"Field '{field}' should be of type {expected_type.__name__}")

        finally:
            Path(parquet_path).unlink(missing_ok=True)

    def test_sugar_audit_fields_json_serialization(self):
        """
        Test that JSON fields are properly serializable and deserializable.
        """
        test_molecule = "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"  # Quercetin

        config = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic', 'audit': True}
        )

        products, qa_stats = enumerate_with_stats(test_molecule, config)
        self.assertGreater(len(products), 0, "Should generate products for testing")

        # Get a product record
        record = products[0]

        # Check that JSON fields can be serialized and deserialized
        json_fields = ['sugar_mask_atoms', 'sugar_rings']
        for field in json_fields:
            self.assertIn(field, record, f"Field '{field}' should be present")

            # Test JSON serialization
            try:
                json_str = json.dumps(record[field])
                self.assertIsInstance(json_str, str, f"Field '{field}' should be JSON serializable")

                # Test JSON deserialization
                deserialized = json.loads(json_str)
                self.assertEqual(deserialized, record[field],
                               f"Field '{field}' should round-trip through JSON correctly")
            except Exception as e:
                self.fail(f"Field '{field}' JSON serialization failed: {e}")

    def test_sugar_audit_fields_with_audit_disabled(self):
        """
        Test that audit fields are present but empty when audit is disabled.
        """
        test_molecule = "c1ccc(O)cc1"  # Simple phenol

        config = EnumConfig(
            k_max=1,
            halogens=('F',),
            rules=('R1',),
            sugar_cfg={'mode': 'heuristic', 'audit': False}  # Audit disabled
        )

        products, qa_stats = enumerate_with_stats(test_molecule, config)
        self.assertGreater(len(products), 0, "Should generate products for testing")

        # Check that audit fields are present but empty/zero
        record = products[0]

        self.assertIn('sugar_mask_atoms', record)
        self.assertIn('sugar_rings', record)
        self.assertIn('masked_oh_count', record)
        self.assertIn('masked_bridge_o_count', record)

        # Should be empty/zero when audit is disabled
        self.assertEqual(record['sugar_mask_atoms'], [])
        self.assertEqual(record['sugar_rings'], [])
        self.assertEqual(record['masked_oh_count'], 0)
        self.assertEqual(record['masked_bridge_o_count'], 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)