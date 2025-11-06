# -*- coding: ascii -*-
"""Pre-commit schema validation check."""

from src.halogenator.schema import (
    validate_products_schema,
    get_all_known_columns,
    validate_products_records,
    empty_qa_paths,
    SUGAR_AUDIT_KEYS,
    CORE_SUGAR_AUDIT_KEYS,
    FALLBACK_SUGAR_AUDIT_KEYS,
)


def main():
    """Run schema validation checks for pre-commit hook."""
    # Validate all known columns
    validate_products_schema(get_all_known_columns())

    # Test minimal record validation with both P0 and P1 fields
    rec = {
        'inchikey': 'TEST-KEY',
        'smiles': 'CCF',
        'product_smiles': 'CCF',
        'parent_inchikey': 'PARENT-KEY',
        'parent_smiles': 'CC',
        'parent_type': 'manual',
        'rule': 'R1',
        'rule_family': 'R1',
        'halogen': 'F',
        'k': 1,
        'depth': 1,
        'k_ops': 1,
        'constraints_ok': True,
    }
    validate_products_records([rec])

    # Test QA paths factory
    qa = empty_qa_paths()
    assert isinstance(qa, dict) and len(qa) > 0

    # Test sugar audit schema constants
    assert len(CORE_SUGAR_AUDIT_KEYS) > 0, "Core sugar audit keys should not be empty"
    assert len(FALLBACK_SUGAR_AUDIT_KEYS) > 0, "Fallback sugar audit keys should not be empty"
    assert len(SUGAR_AUDIT_KEYS) > 0, "Combined sugar audit keys should not be empty"

    # Check no duplicates in sugar audit keys
    assert len(SUGAR_AUDIT_KEYS) == len(set(SUGAR_AUDIT_KEYS)), "Sugar audit keys should not contain duplicates"

    # Check core keys are subset of combined
    core_set = set(CORE_SUGAR_AUDIT_KEYS)
    combined_set = set(SUGAR_AUDIT_KEYS)
    assert core_set.issubset(combined_set), "Core sugar audit keys should be subset of combined keys"

    print('Schema validation OK')


if __name__ == '__main__':
    main()