"""
测试io_utils的JSON字段类型安全修复

验证：
1. substitutions等list字段缺失时默认为[]而非{}
2. constraints_violations等dict字段缺失时默认为{}
3. 读写往返一致性
"""

import json
import os
import tempfile
from halogenator.io_utils import (
    _prepare_records_for_table,
    _get_json_default,
    _JSON_LIST_FIELDS,
    _JSON_DICT_FIELDS,
    write_table,
    read_table
)


def test_json_default_types():
    """测试_get_json_default为不同字段类型返回正确的默认值"""

    # List字段应返回[]
    for field in _JSON_LIST_FIELDS:
        default = _get_json_default(field)
        assert isinstance(default, list), f"{field} should default to list, got {type(default)}"
        assert default == [], f"{field} should default to empty list"

    # Dict字段应返回{}
    for field in _JSON_DICT_FIELDS:
        default = _get_json_default(field)
        assert isinstance(default, dict), f"{field} should default to dict, got {type(default)}"
        assert default == {}, f"{field} should default to empty dict"

    print("[PASS] All JSON fields have correct default types")


def test_prepare_records_missing_fields():
    """测试_prepare_records_for_table在字段缺失时使用正确的默认值"""

    # 模拟一个缺少所有JSON字段的record
    records = [
        {
            'smiles': 'CCO',
            'inchikey': 'TEST-INCHIKEY',
            'k': 1,
            # 故意不包含substitutions等字段
        }
    ]

    prepared = _prepare_records_for_table(records)

    assert len(prepared) == 1
    rec = prepared[0]

    # 验证list字段的默认值
    assert 'substitutions_json' in rec
    subs_val = json.loads(rec['substitutions_json'])
    assert isinstance(subs_val, list), f"substitutions should serialize to list, got {type(subs_val)}"
    assert subs_val == [], f"substitutions should serialize to [], got {subs_val}"

    assert 'sugar_mask_atoms_json' in rec
    mask_val = json.loads(rec['sugar_mask_atoms_json'])
    assert isinstance(mask_val, list), f"sugar_mask_atoms should serialize to list"
    assert mask_val == []

    assert 'sugar_rings_json' in rec
    rings_val = json.loads(rec['sugar_rings_json'])
    assert isinstance(rings_val, list), f"sugar_rings should serialize to list"
    assert rings_val == []

    # 验证dict字段的默认值
    assert 'constraints_violations_json' in rec
    cons_val = json.loads(rec['constraints_violations_json'])
    assert isinstance(cons_val, dict), f"constraints_violations should serialize to dict"
    assert cons_val == {}

    print("[PASS] Missing fields get correct default types")


def test_prepare_records_with_values():
    """测试_prepare_records_for_table正确序列化有值的字段"""

    records = [
        {
            'smiles': 'CCO',
            'substitutions': [
                {'rule': 'R1', 'site': 1, 'halogen': 'F', 'atom_cost': 1},
                {'rule': 'R3', 'site': 5, 'halogen': 'Cl', 'atom_cost': 1}
            ],
            'constraints_violations': {'max_f': 1},
            'sugar_mask_atoms': [10, 11, 12],
            'sugar_rings': [[10, 11, 12, 13, 14]]
        }
    ]

    prepared = _prepare_records_for_table(records)
    rec = prepared[0]

    # 验证list字段正确序列化
    subs = json.loads(rec['substitutions_json'])
    assert isinstance(subs, list)
    assert len(subs) == 2
    assert subs[0]['rule'] == 'R1'

    # 验证dict字段正确序列化
    cons = json.loads(rec['constraints_violations_json'])
    assert isinstance(cons, dict)
    assert cons == {'max_f': 1}

    # 验证其他list字段
    mask_atoms = json.loads(rec['sugar_mask_atoms_json'])
    assert mask_atoms == [10, 11, 12]

    rings = json.loads(rec['sugar_rings_json'])
    assert rings == [[10, 11, 12, 13, 14]]

    print("[PASS] Fields with values serialize correctly")


def test_write_read_roundtrip():
    """测试写入和读取的往返一致性"""

    records = [
        {
            'smiles': 'CCO',
            'inchikey': 'TEST1',
            'k': 1,
            'substitutions': [{'rule': 'R1', 'site': 1}],
            'constraints_violations': {},
            'sugar_mask_atoms': [],
            'sugar_rings': []
        },
        {
            'smiles': 'CCCO',
            'inchikey': 'TEST2',
            'k': 2,
            # 缺少JSON字段，测试默认值
        }
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        path = os.path.join(tmpdir, 'test.parquet')

        # 写入
        write_table(records, path)

        # 读取
        read_records = read_table(path)

        assert len(read_records) == 2

        # 第一个record：有值的字段应保持一致
        rec1 = read_records[0]
        assert rec1['inchikey'] == 'TEST1'
        assert isinstance(rec1['substitutions'], list)
        assert len(rec1['substitutions']) == 1
        assert rec1['substitutions'][0]['rule'] == 'R1'

        # 第二个record：缺失的字段应有正确的默认类型
        rec2 = read_records[1]
        assert rec2['inchikey'] == 'TEST2'
        assert isinstance(rec2['substitutions'], list), "substitutions should be list"
        assert rec2['substitutions'] == [], "substitutions should default to []"
        assert isinstance(rec2['constraints_violations'], dict), "constraints_violations should be dict"
        assert rec2['constraints_violations'] == {}, "constraints_violations should default to {}"

    print("[PASS] Write-read roundtrip maintains correct types")


def test_regression_prevent_empty_object_for_list_fields():
    """
    回归测试：防止list字段被序列化为'{}'

    这是之前R6_methyl k-level bug的根本原因：
    - substitutions字段语义上是list
    - 但io_utils默认用{}序列化缺失字段
    - 导致'substitutions_json': '{}'而非'[]'
    - 后续解析和验证失败
    """

    # 模拟emit_product忘记设置substitutions的情况
    record = {
        'smiles': 'CCO',
        'inchikey': 'TEST',
        'k_ops': 1,
        # 故意不设置substitutions
    }

    prepared = _prepare_records_for_table([record])

    # 验证substitutions_json是'[]'而非'{}'
    subs_json_str = prepared[0]['substitutions_json']
    assert subs_json_str == '[]', f"Expected '[]', got '{subs_json_str}'"

    # 验证解析后是list
    subs = json.loads(subs_json_str)
    assert isinstance(subs, list), f"Expected list, got {type(subs)}"
    assert subs == [], f"Expected empty list, got {subs}"

    print("[PASS] Regression prevented: list fields no longer default to '{}'")


def test_all_json_fields_covered():
    """确保所有_JSON_FIELDS都在类型表中"""
    from halogenator.io_utils import _JSON_FIELDS

    covered_fields = _JSON_LIST_FIELDS | _JSON_DICT_FIELDS

    for field in _JSON_FIELDS:
        assert field in covered_fields, f"Field '{field}' not in type table"

    print("[PASS] All JSON fields are covered in type table")


if __name__ == '__main__':
    print("=" * 80)
    print("io_utils JSON Type Safety Tests")
    print("=" * 80)

    try:
        test_json_default_types()
        test_prepare_records_missing_fields()
        test_prepare_records_with_values()
        test_write_read_roundtrip()
        test_regression_prevent_empty_object_for_list_fields()
        test_all_json_fields_covered()

        print("\n" + "=" * 80)
        print("[ALL TESTS PASSED]")
        print("=" * 80)

    except AssertionError as e:
        print(f"\n[TEST FAILED] {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    except Exception as e:
        print(f"\n[ERROR] {e}")
        import traceback
        traceback.print_exc()
        exit(1)
