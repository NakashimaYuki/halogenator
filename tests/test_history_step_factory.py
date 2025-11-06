"""
测试make_history_step工厂函数

验证：
1. 所有必需字段都存在
2. 可选字段正确处理
3. 不同规则的使用场景
4. 防止字段遗漏
"""

from halogenator.enumerate_k import make_history_step


def test_basic_r1_r2_r3_step():
    """测试R1/R2/R3规则的基本step"""

    step = make_history_step(
        rule='R1',
        site=5,
        halogen='F',
        atom_cost=1,
        depth=1,
        sym=2,
        ring_tag='A'
    )

    # 验证必需字段
    assert step['rule'] == 'R1'
    assert step['site'] == 5
    assert step['halogen'] == 'F'
    assert step['atom_cost'] == 1
    assert step['depth'] == 1

    # 验证R1/R2/R3可选字段
    assert step['sym'] == 2
    assert step['ring_tag'] == 'A'

    # 验证R6字段不存在
    assert 'type' not in step
    assert 'label' not in step

    print("[PASS] R1/R2/R3 basic step")


def test_r1_r2_r3_without_optional_fields():
    """测试R1/R2/R3规则在缺少可选字段时的默认值"""

    step = make_history_step(
        rule='R2a',
        site=10,
        halogen='Cl',
        atom_cost=1,
        depth=2
        # 不提供sym和ring_tag
    )

    assert step['rule'] == 'R2a'
    assert step['site'] == 10
    assert step['halogen'] == 'Cl'
    assert step['atom_cost'] == 1
    assert step['depth'] == 2

    # 验证默认值
    assert step['sym'] is None
    assert step['ring_tag'] == ''

    print("[PASS] R1/R2/R3 without optional fields")


def test_r6_step_mode():
    """测试R6_methyl step mode"""

    step = make_history_step(
        rule='R6_methyl',
        site=20,
        halogen='F',
        atom_cost=1,
        depth=2,
        step_type='step',
        k_ops=2,
        k_atoms=2,
        budget_mode='ops'
    )

    assert step['rule'] == 'R6_methyl'
    assert step['site'] == 20
    assert step['halogen'] == 'F'
    assert step['atom_cost'] == 1
    assert step['depth'] == 2

    # 验证R6特定字段
    assert step['type'] == 'step'
    assert 'label' not in step  # step mode没有label

    # 验证legacy字段
    assert step['k_ops'] == 2
    assert step['k_atoms'] == 2
    assert step['budget_mode'] == 'ops'

    print("[PASS] R6_methyl step mode")


def test_r6_macro_mode():
    """测试R6_methyl macro mode"""

    step = make_history_step(
        rule='R6_methyl',
        site=25,
        halogen='F',
        atom_cost=3,  # CF3 = 3 atoms
        depth=2,
        step_type='macro',
        macro_label='CF3',
        k_ops=2,
        k_atoms=5,
        budget_mode='atoms'
    )

    assert step['rule'] == 'R6_methyl'
    assert step['site'] == 25
    assert step['halogen'] == 'F'
    assert step['atom_cost'] == 3
    assert step['depth'] == 2

    # 验证R6特定字段
    assert step['type'] == 'macro'
    assert step['label'] == 'CF3'

    # 验证legacy字段
    assert step['k_ops'] == 2
    assert step['k_atoms'] == 5
    assert step['budget_mode'] == 'atoms'

    print("[PASS] R6_methyl macro mode")


def test_all_required_fields_present():
    """验证所有必需字段都存在"""

    # 最简单的step（只有必需字段）
    step = make_history_step(
        rule='R3',
        site=1,
        halogen='Br',
        atom_cost=1,
        depth=1
    )

    required_fields = ['rule', 'site', 'halogen', 'atom_cost', 'depth']

    for field in required_fields:
        assert field in step, f"Required field '{field}' missing"

    # 即使不提供可选字段，sym和ring_tag也应该有默认值
    assert 'sym' in step
    assert 'ring_tag' in step

    print("[PASS] All required fields present")


def test_atom_cost_correctness():
    """验证atom_cost字段的正确性"""

    # Single atom (R1/R2/R3, R6 step)
    single_step = make_history_step(
        rule='R1',
        site=1,
        halogen='F',
        atom_cost=1,
        depth=1
    )
    assert single_step['atom_cost'] == 1

    # Macro (R6 CF3/CCl3)
    macro_step = make_history_step(
        rule='R6_methyl',
        site=1,
        halogen='F',
        atom_cost=3,
        depth=1,
        step_type='macro',
        macro_label='CF3'
    )
    assert macro_step['atom_cost'] == 3

    print("[PASS] atom_cost correctness")


def test_legacy_fields_optional():
    """验证legacy字段是可选的"""

    # 不提供legacy字段
    step = make_history_step(
        rule='R1',
        site=1,
        halogen='F',
        atom_cost=1,
        depth=1
    )

    # Legacy字段不应存在
    assert 'k_ops' not in step
    assert 'k_atoms' not in step
    assert 'budget_mode' not in step

    print("[PASS] Legacy fields optional")


def test_regression_prevent_missing_atom_cost():
    """
    回归测试：防止遗漏atom_cost字段

    这是之前R6_methyl bug的根本原因之一：
    - 某些规则的history item缺少atom_cost
    - 导致k_atoms计算错误
    - 最终导致k-level错标
    """

    # 使用工厂函数，atom_cost是必需参数，无法遗漏
    try:
        # 尝试不提供atom_cost（应该失败）
        step = make_history_step(
            rule='R1',
            site=1,
            halogen='F',
            depth=1
            # 故意不提供atom_cost
        )
        # 如果能执行到这里，说明函数没有正确要求atom_cost
        assert False, "Should have raised TypeError for missing atom_cost"
    except TypeError as e:
        # 预期行为：缺少必需参数
        assert 'atom_cost' in str(e)

    print("[PASS] Regression: atom_cost is required")


def test_history_composition():
    """测试构建完整的history链"""

    # 模拟k=2的history构建
    history = []

    # k=1: R1规则
    step1 = make_history_step(
        rule='R1',
        site=5,
        halogen='F',
        atom_cost=1,
        depth=1,
        sym=2,
        ring_tag='A'
    )
    history.append(step1)

    # k=2: R3规则
    step2 = make_history_step(
        rule='R3',
        site=10,
        halogen='Cl',
        atom_cost=1,
        depth=2,
        sym=None,
        ring_tag=''
    )
    history.append(step2)

    # 验证history链
    assert len(history) == 2

    # 验证k_ops和k_atoms可以从history计算
    k_ops = len(history)
    k_atoms = sum(step.get('atom_cost', 0) for step in history)

    assert k_ops == 2
    assert k_atoms == 2

    print("[PASS] History composition")


if __name__ == '__main__':
    print("=" * 80)
    print("make_history_step Factory Function Tests")
    print("=" * 80)

    try:
        test_basic_r1_r2_r3_step()
        test_r1_r2_r3_without_optional_fields()
        test_r6_step_mode()
        test_r6_macro_mode()
        test_all_required_fields_present()
        test_atom_cost_correctness()
        test_legacy_fields_optional()
        test_regression_prevent_missing_atom_cost()
        test_history_composition()

        print("\n" + "=" * 80)
        print("[ALL TESTS PASSED]")
        print("=" * 80)
        print("\nUsage Note:")
        print("This factory function provides a single, consistent way to create")
        print("history items. Use it in new code to avoid field omission bugs.")
        print("Existing code can continue to use inline dict construction,")
        print("but consider migrating to the factory for better maintainability.")

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
