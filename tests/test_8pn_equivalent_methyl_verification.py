"""
验证8-PN场景：异戊烯基两个甲基产生的产物是否为同一结构
这是raw模式（无折叠、无去重）的预期行为，不是bug
"""

import pandas as pd
import json
import sys


def verify_8pn_equivalent_methyl():
    """验证8-PN异戊烯基两个甲基产生的产物是否真的相同"""

    print("=" * 80)
    print("8-PN 异戊烯基等价甲基验证")
    print("=" * 80)

    # 读取数据
    df = pd.read_parquet('data/output/8pn_raw_v4/products_k2.parquet')

    print("\n【场景配置】")
    print("- 分子: 8-prenylnaringenin (异戊烯基黄酮)")
    print("- 模式: Raw (--no-constraints --no-sugar-mask --no-sym-fold --no-dedup --r2-fallback)")
    print("- 异戊烯基: 含有两个等价的甲基（gem-dimethyl）")

    # 分析k=1的R6_methyl产物
    k1_df = df[df['k'] == 1]
    r6_k1 = k1_df[k1_df['rule'] == 'R6_methyl']

    print("\n【K=1 R6_methyl产物分析】")
    print(f"总数: {len(r6_k1)}")

    # 按卤素分组，查看是否有重复的InChIKey
    for halogen in ['F', 'Cl', 'Br', 'I']:
        hal_subset = r6_k1[r6_k1['halogen'] == halogen]
        print(f"\n{halogen}:")
        print(f"  产物数: {len(hal_subset)}")

        # 检查InChIKey重复情况
        inchikey_counts = hal_subset['inchikey'].value_counts()
        if len(inchikey_counts) < len(hal_subset):
            print(f"  唯一InChIKey: {len(inchikey_counts)} （存在重复！）")
            duplicates = inchikey_counts[inchikey_counts > 1]
            for ik, count in duplicates.items():
                print(f"    {ik}: {count} 个记录")
                # 显示这些重复记录的site信息
                dup_records = hal_subset[hal_subset['inchikey'] == ik]
                for idx, rec in dup_records.iterrows():
                    subs = json.loads(rec['substitutions_json']) if rec['substitutions_json'] else []
                    sites = [s.get('site') for s in subs]
                    print(f"      - site={sites}, smiles={rec['smiles'][:60]}...")
        else:
            print(f"  唯一InChIKey: {len(inchikey_counts)} （无重复）")

    # K=2也做类似分析
    print("\n【K=2 R6_methyl相关产物分析】")
    k2_df = df[df['k'] == 2]

    # 找出第二步是R6_methyl的产物
    r6_k2_records = []
    for idx, row in k2_df.iterrows():
        subs = json.loads(row['substitutions_json']) if row['substitutions_json'] else []
        if len(subs) >= 2 and subs[1].get('rule') == 'R6_methyl':
            r6_k2_records.append(row)

    if r6_k2_records:
        print(f"K=2中第二步为R6_methyl的产物数: {len(r6_k2_records)}")
        # 检查这些产物是否也有InChIKey重复
        r6_k2_df = pd.DataFrame(r6_k2_records)
        inchikey_counts = r6_k2_df['inchikey'].value_counts()
        duplicates = inchikey_counts[inchikey_counts > 1]
        if len(duplicates) > 0:
            print(f"  存在{len(duplicates)}个重复的InChIKey:")
            for ik, count in list(duplicates.items())[:3]:  # 只显示前3个
                print(f"    {ik}: {count} 个记录")

    # 结论
    print("\n【验证结论】")

    k1_duplicates_exist = False
    for halogen in ['F', 'Cl', 'Br', 'I']:
        hal_subset = r6_k1[r6_k1['halogen'] == halogen]
        if len(hal_subset['inchikey'].unique()) < len(hal_subset):
            k1_duplicates_exist = True
            break

    if k1_duplicates_exist:
        print("[PASS] 确认：异戊烯基的两个甲基在raw模式下产生了**结构相同**的产物")
        print("   - 这些产物有不同的site（原子位点编号）")
        print("   - 但有相同的InChIKey（化学结构一致）")
        print("   - 在PyMOL中显示外观完全一致是正常的")
    else:
        print("[FAIL] 未发现重复的InChIKey，可能异戊烯基的两个甲基不等价")

    print("\n【这是Bug吗？】")
    print("[NO] 不是Bug。这是raw模式的预期行为：")
    print("   1. Raw模式下关闭了对称折叠（--no-sym-fold）")
    print("   2. 也关闭了去重（--no-dedup）")
    print("   3. 因此等价位点会被独立枚举，产生'路径不同但结构相同'的产物")
    print("   4. 这保留了路径的多样性，适合某些分析场景")

    print("\n【如何避免？】")
    print("如果只想要结构唯一的产物，有两种方法：")
    print("   1. 启用对称折叠（去掉 --no-sym-fold）→ 在生成前减少等价位点枚举")
    print("   2. 启用去重（去掉 --no-dedup）→ 在生成后合并相同结构")
    print("   3. 或实现'规则级微折叠'（见可选优化A）")

    return k1_duplicates_exist


if __name__ == '__main__':
    try:
        result = verify_8pn_equivalent_methyl()
        sys.exit(0 if result else 1)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
