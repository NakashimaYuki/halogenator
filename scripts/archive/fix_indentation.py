"""Fix indentation in enumerate_k1.py after the edit"""
import re

file_path = r'E:\Projects\halogenator\src\halogenator\enumerate_k1.py'

with open(file_path, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find the line "for rule in reaction_based_rules:" and fix indentation for the block
in_block = False
block_start = -1
for i, line in enumerate(lines):
    if 'for rule in reaction_based_rules:' in line:
        in_block = True
        block_start = i
        continue

    if in_block:
        # Check if we've exited the block (line starts with less indentation than expected)
        if line.strip() and not line.startswith('        ') and not line.startswith('\t'):
            # Reached end of block at first non-indented line
            if i > block_start + 10:  # Only break if we've processed some lines
                break

        # Fix indentation: change 16 spaces (or more) at start to 12 spaces
        if line.startswith('                '):  # 16+ spaces
            # Count leading spaces
            stripped = line.lstrip(' ')
            leading_spaces = len(line) - len(stripped)
            if leading_spaces >= 16:
                # Reduce by 4 spaces
                lines[i] = ' ' * (leading_spaces - 4) + stripped

with open(file_path, 'w', encoding='utf-8') as f:
    f.writelines(lines)

print("Indentation fixed!")
