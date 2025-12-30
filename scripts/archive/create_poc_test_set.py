"""Create a small POC test set of 10 terpenoids"""
import pandas as pd
from pathlib import Path

# Read terpenoid base library
input_file = Path('E:/Projects/halogenator/data/output/nplike/terpenoid/base_clean.parquet')
df = pd.read_parquet(input_file)

# Take first 10 molecules
poc_df = df.head(10)

# Save to new file
output_file = Path('E:/Projects/halogenator/data/output/nplike/terpenoid_poc10.parquet')
poc_df.to_parquet(output_file, index=False)

print(f"Created POC test set with {len(poc_df)} molecules")
print(f"Output: {output_file}")
