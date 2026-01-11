import pandas as pd
import sys
import re

# 1. Load data
raw_df = pd.read_csv(sys.argv[1], sep='\t')
metab_df = pd.read_csv(sys.argv[2], index_col=0)

# 2. Stronger ID Normalization: Extract just the numeric part of the GCF
def get_numeric_gcf(id_str):
    # Regex to find 'GCF' followed by numbers
    match = re.search(r'GCF[._]?(\d+)', str(id_str))
    if match:
        return f"GCF_{match.group(1)}"
    return str(id_str).strip()

raw_df['ID_Match'] = raw_df['Assembly'].apply(get_numeric_gcf)
metab_df['ID_Match'] = metab_df.index.to_series().apply(get_numeric_gcf)

# 3. Create Matrix using these robust IDs
defense_matrix = pd.crosstab(raw_df['ID_Match'], raw_df['System'])
defense_matrix = defense_matrix.clip(upper=1)

# 4. Use ID_Match as the index for metabolic data
metab_df = metab_df.set_index('ID_Match')
# Keep only the first if duplicates arise from stripping versions
metab_df = metab_df[~metab_df.index.duplicated(keep='first')]

# 5. Intersect
common = defense_matrix.index.intersection(metab_df.index)

if len(common) < 5:
    print(f"FAILED: Only found {len(common)} matches.")
    print(f"Defense IDs look like: {list(defense_matrix.index[:3])}")
    print(f"Metabolic IDs look like: {list(metab_df.index[:3])}")
else:
    defense_matrix.loc[common].to_csv('defense_presence_matrix.csv')
    metab_df.loc[common].to_csv('metabolic_scorecard_aligned.csv')
    print(f"SUCCESS: Matrix created with {len(common)} GCFs.")

