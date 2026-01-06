import pandas as pd

# Load the 820-row scorecard
df = pd.read_csv('metabolic_scorecard.tsv', sep='\t', index_col=0)

# 1. Keep columns where at least ONE strain has a score > 0
# This removes "Methanogenesis" but KEEPS anything that is present in even 1 strain
df_variable = df.loc[:, (df > 0).any(axis=0)]

# 2. Identify "Loss" events
# If a pathway is 1.0 in most but 0.0 in high-defense strains, we want to see it.
print(f"Kept {df_variable.shape[1]} informative pathways.")

# 3. Save for correlation
df_variable.to_csv('informative_metabolism.csv')
