Last login: Sun Jan 11 13:07:28 on ttys000
(base) giusym@Giuseppinas-MacBook-Pro MetabolismCorrel % nano defense-matrix-gen.py 
(base) giusym@Giuseppinas-MacBook-Pro MetabolismCorrel % nano calculate-epistasis.py 
(base) giusym@Giuseppinas-MacBook-Pro MetabolismCorrel % nano plot1.py 
(base) giusym@Giuseppinas-MacBook-Pro MetabolismCorrel % nano calculate-epistasis.py
(base) giusym@Giuseppinas-MacBook-Pro MetabolismCorrel % nano plot3.py              



























































  UW PICO 5.09                                                                                            File: plot3.py                                                                                              

import pandas as pd
import numpy as np
import statsmodels.api as sm
from itertools import combinations

# 1. Load the aligned files
df_defense = pd.read_csv('defense_presence_matrix.csv', index_col=0)
df_metab = pd.read_csv('metabolic_scorecard_aligned.csv', index_col=0)

results = []
systems = df_defense.columns
pathways = df_metab.columns

print(f"Hunting for Metabolic Synergy across {len(systems)} systems...")

# 2. Iterate through single systems first (Main Effects)
for sys in systems:
    if df_defense[sys].sum() < 5: continue # Skip extremely rare systems

    for path in pathways:
        y = df_metab[path]
        X = df_defense[[sys]]
        X = sm.add_constant(X)

        try:
            model = sm.OLS(y, X).fit()
            coef = model.params[sys]
            p_val = model.pvalues[sys]

            # Look for POSITIVE coefficients (Presence = Higher Metabolic Score)
            if p_val < 0.05 and coef > 0:
                results.append({
                    'System': sys,
                    'Pathway': path,
                    'Type': 'Positive Association',
                    'Benefit_Score': coef,
                    'P_Value': p_val
                })
        except:
            continue

# 3. Save the "Synergy" results
synergy_df = pd.DataFrame(results).sort_values('Benefit_Score', ascending=False)
synergy_df.to_csv('metabolic_synergy_results.csv', index=False)
print("Synergy analysis complete. Saved to 'metabolic_synergy_results.csv'.")


