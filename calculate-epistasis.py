Last login: Sun Jan 11 13:07:28 on ttys000
(base) giusym@Giuseppinas-MacBook-Pro MetabolismCorrel % nano defense-matrix-gen.py 
(base) giusym@Giuseppinas-MacBook-Pro MetabolismCorrel % nano calculate-epistasis.py 






























































  UW PICO 5.09                                                                                     File: calculate-epistasis.py                                                                                       

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

print(f"Analyzing {len(systems)} systems across {len(pathways)} pathways...")

# 2. Iterate through system pairs
for sysA, sysB in combinations(systems, 2):
    # Check co-occurrence frequency
    co_occur = df_defense[(df_defense[sysA] == 1) & (df_defense[sysB] == 1)]

    # Identify "Statistically Forbidden" (Rarity)
    if len(co_occur) < 3:
        results.append({
            'Pair': f"{sysA} x {sysB}",
            'Type': 'Forbidden (Too Rare)',
            'Interaction_Score': -1.0,
            'P_Value': 0.001
        })
        continue

    # 3. Model Interaction for each pathway: Metabolic_Score ~ A + B + (A*B)
    for path in pathways:
        y = df_metab[path]
        X = df_defense[[sysA, sysB]].copy()
        X['Interaction'] = X[sysA] * X[sysB]
        X = sm.add_constant(X)

        try:
            model = sm.OLS(y, X).fit()
            inter_coef = model.params['Interaction']
            p_val = model.pvalues['Interaction']

            # Focus on Significant Synergistic Costs
            if p_val < 0.05 and inter_coef < 0:
                results.append({
                    'Pair': f"{sysA} x {sysB}",
                    'Pathway': path,
                    'Type': 'Synergistic Cost',
                    'Interaction_Score': inter_coef,
                    'P_Value': p_val
                })
        except:
            continue

# 4. Save the results
results_df = pd.DataFrame(results).sort_values('P_Value')
results_df.to_csv('forbidden_combinations_results.csv', index=False)
print("Analysis complete. Results saved to 'forbidden_combinations_results.csv'.")
