import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load Data
total_load = pd.read_csv('total_load_costs.csv')
spec_costs = pd.read_csv('specific_system_costs.csv')

# --- Plot 1: Global Correlations (Top 10 each side) ---
top_neg = total_load.nsmallest(10, 'Correlation')
top_pos = total_load.nlargest(10, 'Correlation')
plot_data = pd.concat([top_neg, top_pos]).sort_values('Correlation')

plt.figure(figsize=(10, 8))
colors = ['#d7191c' if x < 0 else '#2c7bb6' for x in plot_data['Correlation']]
sns.barplot(x='Correlation', y='Pathway', data=plot_data, palette=colors)
plt.title('Global Metabolic Costs and Synergies (Excluding PD-T4-6)', fontsize=14)
plt.axvline(0, color='black', linewidth=0.8)
plt.tight_layout()
plt.savefig('global_metabolic_costs.png')

# --- Plot 2: Specific System Trade-offs (Bubble Plot) ---
sig_spec = spec_costs[spec_costs['p_adj'] < 0.05].nsmallest(30, 'p_adj').copy()
sig_spec['neg_log_p'] = -np.log10(sig_spec['p_adj'])

plt.figure(figsize=(12, 10))
sns.scatterplot(
    data=sig_spec, x='Metabolic_Pathway', y='Defense_System',
    size='neg_log_p', hue='Effect_Size_MeanDiff', palette='RdBu_r',
    sizes=(50, 600), hue_norm=(-1, 1), edgecolor='black', alpha=0.8
)
plt.title('Top 30 Specific System-Pathway Interactions', fontsize=16)
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('specific_system_tradeoffs.png')

