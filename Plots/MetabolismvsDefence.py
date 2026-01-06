import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load the data
total_load = pd.read_csv('total_load_costs.csv')
specific_costs = pd.read_csv('specific_system_costs.csv')

# --- Plot 1: Global Correlations (Top 10 Pos/Neg) ---
top_neg = total_load.nsmallest(10, 'Correlation')
top_pos = total_load.nlargest(10, 'Correlation')
plot_data = pd.concat([top_neg, top_pos]).sort_values('Correlation')

plt.figure(figsize=(10, 8))
colors = ['#d7191c' if x < 0 else '#2c7bb6' for x in plot_data['Correlation']]
sns.barplot(x='Correlation', y='Pathway', data=plot_data, palette=colors)
plt.title('Global Metabolic Costs: Correlation with Total Defense Count', fontsize=14)
plt.grid(False) # Remove internal grid
plt.axvline(0, color='black', linewidth=0.8)
plt.tight_layout()
plt.savefig('global_metabolic_costs.png')

# --- Plot 2: Specific System Trade-offs (Bubble Plot) ---
# Filter for top 30 most significant (p_adj) hits
top_specific = specific_costs[specific_costs['p_adj'] < 0.05].nsmallest(30, 'p_adj').copy()
top_specific['neg_log_p'] = -np.log10(top_specific['p_adj'])

plt.figure(figsize=(14, 10))
ax = sns.scatterplot(
    data=top_specific, x='Metabolic_Pathway', y='Defense_System',
    size='neg_log_p', hue='Effect_Size_MeanDiff', palette='RdBu_r',
    sizes=(100, 1000), hue_norm=(-1, 1), edgecolor='black', alpha=0.8
)

plt.title('Top 30 Specific Defense-Metabolism Interactions', fontsize=16)
plt.xticks(rotation=45, ha='right')

# Move legend outside to prevent overlap
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0, title='-log10(p_adj) & Effect Size')

# Remove grid lines
ax.grid(False)

plt.tight_layout()
plt.savefig('specific_system_tradeoffs.png')

