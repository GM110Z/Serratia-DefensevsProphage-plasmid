import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 1. Plotting Counts per Genome
genome_df = pd.read_csv('genome-filt.txt', sep='\t')
counts = genome_df.groupby('Assembly')['System'].count().reset_index()

plt.figure(figsize=(10, 6))
sns.countplot(data=counts, x='System', color='steelblue')
plt.title('Distribution of Defense System Counts per Genome')
plt.xlabel('Number of Systems')
plt.ylabel('Genome Count')
plt.grid(False)
plt.tight_layout()
plt.savefig('defense_counts_per_genome.pdf')

# 2. Plotting Specific Correlations
spec_costs = pd.read_csv('specific_system_costs.csv')
top_spec = spec_costs[spec_costs['p_adj'] < 0.05].nsmallest(30, 'p_adj').copy()
top_spec['neg_log_p'] = -np.log10(top_spec['p_adj'])

plt.figure(figsize=(15, 11))
ax = sns.scatterplot(
    data=top_spec, x='Metabolic_Pathway', y='Defense_System',
    size='neg_log_p', hue='Effect_Size_MeanDiff', palette='RdBu_r',
    sizes=(100, 1000), hue_norm=(-1, 1), edgecolor='black', alpha=0.8
)
plt.xticks(rotation=40, ha='right')
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', title='-log10(p_adj)')
ax.grid(False)
plt.tight_layout()
plt.savefig('detailed_defense_tradeoffs.pdf')
