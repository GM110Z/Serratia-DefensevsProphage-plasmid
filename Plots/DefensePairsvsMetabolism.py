import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 1. Load the two "Pair" datasets
try:
    df_neg = pd.read_csv('forbidden_combinations_results.csv')
    df_pos = pd.read_csv('positive_pair_results.csv')
except FileNotFoundError:
    print("Error: Ensure results files exist.")
    exit()

def prepare_data(df, threshold=1e-10, top_n=10):
    df = df.dropna(subset=['Pair', 'Pathway']).copy()
    df = df[df['P_Value'] < threshold].copy()
    df['log_p'] = -np.log10(df['P_Value']).clip(upper=300)
    top_pairs = df.groupby('Pair')['log_p'].max().nlargest(top_n).index
    return df[df['Pair'].isin(top_pairs)]

# Prepare datasets
df_neg_clean = prepare_data(df_neg, threshold=1e-50, top_n=10)
df_pos_clean = prepare_data(df_pos, threshold=1e-5, top_n=10)

# --- FIGURE 1: COLORED BY SIGNIFICANCE ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
sns.set_style("whitegrid")
# Red for Conflicts, Green for Synergy
sc1 = ax1.scatter(x=df_neg_clean['Pathway'], y=df_neg_clean['Pair'], s=df_neg_clean['log_p']*8, 
                  c=df_neg_clean['log_p'], cmap='Reds', alpha=0.9, edgecolors='black')
sc2 = ax2.scatter(x=df_pos_clean['Pathway'], y=df_pos_clean['Pair'], s=df_pos_clean['log_p']*8, 
                  c=df_pos_clean['log_p'], cmap='Greens', alpha=0.9, edgecolors='black')
plt.savefig('Interaction_Map_Significance.png', dpi=400)
plt.close()

# --- FIGURE 2: COLORED BY INTERACTION SCORE ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
sns.set_style("whitegrid")
# Reds_r for negative scores, Greens for positive
sc1 = ax1.scatter(x=df_neg_clean['Pathway'], y=df_neg_clean['Pair'], s=df_neg_clean['log_p']*8, 
                  c=df_neg_clean['Interaction_Score'], cmap='Reds_r', alpha=0.9, edgecolors='black')
sc2 = ax2.scatter(x=df_pos_clean['Pathway'], y=df_pos_clean['Pair'], s=df_pos_clean['log_p']*8, 
                  c=df_pos_clean['Interaction_Score'], cmap='Greens', alpha=0.9, edgecolors='black')
plt.savefig('Interaction_Map_Scores.png', dpi=400)
plt.close()
