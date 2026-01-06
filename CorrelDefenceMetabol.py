import pandas as pd
import numpy as np
import sys
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Use the files provided in the command line or defaults
DEFENSE_FILE = sys.argv[1] if len(sys.argv) > 1 else 'genome.txt'
METABOLIC_FILE = sys.argv[2] if len(sys.argv) > 2 else 'informative_metabolism.csv'

def run_cost_analysis():
    print(f"--- 1. Loading Datasets ---")
    print(f"Metabolic File: {METABOLIC_FILE}")
    print(f"Defense File: {DEFENSE_FILE}")
    
    # 1. Load Metabolism with space-stripping
    # sep=None, engine='python' detects if it is Tab or Comma automatically
    met_df = pd.read_csv(METABOLIC_FILE, sep=None, engine='python', index_col=0)
    
    # Clean whitespace from column names and index
    met_df.index = [str(i).strip().replace('.', '_') for i in met_df.index]
    met_df.columns = [c.strip() for c in met_df.columns]
    
    # CRITICAL: Strip spaces from all cells and convert to numeric
    # This fixes the "TypeError: '>' not supported between instances of 'str' and 'int'"
    met_df = met_df.apply(lambda x: pd.to_numeric(x.astype(str).str.strip(), errors='coerce'))
    
    # Filter out pathways that are 0 or NaN everywhere
    met_df = met_df.loc[:, (met_df > 0).any(axis=0)].dropna(axis=1, how='all')
    print(f"Loaded {len(met_df)} genomes and {met_df.shape[1]} active pathways.")

    # 2. Load Defense Finder
    def_raw = pd.read_csv(DEFENSE_FILE, sep=None, engine='python')
    
    # Ensure Assembly IDs match the metabolic IDs
    def_raw['Assembly'] = def_raw['Assembly'].astype(str).str.strip().str.replace('.', '_', regex=False)

    print("--- 2. Processing Defense Metrics ---")
    
    # Analysis A: Total Load (Count of systems per Assembly)
    def_counts = def_raw.groupby('Assembly').size().rename('total_defense_systems')
    
    # Analysis B: Specific Systems (Presence/Absence matrix)
    def_binary = pd.crosstab(def_raw['Assembly'], def_raw['System'])
    def_binary = (def_binary > 0).astype(int)

    # 3. MERGE
    full_df = met_df.join(def_counts, how='inner').join(def_binary, how='inner')
    print(f"Genomes matched for analysis: {len(full_df)}")
    
    if len(full_df) == 0:
        print("ERROR: No matching GCF IDs found between the two files!")
        print(f"Example Metabolic ID: {met_df.index[0]}")
        print(f"Example Defense ID: {def_raw['Assembly'].iloc[0]}")
        return

    # ---------------------------------------------------------
    # ANALYSIS 1: Total Count (Spearman Correlation)
    # ---------------------------------------------------------
    print("--- 3. Calculating Total Load Costs ---")
    spearman_results = []
    for pathway in met_df.columns:
        coeff, pval = spearmanr(full_df['total_defense_systems'], full_df[pathway], nan_policy='omit')
        spearman_results.append({'Pathway': pathway, 'Correlation': coeff, 'p_value': pval})
    
    spearman_df = pd.DataFrame(spearman_results)
    spearman_df = spearman_df.dropna()
    if not spearman_df.empty:
        spearman_df['p_adj'] = multipletests(spearman_df['p_value'], method='fdr_bh')[1]
    
    # ---------------------------------------------------------
    # ANALYSIS 2: Specific Systems (Mann-Whitney U)
    # ---------------------------------------------------------
    print("--- 4. Calculating Specific System Costs ---")
    mw_results = []
    systems_to_test = [s for s in def_binary.columns if s in full_df.columns and full_df[s].sum() >= 5]
    
    for system in systems_to_test:
        for pathway in met_df.columns:
            has_sys = full_df[full_df[system] == 1][pathway]
            no_sys = full_df[full_df[system] == 0][pathway]
            
            if len(has_sys) < 1 or len(no_sys) < 1: continue
            
            stat, pval = mannwhitneyu(has_sys, no_sys, alternative='two-sided')
            mean_diff = has_sys.mean() - no_sys.mean()
            
            mw_results.append({
                'Defense_System': system, 
                'Metabolic_Pathway': pathway, 
                'Effect_Size_MeanDiff': mean_diff, 
                'p_value': pval
            })
            
    mw_df = pd.DataFrame(mw_results)
    if not mw_df.empty:
        mw_df['p_adj'] = multipletests(mw_df['p_value'], method='fdr_bh')[1]

    # --- 5. SAVE RESULTS ---
    spearman_df.sort_values('Correlation').to_csv('total_load_costs.csv', index=False)
    mw_df.sort_values('p_adj').to_csv('specific_system_costs.csv', index=False)
    
    print("\nSUCCESS!")
    print("Files created: 'total_load_costs.csv' and 'specific_system_costs.csv'")

if __name__ == "__main__":
    run_cost_analysis()
