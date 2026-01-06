import sys
import pandas as pd

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <padloc.tsv> <defensefinder.tsv> [output.tsv]")
        sys.exit(1)

    padloc_file = sys.argv[1]
    df_file = sys.argv[2]
    out_file = sys.argv[3] if len(sys.argv) >= 4 else "nonredundant_systems.tsv"

    # ---------- READ INPUT ----------
    padloc = pd.read_csv(padloc_file, sep="\t")
    df = pd.read_csv(df_file, sep="\t")

    # ---------- PADLOC: collect all protein IDs ----------
    # assumes PADLOC protein column is 'target.name' with WP_...
    padloc_proteins = (
        padloc["target.name"]
        .astype(str)
        .str.strip()
        .dropna()
        .unique()
    )
    padloc_protein_set = set(padloc_proteins)

    # ---------- DEFENSEFINDER: parse protein_in_syst ----------
    df_expanded = df.assign(
        protein=df["protein_in_syst"].str.split(",")
    ).explode("protein")

    df_expanded["protein"] = (
        df_expanded["protein"]
        .astype(str)
        .str.strip()
    )

    # For each DF system (sys_id), check if it has ANY protein not in PADLOC
    def has_any_new_protein(protein_series: pd.Series) -> bool:
        for p in protein_series:
            if p and p not in padloc_protein_set:
                # at least one new protein -> keep this DF system
                return True
        # all proteins already in PADLOC -> drop
        return False

    # This gives a Series indexed by sys_id with boolean keep/drop
    df_sys_flag = df_expanded.groupby("sys_id")["protein"].apply(has_any_new_protein)

    # sys_ids to keep
    keep_ids = set(df_sys_flag[df_sys_flag].index)

    # Filter original DF table to only those sys_ids
    df_dedup = df[df["sys_id"].isin(keep_ids)].copy()

    # ---------- COMBINE PADLOC + filtered DF ----------
    padloc_out = padloc.copy()
    padloc_out["source"] = "PADLOC"

    df_out = df_dedup.copy()
    df_out["source"] = "DefenseFinder"

    combined = pd.concat([padloc_out, df_out], ignore_index=True)
    combined.to_csv(out_file, sep="\t", index=False)

    print(f"Done! Wrote: {out_file}")
    print(f"PADLOC rows kept: {len(padloc_out)}")
    total_df_sys = df['sys_id'].nunique()
    kept_df_sys = df_out['sys_id'].nunique()
    print(f"DefenseFinder systems removed as duplicates (all proteins already in PADLOC): "
          f"{total_df_sys - kept_df_sys}")
    print(f"DefenseFinder systems kept (≥1 new protein): {kept_df_sys}")

if __name__ == "__main__":
    main()
