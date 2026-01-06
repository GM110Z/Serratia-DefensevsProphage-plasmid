import sys
import time
import pandas as pd
from Bio import Entrez

# -------------- CONFIG --------------
Entrez.email = "your.email@somewhere"   # <-- put your email here
Entrez.api_key = None                   # optional, if you have an NCBI key
SLEEP = 0.34                            # ~3 requests/sec without API key
# ------------------------------------


def get_assembly_for_nuccore(acc):
    """
    Given a nuccore accession (e.g. NZ_GG753567.1),
    try to find the linked assembly accession (e.g. GCF_000163595.1).

    Returns a string (assembly accession) or None if not found.
    """
    from Bio import Entrez

    # Step 1: elink from nuccore to assembly
    # Entrez.elink can take accession strings directly in 'id'
    try:
        handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=acc, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"ELINK error for {acc}: {e}")
        return None

    if not records or "LinkSetDb" not in records[0] or len(records[0]["LinkSetDb"]) == 0:
        return None

    # there may be multiple assemblies; we just take the first for now
    links = records[0]["LinkSetDb"][0]["Link"]
    if not links:
        return None

    assembly_uid = links[0]["Id"]

    # Step 2: esummary on assembly uid to get the accession
    try:
        handle = Entrez.esummary(db="assembly", id=assembly_uid, retmode="xml")
        summary = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"ESUMMARY error for assembly UID {assembly_uid} (from {acc}): {e}")
        return None

    # The structure is DocumentSummarySet / DocumentSummary list
    try:
        docsum = summary["DocumentSummarySet"]["DocumentSummary"][0]
        # The key name is typically 'AssemblyAccession' or similar;
        # if unsure, print(docsum.keys()) once to inspect.
        asm_acc = docsum.get("AssemblyAccession", None)
        return asm_acc
    except Exception as e:
        print(f"Could not parse assembly accession for nuccore {acc}: {e}")
        return None


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.xlsx> [output.xlsx]")
        sys.exit(1)

    in_file = sys.argv[1]
    out_file = sys.argv[2] if len(sys.argv) >= 3 else "Combo-pad-def.fin.ALL.withAssembly.xlsx"

    # Load your merged PADLOC+DefenseFinder file
    df = pd.read_excel(in_file)

    # Expecting columns: 'Nuccore', 'System', 'Assembly'
    if "Nuccore" not in df.columns or "Assembly" not in df.columns:
        raise ValueError("Input file must have at least 'Nuccore' and 'Assembly' columns.")

    # Identify rows needing assembly info
    mask_missing = df["Assembly"].isna() | (df["Assembly"].astype(str).str.strip() == "")
    nuccore_to_fill = df.loc[mask_missing, "Nuccore"].dropna().astype(str).unique()

    print(f"Unique Nuccore accessions needing Assembly: {len(nuccore_to_fill)}")

    nuccore_to_assembly = {}

    for i, acc in enumerate(nuccore_to_fill, start=1):
        print(f"[{i}/{len(nuccore_to_fill)}] Querying assembly for {acc}...")
        asm = get_assembly_for_nuccore(acc)
        if asm is None:
            print(f"  -> No assembly found for {acc}")
        else:
            print(f"  -> {acc} -> {asm}")
        nuccore_to_assembly[acc] = asm
        time.sleep(SLEEP)

    # Fill Assembly column where missing
    def fill_asm(row):
        if pd.isna(row["Assembly"]) or str(row["Assembly"]).strip() == "":
            acc = str(row["Nuccore"])
            return nuccore_to_assembly.get(acc, row["Assembly"])
        else:
            return row["Assembly"]

    df["Assembly"] = df.apply(fill_asm, axis=1)

    # Save result
    df.to_excel(out_file, index=False)
    print(f"Written output with filled Assembly column to: {out_file}")


if __name__ == "__main__":
    main()
