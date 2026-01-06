# 1. Define paths to KoFamScan database files
KO_LIST="/home/giusym/Bioinfo/DB/kofam/ko_list" #change to your paths
PROFILES="/home/giusym/Bioinfo/DB/kofam/profiles" #change to your paths

# 2. Create an output directory to keep things tidy
mkdir -p kofam_outputs

# 3. Start the loop
for faa in *.faa; do
    # Get the filename without the extension (e.g., "GFC_01" from "GFC_01.faa")
    name=$(basename "$faa" .faa)
    
    echo "Processing $name..."

    # Step A: Remove duplicate headers on the fly so KOfamScan doesn't crash
    # We save a temporary clean file
    awk '/^>/ { f = !a[$0]++ } f' "$faa" > "${name}_clean.tmp"

    # Step B: Run KOfamScan
    exec_annotation -o "kofam_outputs/${name}.txt" \
                    --ko-list "$KO_LIST" \
                    --profile "$PROFILES" \
                    --cpu 32 \ #change to CPUs for your machine
                    --format mapper \
                    "${name}_clean.tmp"

    # Step C: Remove the temporary file
    rm "${name}_clean.tmp"
    
    echo "Done with $name. Output saved to kofam_outputs/${name}.txt"
done
