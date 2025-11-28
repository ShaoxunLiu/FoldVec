#!/bin/bash

# Configuration
export PYTHONNOUSERSITE=1
export TORCH_HOME="$(pwd)"

# Defaults
MODE="structure"
INPUT_FILE=""
OUTPUT_FILE="embeddings_table.csv"
TEMP_DIR="temp_pdbs"

# Parse Flags
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -m|--mode) MODE="$2"; shift ;;
        -i|--input) INPUT_FILE="$2"; shift ;;
        -o|--output) OUTPUT_FILE="$2"; shift ;;
        *) 
            # Assume positional arguments if no flags used: ./script.sh input output
            if [ -z "$INPUT_FILE" ]; then INPUT_FILE="$1";
            elif [ -z "$OUTPUT_FILE" ]; then OUTPUT_FILE="$1"; fi
            ;;
    esac
    shift
done

# Help / Validation
if [ -z "$INPUT_FILE" ]; then
    echo "Usage: ./generate_embeddings.sh -i <input_tsv> [-o <output_csv>] [-m structure|sequence]"
    echo "  --mode structure : Download AF PDBs and encode structure (Default, 512 dim)"
    echo "  --mode sequence  : Read sequences from TSV and encode (Fast, 320 dim)"
    exit 1
fi

echo "=========================================="
echo "Starting Batch Processing"
echo "Mode:   $MODE"
echo "Input:  $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo "=========================================="

# ==============================================================================
# MODE: SEQUENCE (ESM-2)
# ==============================================================================
if [ "$MODE" == "sequence" ]; then
    echo "[INFO] Using ESM-2 (Sequence only)..."
    python3 seq_to_vec.py "$INPUT_FILE" "$OUTPUT_FILE"
    
    if [ $? -eq 0 ]; then
        echo "=========================================="
        echo "Sequence processing complete."
    else
        echo "[FAIL] Sequence processing failed."
        exit 1
    fi

# ==============================================================================
# MODE: STRUCTURE (ESM-IF1 + AlphaFold)
# ==============================================================================
elif [ "$MODE" == "structure" ]; then
    echo "[INFO] Using ESM-IF1 (Structure / Inverse Folding)..."
    
    # Check for jq
    if ! command -v jq &> /dev/null; then
        echo "Error: 'jq' is not installed. Please install it."
        exit 1
    fi

    # Prepare Headers and Temp Dir
    mkdir -p "$TEMP_DIR"
    echo "Generating CSV header (512 dimensions)..."
    HEADER="UniProt_ID"
    for i in {0..511}; do HEADER="$HEADER,Dim_$i"; done
    echo "$HEADER" > "$OUTPUT_FILE"

    # Loop through TSV
    while read -r line || [ -n "$line" ]; do
        # Skip empty or header lines (naive check for 'Entry')
        if [ -z "$line" ] || [[ "$line" == "Entry"* ]]; then continue; fi
        
        # Extract first column (UniProt ID)
        uniprot_id=$(echo "$line" | awk '{print $1}')
        echo "Processing: $uniprot_id"
        
        # 1. Download PDB
        api_url="https://alphafold.ebi.ac.uk/api/prediction/$uniprot_id"
        pdb_url=$(curl -s "$api_url" | jq -r '.[0].pdbUrl')
        
        if [ "$pdb_url" == "null" ] || [ -z "$pdb_url" ]; then
            echo "  [WARN] No structure found for $uniprot_id"
            continue
        fi
        
        pdb_filename="$TEMP_DIR/${uniprot_id}.pdb"
        curl -s -o "$pdb_filename" "$pdb_url"
        
        # 2. Encode
        npy_filename="$TEMP_DIR/${uniprot_id}.npy"
        python3 pdb_to_vec.py "$pdb_filename" --out "$npy_filename" > /dev/null 2>&1
        
        if [ -f "$npy_filename" ]; then
            vector_string=$(python3 -c "import numpy as np; v = np.load('$npy_filename'); print(','.join(map(str, v)))")
            echo "$uniprot_id,$vector_string" >> "$OUTPUT_FILE"
            echo "  [OK] Embedded."
            rm "$npy_filename"
        else
            echo "  [FAIL] Encoding failed."
        fi
        rm "$pdb_filename"

    done < "$INPUT_FILE"
    
    rmdir "$TEMP_DIR"
    echo "=========================================="
    echo "Structure processing complete."

else
    echo "Error: Unknown mode '$MODE'. Use 'structure' or 'sequence'."
    exit 1
fi