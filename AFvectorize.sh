#!/bin/bash

# Configuration
INPUT_FILE="$1"
OUTPUT_FILE="${2:-embeddings_table.csv}"
TEMP_DIR="temp_pdbs"

# 1. Environment Safety Checks
export PYTHONNOUSERSITE=1
# Download models to the current folder instead of ~/.cache
export TORCH_HOME="$(pwd)"

if [ -z "$INPUT_FILE" ]; then
    echo "Usage: ./generate_embeddings.sh <input_tsv> [output_csv]"
    echo "Example: ./generate_embeddings.sh uniprot_ids.tsv results.csv"
    exit 1
fi

if ! command -v jq &> /dev/null; then
    echo "Error: 'jq' is not installed. Please update your environment."
    exit 1
fi

# 2. Prepare Output & Temp Directory
mkdir -p "$TEMP_DIR"
# Create CSV Header: UniProt_ID followed by Dim_0 to Dim_511
echo "Generating header..."
HEADER="UniProt_ID"
for i in {0..511}; do HEADER="$HEADER,Dim_$i"; done
echo "$HEADER" > "$OUTPUT_FILE"

echo "=========================================="
echo "Starting Batch Processing"
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo "=========================================="

# 3. Process each line in the TSV
while read -r line || [ -n "$line" ]; do
    # Skip empty lines
    if [ -z "$line" ]; then continue; fi
    
    # Extract first column (UniProt ID)
    uniprot_id=$(echo "$line" | awk '{print $1}')
    
    echo "Processing: $uniprot_id"
    
    # --- STEP A: Download from AlphaFold DB ---
    api_url="https://alphafold.ebi.ac.uk/api/prediction/$uniprot_id"
    pdb_url=$(curl -s "$api_url" | jq -r '.[0].pdbUrl')
    
    if [ "$pdb_url" == "null" ] || [ -z "$pdb_url" ]; then
        echo "  [WARN] No structure found for $uniprot_id on AlphaFold DB."
        continue
    fi
    
    pdb_filename="$TEMP_DIR/${uniprot_id}.pdb"
    curl -s -o "$pdb_filename" "$pdb_url"
    
    # --- STEP B: Encode using pdb_to_vec.py ---
    # We use a temporary .npy file for robustness, then clean it up
    npy_filename="$TEMP_DIR/${uniprot_id}.npy"
    
    # Run the python script, sending logs to stderr and saving output to file
    python3 pdb_to_vec.py "$pdb_filename" --out "$npy_filename" > /dev/null 2>&1
    
    # Check if the .npy file was created successfully
    if [ -f "$npy_filename" ]; then
        # Parse the NPY file into a CSV string using a python one-liner
        vector_string=$(python3 -c "import numpy as np; v = np.load('$npy_filename'); print(','.join(map(str, v)))")
        
        echo "$uniprot_id,$vector_string" >> "$OUTPUT_FILE"
        echo "  [OK] Embedded."
        
        # --- CLEANUP NPY FILE ---
        rm "$npy_filename"
    else
        echo "  [FAIL] Encoding failed."
    fi
    
    # --- CLEANUP PDB FILE ---
    rm "$pdb_filename"

done < "$INPUT_FILE"

# 4. Cleanup Temp Directory
rmdir "$TEMP_DIR"
echo "=========================================="
echo "Batch processing complete."
echo "Results saved to $OUTPUT_FILE"