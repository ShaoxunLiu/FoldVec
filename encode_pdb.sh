#!/bin/bash

# CRITICAL FIX: Prevent Python from looking at your global/local libraries
export PYTHONNOUSERSITE=1

# CRITICAL FIX: Download models to the current folder instead of ~/.cache
# This forces PyTorch/ESM to store weights in ./hub/checkpoints/ relative to this script
export TORCH_HOME="$(pwd)"

# Configuration
PDB_FILE="$1"
CHAIN_ID="$2"
OUTPUT_FILE="structure_vector.npy"

# Help Message
if [ -z "$PDB_FILE" ]; then
    echo "Usage: ./encode_pdb.sh <pdb_file> [chain_id]"
    echo "Example: ./encode_pdb.sh my_protein.pdb A"
    exit 1
fi

# 1. Check for Python
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 is not installed."
    exit 1
fi

# 2. Check for dependencies (simple check)
if ! python3 -c "import torch, esm" 2>/dev/null; then
    echo "Dependencies missing. Please activate your conda environment: conda activate esm_env"
    exit 1
fi

# 3. Run the Encoder
echo "=========================================="
echo "Encoding $PDB_FILE..."
echo "=========================================="

if [ -z "$CHAIN_ID" ]; then
    # No chain specified
    python3 pdb_to_vec.py "$PDB_FILE" --out "$OUTPUT_FILE"
else
    # Chain specified
    python3 pdb_to_vec.py "$PDB_FILE" --chain "$CHAIN_ID" --out "$OUTPUT_FILE"
fi

# 4. Report Result
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "Done! The fixed-length vector is saved in: $OUTPUT_FILE"
    echo "Vector Dimension: 512"
    echo "=========================================="
    
    # Optional: Preview the first 5 values
    echo "Preview of vector (first 5 values):"
    python3 -c "import numpy as np; v = np.load('$OUTPUT_FILE'); print(v[:5])"
else
    echo "Encoding failed."
    exit 1
fi