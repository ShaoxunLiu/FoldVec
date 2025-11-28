#!/bin/bash

# Check if jq is installed
if ! command -v jq &> /dev/null; then
    echo "Error: 'jq' is not installed. Please install it (e.g., sudo apt install jq / brew install jq)."
    exit 1
fi

# Check if at least one ID is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <UniProtID_1> <UniProtID_2> ..."
    echo "Example: $0 Q04228 P53"
    exit 1
fi

# Base API URL
API_URL="https://alphafold.ebi.ac.uk/api/prediction"

# Loop through all provided UniProt IDs
for UNIPROT_ID in "$@"; do
    echo "Processing $UNIPROT_ID..."

    # Fetch metadata from AlphaFold API
    RESPONSE=$(curl -s "${API_URL}/${UNIPROT_ID}")

    # Check if the response is an empty list (ID not found)
    if [[ "$RESPONSE" == "[]" ]]; then
        echo "  Warning: No AlphaFold structure found for ID: $UNIPROT_ID"
        continue
    fi

    # Extract PDB URLs using jq
    # The API returns an array because large proteins are split into fragments. 
    # This loop handles single files AND multi-fragment files.
    PDB_URLS=$(echo "$RESPONSE" | jq -r '.[].pdbUrl')

    if [ -z "$PDB_URLS" ] || [ "$PDB_URLS" == "null" ]; then
        echo "  Error: Could not extract URL for $UNIPROT_ID"
        continue
    fi

    # Download each URL found
    for URL in $PDB_URLS; do
        FILENAME=$(basename "$URL")
        echo "  Downloading $FILENAME..."
        curl -s -O "$URL"
    done
done

echo "Done!"