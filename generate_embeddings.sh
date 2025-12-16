#!/bin/bash

# Configuration
export PYTHONNOUSERSITE=1
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
            if [ -z "$INPUT_FILE" ]; then INPUT_FILE="$1";
            elif [ -z "$OUTPUT_FILE" ]; then OUTPUT_FILE="$1"; fi
            ;;
    esac
    shift
done

# Help / Validation
if [ -z "$INPUT_FILE" ]; then
    echo "Usage: ./generate_embeddings.sh -i <input_tsv> [-o <output_csv>] [-m <mode>]"
    echo "Modes:"
    echo "  --mode structure : ESM-IF (Slow, 512 dim)"
    echo "  --mode angles    : Phi/Psi Torsion (Fast, 400 dim, starts at Psi_0)"
    echo "  --mode sequence  : ESM-2 (Fast, 320 dim)"
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

# ==============================================================================
# MODE: STRUCTURE (Original ESM-IF)
# ==============================================================================
elif [ "$MODE" == "structure" ]; then
    echo "[INFO] Using ESM-IF1 (Structure)..."
    if ! command -v jq &> /dev/null; then echo "Error: 'jq' missing."; exit 1; fi
    mkdir -p "$TEMP_DIR"
    
    HEADER="UniProt_ID"
    for i in {0..511}; do HEADER="$HEADER,Dim_$i"; done
    echo "$HEADER" > "$OUTPUT_FILE"

    while read -r line || [ -n "$line" ]; do
        if [ -z "$line" ] || [[ "$line" == "Entry"* ]]; then continue; fi
        uniprot_id=$(echo "$line" | awk '{print $1}')
        
        # Download
        api_url="https://alphafold.ebi.ac.uk/api/prediction/$uniprot_id"
        pdb_url=$(curl -s "$api_url" | jq -r '.[0].pdbUrl')
        if [ "$pdb_url" == "null" ] || [ -z "$pdb_url" ]; then continue; fi
        
        pdb_file="$TEMP_DIR/${uniprot_id}.pdb"
        curl -s -o "$pdb_file" "$pdb_url"
        
        # Encode
        npy_file="$TEMP_DIR/${uniprot_id}.npy"
        python3 pdb_to_vec.py "$pdb_file" --out "$npy_file" > /dev/null 2>&1
        
        if [ -f "$npy_file" ]; then
            vec=$(python3 -c "import numpy as np; v = np.load('$npy_file'); print(','.join(map(str, v)))")
            echo "$uniprot_id,$vec" >> "$OUTPUT_FILE"
            echo "  [OK] $uniprot_id"
            rm "$npy_file"
        fi
        rm "$pdb_file"
    done < "$INPUT_FILE"
    rmdir "$TEMP_DIR" 2>/dev/null

# ==============================================================================
# MODE: ANGLES (Phi/Psi Sequence, 400 dim)
# ==============================================================================
elif [ "$MODE" == "angles" ]; then
    echo "[INFO] Using Phi/Psi Torsion Angles (400 dimensions)..."
    echo "[NOTE] Removing Phi_0 column (undefined). Starting at Psi_0."
    
    if ! command -v jq &> /dev/null; then echo "Error: 'jq' required."; exit 1; fi
    mkdir -p "$TEMP_DIR"

    # 1. Generate Header: Psi_0, Phi_1, Psi_1 ... 
    # Logic: We generate 400 labels alternating.
    # The first valid angle is Psi_0. Then Phi_1, Psi_1, Phi_2...
    HEADER=$(python3 -c '
header = "UniProt_ID"
# We want 400 dimensions. 
# Sequence of valid angles: Psi0, Phi1, Psi1, Phi2, Psi2...
for i in range(400):
    residue_idx = (i + 1) // 2
    if i % 2 == 0:
        label = f"Psi_{residue_idx}"  # Evens (0, 2..): Psi0, Psi1...
    else:
        label = f"Phi_{residue_idx}"  # Odds (1, 3..): Phi1, Phi2...
    header += "," + label
print(header)
    ')
    echo "$HEADER" > "$OUTPUT_FILE"

    while read -r line || [ -n "$line" ]; do
        if [ -z "$line" ] || [[ "$line" == "Entry"* ]]; then continue; fi
        uniprot_id=$(echo "$line" | awk '{print $1}')
        
        # Download PDB
        pdb_url="https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-model_v4.pdb"
        pdb_file="$TEMP_DIR/${uniprot_id}.pdb"
        
        if ! curl -s -f -o "$pdb_file" "$pdb_url"; then
            api_url="https://alphafold.ebi.ac.uk/api/prediction/$uniprot_id"
            real_pdb_url=$(curl -s "$api_url" | jq -r '.[0].pdbUrl')
            if [ -z "$real_pdb_url" ] || [ "$real_pdb_url" == "null" ]; then
                echo "  [WARN] No structure found for $uniprot_id"
                continue
            fi
            curl -s -o "$pdb_file" "$real_pdb_url"
        fi
        
        # 2. CALC ANGLES (Pure AWK)
        vector_csv=$(awk '
        function torsion(ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz) {
            v1x=bx-ax; v1y=by-ay; v1z=bz-az;
            v2x=cx-bx; v2y=cy-by; v2z=cz-bz;
            v3x=dx-cx; v3y=dy-cy; v3z=dz-cz;
            
            n1x=v1y*v2z-v1z*v2y; n1y=v1z*v2x-v1x*v2z; n1z=v1x*v2y-v1y*v2x;
            n2x=v2y*v3z-v2z*v3y; n2y=v2z*v3x-v2x*v3z; n2z=v2x*v3y-v2y*v3x;
            
            ln1=sqrt(n1x^2+n1y^2+n1z^2); ln2=sqrt(n2x^2+n2y^2+n2z^2);
            if(ln1==0||ln2==0) return 0; 

            x = n1x*n2x + n1y*n2y + n1z*n2z;
            c1x=n1y*n2z-n1z*n2y; c1y=n1z*n2x-n1x*n2z; c1z=n1x*n2y-n1y*n2x;
            
            lv2=sqrt(v2x^2+v2y^2+v2z^2);
            if(lv2==0) return 0;
            
            y = c1x*v2x + c1y*v2y + c1z*v2z;
            
            rads = atan2(y/lv2, x);
            deg = rads * 180 / 3.14159265;
            return deg;
        }

        BEGIN {
            prev_res_idx = -1;
            TARGET_DIMS = 400; 
        }

        /^ATOM/ && ($3 == "N" || $3 == "CA" || $3 == "C") {
            res_id = substr($0, 23, 4) + 0;
            atom_type = $3;
            x = substr($0, 31, 8) + 0; y = substr($0, 39, 8) + 0; z = substr($0, 47, 8) + 0;
            
            atoms[res_id, atom_type, "x"] = x;
            atoms[res_id, atom_type, "y"] = y;
            atoms[res_id, atom_type, "z"] = z;
            atoms[res_id, atom_type, "seen"] = 1;
            if (res_id > max_res_id) max_res_id = res_id;
            count++;
        }
        
        END {
            if (count == 0) { print "ERR"; exit; }

            # 1. Compute Raw Angles List
            val_count = 0;
            
            for (i=1; i<=max_res_id; i++) {
                if (!atoms[i, "CA", "seen"]) continue;

                # PHI (Skip for i=1, as it is always 0)
                # We calculate it, but we only store it if i > 1
                if (i > 1) {
                    phi = 0.0;
                    if (atoms[i-1, "C", "seen"] && atoms[i, "N", "seen"] && atoms[i, "CA", "seen"] && atoms[i, "C", "seen"]) {
                        phi = torsion( \
                            atoms[i-1,"C","x"], atoms[i-1,"C","y"], atoms[i-1,"C","z"], \
                            atoms[i,"N","x"],   atoms[i,"N","y"],   atoms[i,"N","z"],   \
                            atoms[i,"CA","x"],  atoms[i,"CA","y"],  atoms[i,"CA","z"],  \
                            atoms[i,"C","x"],   atoms[i,"C","y"],   atoms[i,"C","z"] );
                    }
                    vals[val_count++] = phi;
                }

                # PSI
                psi = 0.0;
                if (atoms[i, "N", "seen"] && atoms[i, "CA", "seen"] && atoms[i, "C", "seen"] && atoms[i+1, "N", "seen"]) {
                    psi = torsion( \
                        atoms[i,"N","x"],   atoms[i,"N","y"],   atoms[i,"N","z"],   \
                        atoms[i,"CA","x"],  atoms[i,"CA","y"],  atoms[i,"CA","z"],  \
                        atoms[i,"C","x"],   atoms[i,"C","y"],   atoms[i,"C","z"],   \
                        atoms[i+1,"N","x"], atoms[i+1,"N","y"], atoms[i+1,"N","z"] );
                }
                vals[val_count++] = psi;
            }

            # 2. Pad (Concatenate self) or Truncate
            out_str = "";
            for (k=0; k<TARGET_DIMS; k++) {
                if (val_count == 0) val = 0; 
                else val = vals[k % val_count];
                
                if (k == 0) out_str = sprintf("%.4f", val);
                else out_str = out_str sprintf(",%.4f", val);
            }
            
            print out_str;
        }
        ' "$pdb_file")
        
        if [[ "$vector_csv" == "ERR" ]]; then
             echo "  [WARN] Failed $uniprot_id"
        else
             echo "$uniprot_id,$vector_csv" >> "$OUTPUT_FILE"
             echo "  [OK] $uniprot_id"
        fi
        rm "$pdb_file"

    done < "$INPUT_FILE"
    rmdir "$TEMP_DIR" 2>/dev/null
fi