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
            # Assume positional arguments if no flags used
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
    echo "  --mode structure : (Original) Download AF PDBs and use ESM-IF (Slow, 512 dim)"
    echo "  --mode 2resi     : (New) Download AF PDBs and calc Di-peptides + Geom (Fast, 405 dim)"
    echo "  --mode sequence  : Read sequences from TSV and use ESM-2 (Fast, 320 dim)"
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
        echo "Sequence processing complete."
    else
        echo "[FAIL] Sequence processing failed."
        exit 1
    fi

# ==============================================================================
# MODE: STRUCTURE (Original ESM-IF)
# ==============================================================================
elif [ "$MODE" == "structure" ]; then
    echo "[INFO] Using ESM-IF1 (Structure / Inverse Folding)..."
    
    if ! command -v jq &> /dev/null; then echo "Error: 'jq' is not installed."; exit 1; fi

    mkdir -p "$TEMP_DIR"
    
    # Header for ESM-IF (512 dims)
    HEADER="UniProt_ID"
    for i in {0..511}; do HEADER="$HEADER,Dim_$i"; done
    echo "$HEADER" > "$OUTPUT_FILE"

    while read -r line || [ -n "$line" ]; do
        if [ -z "$line" ] || [[ "$line" == "Entry"* ]]; then continue; fi
        
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
        
        # 2. Encode with ESM-IF Python script
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
    echo "Structure processing complete."

# ==============================================================================
# MODE: 2RESI (Fast Di-peptide + Geometry)
# ==============================================================================
elif [ "$MODE" == "2resi" ]; then
    echo "[INFO] Using 2-Residue Fingerprint (405 dimensions)..."
    
    if ! command -v jq &> /dev/null; then echo "Error: 'jq' is not installed."; exit 1; fi

    mkdir -p "$TEMP_DIR"

    # 1. Generate Header (400 Di-peptides + 5 Geom)
    # Using python purely to generate the CSV header string cleanly once
    HEADER=$(python3 -c '
aas = "ACDEFGHIKLMNPQRSTVWY"
header = "UniProt_ID"
for i in aas:
    for j in aas:
        header += f",{i}{j}"
header += ",Len,Rg,BoxDiag,Hydro,Charge"
print(header)
    ')
    echo "$HEADER" > "$OUTPUT_FILE"

    while read -r line || [ -n "$line" ]; do
        if [ -z "$line" ] || [[ "$line" == "Entry"* ]]; then continue; fi
        
        uniprot_id=$(echo "$line" | awk '{print $1}')
        
        # Download PDB (Fast v4 direct link first, then API fallback)
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
        
        # 2. CALC 405-DIM FINGERPRINT (Pure AWK)
        vector_csv=$(awk '
        BEGIN {
            # Hydrophobicity map
            split("1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3", h_vals, ",");
            split("A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y", aa_order, ",");
            for(i=1;i<=20;i++) { h_map[aa_order[i]] = h_vals[i]; aa_idx[aa_order[i]] = i; }
            
            # Charge map
            charge_map["K"]=1; charge_map["R"]=1; charge_map["D"]=-1; charge_map["E"]=-1; charge_map["H"]=0.1;

            min_x=9999; max_x=-9999; min_y=9999; max_y=-9999; min_z=9999; max_z=-9999;
            prev_res = "";
        }
        /^ATOM/ && $3 == "CA" { 
            # Get AA
            res = substr($0, 18, 3);
            if (res == "ALA") c="A"; else if (res == "CYS") c="C"; else if (res == "ASP") c="D";
            else if (res == "GLU") c="E"; else if (res == "PHE") c="F"; else if (res == "GLY") c="G";
            else if (res == "HIS") c="H"; else if (res == "ILE") c="I"; else if (res == "LYS") c="K";
            else if (res == "LEU") c="L"; else if (res == "MET") c="M"; else if (res == "ASN") c="N";
            else if (res == "PRO") c="P"; else if (res == "GLN") c="Q"; else if (res == "ARG") c="R";
            else if (res == "SER") c="S"; else if (res == "THR") c="T"; else if (res == "VAL") c="V";
            else if (res == "TRP") c="W"; else if (res == "TYR") c="Y"; else c="X";

            if (c == "X") next; # Skip unknown

            # Geometry
            x = substr($0, 31, 8) + 0; y = substr($0, 39, 8) + 0; z = substr($0, 47, 8) + 0;
            sum_x += x; sum_y += y; sum_z += z;
            xs[count+1]=x; ys[count+1]=y; zs[count+1]=z;
            if(x<min_x) min_x=x; if(x>max_x) max_x=x;
            if(y<min_y) min_y=y; if(y>max_y) max_y=y;
            if(z<min_z) min_z=z; if(z>max_z) max_z=z;

            # Properties
            total_hydro += h_map[c];
            total_charge += charge_map[c];

            # Di-peptide counting
            if (count > 0) {
                pair = prev_res c;
                pair_counts[pair]++;
            }
            prev_res = c;
            count++;
        }
        END {
            if (count == 0) { print "ERR"; exit; }
            
            # 1. Output 400 Di-peptide frequencies (row-major order A-Y)
            out = "";
            # Normalized by (count - 1) pairs
            norm = (count > 1) ? (count - 1) : 1;
            
            for(i=1; i<=20; i++) {
                for(j=1; j<=20; j++) {
                    key = aa_order[i] aa_order[j];
                    val = (pair_counts[key] / norm);
                    out = out sprintf("%.4f,", val);
                }
            }
            
            # 2. Geometry
            cx = sum_x / count; cy = sum_y / count; cz = sum_z / count;
            sq_dist_sum = 0;
            for(i=1; i<=count; i++) {
                sq_dist_sum += (xs[i]-cx)^2 + (ys[i]-cy)^2 + (zs[i]-cz)^2;
            }
            rg = sqrt(sq_dist_sum / count);
            box_diag = sqrt((max_x-min_x)^2 + (max_y-min_y)^2 + (max_z-min_z)^2);
            hydro = total_hydro / count;
            
            # Output final string
            printf "%s%d,%.4f,%.4f,%.4f,%.4f", out, count, rg, box_diag, hydro, total_charge;
        }
        ' "$pdb_file")

        if [[ "$vector_csv" == "ERR" ]]; then
             echo "  [WARN] Failed to parse structure for $uniprot_id"
        else
             echo "$uniprot_id,$vector_csv" >> "$OUTPUT_FILE"
             echo "  [OK] $uniprot_id"
        fi
        
        rm "$pdb_file"

    done < "$INPUT_FILE"
    
    rmdir "$TEMP_DIR" 2>/dev/null
    echo "2resi processing complete."

else
    echo "Error: Unknown mode '$MODE'. Use 'structure', '2resi', or 'sequence'."
    exit 1
fi