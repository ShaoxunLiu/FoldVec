import argparse
import torch
import esm
import esm.inverse_folding.util as if_util
import sys
import numpy as np
import warnings

# Filter out the specific ESM warning about regression weights.
# This is safe because we are generating embeddings, not predicting contacts.
warnings.filterwarnings("ignore", message="Regression weights not found")

def get_structure_embedding(pdb_path, chain_id=None):
    """
    Loads a PDB, extracts coordinates, and runs them through the ESM-IF
    encoder to get a fixed-length structural embedding.
    """
    # 1. Load the ESM-IF model (Structure-aware)
    # This downloads the model to ~/.cache/torch/hub/ if not present
    print(f"[INFO] Loading ESM-IF1 model...", file=sys.stderr)
    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model.eval()
    
    if torch.cuda.is_available():
        model = model.cuda()
        print("[INFO] Using GPU.", file=sys.stderr)
    else:
        print("[INFO] Using CPU (this might be slow for large proteins).", file=sys.stderr)

    # 2. Parse the PDB file structure
    print(f"[INFO] Parsing PDB: {pdb_path}", file=sys.stderr)
    structure = if_util.load_structure(pdb_path, chain_id)
    
    # Extract coordinates (N, Ca, C) and sequence
    # coords shape: [Length, 3, 3]
    coords, native_seq = if_util.extract_coords_from_structure(structure)
    
    # 3. Prepare inputs for the model
    # The model expects a batch dimension, so we unsqueeze
    coords = torch.tensor(coords)
    if torch.cuda.is_available():
        coords = coords.cuda()
        
    # ESM-IF expects: coords, padding_mask, confidence
    # We treat it as a single batch of size 1
    # padding_mask is all False (no padding)
    # confidence is 1.0 (we trust the crystal structure)
    batch_coords = coords.unsqueeze(0) # Shape: [1, L, 3, 3]
    padding_mask = torch.zeros(1, len(native_seq), dtype=torch.bool)
    confidence = torch.ones(1, len(native_seq), dtype=torch.float)
    
    if torch.cuda.is_available():
        padding_mask = padding_mask.cuda()
        confidence = confidence.cuda()

    # 4. Run the Encoder
    # We only run the 'encoder', not the decoder (which predicts sequence)
    # The encoder outputs the structural representation.
    print("[INFO] Encoding structure...", file=sys.stderr)
    with torch.no_grad():
        encoder_out = model.encoder(batch_coords, padding_mask, confidence)
    
    # encoder_out['encoder_out'][0] has shape [Length, Batch, HiddenDim]
    # ESM-IF Hidden Dimension is typically 512
    residue_embeddings = encoder_out['encoder_out'][0].squeeze(1) # [Length, 512]
    
    # 5. Pooling (Convert L x 512 -> 1 x 512)
    # We take the mean across the sequence length to get a fixed vector
    fixed_vector = residue_embeddings.mean(dim=0).cpu().numpy()
    
    return fixed_vector

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Encode PDB structure to fixed vector")
    parser.add_argument("pdb_file", help="Path to the .pdb file")
    parser.add_argument("--chain", help="Chain ID to extract (optional, defaults to first chain)", default=None)
    parser.add_argument("--out", help="Output file for the vector (npy format). If not set, prints to stdout.", default=None)
    
    args = parser.parse_args()
    
    try:
        vector = get_structure_embedding(args.pdb_file, args.chain)
        
        if args.out:
            np.save(args.out, vector)
            print(f"[SUCCESS] Vector saved to {args.out}", file=sys.stderr)
        else:
            # Print as comma-separated values for easy bash parsing
            print(",".join(map(str, vector)))
            
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)