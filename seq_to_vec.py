import argparse
import torch
import esm
import sys
import csv
import os

def process_tsv(input_file, output_file, batch_size=1):
    # 1. Load ESM-2 Model (8M parameters - Smallest)
    # Output dimension: 320
    model_name = "esm2_t6_8M_UR50D"
    print(f"[INFO] Loading {model_name}...", file=sys.stderr)
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    if torch.cuda.is_available():
        model = model.cuda()
        print("[INFO] Using GPU.", file=sys.stderr)
    else:
        print("[INFO] Using CPU.", file=sys.stderr)

    # 2. Open Files
    print(f"[INFO] Processing {input_file}...", file=sys.stderr)
    
    with open(input_file, 'r') as f_in, open(output_file, 'w', newline='') as f_out:
        # Detect TSV format
        reader = csv.DictReader(f_in, delimiter='\t')
        
        # specific handling for UniProt TSV headers
        # We look for 'Sequence' and 'Entry' (or 'Entry name')
        headers = reader.fieldnames
        seq_col = next((h for h in headers if 'Sequence' in h), None)
        id_col = next((h for h in headers if 'Entry' in h), None)

        if not seq_col or not id_col:
            print(f"[ERROR] Could not find 'Entry' or 'Sequence' columns in TSV. Found: {headers}", file=sys.stderr)
            sys.exit(1)

        # Write CSV Header
        # ESM-2 8M has 320 dimensions
        embed_dim = 320
        csv_header = ["UniProt_ID"] + [f"Dim_{i}" for i in range(embed_dim)]
        writer = csv.writer(f_out)
        writer.writerow(csv_header)

        # 3. Iterate and Process
        data = []
        
        for i, row in enumerate(reader):
            seq = row[seq_col]
            entry_id = row[id_col]
            
            # Simple validation
            if len(seq) == 0: continue
            
            # Truncate to 1022 to avoid OOM on standard GPUs (1024 limit with BOS/EOS)
            if len(seq) > 1022:
                seq = seq[:1022]

            data.append((entry_id, seq))

            # Process batch (currently batch_size=1 for safety, can increase)
            if len(data) >= batch_size:
                process_batch(model, batch_converter, data, writer)
                data = []
                
                if (i + 1) % 100 == 0:
                    print(f"[INFO] Processed {i + 1} sequences...", file=sys.stderr)

        # Process remaining
        if data:
            process_batch(model, batch_converter, data, writer)

    print(f"[SUCCESS] Saved sequence embeddings to {output_file}", file=sys.stderr)

def process_batch(model, batch_converter, data, writer):
    # data is list of (id, seq)
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    
    if torch.cuda.is_available():
        batch_tokens = batch_tokens.cuda()

    # The 8M parameter model has 6 layers. We extract the last layer (index 6).
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[6], return_contacts=False)
    
    # Extract the representation from the last layer (6 for this model)
    token_representations = results["representations"][6]

    # Generate per-sequence representations via averaging
    # NOTE: token 0 is BOS, so start at 1
    for i, (_, seq) in enumerate(data):
        # We slice from 1 : len(seq) + 1 to exclude BOS and EOS/Padding
        # Use batch_lens logic if batching multiple lengths
        seq_len = len(seq)
        embedding = token_representations[i, 1 : seq_len + 1].mean(0)
        
        # Write to CSV
        vector_list = embedding.cpu().numpy().tolist()
        writer.writerow([data[i][0]] + vector_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Encode Protein Sequences to Vectors using ESM-2")
    parser.add_argument("input_tsv", help="Input TSV file")
    parser.add_argument("output_csv", help="Output CSV file")
    
    args = parser.parse_args()
    process_tsv(args.input_tsv, args.output_csv)