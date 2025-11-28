FoldVec

FoldVec is a lightweight tool that encodes protein 3D structures into fixed-length vectors (embeddings) using the ESM-IF1 (Inverse Folding) model.

This repository is designed to build a database for machine learning training that combines gene names, functional annotations, interactomes, and predicted protein structure data into a unified format.

🚀 Features

Structure-Aware: Uses geometric deep learning (GVP-GNN) to encode spatial coordinates, not just sequence.

Batch Processing: Handles lists of UniProt IDs automatically.

AlphaFold Integration: Auto-fetches the latest models from AlphaFold DB.

Fixed Output: Maps proteins of any length to a consistent $1 \times 512$ vector.

🛠️ Installation

1. Clone the repository

git clone [https://github.com/ShaoxunLiu/FoldVec.git](https://github.com/ShaoxunLiu/FoldVec.git)
cd FoldVec


2. Set up the Environment

We strongly recommend using Conda to manage the dependencies (PyTorch, PyG, ESM, and jq). The provided environment.yml handles version mismatches (specifically for CUDA and Numpy 2.0 conflicts).

# Create the environment
conda env create -f environment.yml

# Activate it
conda activate esm_env


📊 Data Preparation

This tool is optimized for large datasets. You can download the Mammalia protein entries from UniProt (~7M entries) using the command below. This creates the uniprot_data.tsv input file required for processing.

wget -O uniprot_data.tsv "[https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cgene_names%2Csequence%2Cprotein_name%2Cgo%2Cxref_alphafolddb%2Ccc_subcellular_location%2Ckeyword%2Ccc_interaction&format=tsv&query=%28%28taxonomy_id%3A40674%29%29](https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cgene_names%2Csequence%2Cprotein_name%2Cgo%2Cxref_alphafolddb%2Ccc_subcellular_location%2Ckeyword%2Ccc_interaction&format=tsv&query=%28%28taxonomy_id%3A40674%29%29)"


Input Data Schema

The downloaded uniprot_data.tsv contains 9 fields:

Entry (str): UniProt ID (Accession).

Gene Names (StrList): Gene names as a list of strings.

Sequence (Str): Protein sequence as a long string.

Protein Names (StrList): Protein names as a list of strings.

Gene Ontology (StrList): GO accession number and name as a list of strings.

AlphaFoldDB cross-reference (Str): Identical to Entry if the protein has an AlphaFold prediction.

Subcellular Location (Str): List of subcellular locations and Note as long string.

Keywords (StrList): Tags describing the protein as a list of strings.

Interacts with (StrList): Interactome of the protein as a list of strings.

📖 Usage

Batch Processing

To generate the 512-dimension structure embeddings from the downloaded TSV, run the batch processing script (referred to as AFvectorize.sh or generate_embeddings.sh depending on your local filename).

# Make sure the script is executable
chmod +x generate_embeddings.sh

# Run the pipeline
# This reads 'uniprot_data.tsv' and outputs 'embeddings_table.csv'
./generate_embeddings.sh uniprot_data.tsv embeddings_table.csv


Output:
The script generates a CSV file where the first column is the UniProt ID, followed by the 512 dimensions of the vector (Dim_0 ... Dim_511).

Single PDB Processing

If you have a local .pdb file and want to process it individually:

python3 pdb_to_vec.py my_structure.pdb --out vector.npy


🧠 Methodology

This tool uses the encoder from ESM-IF1 (Evolutionary Scale Modeling - Inverse Folding).

Input: The tool reads the N, C$\alpha$, and C backbone coordinates from the PDB file.

Encoding: It passes these coordinates through Geometric Vector Perceptrons (GVP) layers, which are invariant to rotation and translation.

Pooling: The model outputs a vector for every residue ($L \times 512$). We perform Mean Pooling across the sequence length to derive a single global structure fingerprint ($1 \times 512$).

⚠️ Notes

Model Weights: The first time you run the script, it will download the ESM-IF1 model weights (~142MB) to a local hub/ directory inside the project folder.

Numpy Compatibility: This tool explicitly requires numpy<2 due to binary incompatibilities with the biotite library used for PDB parsing. The environment.yml handles this automatically.

📚 References

ESM-IF1 (Inverse Folding Model):

Hsu, C., Verkuil, R., Liu, J., Lin, Z., Hie, B., Sercu, T., ... & Rives, A. (2022). Learning inverse folding from millions of predicted structures. International Conference on Machine Learning (pp. 8946-8970). PMLR.

AlphaFold (Structure Predictions):

Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., ... & Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596(7873), 583-589.