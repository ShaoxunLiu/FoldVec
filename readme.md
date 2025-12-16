# FoldVec

This repository is designed to build a database for machine learning training that combines gene names, functional annotations, interactomes, and protein embeddings into a unified table.

**FoldVec** contains a lightweight tool that encodes proteins into fixed-length vectors (embeddings). It supports two distinct modes:

1.  **Structure Mode (Default):** Uses **ESM-IF1** to encode 3D backbone coordinates from AlphaFold predictions.
2.  **Sequence Mode:** Uses **ESM-2** to encode amino acid sequences directly.
3.  **Psi-Phi Angle Mode:** Uses Psi-Phi angles as structural property to encode AlphaFold predictions.


## 🚀 Features

* **Structure-Aware:** Uses geometric deep learning (GVP-GNN) to encode spatial coordinates.
* **Sequence-Aware:** Uses the efficient ESM-2 (8M) model for fast sequence encoding.
* **Batch Processing:** Handles lists of UniProt IDs automatically.
* **AlphaFold Integration:** Auto-fetches the latest models from AlphaFold DB.
* **Fixed Output:**
    * Structure Mode: $1 \times 512$ vector.
    * Sequence Mode: $1 \times 320$ vector.

---

## 🛠️ Installation

### 1. Clone the repository

```bash
git clone https://github.com/ShaoxunLiu/FoldVec.git
cd FoldVec
```

### 2. Set up the Environment

We strongly recommend using **Conda** to manage the dependencies (PyTorch, PyG, ESM, and jq). The provided `environment.yml` handles version mismatches (specifically for CUDA and Numpy 2.0 conflicts).

```bash
# Create the environment
conda env create -f environment.yml
# For CPU-only environment
conda env create -f environment_cpu.yml

# Activate it
conda activate esm_env
```

---

## 📊 Data Preparation

This tool is optimized for large datasets. You can download the Mammalia protein entries from UniProt (~7M entries) using the command below. This creates the `uniprot_data.tsv` input file required for processing.

```bash
wget -O uniprot_data.tsv "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cgene_names%2Csequence%2Cprotein_name%2Cgo%2Cxref_alphafolddb%2Ccc_subcellular_location%2Ckeyword%2Ccc_interaction%2Ccc_function&format=tsv&query=%28%28taxonomy_id%3A40674%29%29"
```

For a smaller Toy example, you can instead download the Mammalia disease-related protein entries from UniProt (~5k entries) using the command below. 

```bash
wget -O uniprot_data.tsv "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cgene_names%2Csequence%2Cprotein_name%2Cgo%2Cxref_alphafolddb%2Ccc_subcellular_location%2Ckeyword%2Ccc_interaction%2Ccc_function&format=tsv&query=%28%28taxonomy_id%3A40674%29%29+AND+%28proteins_with%3A20%29"
```


### Input Data Schema

The downloaded `uniprot_data.tsv` contains 9 fields:

| Field | Type | Description |
| :--- | :--- | :--- |
| **Entry** | `str` | UniProt ID (Accession). |
| **Gene Names** | `list` | Gene names as a list of strings. |
| **Sequence** | `str` | Protein sequence as a long string. |
| **Protein Names** | `list` | Protein names as a list of strings. |
| **Gene Ontology** | `list` | GO accession number and name. |
| **AlphaFoldDB Ref** | `str` | Identical to Entry if the protein has an AlphaFold prediction. |
| **Subcellular Location** | `str` | List of subcellular locations and Note as a long string. |
| **Keywords** | `list` | Tags describing the protein. |
| **Interacts with** | `list` | Interactome of the protein. |
| **Function** | `str` | Protein function as a long string. |

---

## 📖 Usage

You can use the wrapper script `generate_embeddings.sh` for all batch operations.

### 1. Structure Mode (Default)
Downloads PDBs from AlphaFold and encodes the 3D structure using ESM-IF1.

* **Flag:** `--mode structure` (or omit, as it is default)
* **Output Dimension:** 512

```bash
bash generate_embeddings.sh -i uniprot_data.tsv -o struct_results.csv --mode structure
```

### 2. Sequence Mode
Reads sequences directly from the TSV file and encodes them using ESM-2 (8M). This is significantly faster as it skips PDB downloading.

* **Flag:** `--mode sequence`
* **Output Dimension:** 320

```bash
bash generate_embeddings.sh -i uniprot_data.tsv -o seq_results.csv --mode sequence
```

### 2. Sequence Mode
Downloads PDBs from AlphaFold and encode structures using psi-phi angles (significantly faster than ESM-IF)

* **Flag:** `--mode angles`
* **Output Dimension:** 400

```bash
bash generate_embeddings.sh -i uniprot_data.tsv -o seq_results.csv --mode angles
```

### Single File Processing (Manual)
If you want to run the Python scripts directly on single files:

```bash
# Structure (PDB to Vector)
python3 pdb_to_vec.py my_structure.pdb --out vector.npy

# Sequence (TSV to CSV)
python3 seq_to_vec.py input_sequences.tsv output_embeddings.csv
```

---

## 🧠 Methodology

| Feature | Structure Mode (ESM-IF1) | Sequence Mode (ESM-2) |
| :--- | :--- | :--- |
| **Input** | N, C$\alpha$, and C backbone coordinates (PDB). | Amino acid sequence strings. |
| **Model** | Geometric Vector Perceptrons (GVP); invariant to rotation. | ESM-2 (8M Parameters) Transformer. |
| **Pooling** | Mean pooling across sequence length. | Mean pooling of the last hidden layer. |
| **Dimension** | $1 \times 512$ | $1 \times 320$ |

---

## ⚠️ Notes

> **Model Weights:** The first time you run the script, it will download model weights to a local `hub/` directory inside the project folder.

> **Numpy Compatibility:** This tool explicitly requires `numpy<2` due to binary incompatibilities with the `biotite` library used for PDB parsing.

---

## 📚 References

If you use this tool or data, please cite the following:

**ESM-IF1 (Structure) & ESM-2 (Sequence):**
* Hsu, C., Verkuil, R., Liu, J., Lin, Z., Hie, B., Sercu, T., ... & Rives, A. (2022). *Learning inverse folding from millions of predicted structures*. International Conference on Machine Learning.
* Lin, Z., Akin, H., Rao, R., Hie, B., Zhu, Z., Lu, W., ... & Rives, A. (2023). *Evolutionary-scale prediction of atomic-level protein structure with a language model*. Science, 379(6637).

**AlphaFold (Structure Predictions):**
* Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., ... & Hassabis, D. (2021). *Highly accurate protein structure prediction with AlphaFold*. Nature, 596(7873).