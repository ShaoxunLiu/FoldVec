This database is used to for machine learning training with mixed fields of gene names, functions, interactomes, and predicted protein structure data. 

Download Mammalia protein protien entries from UniProt (7M entries): 
wget -O uniprot_data.tsv "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cgene_names%2Csequence%2Cprotein_name%2Cgo%2Cxref_alphafolddb%2Ccc_subcellular_location%2Ckeyword%2Ccc_interaction&format=tsv&query=%28%28taxonomy_id%3A40674%29%29"

Downloaded data has 9 fields
    Entry(str):                         Uniprot ID
    Gene Names(StrList):                Gene names as list of strings
    Sequence(Str):                      Protein sequence as a long string
    Protein Names(StrList):             Protein names as list of strings
    Gene Ontology(StrList):             GO assession number and name as a list of strings
    AlphaFoldDB cross-reference(Str):   Identicel to Entry if the protein has AlphaFold prediction
	Subcellular Location(Str):          List of subcellular locations and Note as long string
    Keywords(StrList):                  Tags describing the protien as list of strings
    Interacts with(StrList):            Interactome of the protein as list of strings

Run 
bash AFvectorize.sh uniprot_data.tsv
To generate 512-dimention embbedding from AlphaFold predicted structure using ESM-IF1


Reference:
Hsu, C., Verkuil, R., Liu, J., Lin, Z., Hie, B., Sercu, T., ... & Rives, A. (2022, June). Learning inverse folding from millions of predicted structures. In International conference on machine learning (pp. 8946-8970). PMLR.
Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., ... & Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. nature, 596(7873), 583-589.


