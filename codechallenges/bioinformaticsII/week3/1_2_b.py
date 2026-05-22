from pathlib import Path
from bioalgs.translation import findSubstringsThatEncodeAA, load_codon_table

# point to the actual dataset file
p = Path(__file__).parent / "datasets" / "dataset_30213_7.txt"

dna, peptide = p.read_text().splitlines()
table = load_codon_table()
results = findSubstringsThatEncodeAA(dna, peptide, table)

with open("output.txt", "w") as f:
    f.writelines(s + "\n" for s in results)
