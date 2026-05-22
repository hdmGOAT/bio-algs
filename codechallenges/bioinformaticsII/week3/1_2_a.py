from pathlib import Path

from bioalgs.translation import load_codon_table, translateRNAtoAminoAcidstr


p = Path(__file__).parent / "datasets" / "dataset_30213_4.txt"

table = load_codon_table
s = p.read_text()
with open("output.txt", "w") as f:
    f.write(translateRNAtoAminoAcidstr(s, table))
