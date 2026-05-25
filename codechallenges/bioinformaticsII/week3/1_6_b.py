from pathlib import Path

from bioalgs.translation import cyclopeptideSequencing, load_integer_mass_table

p = Path(__file__).parent / "datasets" / "dataset_30217_6.txt"

spectrum = list(map(int, p.read_text().split()))
table = load_integer_mass_table()

peptides = cyclopeptideSequencing(spectrum, table)

out = []
for peptide in peptides:
    masses = [str(table[aa]) for aa in peptide]
    out.append("-".join(masses))

result = " ".join(out)

with open("output.txt", "w") as f:
    f.write(result)
