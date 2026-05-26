from pathlib import Path

from bioalgs.translation import cyclopeptideSequencingMass, load_integer_mass_table

p = Path(__file__).parent / "datasets" / "dataset_30217_6.txt"

spectrum = list(map(int, p.read_text().split()))
table = load_integer_mass_table()

peptides = cyclopeptideSequencingMass(spectrum, table)

out = []
for peptide in peptides:
    masses = [str(m) for m in peptide]
    out.append("-".join(masses))

result = " ".join(out)

with open("output.txt", "w") as f:
    f.write(result)
