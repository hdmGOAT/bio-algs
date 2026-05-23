from pathlib import Path

from bioalgs.translation import linearSpectrum, load_integer_mass_table

p = Path(__file__).parent / "datasets" / "dataset_30248_2.txt"

pep = p.read_text().strip()
table = load_integer_mass_table()
with open("output.txt", "w") as f:
    f.write(" ".join(map(str, linearSpectrum(pep, table))))
