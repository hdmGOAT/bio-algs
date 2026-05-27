from pathlib import Path

from bioalgs.translation import load_integer_mass_table, trim


p = Path(__file__).parent / "datasets" / "dataset_30249_3.txt"

peptides_str, spec_str, n_str = p.read_text().splitlines()

peptides = peptides_str.split()
spectrum = list(map(int, spec_str.split()))
n = int(n_str)
table = load_integer_mass_table()
with open("output.txt", "w") as f:
    f.write(" ".join(trim(peptides, spectrum, n, table)))
