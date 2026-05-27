from pathlib import Path

from bioalgs.translation import load_integer_mass_table, scoreCyclopeptide


p = Path(__file__).parent / "datasets" / "dataset_30244_3.txt"

pep, spec_str = p.read_text().splitlines()

spectrum = list(map(int, spec_str.split()))
table = load_integer_mass_table()
print(scoreCyclopeptide(pep, spectrum, table))
