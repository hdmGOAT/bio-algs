from pathlib import Path

from bioalgs.translation import linearScoreCyclopeptide, load_integer_mass_table

p = Path(__file__).parent / "datasets" / "dataset_30249_1.txt"

pep, spec = p.read_text().splitlines()
table = load_integer_mass_table()
spec_split = list(map(int, spec.split()))
print(linearScoreCyclopeptide(pep, spec_split, table))
