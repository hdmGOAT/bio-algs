from pathlib import Path

from bioalgs.translation import BFCountPeptides, load_integer_mass_table


p = Path(__file__).parent / "datasets" / "dataset_30216_2.txt"

table = load_integer_mass_table()
n = int(p.read_text())
print(BFCountPeptides(n, table))
