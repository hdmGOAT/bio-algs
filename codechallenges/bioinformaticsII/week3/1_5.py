from pathlib import Path

from bioalgs.translation import BFCountPeptides, load_integer_mass_table


p = Path(__file__).parent / "datasets" / "dataset_30216_2.txt"

table = load_integer_mass_table()
n = int(p.read_text())
nCount = BFCountPeptides(n, table)
print(nCount)


# estimate C
Cs = []
for m in range(1200, 1300):
    Cs.append(BFCountPeptides(m+1, table) / BFCountPeptides(m, table))

C = sum(Cs) / len(Cs)

print(C)
