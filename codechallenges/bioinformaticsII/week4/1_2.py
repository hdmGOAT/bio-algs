from pathlib import Path
from bioalgs.translation import extendedLeaderboardCyclopeptideSequencing


p = Path(__file__).parent / "datasets" / "dataset_30245_2.txt"

n, spec_str = p.read_text().splitlines()

n = int(n)
spectrum = list(map(int, spec_str.split()))
masses = [n for n in range(57, 201)]
peptides = extendedLeaderboardCyclopeptideSequencing(spectrum, n, masses)
result = " ".join("-".join(str(m) for m in peptide) for peptide in peptides)

with open("output.txt", "w") as f:
    f.write(result)

