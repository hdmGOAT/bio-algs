from pathlib import Path

from bioalgs.translation import convolutionCyclopeptideSequencing


dataset_path = Path(__file__).parent / "datasets" / "dataset_30246_8.txt"

raw = dataset_path.read_text().strip()

lines = [line.strip() for line in raw.splitlines() if line.strip()]
m = int(lines[0])
n = int(lines[1])
spectrum = list(map(int, " ".join(lines[2:]).split()))

peptides = convolutionCyclopeptideSequencing(spectrum, m, n)
result = " ".join("-".join(str(mass) for mass in peptide) for peptide in sorted(peptides))
print(len(peptides))
with open("output.txt", "w") as f:
	f.write(result)
