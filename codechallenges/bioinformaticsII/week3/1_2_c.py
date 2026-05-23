from pathlib import Path

from bioalgs.translation import findSubstringsThatEncodeAA, load_codon_table

p = Path(__file__).parent / "datasets" / "Bacillus_brevis.txt"

lines = p.read_text().splitlines()

s = ""
for line in lines:
    s += line.strip()
peptide = "VKLFPWFNQY"
table = load_codon_table()
locs = findSubstringsThatEncodeAA(s, peptide ,table)
print(len(locs))
