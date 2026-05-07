from pathlib import Path

from bioalgs.motifs_deterministic import MotifEnumeration
p = Path(__file__).parent / "datasets" / "dataset_30302_8.txt"

nums, dna = p.read_text().splitlines()
k, d = nums.split()

with open("output.txt", "w") as f:
    f.write(" ".join(MotifEnumeration(dna, int(k), int(d))))
