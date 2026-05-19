from pathlib import Path

from bioalgs.genome_assembly import stringReconstruction

p = Path(__file__).parent / "datasets" / "dataset_30187_7.txt"

k, pat = p.read_text().splitlines()
patterns = pat.strip().split()

with open("output.txt", "w") as f:
    f.write(stringReconstruction(patterns))
