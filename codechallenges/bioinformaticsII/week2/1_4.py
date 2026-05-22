from pathlib import Path
from bioalgs.genome_assembly import generateContigs

p = Path(__file__).parent / "datasets" / "dataset_30189_5.txt"
s = p.read_text().split()

with open("output.txt", "w") as f:
    f.write(" ".join(generateContigs(s)))

