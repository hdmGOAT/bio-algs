from pathlib import Path

from bioalgs.genome_assembly import kUniversalString

p = Path(__file__).parent / "datasets" / "dataset_30187_11.txt"

k = p.read_text().strip()

with open("output.txt", "w") as f:
    f.write(kUniversalString(int(k)))
