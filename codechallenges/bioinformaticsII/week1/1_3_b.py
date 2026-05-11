from pathlib import Path

from bioalgs.genome_assembly import overlap

p = Path(__file__).parent / "datasets" / "dataset_30182_10.txt"

patterns = p.read_text().strip().split()

d = overlap(patterns)

with open("output.txt", "w") as f:
    for key, values in d.items():
        line = f"{key}: {' '.join(values)}\n"
        f.write(line)
