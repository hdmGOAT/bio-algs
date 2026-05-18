from pathlib import Path

from bioalgs.genome_assembly import debrujinFromKmers

p = Path(__file__).parent / "datasets" / "dataset_30184_8.txt"
patterns = p.read_text().strip().split()
d = debrujinFromKmers(patterns)
with open("output.txt", "w") as f:
    for key, values in d.items():
        line = f"{key}: {' '.join(values)}\n"
        f.write(line)
