from pathlib import Path
from bioalgs.genome_assembly import debrujinFromString
p = Path(__file__).parent / "datasets" / "dataset_30183_6.txt"

k, text = p.read_text().splitlines()
d = debrujinFromString(text, int(k))
with open("output.txt", "w") as f:
    for key, values in d.items():
        line = f"{key}: {' '.join(values)}\n"
        f.write(line)
