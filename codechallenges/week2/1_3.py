from bioalgs.genome_metrics import MinimumSkew
from pathlib import Path
path = Path(__file__).parent / "datasets" / "dataset_30277_10.txt"

inp = path.read_text('UTF-8').strip()

items = MinimumSkew (inp)

with open("output.txt", "w") as f:
    f.write(" ".join(map(str, items)))
