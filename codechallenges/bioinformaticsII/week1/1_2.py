from pathlib import Path

from bioalgs.genome_assembly import composition

path = Path(__file__).parent / "datasets" / "dataset_30153_3.txt"

k, text = path.read_text().splitlines()
with open("output.txt", "w") as f:
    f.write(" ".join(composition(text, int(k))))

