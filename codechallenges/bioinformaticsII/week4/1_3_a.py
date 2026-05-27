from pathlib import Path

from bioalgs.translation import convolution


p = Path(__file__).parent / "datasets" / "dataset_30246_4.txt"

spectrum = list(map(int,p.read_text().split()))

with open("output.txt", "w") as f:
    f.write(" ".join(map(str, convolution(spectrum))))
