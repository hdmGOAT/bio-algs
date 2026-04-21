from pathlib import Path

from bioalgs.sequence_patterns import Neighbors
p = Path(__file__).parent / "datasets" / "dataset_30282_4.txt"

pattern, d = p.read_text().splitlines()

with open("output.txt", "w") as f:
    f.write(" ".join(Neighbors(pattern, int(d))))


