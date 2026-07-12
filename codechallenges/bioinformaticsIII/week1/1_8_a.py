from bioalgs import AlignmentGraph, OutputLCS

from pathlib import Path

p = Path(__file__).parent / "datasets" / "dataset_30197_5.txt"

v, w = p.read_text().splitlines()

_, bt =  AlignmentGraph(v, w)

with open("output.txt", "w") as f:
    f.write(OutputLCS(bt, v, len(v), len(w)))
