from pathlib import Path
from bioalgs.sequence_patterns import PatternMatching

path = Path(__file__).parent / "datasets" / "dataset_30273_5.txt"
pattern, genome = path.read_text('UTF-8').splitlines()

positions = PatternMatching(pattern, genome)

with open("output.txt", "w") as f:
    f.write(" ".join(map(str, positions)))
