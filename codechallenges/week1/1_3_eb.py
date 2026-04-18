from pathlib import Path
from bioalgs.sequence_patterns import PatternMatching

path = Path(__file__).parent / "datasets" / """Vibrio_cholerae.txt"""

genome = path.read_text('UTF-8')
pattern = "CTTGATCAT"
positions = PatternMatching(pattern, genome)

with open("output.txt", "w") as f:
    f.write(" ".join(map(str, positions)))
