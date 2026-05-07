
from pathlib import Path
from bioalgs.sequence_patterns import FindClumps

p = Path(__file__).parent / "datasets" / "E_coli.txt"

text = p.read_text('UTF-8')

c = FindClumps(text, 9, 500, 3)
print(len(c))
