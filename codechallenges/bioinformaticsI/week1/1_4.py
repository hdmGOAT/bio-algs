from pathlib import Path
from bioalgs.sequence_patterns import FindClumps

p = Path(__file__).parent / "datasets" / "dataset_30274_5.txt"

text, nums = p.read_text('UTF-8').splitlines()
k, l, t = nums.split()

print(FindClumps(text, int(k), int(l), int(t)))
