
from pathlib import Path

from bioalgs.motifs_deterministic import GreedyMotifSearch

p = Path(__file__).parent / "datasets" / "dataset_30305_5.txt"

nums, dna = p.read_text().splitlines()
k, t = nums.split()

dna_sp = dna.split()

with open("output.txt", "w") as f:
    f.write(' '.join(GreedyMotifSearch(dna_sp, int(k), int(t))))
