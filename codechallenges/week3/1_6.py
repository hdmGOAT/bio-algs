
from pathlib import Path

from bioalgs.motifs_stochastic import GreedyMotifSearchWithPseudocounts

p = Path(__file__).parent / "datasets" / "dataset_30306_9.txt"

nums, dna = p.read_text().splitlines()
k, t = nums.split()

dna_sp = dna.split()

with open("output.txt", "w") as f:
    f.write(' '.join(GreedyMotifSearchWithPseudocounts(dna_sp, int(k), int(t))))
