from pathlib import Path

from bioalgs.motifs_stochastic import RandomizedMotifSearch


p = Path(__file__).parent / "datasets" / "dataset_30307_5.txt"

nums, text = p.read_text().splitlines()
k, t = nums.split()

dna = text.split()
with open("output.txt", "w") as f:
    f.write(" ".join(RandomizedMotifSearch(dna, int(k), int(t))))
