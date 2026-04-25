from pathlib import Path

from bioalgs.motifs_stochastic import RandomizedMotifSearchNTimes


p = Path(__file__).parent / "datasets" / "dataset_30307_5.txt"

nums, text = p.read_text().splitlines()
k, t = nums.split()

dna = text.split()

t = " ".join(RandomizedMotifSearchNTimes(dna, int(k), int(t), 1000))


with open("output.txt", "w") as f:
    f.write(t)
