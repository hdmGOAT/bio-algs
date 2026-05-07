from pathlib import Path

from bioalgs.motifs_deterministic import MedianString

p = Path(__file__).parent / "datasets" / "dataset_30304_9.txt"

k, dna = p.read_text().splitlines()

print(MedianString(dna, int(k)))
