from pathlib import Path

from bioalgs.sequence_patterns import FrequentWordsWithMismatchesAndRC


p = Path(__file__).parent / "datasets" / "dataset_30278_10.txt"


text, nums = p.read_text().splitlines()
d, k = nums.split()

with open("output.txt", "w") as f:
    f.write(" ".join(FrequentWordsWithMismatchesAndRC(text, int(d), int(k))))
