from pathlib import Path

from bioalgs.sequence_patterns import FrequentWordsWithMismatches


p = Path(__file__).parent / "datasets" /  "dataset_30278_9.txt"

text, nums = p.read_text().splitlines()
k, d = nums.split();

with open("output.txt", "w") as f:
    f.write(" ".join(FrequentWordsWithMismatches(text, int(k), int(d))))
