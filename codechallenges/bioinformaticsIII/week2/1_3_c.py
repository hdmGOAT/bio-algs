from pathlib import Path

from bioalgs.sequence_alignment import OutputOA, overlapAlignment

p = Path(__file__).parent / "datasets" / "dataset_30200_7.txt"
nums_str, a, b = p.read_text().splitlines()

match, mismatch, indel = map(int, nums_str.split())

score, bt, max_score, max_j = overlapAlignment(a, b, match, mismatch, indel)
aout, bout = OutputOA(bt, a, b, len(a), max_j)

with open("output.txt", "w") as f:
    f.write(f"{max_score}\n")
    f.write(f"{aout}\n")
    f.write(f"{bout}\n")
