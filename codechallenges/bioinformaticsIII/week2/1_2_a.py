from bioalgs import gaGraph, OutputGA

from pathlib import Path

p = Path(__file__).parent / "datasets" / "dataset_30199_3.txt"

nums, a, b = p.read_text().splitlines()
match, mis, indel = map(int,nums.split())

score, bt = gaGraph(a, b, match, mis, indel)

final_score = score[len(b)][len(a)]
aligned_a, aligned_b = OutputGA(bt, a, b)

with open("output.txt", "w") as f:
    f.write(f"{final_score}\n")
    f.write(f"{aligned_a}\n")
    f.write(f"{aligned_b}\n")
