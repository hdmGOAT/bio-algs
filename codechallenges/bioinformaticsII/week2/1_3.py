from pathlib import Path
from bioalgs.genome_assembly import (
    eulerianPath,
    pairedDeburjinFromKDpairs,
    stringSpelledByGappedPatterns,
)

p = Path(__file__).parent / "datasets" / "dataset_30188_16.txt"

nums, pairs_text = p.read_text().splitlines()
k, d = map(int, nums.split())

pairs = [tuple(p.split('|')) for p in pairs_text.split()]
graph = pairedDeburjinFromKDpairs(pairs)
path = eulerianPath(graph)

result = stringSpelledByGappedPatterns(path, k, d)

with open("output.txt", "w") as f:
    if result is not None:
        f.write(result)
