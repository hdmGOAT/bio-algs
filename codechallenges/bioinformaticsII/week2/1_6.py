
from pathlib import Path

from bioalgs.genome_assembly import maximalNonBranchingPaths

p = Path(__file__).parent / "datasets" / ""

#lines = p.read_text().splitlines()
lines = """1: 2
2: 3
3: 4 5
6: 7
7: 6""".splitlines()

graph = {}
for line in lines:
    key, vals = line.split(":")
    key = int(key.strip())
    if vals.strip():
        graph[key] = list(map(int, vals.split()))
    else:
        graph[key] = []



print(maximalNonBranchingPaths(graph))
