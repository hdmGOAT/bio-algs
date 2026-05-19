
from pathlib import Path

from bioalgs.genome_assembly import eulerianPath

p = Path(__file__).parent / "datasets" / "dataset_30187_6.txt"

lines = p.read_text().splitlines()
graph = {}
for line in lines:
    key, vals = line.split(":")
    key = int(key.strip())
    if vals.strip():
        graph[key] = list(map(int, vals.split()))
    else:
        graph[key] = []

cycle = eulerianPath(graph)
with open("output.txt", "w") as f:
    f.write(" ".join(map(str, cycle)))
