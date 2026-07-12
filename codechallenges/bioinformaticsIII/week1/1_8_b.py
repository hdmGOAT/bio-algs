from pathlib import Path

from bioalgs import LongestPathInDAG

p = Path(__file__).parent / "datasets" / "dataset_30197_7.txt"

lines = p.read_text().splitlines()
start, end = map(int, lines[0].split())
mat = lines[1:]

edges = [list(map(int, row.split())) for row in mat]

num_nodes = max(
    max(u, v)
    for u, v, _ in edges
) + 1

adj = [[] for _ in range(num_nodes)]

for u, v, w in edges:
    adj[u].append((v, w))

with open("output.txt", "w") as f:
    dist, path = LongestPathInDAG(start, end, adj)

    f.write(f"{dist}\n")
    f.write(" ".join(map(str, path)))
