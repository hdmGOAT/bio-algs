from pathlib import Path

from bioalgs.sequence_alignment import gaGraph
p = Path(__file__).parent / "datasets" / "dataset_30200_3.txt"

a,b = p.read_text().splitlines()
sc, bt = gaGraph(a, b, 0, 1, 1)

final_score = sc[len(b)][len(a)]

print(-final_score)
