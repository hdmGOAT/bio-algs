from pathlib import Path

from bioalgs.motifs_deterministic import Score
from bioalgs.motifs_stochastic import GibbsSampler
p = Path(__file__).parent / "datasets" / "dataset_30309_11.txt"
nums, text = p.read_text().splitlines()

k, t, N = nums.split()
dna = text.split()
best_score = float('inf')
for i in range(20):
    m = GibbsSampler(dna, int(k), int(t), int(N))
    score = Score(m)
    if score < best_score:
        best_score = score
        motifs = m

with open("output.txt", "w") as f:
    f.write(" ".join(motifs))
