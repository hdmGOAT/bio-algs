from pathlib import Path

from bioalgs.sequence_alignment import load_scoring_matrix, localAllignment, OutputLA

p = Path(__file__).parent / "datasets" / "dataset_30199_10.txt"
a, b = p.read_text().splitlines()
pam = load_scoring_matrix()
score, bt, ms, i, j = localAllignment(a, b, pam)

ala, alb = OutputLA(bt, a, b, i, j)

with open("output.txt", "w") as f:
    f.write(f"{score[i][j]}\n")
    f.write(f"{ala}\n")
    f.write(f"{alb}\n")
