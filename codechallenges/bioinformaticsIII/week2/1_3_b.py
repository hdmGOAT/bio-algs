from pathlib import Path

from bioalgs.sequence_alignment import OutputLA, load_scoring_matrix, localAllignment

# p = Path(__file__).parent / "datasets" / ""

# a, b = p.read_text().splitlines()

a = "DISCREPANTLY"
b = "PATENT"
blosum_p = Path(__file__).parent / "datasets" / "BLOSUM62.txt"
blosum = load_scoring_matrix(blosum_p)

out, bt, max_score, i, j =localAllignment(a, b, blosum, 1)

out_a, out_b = OutputLA(bt, a, b, i, j)
best = out[i][j]

with open("output.txt", "w") as f:
    f.write(f"{best}\n")
    f.write(f"{out_a}\n")
    f.write(f"{out_b}\n")

