def localAllignment(a, b, pam, indel=5):
    a_len = len(a)
    b_len = len(b)
    out = [[0] * (a_len + 1) for _ in range(b_len + 1)]
    bt = [[""] * (a_len + 1) for _ in range(b_len + 1)]
    
    max_score = 0
    max_i, max_j = 0, 0

    for i in range(1, b_len + 1):
        for j in range(1, a_len + 1):
            a_let = a[j-1]
            b_let = b[i-1]
            # diag
            diag = out[i-1][j-1] + pam[(a_let, b_let)]
            # left
            left = out[i][j-1] - indel
            # up
            up = out[i-1][j] - indel

            best = max(diag, left, up, 0)
            out[i][j] = best
            
            if best > max_score:
                max_score = best
                max_i = i
                max_j = j

            if best == 0:
                bt[i][j] = ""
            elif best == diag:
                bt[i][j] = "d"
            elif best == left:
                bt[i][j] = "l"
            elif best == up:
                bt[i][j] = "u"

    return out, bt, max_score, max_i, max_j

def OutputLA(bt, a, b, i, j):
    aligned_a = []
    aligned_b = []

    while i > 0 or j > 0:
        if i > 0 and j > 0 and bt[i][j] == "d":
            aligned_a.append(a[j-1])
            aligned_b.append(b[i-1])
            i -= 1
            j -= 1
        elif j > 0 and bt[i][j] == "l":
            aligned_a.append(a[j-1])
            aligned_b.append("-")
            j -= 1
        elif i > 0 and bt[i][j] == "u":
            aligned_a.append("-")
            aligned_b.append(b[i-1])
            i -= 1
        elif bt[i][j] == "":
            break

    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b))

from pathlib import Path

from bioalgs.sequence_alignment import load_scoring_matrix

p = Path(__file__).parent / "datasets" / "dataset_30199_10.txt"
a, b = p.read_text().splitlines()
pam = load_scoring_matrix()
score, bt, ms, i, j = localAllignment(a, b, pam)

ala, alb = OutputLA(bt, a, b, i, j)

with open("output.txt", "w") as f:
    f.write(f"{score[i][j]}\n")
    f.write(f"{ala}\n")
    f.write(f"{alb}\n")
