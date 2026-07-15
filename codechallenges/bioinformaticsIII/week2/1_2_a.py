def gaGraph(a, b, match, mismatch, indel):
    a_len = len(a)
    b_len = len(b)

    # rows = b, cols = a
    score = [[0] * (a_len + 1) for _ in range(b_len + 1)]
    bt = [[""] * (a_len + 1) for _ in range(b_len + 1)]

    # initialize first row
    for j in range(1, a_len + 1):
        score[0][j] = score[0][j - 1] - indel
        bt[0][j] = "l"

    # initialize first column
    for i in range(1, b_len + 1):
        score[i][0] = score[i - 1][0] - indel
        bt[i][0] = "u"

    # fill DP table
    for i in range(1, b_len + 1):
        for j in range(1, a_len + 1):
            b_let = b[i - 1]
            a_let = a[j - 1]

            # diagonal (match/mismatch)
            if a_let == b_let:
                across = score[i - 1][j - 1] + match
            else:
                across = score[i - 1][j - 1] - mismatch

            # gap penalties
            left = score[i][j - 1] - indel
            up = score[i - 1][j] - indel

            best = max(across, left, up)
            score[i][j] = best

            if best == across:
                bt[i][j] = "d"
            elif best == left:
                bt[i][j] = "l"
            else:
                bt[i][j] = "u"

    return score, bt

def OutputGA(bt, a, b):
    i = len(b)
    j = len(a)
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

    return "".join(reversed(aligned_a)), "".join(reversed(aligned_b))

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
