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
