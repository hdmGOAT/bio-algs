def AlignmentGraph(v, w):
    n = len(v)
    m = len(w)

    out = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack = [[""] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            up = out[i - 1][j]
            left = out[i][j - 1]

            if up > left:
                out[i][j] = up
                backtrack[i][j] = "u"
            else:
                out[i][j] = left
                backtrack[i][j] = "l"
            if v[i - 1] == w[j - 1]:
                diag = out[i - 1][j - 1] + 1

                if diag > out[i][j]:
                    out[i][j] = diag
                    backtrack[i][j] = "tl"

    return out, backtrack

def OutputLCS(backtrack, v, i, j):
    lcs = []

    while i > 0 and j > 0:
        if backtrack[i][j] == "tl":
            lcs.append(v[i - 1])
            i -= 1
            j -= 1
        elif backtrack[i][j] == "u":
            i -= 1
        else:
            j -= 1

    return "".join(reversed(lcs))

def LongestPathInDAG(start, end, adj):
    n = len(adj)

    dist = [-float("inf")] * n
    parent: list[int | None] = [None] * n
    dist[start] = 0

    for u in range(start, n):
        for v, w in adj[u]:
            if dist[u] + w > dist[v]:
                dist[v] = dist[u] + w
                parent[v] = u

    path = []
    curr = end

    while curr is not None:
        path.append(curr)
        curr = parent[curr]

    path.reverse()

    return dist[end], path

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
