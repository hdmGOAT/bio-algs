'''
sample input
4 4
1 0 2 4 3   
4 6 5 2 1
4 4 5 2 1
5 6 8 5 3
-
3 2 4 0
3 2 4 2
0 7 3 3
3 3 0 2
1 3 2 2
'''
from pathlib import Path
def manhattanTourist(n, m, down, right):
    out = [[0] * (m + 1) for _ in range(n + 1)]

    # initialize top row
    for i in range(1, len(out[0])): out[0][i] = out[0][i-1] + right[0][i-1]

    # initialize 1st col
    for j in range(1, len(out)): out[j][0] = out[j-1][0] + down[j-1][0]

    for i in range(1, len(out)):
        for j in range(1, len(out[0])):
            # down check
            dc = out[i-1][j] + down[i-1][j]

            # right check
            rc = out[i][j-1] + right[i][j-1]
            out[i][j] = max(dc,rc)

    return out[n][m]

p = Path(__file__).parent / "datasets" / "dataset_30205_10.txt"

lines = p.read_text().splitlines()

nm = lines[0]
matrices = lines[1:]

sep = matrices.index("-")

down = [list(map(int, row.split())) for row in matrices[:sep]]
right = [list(map(int, row.split())) for row in matrices[sep+1:]]

n, m = map(int, nm.split())

print(manhattanTourist(n,m,down,right))
