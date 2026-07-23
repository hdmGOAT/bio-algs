from io import UnsupportedOperation


def alignmentWithAffineGapsGraph(a, b, match, mis, go, ge):
    upper = [[0] * len(a) for _ in range(len(b))]
    middle = [[0] * len(a) for _ in range(len(b))]
    lower = [[0] * len(a) for _ in range(len(b))]
