def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i : i + len(Pattern)] == Pattern:
            count += 1
    return count


def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    extended = Genome + Genome[0 : n // 2]
    for i in range(n):
        array[i] = PatternCount(extended[i : i + (n // 2)], symbol)
    return array


def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    extended = Genome + Genome[0 : n // 2]

    array[0] = PatternCount(Genome[0 : n // 2], symbol)
    for i in range(1, n):
        array[i] = array[i - 1]
        if extended[i - 1] == symbol:
            array[i] -= 1
        if extended[i + (n // 2) - 1] == symbol:
            array[i] += 1
    return array


def SkewArray(Genome):
    out = [0]
    for char in Genome:
        if char == "C":
            out.append(out[-1] - 1)
        elif char == "G":
            out.append(out[-1] + 1)
        else:
            out.append(out[-1])
    return out


def MinimumSkew(Genome):
    skew = SkewArray(Genome)
    minimum = min(skew)
    return [i for i, value in enumerate(skew) if value == minimum]


def HammingDistance(p, q):
    count = 0
    for i, char in enumerate(p):
        if char != q[i]:
            count += 1
    return count


def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i : i + len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions


def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Text[i : i + len(Pattern)], Pattern) <= d:
            count += 1
    return count
