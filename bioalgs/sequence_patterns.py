def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i : i + len(Pattern)] == Pattern:
            count += 1
    return count

def FindClumps(genome, k, L, t):
    patterns = set()
    n = len(genome)

    window_map = FrequencyMap(genome[0:L], k)
    patterns.update(key for key, val in window_map.items() if val >= t)

    for i in range(1, n - L + 1):
        old = genome[i-1:i-1+k]
        window_map[old] -= 1
        if window_map[old] == 0:
            del window_map[old]

        new = genome[i+L-k:i+L]
        window_map[new] = window_map.get(new, 0) + 1

        if window_map[new] >= t:
            patterns.add(new)

    return list(patterns)

def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        pattern = Text[i : i + k]
        freq[pattern] = freq.get(pattern, 0) + 1
    return freq


def FrequentWords(Text, k):
    freq = FrequencyMap(Text, k)
    max_count = max(freq.values())
    return [pattern for pattern, count in freq.items() if count == max_count]


def Complement(Pattern):
    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}
    out = ""
    for char in Pattern:
        out += complements[char]
    return out


def ReverseComplement(Pattern):
    return Complement("".join(reversed(Pattern)))


def PatternMatching(Pattern, Genome):
    positions = []
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i : i + len(Pattern)] == Pattern:
            positions.append(i)
    return positions
