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

def Neighbors(Pattern, d):
    if d == 0:
        return {Pattern}
    a = 'ACGT'
    neighborhood = {Pattern}
    frontier = {Pattern}
    for _ in range(d):
        new_frontier = set()
        for t in frontier:
            for i, ch in enumerate(t):
                for c in a:
                    if c != ch:
                        u = t[:i] + c + t[i+1:]
                        if u not in neighborhood:
                            neighborhood.add(u)
                            new_frontier.add(u)
        frontier = new_frontier
    return neighborhood
        

def FrequentWordsWithMismatches(text, k, d):
    freqMap = {}
    m = 0
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i:i+k]
        neighborhood = Neighbors(pattern, d)
        for neighbor in neighborhood:
            c =  freqMap.get(neighbor, 0) + 1
            freqMap[neighbor] = c
            m = max(m, c)
    return [key for key, val in freqMap.items() if val == m] 

    
def FrequentWordsWithMismatchesAndRC(text, k, d):
    freqMap = {}
    m = 0
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i:i+k]
        neighborhood = Neighbors(pattern, d)
        for neighbor in neighborhood:
            rc = ReverseComplement(neighbor)
            freqMap[neighbor] = freqMap.get(neighbor, 0) + 1
            freqMap[rc] = freqMap.get(rc, 0) + 1
            m = max(m, freqMap[neighbor], freqMap[rc])
    return [key for key, val in freqMap.items() if val == m] 


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
