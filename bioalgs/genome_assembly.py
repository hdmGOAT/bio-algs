from collections import defaultdict
from itertools import product
def composition(text, k):
    c = []
    seen = set()

    for i in range(len(text) - k + 1):
        word = text[i:i+k]

        if word not in seen:
            c.append(word)
            seen.add(word)

    return c

def genome_path(kmers):
    return kmers[0] + "".join(kmer[-1] for kmer in kmers[1:])

def overlap(patterns):
    overlap = defaultdict(list)
    for a in patterns:
        suffix = a[1:]
        for b in patterns:
            if b == a:
                continue
            if b.startswith(suffix):
                overlap[a].append(b)
    return overlap

def debrujinFromString(text, k):
    g = defaultdict(list)
    for s in range(1, len(text)-k+2):
        e = s+k-1
        last = text[s-1:e-1]
        now = text[s:e]
        g[last].append(now)

    return g

def debrujinFromKmers(kmers):
    g = defaultdict(list)
    for kmer in kmers:
        g[kmer[:-1]].append(kmer[1:])
    return g

def eulerianCycle(graph):
    curr = None
    for v, adj in graph.items():
        if adj:
            curr = v
            break

    if curr is None:
        return []

    stack = []
    path = []
    while graph[curr] or stack:
        if not graph[curr]:
            path.append(curr)
            curr = stack.pop()
            continue
        next = graph[curr].pop()
        stack.append(curr)
        curr=next
    path.append(curr)

    return list(reversed(path))

def find_start(graph):
    indeg = defaultdict(int)
    outdeg = defaultdict(int)

    for u in graph:
        for v in graph[u]:
            outdeg[u] += 1
            indeg[v] += 1

    start = None

    for node in set(list(indeg.keys()) + list(outdeg.keys())):
        out = outdeg[node]
        inn = indeg[node]

        if out - inn == 1:
            return node

        if out > 0:
            start = node

    return start

def eulerianPath(graph):
    start = find_start(graph)
    if start is None:
        return []

    stack = []
    path = []
    curr = start

    while stack or graph.get(curr, []):
        if not graph.get(curr, []):
            path.append(curr)
            curr = stack.pop()
            continue

        stack.append(curr)
        curr = graph[curr].pop()

    path.append(curr)
    path.reverse()
    return path

def stringReconstruction(patterns):
    db = debrujinFromKmers(patterns)
    path = eulerianPath(db)
    text  = genome_path(path)

    return text

def kUniversalString(k):
    pieces = [''.join(p) for p in product('01', repeat=k)]     

    graph = debrujinFromKmers(pieces)
    cycle = eulerianCycle(graph)

    text = genome_path(cycle)

    return text[:-k+1]
