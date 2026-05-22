from collections import defaultdict, deque
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

def pairedDeburjinFromKDpairs(
    kdpairs
):
    g = defaultdict(list)
    for left, right in kdpairs:
        prefix = (left[:-1], right[:-1])
        suffix = (left[1:], right[1:])
        g[prefix].append(suffix)
    return g

def stringSpelledByGappedPatterns(path, k, d):
    first = [p[0] for p in path]
    second = [p[1] for p in path]

    prefix_string = genome_path(first)
    suffix_string = genome_path(second)

    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return None 

    return prefix_string + suffix_string[-(k + d):]

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

    adj = {node: deque(edges) for node, edges in graph.items()}

    def has_edges(node):
        neighbors = adj.get(node)
        return neighbors is not None and len(neighbors) > 0

    def next_edge(node):
        return adj[node].popleft()

    stack = []
    path = []
    curr = start

    while stack or has_edges(curr):
        if not has_edges(curr):
            path.append(curr)
            curr = stack.pop()
            continue

        stack.append(curr)
        curr = next_edge(curr)

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

def in_out_degrees(graph):
    indeg = defaultdict(int)
    outdeg = defaultdict(int)

    for v, adj in graph.items():
        outdeg[v] = len(adj)
        for w in adj:
            indeg[w] += 1

    for v in graph:
        indeg.setdefault(v, 0)
        outdeg.setdefault(v, 0)

    return indeg, outdeg

def maximalNonBranchingPaths(graph):
    paths = []
    ind, outd = in_out_degrees(graph)
    visited_nodes = set()

    def is_1_in_1_out(v):
        return ind[v] == 1 and outd[v] == 1

    for v in graph:
        if not is_1_in_1_out(v) and outd[v] > 0:
            for w in graph[v]:
                path = [v, w]
                curr = w
                
                while is_1_in_1_out(curr):
                    visited_nodes.add(curr)
                    u = graph[curr][0]
                    path.append(u)
                    curr = u
                paths.append(path)

    for v in graph:
        if is_1_in_1_out(v) and v not in visited_nodes:
            visited_nodes.add(v)
            cycle = [v]
            curr = graph[v][0]
            
            while curr != v:
                cycle.append(curr)
                visited_nodes.add(curr)
                curr = graph[curr][0]
                
            cycle.append(v)
            paths.append(cycle)

    return paths
def generateContigs(patterns):
    g = debrujinFromKmers(patterns)
    paths = maximalNonBranchingPaths(g)
    strings =[genome_path(path) for path in paths] 
    out = sorted(strings)
    return out
