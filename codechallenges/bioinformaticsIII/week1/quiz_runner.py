'''
Question 1

There is a unique longest common subsequence of the strings CTCGAT and TACGTC.  What is it?
'''

from bioalgs.sequence_alignment import AlignmentGraph, OutputLCS


g, bt  = AlignmentGraph("CTCGAT", "TACGTC")
lcs = OutputLCS(bt, "CTCGAT", 5, 5)
print(lcs)
'''
a -> b: 5

a -> c: 6

a -> d: 5

b -> c: 2

b -> f: 9

c -> e: 4

c -> f: 3

c -> g: 7

d -> e: 4

d -> f: 5

e -> g: 2

f -> g: 1

What is the longest path in this graph?  Give your answer as a sequence of nodes separated by spaces.  (Note: a, b, c, d, e, f, g is a topological order for this graph.)
'''
from bioalgs import LongestPathInDAG

edges_str = """a -> b: 5
a -> c: 6
a -> d: 5
b -> c: 2
b -> f: 9
c -> e: 4
c -> f: 3
c -> g: 7
d -> e: 4
d -> f: 5
e -> g: 2
f -> g: 1"""

nodes = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
node_to_idx = {n: i for i, n in enumerate(nodes)}

adj = [[] for _ in range(len(nodes))]
for line in edges_str.splitlines():
    u, rest = line.split(' -> ')
    v, w = rest.split(': ')
    adj[node_to_idx[u]].append((node_to_idx[v], int(w)))

dist, path_idx = LongestPathInDAG(node_to_idx['a'], node_to_idx['g'], adj)
path = [nodes[i] for i in path_idx]
print(" ".join(path))
