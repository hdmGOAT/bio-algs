from bioalgs.genome_assembly import debrujinFromKmers, eulerianPath, genome_path, in_out_degrees, pairedDeburjinFromKDpairs, stringSpelledByGappedPatterns


mers = """AAAT
AATG
ACCC
ACGC
ATAC
ATCA
ATGC
CAAA
CACC
CATA
CATC
CCAG
CCCA
CGCT
CTCA
GCAT
GCTC
TACG
TCAC
TCAT
TGCA""".splitlines()

g1 = debrujinFromKmers(mers)
path = eulerianPath(g1)
print(genome_path(path))

g2 = {
    1: [2, 3, 5],
    2: [4],
    3: [2,5],
    4: [1,2,5],
    5: [3]
}

i, o = in_out_degrees(g2)
edges_needed = 0
for node in range(1, 6):
    if o[node] > i[node]:
        edges_needed += o[node] - i[node]

print(edges_needed)
tups = [("ACC","ATA"), ("ACT","ATT"), ("ATA","TGA"),("ATT","TGA"), ("CAC","GAT"), ("CCG","TAC"), ("CGA","ACT"),("CTG","AGC"), ("CTG","TTC"), ("GAA","CTT"), ("GAT","CTG"), ("GAT","CTG"), ("TAC","GAT"), ("TCT","AAG"), ("TGA","GCT"), ("TGA","TCT"), ("TTC","GAA")]

g3 = pairedDeburjinFromKDpairs(tups)

path = eulerianPath(g3)

print(stringSpelledByGappedPatterns(path, 3, 1))
