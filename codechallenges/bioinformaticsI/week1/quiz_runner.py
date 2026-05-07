from bioalgs.sequence_patterns import FrequencyMap, PatternCount, PatternMatching, ReverseComplement


print(PatternCount("CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC", "CGCG"))

m = FrequencyMap("TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT", 3)
ma = 0
l = None
for key, val in m.items():
    if val > ma:
        ma = val
        l = key
    elif val == ma:
        print("tie")
        print(key)
print(l)


print(ReverseComplement("TTGTGTC"))

print(PatternMatching(Genome="GACGATATACGACGATA", Pattern="ATA"))
