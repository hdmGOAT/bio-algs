'''
Consider the following alignment:

TCGAC--ATT

CC---GAA-T

What is the score of this alignment if the match score is 1, the mismatch penalty is 1, and the indel penalty is 2?
'''

from bioalgs.sequence_alignment import AlignmentGraph, overlapAlignment


def scoreGetter(a, b, match, mis, ind):
    score = 0
    for i in range(len(a)):
        if a[i] == b[i]:
            score += match
        elif a[i] == "-" or b[i] == "-":
            score -= ind
        else:
            score -= mis
    return score

# 1
print(scoreGetter("TCGAC--ATT", "CC---GAA-T", 1, 1, 2))

# 2
print(scoreGetter("ATAGCGACGCCT", "ATA-CGATA-CA", 1, 3, 1))

# 3
print(scoreGetter("GATA-CACT", "GATACCGCT", 1, 1, 1))

# 4
a = "AGTACATCAGAGGAGTT-ACATACTAACG"
b= "AGTTCACAGGCTA-CGTACAGATATTACGACAGGCAGA"

_, _, score, _ = overlapAlignment(a,b,1,0,2)
print(score)
