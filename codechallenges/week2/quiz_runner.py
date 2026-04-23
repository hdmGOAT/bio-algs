from bioalgs.genome_metrics import HammingDistance, MinimumSkew, SkewArray
from bioalgs.sequence_patterns import Neighbors, PatternCount


a = "CTACAGCAATACGATCATATGCGGATCCGCAGTGGCCGGTAGACACACGT"
b = "CTACCCCGCTGCTCAATGACCGGGACTAAAGAGGCGAAGATTATGGTGTG"

print("1. " + str(HammingDistance(a, b)))

c = "CATTCCAGTACTTCATGATGGCGTGAAGA" 
print("2. " + str(SkewArray(c)))


d = "CGTGACAGTGTATGGGCATCTTT"
print("3. " + str(PatternCount(d, "TGT")))

e = "TGCAT"
print("4. " + str(len(Neighbors(e, 2))))
