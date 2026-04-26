from bioalgs.motifs_deterministic import Score
from bioalgs.motifs_stochastic import Motifs, ProfileWithPseudocounts, RandomMotifs


dna = ["AAGCCAAA", "AATCCTGG", "GCTACTTG", "ATGTTTTG"]


motifs = ["CCA", "CCT", "CTT", "TTG"]
best_motifs = motifs[:]

while True:
    profile = ProfileWithPseudocounts(motifs)
    motifs = Motifs(profile, dna)
    if Score(motifs) < Score(best_motifs):
        best_motifs = motifs
    else:
        print(" ".join(best_motifs))
        break
