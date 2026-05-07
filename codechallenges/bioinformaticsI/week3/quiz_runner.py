from bioalgs.motifs_deterministic import MedianStringList, MedianStringList, Pr


m = "CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"

print("1. " + str(MedianStringList(m, 7)))

prof = {
    "A": [0.4, 0.3, 0.0, 0.1, 0.0, 0.9],
    "C": [0.2, 0.3, 0.0, 0.4, 0.0, 0.1],
    "G": [0.1, 0.3, 1.0, 0.1, 0.5, 0.0],
    "T": [0.3, 0.1, 0.0, 0.4, 0.5, 0.0],
}

print("2. " + str(Pr("AAGTTC", prof)))
