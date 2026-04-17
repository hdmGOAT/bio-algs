import random

from .motifs_deterministic import Pr, Score


def CountWithPseudocounts(Motifs):
    k = len(Motifs[0])
    counts = {base: [1] * k for base in "ACGT"}
    for i in range(k):
        for motif in Motifs:
            counts[motif[i]][i] += 1
    return counts


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    profile = CountWithPseudocounts(Motifs)
    for base in profile:
        for i in range(len(profile[base])):
            profile[base][i] = profile[base][i] / (t + 4)
    return profile


def ProfileMostProbablePattern(Text, k, profile):
    max_prob = -1.0
    most_probable = ""
    for i in range(len(Text) - k + 1):
        pattern = Text[i : i + k]
        prob = Pr(pattern, profile)
        if prob > max_prob:
            max_prob = prob
            most_probable = pattern
    return most_probable


def Motifs(Profile, Dna):
    k = len(Profile["A"])
    return [ProfileMostProbablePattern(text, k, Profile) for text in Dna]


def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    n = len(Dna[0])
    best = [Dna[i][0:k] for i in range(t)]
    for i in range(n - k + 1):
        motifs = [Dna[0][i : i + k]]
        for j in range(1, t):
            profile = ProfileWithPseudocounts(motifs)
            motifs.append(ProfileMostProbablePattern(Dna[j], k, profile))
        if Score(motifs) < Score(best):
            best = motifs
    return best


def RandomMotifs(Dna, k, t):
    motifs = []
    for strand in Dna:
        start = random.randint(0, len(strand) - k)
        motifs.append(strand[start : start + k])
    return motifs


def RandomizedMotifSearch(Dna, k, t):
    motifs = RandomMotifs(Dna, k, t)
    best_motifs = motifs[:]

    while True:
        profile = ProfileWithPseudocounts(motifs)
        motifs = Motifs(profile, Dna)
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


def Normalize(Probabilities):
    total = sum(Probabilities.values())
    return {kmer: prob / total for kmer, prob in Probabilities.items()}


def WeightedDie(Probabilities):
    r = random.uniform(0, 1)
    for key, value in Probabilities.items():
        r -= value
        if r <= 0:
            return key
    return next(iter(Probabilities))


def ProfileGeneratedString(Text, profile, k):
    probabilities = {}
    for i in range(len(Text) - k + 1):
        kmer = Text[i : i + k]
        probabilities[kmer] = Pr(kmer, profile)
    return WeightedDie(Normalize(probabilities))


def GibbsSampler(Dna, k, t, N):
    motifs = RandomMotifs(Dna, k, t)
    best_motifs = motifs[:]

    for _ in range(N):
        i = random.randint(0, t - 1)
        excluded = motifs[:i] + motifs[i + 1 :]
        profile = ProfileWithPseudocounts(excluded)
        motifs[i] = ProfileGeneratedString(Dna[i], profile, k)

        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs[:]

    return best_motifs
