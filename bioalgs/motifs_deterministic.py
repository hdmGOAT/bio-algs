from bioalgs.sequence_patterns import DistanceBetweenPatternAndStrings, FrequentWordsWithMismatches, Neighbors

import math

def MotifsEntropy(motifs):
    profile = Profile(motifs)
    k = len(motifs[0])
    entropy = 0.0

    for j in range(k):
        for base in "ACGT":
            p = profile[base][j]
            if p > 0:
                entropy -= p * math.log2(p)

    return entropy

def MedianString(dna, k):
    distance = float('inf')
    patterns = Neighbors('A' * k, k)
    for pattern in patterns:
        d = DistanceBetweenPatternAndStrings(pattern, dna)
        if distance > d:
            distance = d
            median = pattern

    return median

def MedianStringList(dna, k):
    best_distance = float('inf')
    best_patterns = []
    patterns = Neighbors('A' * k, k)

    for pattern in patterns:
        d = DistanceBetweenPatternAndStrings(pattern, dna)

        if d < best_distance:
            best_distance = d
            best_patterns = [pattern]
        elif d == best_distance:
            best_patterns.append(pattern)

    return best_patterns

def MotifEnumeration(Dna, k, d):
    dna_list = Dna.split()
    neighbor_sets = []

    for dna in dna_list:
        motifs = set()
        for i in range(len(dna) - k + 1):
            kmer = dna[i:i+k]
            motifs |= Neighbors(kmer, d)
        neighbor_sets.append(motifs)

    return set.intersection(*neighbor_sets)

def Count(Motifs):
    count = {base: [0] * len(Motifs[0]) for base in "ACGT"}
    for motif in Motifs:
        for j, char in enumerate(motif):
            count[char][j] += 1
    return count


def Profile(Motifs):
    t = len(Motifs)
    profile = Count(Motifs)
    for base in profile:
        for i in range(len(profile[base])):
            profile[base][i] = profile[base][i] / t
    return profile


def Consensus(Motifs):
    count = Count(Motifs)
    n = len(Motifs[0])
    consensus = ""
    bases = list(count.keys())
    for i in range(n):
        best = bases[0]
        for base in bases:
            if count[base][i] > count[best][i]:
                best = base
        consensus += best
    return consensus


def Score(Motifs):
    consensus = Consensus(Motifs)
    score = 0
    for j, base in enumerate(consensus):
        for motif in Motifs:
            if motif[j] != base:
                score += 1
    return score


def Pr(Text, Profile):
    probability = 1.0
    for i, char in enumerate(Text):
        probability *= Profile[char][i]
    return probability


def ProfileMostProbableKmer(text, k, profile):
    max_probability = -1.0
    most_probable = ""

    for i in range(len(text) - k + 1):
        kmer = text[i : i + k]
        probability = Pr(kmer, profile)
        if probability > max_probability:
            max_probability = probability
            most_probable = kmer

    return most_probable


def GreedyMotifSearch(Dna, k, t):
    n = len(Dna[0])
    best_motifs = [Dna[i][0:k] for i in range(t)]
    for i in range(n - k + 1):
        motifs = [Dna[0][i : i + k]]
        for j in range(1, t):
            profile = Profile(motifs[0:j])
            motifs.append(ProfileMostProbableKmer(Dna[j], k, profile))
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
    return best_motifs
