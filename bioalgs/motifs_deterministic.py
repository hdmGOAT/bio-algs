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
