from pathlib import Path
from collections import Counter

from bioalgs.sequence_patterns import ReverseComplement
def load_codon_table(
    path: Path = Path(__file__).parent / "constants" / "RNA_codon_table.txt"
):
    return {
        parts[0]: parts[1] if len(parts) == 2 else "*"
        for parts in (
            line.split()
            for line in path.read_text().splitlines()
            if line.strip()
        )
    }

def load_integer_mass_table(path: Path = Path(__file__).parent / "constants" / "integer_mass_table.txt"):
    return {
        parts[0]: int(parts[1])
        for parts in (
            line.split()
            for line in path.read_text().splitlines()
            if line.strip()
        )
    }

def translateRNAtoAminoAcidstr(rna, table) -> str:
    protein = []

    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        amino = table.get(codon)

        if amino == "*" or amino is None:
            break

        protein.append(amino)

    return "".join(protein)

def findSubstringsThatEncodeAA(dna: str, peptide: str, table):
    result = []
    k = len(peptide)
    window_len = 3 * k

    for i in range(len(dna) - window_len + 1):
        window = dna[i:i + window_len]

        rna_fwd = window.replace("T", "U")
        if translateRNAtoAminoAcidstr(rna_fwd, table) == peptide:
            result.append(window)
            continue

        rev_window = ReverseComplement(window)
        rna_rev = rev_window.replace("T", "U")
        if translateRNAtoAminoAcidstr(rna_rev, table) == peptide:
            result.append(window)

    return result

def linearSpectrum(peptide, aaMass):
    prefixMass = [0]
    for i in range(len(peptide)):
        prefixMass.append(prefixMass[i] + aaMass[peptide[i]])

    lS = [0]
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            lS.append(prefixMass[j] - prefixMass[i])

    return sorted(lS)

def linearSpectrumMass(peptide_masses):
    prefixMass = [0]
    for i in range(len(peptide_masses)):
        prefixMass.append(prefixMass[i] + peptide_masses[i])

    lS = [0]
    for i in range(len(peptide_masses)):
        for j in range(i + 1, len(peptide_masses) + 1):
            lS.append(prefixMass[j] - prefixMass[i])

    return sorted(lS)

def cyclicSpectrum(peptide, aaMass):
    prefixMass = [0]
    for i in range(len(peptide)):
        prefixMass.append(prefixMass[i] + aaMass[peptide[i]])

    peptideMass = prefixMass[-1]
    cS = [0]

    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            sub_mass = prefixMass[j] - prefixMass[i]
            cS.append(sub_mass)

            if i > 0 and j < len(peptide):
                cS.append(peptideMass - sub_mass)

    return sorted(cS)

def cyclicSpectrumMass(peptide_masses):
    prefixMass = [0]
    for i in range(len(peptide_masses)):
        prefixMass.append(prefixMass[i] + peptide_masses[i])

    peptideMass = prefixMass[-1]
    cS = [0]

    for i in range(len(peptide_masses)):
        for j in range(i + 1, len(peptide_masses) + 1):
            sub_mass = prefixMass[j] - prefixMass[i]
            cS.append(sub_mass)

            if i > 0 and j < len(peptide_masses):
                cS.append(peptideMass - sub_mass)

    return sorted(cS)

def BFCountPeptides(mass, aaMass):
    masses = sorted(set(aaMass.values()))

    dp = [0] * (mass + 1)
    dp[0] = 1

    for i in range(1, mass + 1):
        for a in masses:
            if i >= a:
                dp[i] += dp[i - a]

    return dp[mass]


def is_consistent(peptide, spectrum, aaMass):
    # peptide is a tuple of masses
    return not (Counter(linearSpectrum(peptide, aaMass)) - Counter(spectrum))

def is_consistent_mass(peptide_masses, spectrum):
    return not (Counter(linearSpectrumMass(peptide_masses)) - Counter(spectrum))


def expand(peptides, aaMass):
    masses = sorted(set(aaMass.values()))
    expanded = set()
    for peptide in peptides:
        for m in masses:
            expanded.add(peptide + (m,))
    return expanded

def expandExtended(peptides, masses):
    expanded = set()
    for peptide in peptides:
        for m in masses:
            expanded.add(peptide + (m,))
    return expanded

def mass(peptide):
    return sum(peptide)


def canonical(peptide):
    n = len(peptide)
    return min(peptide[i:] + peptide[:i] for i in range(n))


def cyclopeptideSequencing(spectrum, aaMass):
    parentMass = max(spectrum)
    candidates = {()}
    final = set()

    while candidates:
        candidates = expand(candidates, aaMass)

        for peptide in list(candidates):
            m = sum(peptide)

            if m == parentMass:
                if cyclicSpectrum(peptide, aaMass) == spectrum:
                    final.add(canonical(peptide))
                candidates.remove(peptide)

            elif m > parentMass or not is_consistent(peptide, spectrum, aaMass):
                candidates.remove(peptide)

    return final

def cyclopeptideSequencingMass(spectrum, aaMass):
    parentMass = max(spectrum)
    candidates = {()}
    final = set()

    while candidates:
        candidates = expand(candidates, aaMass)

        for peptide in list(candidates):
            m = sum(peptide)

            if m == parentMass:
                if cyclicSpectrumMass(peptide) == spectrum:
                    final.add(peptide)
                candidates.remove(peptide)

            elif m > parentMass or not is_consistent_mass(peptide, spectrum):
                candidates.remove(peptide)

    return final

def scoreCyclopeptide(peptide, spectrum, aaMass):
    pepspec = cyclicSpectrum(peptide, aaMass)
    pCount = Counter(pepspec)
    sCount = Counter(spectrum)

    similarity = sum((pCount & sCount).values())

    return similarity

def scoreCyclopeptideMass(peptide_masses, spectrum):
    pepspec = cyclicSpectrumMass(peptide_masses)
    pCount = Counter(pepspec)
    sCount = Counter(spectrum)
    similarity = sum((pCount & sCount).values())

    return similarity

def linearScoreCyclopeptideMass(peptide_masses, spectrum):
    pepspec = linearSpectrumMass(peptide_masses)
    pCount = Counter(pepspec)
    sCount = Counter(spectrum)
    similarity = sum((pCount & sCount).values())

    return similarity

def linearScoreCyclopeptide(peptide, spectrum, aaMass):
    pepspec = linearSpectrum(peptide, aaMass)
    pCount = Counter(pepspec)
    sCount = Counter(spectrum)
    similarity = sum((pCount & sCount).values())

    return similarity

def trimMass(leaderboard, spectrum, n):
    scored = [
        (peptide, linearScoreCyclopeptideMass(peptide, spectrum))
        for peptide in leaderboard
    ]

    scored.sort(key=lambda x: x[1], reverse=True)

    if len(scored) <= n:
        return [p for p, _ in scored]

    cutoff = scored[n-1][1]

    return {p for p, score in scored if score >= cutoff}

def trim(leaderboard, spectrum, n, aaMass):
    scored = [
        (peptide, linearScoreCyclopeptide(peptide, spectrum, aaMass))
        for peptide in leaderboard
    ]

    scored.sort(key=lambda x: x[1], reverse=True)

    if len(scored) <= n:
        return {p for p, _ in scored}
    cutoff = scored[n-1][1]

    return {p for p, score in scored if score >= cutoff}

def leaderboardCyclopeptideSequencing(spectrum, n, aaMass):
    leaderboard = {()}
    leaderPeptides = set()
    leaderScore = -1
    parentMass = max(spectrum)

    while leaderboard:
        leaderboard = expand(leaderboard, aaMass)

        for peptide in list(leaderboard):
            if mass(peptide) == parentMass:
                score = scoreCyclopeptideMass(peptide, spectrum)
                if score > leaderScore:
                    leaderPeptides = {peptide}
                    leaderScore = score
                elif score == leaderScore:
                    leaderPeptides.add(peptide)
            elif mass(peptide) > parentMass:
                leaderboard.remove(peptide)

        leaderboard = trimMass(leaderboard, spectrum, n)

    return leaderPeptides
def extendedLeaderboardCyclopeptideSequencing(spectrum, n, masses):
    leaderboard = {()}
    leaderPeptides = set()
    leaderScore = -1
    parentMass = max(spectrum)
    
    while leaderboard:
        leaderboard = expandExtended(leaderboard, masses)

        for peptide in list(leaderboard):
            if mass(peptide) == parentMass:
                score = scoreCyclopeptideMass(peptide, spectrum)
                if score > leaderScore:
                    leaderPeptides = {peptide}
                    leaderScore = score
                elif score == leaderScore:
                    leaderPeptides.add(peptide)
            elif mass(peptide) > parentMass:
                leaderboard.remove(peptide)

        leaderboard = trimMass(leaderboard, spectrum, n)

    return leaderPeptides

def convolution(spectrum):
    conv = []
    for i in range(len(spectrum)):
        for j in range(i + 1, len(spectrum)):
            diff = spectrum[j] - spectrum[i]
            if diff >= 57 and diff <= 200:
                conv.append(diff)
    return conv

def top_m_freq(convs, m):
    counts = Counter(convs)

    freqs = sorted(counts.values(), reverse=True)

    threshold = freqs[m-1]  if m <= len(freqs) else freqs[-1]

    return [mass for mass, c in counts.items() if c >= threshold]

def convolutionCyclopeptideSequencing(spectrum, m, n):
    convs = convolution(spectrum)
    frequent_masses = top_m_freq(convs, m)

    masses = sorted(frequent_masses)
    return extendedLeaderboardCyclopeptideSequencing(spectrum, n, masses)

