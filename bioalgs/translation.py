from pathlib import Path

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

def BFCountPeptides(mass, aaMass):
    masses = sorted(set(aaMass.values()))

    dp = [0] * (mass + 1)
    dp[0] = 1

    for i in range(1, mass + 1):
        for a in masses:
            if i >= a:
                dp[i] += dp[i - a]

    return dp[mass]
