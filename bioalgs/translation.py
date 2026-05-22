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
