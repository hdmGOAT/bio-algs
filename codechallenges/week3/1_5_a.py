from pathlib import Path

from bioalgs.motifs_deterministic import ProfileMostProbableKmer


p = Path(__file__).parent / "datasets" / "dataset_30305_3.txt"

lines = p.read_text().splitlines()

text = lines[0]
k = int(lines[1])

profile_rows = [
    list(map(float, line.split()))
    for line in lines[2:]
]

profile = {base: row for base, row in zip("ACGT", profile_rows)}

print(ProfileMostProbableKmer(text, k, profile))
