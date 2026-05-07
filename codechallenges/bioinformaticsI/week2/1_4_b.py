from pathlib import Path

from bioalgs.genome_metrics import ApproximatePatternMatching


p = Path(__file__).parent / "datasets" / "dataset_30278_4.txt"

pattern, text, d = p.read_text('UTF-8').splitlines()

m = ApproximatePatternMatching(text, pattern, int(d))
with open("output.txt", "w") as file:
    file.write(" ".join(map(str, m)))
