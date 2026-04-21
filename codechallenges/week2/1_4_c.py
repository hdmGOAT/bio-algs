from pathlib import Path
from bioalgs.genome_metrics import ApproximatePatternCount

p = Path(__file__).parent / "datasets" / "dataset_30278_6.txt"

pattern, text, d = p.read_text().splitlines()

print(ApproximatePatternCount(d=int(d), Pattern=pattern, Text=text))
