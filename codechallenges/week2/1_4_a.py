from pathlib import Path

from bioalgs.genome_metrics import HammingDistance


p = Path(__file__).parent / "datasets" / "dataset_30278_3.txt"
text1, text2 = p.read_text('UTF-8').splitlines()
print(HammingDistance(text1, text2))
