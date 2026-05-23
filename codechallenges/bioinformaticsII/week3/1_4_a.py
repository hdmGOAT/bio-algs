from pathlib import Path


p = Path(__file__).parent / "datasets" /  "dataset_30215_3.txt"

n = int(p.read_text())
print(n * (n-1))
