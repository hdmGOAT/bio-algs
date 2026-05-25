from pathlib import Path

p = Path(__file__).parent / "datasets" / "dataset_30217_3.txt"

n = int(p.read_text().strip())

x = (n * (n + 1)) // 2 + 1
print(x)

