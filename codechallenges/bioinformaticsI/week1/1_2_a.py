from pathlib import Path

from bioalgs import PatternCount


dataset_path = Path(__file__).parent / "datasets" / "dataset_30272_6.txt"
text, pattern = dataset_path.read_text(encoding="utf-8").splitlines()
print(PatternCount(text.strip(), pattern.strip()))
