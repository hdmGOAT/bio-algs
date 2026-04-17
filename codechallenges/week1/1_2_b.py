
from pathlib import Path

from bioalgs.sequence_patterns import FrequentWords


dataset_path = Path(__file__).parent / "datasets" / "dataset_30272_13.txt"
text, k = dataset_path.read_text(encoding="utf-8").splitlines()
print(FrequentWords(text, int(k)))
