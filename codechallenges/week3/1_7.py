from pathlib import Path

from bioalgs.sequence_patterns import DistanceBetweenPatternAndStrings

dataset = Path(__file__).parent / "datasets" / "dataset_30312_1.txt"

pattern, dna = dataset.read_text().splitlines()

print(DistanceBetweenPatternAndStrings(pattern, dna))
