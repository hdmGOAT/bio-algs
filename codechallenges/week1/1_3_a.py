from pathlib import Path

from bioalgs.sequence_patterns import ReverseComplement


path = Path(__file__).parent / "datasets" / "dataset_30273_2.txt"

input = path.read_text('UTF-8').strip()


r = ReverseComplement(input)

with open("output.txt", "w") as f:
    f.write(r)
