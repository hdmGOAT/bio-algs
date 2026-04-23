from pathlib import Path
from bioalgs.genome_metrics import MinimumSkew, SkewArray
from bioalgs.sequence_patterns import FindClumps, FrequentWordsWithMismatchesAndRC

p = Path(__file__).parent / "datasets" / "Salmonella_enterica.txt"

whole = p.read_text()

sequence = "".join(whole.splitlines()[1:])

mins = MinimumSkew(sequence)

ave = sum(mins) // len(mins)
candidate_win = sequence[ave-250:ave+250]

interesting = FrequentWordsWithMismatchesAndRC(candidate_win, 9, 1)

print(interesting)
