# 3

from bioalgs.translation import linearScoreCyclopeptide, load_integer_mass_table, scoreCyclopeptide

from collections import Counter
from bioalgs.translation import convolution

pep = "MAMA"

spec = [0, 57, 71, 71, 71, 104, 131, 202, 202, 202, 256, 333, 333, 403, 404]
table = load_integer_mass_table()
print(scoreCyclopeptide(pep, spec, table))

pep = "PEEP"

spec =  [0, 97, 97, 97, 100, 129, 194, 226, 226, 226, 258, 323, 323, 355, 393, 452]

print(linearScoreCyclopeptide(pep, spec, table))


spectrum = [0, 57, 118, 179, 236, 240, 301]
convs = convolution(spectrum)
counts = Counter(convs)
max_multiplicity = max(counts.values()) if counts else 0
most_frequent = [m for m, c in counts.items() if c == max_multiplicity]

print(max_multiplicity)
print(sorted(most_frequent))
