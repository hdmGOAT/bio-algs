from bioalgs.translation import cyclicSpectrum, findSubstringsThatEncodeAA, is_consistent, is_consistent_mass, linearSpectrum, load_codon_table, load_integer_mass_table

from collections import Counter
from functools import reduce
from operator import mul

two = "CCAAGUACAGAGAUUAACCCGAGGACCGAAAUCAACCCAAGAACAGAUAUCAAUCCUCGUACAGAAAUCAAC"
codon_table = load_codon_table()
print(findSubstringsThatEncodeAA(two, "PRTEIN", codon_table))



aa_counts = Counter(codon_table.values())

word = "LEADER"
answer = reduce(mul, (aa_counts[aa] for aa in word), 1)

print(answer)

cyclical_peptides = ["TAIM","MAIT","MTAI","TMLA","TLAM","TMIA"]
mass_table = load_integer_mass_table()
compatible = []

target = [0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]
print("target: " + str(target))
for p in cyclical_peptides:
    ps = cyclicSpectrum(p, mass_table)
    if ps == target:
        compatible.append(p)

print(compatible)

linear_peptides = ["QCV", "AQV", "TCE", "CTQ", "ETC", "CTV"]
compatible = []

target = [0, 71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231, 298, 303, 328, 330, 332, 333]
print("target: " + str(target))
for p in linear_peptides:
    if is_consistent(p, target, mass_table):
        compatible.append(p)

print(compatible)

