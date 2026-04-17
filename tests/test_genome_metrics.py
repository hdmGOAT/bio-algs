from bioalgs.genome_metrics import (
    ApproximatePatternCount,
    ApproximatePatternMatching,
    FasterSymbolArray,
    HammingDistance,
    MinimumSkew,
    SkewArray,
    SymbolArray,
)


def test_symbol_array_matches_faster_symbol_array():
    genome = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    assert SymbolArray(genome, "C") == FasterSymbolArray(genome, "C")


def test_skew_array_and_minimum_skew():
    genome = "CAGCCTA"
    assert SkewArray(genome) == [0, -1, -1, 0, -1, -2, -2, -2]
    assert MinimumSkew(genome) == [5, 6, 7]


def test_hamming_distance_known_example():
    assert HammingDistance("GGGCCGTTGGT", "GGACCGTTGAC") == 3


def test_approximate_pattern_matching_and_count():
    text = "ACGTACGTACGT"
    pattern = "ACG"
    d = 1
    expected_positions = [0, 4, 8]
    assert ApproximatePatternMatching(text, pattern, d) == expected_positions
    assert ApproximatePatternCount(pattern, text, d) == len(expected_positions)


def test_approximate_pattern_count_d_zero_equals_exact_count():
    text = "AAAAA"
    pattern = "AAA"
    assert ApproximatePatternCount(pattern, text, 0) == 3
