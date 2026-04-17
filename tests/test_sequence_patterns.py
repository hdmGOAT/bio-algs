from bioalgs.sequence_patterns import (
    Complement,
    FrequencyMap,
    FrequentWords,
    PatternCount,
    PatternMatching,
    ReverseComplement,
)


def test_pattern_count_overlapping_matches():
    assert PatternCount("ATATAT", "ATA") == 2


def test_frequency_map_counts_all_windows():
    text = "AACAAGCTGATAAACATTTAAAGAG"
    k = 2
    freq = FrequencyMap(text, k)
    assert sum(freq.values()) == len(text) - k + 1
    assert freq["AA"] == 6


def test_frequent_words_returns_all_max_patterns():
    text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    assert set(FrequentWords(text, 4)) == {"CATG", "GCAT"}


def test_complement_and_reverse_complement():
    assert Complement("ATGC") == "TACG"
    assert ReverseComplement("AAAACCCGGT") == "ACCGGGTTTT"


def test_pattern_matching_all_positions():
    genome = "GATATATGCATATACTT"
    assert PatternMatching("ATAT", genome) == [1, 3, 9]
