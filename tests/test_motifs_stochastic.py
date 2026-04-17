import pytest

from bioalgs.motifs_stochastic import (
    CountWithPseudocounts,
    GibbsSampler,
    GreedyMotifSearchWithPseudocounts,
    Motifs,
    Normalize,
    ProfileGeneratedString,
    ProfileMostProbablePattern,
    ProfileWithPseudocounts,
    RandomizedMotifSearch,
    RandomMotifs,
    WeightedDie,
)
from bioalgs.motifs_deterministic import Score


@pytest.fixture
def motifs():
    return ["ATGCA", "AAGCA", "ATGGA", "ATGCT"]


@pytest.fixture
def fake_dna_set():
    return [
        "TTACCTTAAC",
        "GATATCTGTC",
        "ACGGCGTTAG",
        "CCCTAACGAG",
        "CGTCAGAGGT",
    ]


def test_count_with_pseudocounts_and_profile_sums(motifs):
    counts = CountWithPseudocounts(motifs)
    assert counts["A"] == [5, 2, 1, 1, 4]
    assert counts["C"] == [1, 1, 1, 4, 1]
    assert counts["G"] == [1, 1, 5, 2, 1]
    assert counts["T"] == [1, 4, 1, 1, 2]

    profile = ProfileWithPseudocounts(motifs)
    for i in range(len(motifs[0])):
        col_sum = sum(profile[base][i] for base in "ACGT")
        assert col_sum == pytest.approx(1.0)


def test_profile_most_probable_pattern_and_motifs():
    profile = {
        "A": [0.2, 0.2, 0.3],
        "C": [0.3, 0.3, 0.2],
        "G": [0.3, 0.2, 0.3],
        "T": [0.2, 0.3, 0.2],
    }
    assert ProfileMostProbablePattern("ACGTAC", 3, profile) == "GTA"
    assert Motifs(profile, ["ACGTAC", "GTAACG"]) == ["GTA", "GTA"]


def test_greedy_motif_search_with_pseudocounts_valid_solution(fake_dna_set):
    k = 3
    t = len(fake_dna_set)
    initial = [strand[:k] for strand in fake_dna_set]
    result = GreedyMotifSearchWithPseudocounts(fake_dna_set, k, t)

    assert len(result) == t
    assert all(len(motif) == k for motif in result)
    assert Score(result) <= Score(initial)


def test_random_motifs_with_monkeypatched_randint(monkeypatch, fake_dna_set):
    import bioalgs.motifs_stochastic as ms

    monkeypatch.setattr(ms.random, "randint", lambda a, b: 0)
    result = RandomMotifs(fake_dna_set, 3, len(fake_dna_set))
    assert result == [strand[:3] for strand in fake_dna_set]


def test_normalize_and_weighted_die(monkeypatch):
    import bioalgs.motifs_stochastic as ms

    probs = Normalize({"AAA": 2.0, "CCC": 6.0, "GGG": 2.0})
    assert sum(probs.values()) == pytest.approx(1.0)
    assert probs == {"AAA": 0.2, "CCC": 0.6, "GGG": 0.2}

    monkeypatch.setattr(ms.random, "uniform", lambda a, b: 0.05)
    assert WeightedDie(probs) == "AAA"

    monkeypatch.setattr(ms.random, "uniform", lambda a, b: 0.7)
    assert WeightedDie(probs) == "CCC"


def test_profile_generated_string_returns_kmer_in_text(monkeypatch):
    import bioalgs.motifs_stochastic as ms

    profile = {
        "A": [0.25, 0.25, 0.25],
        "C": [0.25, 0.25, 0.25],
        "G": [0.25, 0.25, 0.25],
        "T": [0.25, 0.25, 0.25],
    }
    text = "ACGTAC"
    monkeypatch.setattr(ms.random, "uniform", lambda a, b: 0.0)

    sampled = ProfileGeneratedString(text, profile, 3)
    assert sampled in {"ACG", "CGT", "GTA", "TAC"}


def test_randomized_motif_search_and_gibbs_sampler_shapes(fake_dna_set):
    import random

    random.seed(7)
    k = 3
    t = len(fake_dna_set)

    randomized = RandomizedMotifSearch(fake_dna_set, k, t)
    gibbs = GibbsSampler(fake_dna_set, k, t, N=20)

    assert len(randomized) == t
    assert len(gibbs) == t
    assert all(len(motif) == k for motif in randomized)
    assert all(len(motif) == k for motif in gibbs)
