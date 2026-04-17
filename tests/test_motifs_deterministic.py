import pytest

from bioalgs.motifs_deterministic import (
    Consensus,
    Count,
    GreedyMotifSearch,
    Pr,
    Profile,
    ProfileMostProbableKmer,
    Score,
)


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


def test_count_profile_consensus_score(motifs):
    counts = Count(motifs)
    assert counts["A"] == [4, 1, 0, 0, 3]
    assert counts["C"] == [0, 0, 0, 3, 0]
    assert counts["G"] == [0, 0, 4, 1, 0]
    assert counts["T"] == [0, 3, 0, 0, 1]

    profile = Profile(motifs)
    assert profile["A"] == [1.0, 0.25, 0.0, 0.0, 0.75]
    assert profile["T"] == [0.0, 0.75, 0.0, 0.0, 0.25]

    assert Consensus(motifs) == "ATGCA"
    assert Score(motifs) == 3


def test_pr_multiplies_profile_column_probabilities(motifs):
    profile = Profile(motifs)
    assert Pr("ATGCA", profile) == pytest.approx(0.421875)


def test_profile_most_probable_kmer():
    profile = {
        "A": [0.2, 0.2, 0.3],
        "C": [0.3, 0.3, 0.2],
        "G": [0.3, 0.2, 0.3],
        "T": [0.2, 0.3, 0.2],
    }
    assert ProfileMostProbableKmer("ACGTAC", 3, profile) == "GTA"


def test_greedy_motif_search_returns_valid_and_improving_solution(fake_dna_set):
    k = 3
    t = len(fake_dna_set)
    initial = [strand[:k] for strand in fake_dna_set]
    result = GreedyMotifSearch(fake_dna_set, k, t)

    assert len(result) == t
    assert all(len(motif) == k for motif in result)
    assert Score(result) <= Score(initial)
