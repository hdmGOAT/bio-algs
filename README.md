This repository contains code I wrote for UC San Diego's Bioinformatics coursework.

## Functional Structure

- `bioalgs/sequence_patterns.py`: pattern counting, frequent words, complements, reverse complements, exact pattern matching.
- `bioalgs/genome_metrics.py`: symbol arrays, skew analysis, Hamming distance, approximate pattern matching/counting.
- `bioalgs/motifs_deterministic.py`: count/profile/consensus/score functions and greedy motif search.
- `bioalgs/motifs_stochastic.py`: pseudocount profiles, randomized motif search, Gibbs sampler, and weighted sampling helpers.
- `bioalgs/genome_assembly.py`: k-mers, de Bruijn graphs, Eulerian paths and cycles, string reconstruction, and contig generation.
- `bioalgs/translation.py`: RNA translation, peptide mass spectra, cyclopeptide sequencing, and spectral convolution.
- `bioalgs/__init__.py`: package exports for convenient imports.