# rsomics-align-score

Pairwise sequence alignment scoring — Needleman-Wunsch (global) or
Smith-Waterman (local) with affine gaps. Reads a FASTA of sequences and emits,
for every pair, the optimal alignment score and percent identity.

```
rsomics-align-score seqs.fa [--local] [--match-score 1] [--mismatch -1] [--gap-open -2] [--gap-extend -1]
```

Output is tab-separated `seq1<TAB>seq2<TAB>score<TAB>identity`, one row per
unordered pair. `identity` is matching columns divided by alignment length.

## Origin

This crate scores alignments with the standard affine gap penalty, where a gap
of length `k` costs `gap_open + (k-1)*gap_extend` — the first gap position pays
only `gap_open`. Scores are value-exact to Biopython's
`Bio.Align.PairwiseAligner` (`mode="global"` for Needleman-Wunsch, `"local"`
for Smith-Waterman) with matching `match_score`, `mismatch_score`,
`open_gap_score`, and `extend_gap_score`, and to EMBOSS `needle`/`water` under
the same convention. Verified over hundreds of random global and local pairs;
see `tests/compat.rs`.

License: MIT OR Apache-2.0.
Upstream credit: Biopython (https://biopython.org/, BSD-3-Clause).
