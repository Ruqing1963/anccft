# Ancillary files for *The reverse-complement quotient is a signed graph*

Ruqing Chen, 20 September 2026.

## Contents

| file | what it is |
|---|---|
| `rc_quotient_signed.tex` / `.pdf` | the note (9 pp.) |
| `bidirected_sandpile.py` | the verification battery, 8 checks |
| `bidirected_sandpile_output.txt` | the transcript of one full run |
| `phix174_NC_001422.1.txt` | bacteriophage phiX174, RefSeq NC_001422.1, 5386 bp, plain sequence |

The *E. coli* K-12 MG1655 genome (`NC_000913.3.fasta`) and its feature table
(`NC_000913.3.ft.txt`) are **not** duplicated here; the script reads them through
`ecoli_sandpile.load_genome()` / `load_features()` from the Bio 5 folder, which caches them
on first use.

## Running it

```
python bidirected_sandpile.py
```

About twelve minutes and roughly 3 GB of memory, almost all of it in check (E), which builds
the reverse-complement double cover of a 4.6 Mbp genome at three word lengths. The phiX174
checks (P), (W), (B), (R), (N) and the synthetic check (C) together take under a minute, so
commenting out `check_E()` gives a fast smoke test.

Dependencies: `numpy`, `sympy`. The script imports

* `rc_double_cover.py` from `../Bio 6` — the double cover `D_k`, `rc`, `compacted_multi`,
  `n_components`, `inverted_repeat`;
* `ecoli_sandpile.py` from `../Bio 5` — genome and feature loading, `compacted_graph`,
  `reduced_laplacian`, `elementary`, `gstr`;
* `pgl3_building.py` from `../Paper 10` — `coker_structure` (p-adic Smith form) and
  `bareiss_det`.

## The seven checks

| key | claim | mode |
|---|---|---|
| P | `rho` is free on `V(D_k)` exactly when no repeated k-mer is its own reverse complement: automatic for odd k, and false for phiX174 at k = 10 and k = 12, where the quotient adjacency is not symmetric | exact |
| W | every bundle of `Sigma_k` is monochromatic, and `K(Sigma_k)` is unchanged by switching | exact |
| B | `Sigma_k` is unbalanced iff `D_k` is connected iff S has an inverted repeat of length >= k | exact |
| R | the Reiner–Tseng sequence `0 -> K(Sigma_k) -> K(|D_k|) -> K(|Sigma_k|) -> 0` holds, orders multiply, and the 2-primary part is a non-split extension | exact |
| N | `P Delta(D_k) P = Delta^t` and never `Delta`, so the directed sandpile group sits in no covering sequence | exact |
| D | what `rho` gives instead is a **symmetric** Bowen–Franks linking form on `K(D_k)`; the relevant notion is a metabolizer, not a Lagrangian, and for phiX174 the order is not a square so none exists | exact |
| C | a signed r-bundle chain raises the r-rank by exactly s but does **not** split off `(Z/r)^s`: the series law | exact |
| E | *E. coli*'s seven rRNA operons become genuine 7-bundles of `Sigma_k`, longest chain 728 bp at k = 31, 51, 75 | exact |

"Exact" means exact integer arithmetic throughout — Smith normal form over **Z**, p-adic
valuations, Bareiss determinants, and sparse Gaussian elimination over **F**<sub>7</sub>.
Nothing in this note is numerical or sampled.

## Why this file uses `article` and not `amsart`

All seven Bio notes (Bio 1 through Bio 7) are set in `article`, whereas the main ANCCFT
series (Papers XV through XX) uses `amsart`. That split is deliberate and kept: the Bio
notes are a self-contained applied strand with its own citation and section conventions, and
switching one of them would break the uniformity of the strand rather than restore it. If
the Bio strand is ever reset as its own monograph, the whole strand moves together.

## Two conventions that matter

**The critical group of a signed graph.** Reiner and Tseng set
`K(Sigma) = im(d) / im(Delta)` where `im(d) = {x : sum x_v even}` is an index-two sublattice
of **Z**<sup>V</sup> for a connected unbalanced signed graph. This is **not**
`coker(Delta)`; using the cokernel doubles every order and changes the 2-primary part. The
function `K_signed` implements the Reiner–Tseng group from the explicit basis
`2e_0, e_0+e_1, ..., e_0+e_{n-1}`. Cho, Dochtermann, Inagaki, Oh, Snustad and Zacovic use a
third convention (delete a sink row and column); that group is not computed here.

**Which Laplacian.** `Delta(D_k)` in check (N) is the *out-degree* Laplacian of the
multidigraph, whose transposed reduced cokernel is the sandpile group. `Delta(|D_k|)` and
`Delta(Sigma_k)` are the ordinary and signed Laplacians of undirected multigraphs. The
whole point of check (N) is that the first is not obtainable from the other two.

## Cross-references to the rest of the series

* the double cover `D_k`, its balance, its connectivity criterion, and the annotation of the
  seven operons: Bio 6, *Double strands without a new lemma*;
* the bundle-chain lemma and the splitting theorem for digraphs, which
  Theorem 6 of this note shows has no signed analogue: Bio 3;
* the *E. coli* single-strand decomposition at k = 31, 51, 75, 100, 150: Bio 5;
* Bowen–Franks duality, invoked in Remark 5: Paper VIII, Lemma 1.5.
