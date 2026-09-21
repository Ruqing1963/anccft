# anccft — Algorithmic Non-commutative Class Field Theory

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22883908.svg)](https://doi.org/10.5281/zenodo.22883908)
[![Text: CC BY 4.0](https://img.shields.io/badge/text-CC%20BY%204.0-lightgrey.svg)](LICENSE-TEXT)
[![Code: MIT](https://img.shields.io/badge/code-MIT-lightgrey.svg)](LICENSE)

**Volume I: Trees, Buildings, and Operator Algebras.** Hecke operators, sandpile groups and
operator algebras on Bruhat–Tits trees and Ã₂ / Ã_d buildings, with applications to condensed
matter, p-adic holography, network science and positive systems.

**Author:** Ruqing Chen (GUT Geoservice Inc., Rivière-Beaudette, Québec, Canada) — ruqing@hotmail.com
**Archived book:** Zenodo, DOI [10.5281/zenodo.22883908](https://doi.org/10.5281/zenodo.22883908)
**Status:** unrefereed; every numerical claim in the book is reproduced by a script in this repository.

## What is here

| path | content |
|------|---------|
| `monograph/` | the book: `main.tex` + `chapters/` (one chapter per paper, re-ordered by subject) + `frontmatter/` (preface, introduction, concordance, closing chapter "Solved, open, impossible") + merged `bibliography.tex` (130 external references, no self-citations); compiled `main.pdf` (224 pp, 26 chapters in six parts) |
| `papers/01` … `papers/20` | the twenty papers of the series: `.tex`, `.pdf`, verification script(s), archived output, ancillary README |
| `applications/` | the four applied notes of Part V — `phys-bloch-mass` (Phys 1), `holo-boundary-blind` (Holo 1), `net-tower-lambda` (Net 1), `ctrl-hidden-eigenvalues` (Ctrl 1) — same layout |
| `code/` | `build_monograph.py` (assembles the book from the papers, in this layout or the author's), `make_data_figures.py` (writes `data/` and `figures/`), and a README mapping every chapter to its script |
| `data/` | CSV tables of the verified numbers: the K₄ Iwahori and voltage towers, the composite-q towers, the Ã₂ complexes and their edge flow, the PG(d,q) counts behind the rank-d theorems, the cyclic-class witnesses, the Holo 1 and Net 1 tables |
| `figures/` | PDF figures drawn from `data/` |

The biological strand of the programme — the sandpile groups of de Bruijn graphs, double-stranded
and polyploid assembly, pan-genome graphs and transcript quantification (Bio 1–11) — is
**Volume II, *Algebraic Genomics: sandpile groups, cellular sheaves, and the topology of sequence
assembly***, published separately; it is not in this repository. The closing chapter of Volume I
points to the two questions Volume II answers.

## The mathematics in one paragraph

For a finite (q+1)-regular graph — an arithmetic quotient of the Bruhat–Tits tree of PGL₂ — the
sandpile (critical) group of the Hecke–Laplacian has order ∏(p+1−a_p(f))/n, a product of inverse
local L-values: Frobenius traces precipitate as K-theoretic torsion (Paper I). The boundary crossed
product is arithmetically blind (K₀ = Zʳ ⊕ Z/(r−1)); the Ihara companion cannot be realized by any
non-negative matrix (trace obstruction, Paper II), but a 2n×2n "shadow" graph realizes Z ⊕ Jac(Y)
as K₀ of a Kirchberg algebra (Paper VI), and the tail map of the Iwahori tower assembles these into
a Kirchberg inductive limit whose torsion K₀ is dual to the sandpile pro-module (Paper VIII), with
the kernel of the tail norm equal to the q-torsion for every q, prime or not (Paper XVII), and the
head map giving a C*-correspondence whose K₀-index is the head pushforward (Paper XVI). Abelian
towers obey exact Iwasawa control with a second-order exceptional zero whose L-invariant is
−‖h_α‖² (Papers IV, IX, XIV); no odd cyclic class computes it by a linear construction (Paper XIX);
the Iwahori tower has λ = −1 and is governed by a degeneracy calculus instead of Z_p[[T]]
(Papers III, V, VII, XIII). In rank two the Z/3 grading kills the trace obstruction, thick quotients
have b₁ = 0, and the geodesic edge flow's complement is √q·orthogonal, so the parahoric zeros of the
Ã₂ edge zeta are always critical (Papers X–XII, XV); in rank d the same holds with exponent
(d−1)/2, on an invariant subspace produced by a single flag count in PG(d,q) whose companion matrix
has the local Hecke polynomial of PGL_{d+1} as characteristic polynomial (Papers XVIII, XX). The
applied notes carry one theorem each into Bloch effective mass, p-adic AdS/CFT, tower normality
certificates and positive realizations.

## Reproducing

```
pip install -r requirements.txt
cd papers/10 && python pgl3_building.py        # Ã₂ complexes (10–20 min; imported by 11–13, 15–18, 20 and the applied notes)
cd papers/17 && python composite_tail_norm.py  # ~5 s, six checks on eight towers
cd papers/20 && python ad_invariant_subspace.py   # ~5 s, eight checks
cd code && python make_data_figures.py         # regenerates data/*.csv and figures/*.pdf
python code/build_monograph.py build monograph # regenerates the book's chapters, macros, bibliography, concordance
cd monograph && pdflatex main && pdflatex main && pdflatex main
```

Scripts locate their imports by relative path (`../10/pgl3_building.py`), so keep the folder
layout. Each paper folder's `README_ancillary.md` keys the script's checks to the paper's
verification table and states the runtime. Every check is labelled *exact* (integer/symbolic),
*sampled* (with a p-value) or *numeric* (with tolerance). The book compiles with 0 errors and
0 warnings (pdfTeX, MiKTeX 24.1).

## Reading order

The monograph's introduction and concordance give the logical order; the papers were written in
numerical order. Part I (Papers I, II, VI, III) is the foundation; Part II (IV, IX, XIV, XIX, V,
VII, XIII) the towers; Part III (VIII, XVII, XVI) the operator algebras of the tower; Part IV
(X, XI, XII, XV, XVIII, XX) higher rank; Part V the four applied notes; Part VI the closing
chapter, which lists every result as proved / verified / open / impossible, with the corrections
later papers made to earlier ones and the pointers to Volume II.

## Licence

Text of the papers and the monograph: CC BY 4.0 (`LICENSE-TEXT`). Code: MIT (`LICENSE`).

## Citing

Cite the book by its DOI:

> R. Chen, *Algorithmic Non-commutative Class Field Theory: Hecke operators, sandpile groups and
> operator algebras on trees and buildings, with applications*, 2026. Zenodo.
> https://doi.org/10.5281/zenodo.22883908

`CITATION.cff` carries the same record in machine-readable form. Inside the book the papers cite
one another as chapters; none of them has been published elsewhere, and the bibliography lists only
the external literature.
