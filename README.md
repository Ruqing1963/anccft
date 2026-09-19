# anccft — Algorithmic Non-commutative Class Field Theory

Hecke operators, sandpile groups and operator algebras on Bruhat–Tits trees and Ã₂ buildings, with applications to condensed matter, p-adic holography, network science, positive systems and genome assembly.

**Author:** Ruqing Chen (GUT Geoservice Inc., Rivière-Beaudette, Québec, Canada) — ruqing@hotmail.com
**Status:** unrefereed preprints and a compiled monograph; every numerical claim is reproduced by a script in this repository.

## What is here

| path | content |
|------|---------|
| `monograph/` | the book: `main.tex` + `chapters/` (one chapter per paper, re-ordered by subject) + `frontmatter/` (preface, introduction, concordance, closing chapter "Solved, open, impossible") + merged `bibliography.tex`; compiled `main.pdf` (226 pp) |
| `papers/01` … `papers/16` | the sixteen papers of the series: `.tex`, `.pdf`, verification script(s), archived output, ancillary README |
| `applications/` | the eight applied notes (Phys 1, Ctrl 1, Holo 1, Net 1, Bio 1–4), same layout; `phix174_NC_001422.1.txt` is the φX174 genome (NCBI NC_001422.1, 5386 bp) used by Bio 1–4 |
| `code/` | `build_monograph.py` (assembles the book from the papers), `make_data_figures.py` (writes `data/` and `figures/`), and a README mapping every script to its chapter |
| `data/` | CSV tables of the verified numbers (K₄ Iwahori tower, voltage towers, Ã₂ complexes, φX174 sandpile spectrum and skeleton profile, multi-copy census, …) |
| `figures/` | PDF figures generated from `data/` |

## The mathematics in one paragraph

For a finite (q+1)-regular graph — an arithmetic quotient of the Bruhat–Tits tree of PGL₂ — the sandpile (critical) group of the Hecke–Laplacian has order ∏(p+1−a_p(f))/n, a product of inverse local L-values: Frobenius traces precipitate as K-theoretic torsion (Paper I). The boundary crossed product is arithmetically blind (K₀ = Zʳ ⊕ Z/(r−1)); the Ihara companion cannot be realized by any non-negative matrix (trace obstruction, Paper II), but a 2n×2n "shadow" graph realizes Z ⊕ Jac(Y) as K₀ of a Kirchberg algebra (Paper VI) and the tail map of the Iwahori tower assembles these into a Kirchberg inductive limit whose torsion K₀ is dual to the sandpile pro-module (Paper VIII). Abelian towers obey exact Iwasawa control with a second-order exceptional zero whose L-invariant is −‖h_α‖² (Papers IV, IX, XIV); the Iwahori tower has λ = −1 and is governed by a degeneracy calculus instead of Z_p[[T]] (Papers III, V, VII, XIII). In rank two the Z/3 grading kills the trace obstruction, thick quotients have b₁ = 0, and the geodesic edge flow's complement is √q·orthogonal, so the parahoric zeros of the Ã₂ edge zeta are always critical (Papers X–XII, XV). The applied notes transport the sandpile/tower machinery to Bloch effective mass, p-adic AdS/CFT, tower normality certificates, positive realizations, and — most substantially — genome assembly, where the critical group of a de Bruijn graph splits exactly along repeats and, for two-copy repeats, is the critical group of a one-vertex map.

## Reproducing

```
pip install numpy sympy mpmath scipy matplotlib
cd papers/10 && python pgl3_building.py        # Ã₂ complexes (10–20 min; imported by 11, 12, 13, 15, 16 and by Bio 1–4)
cd applications/bio3-repeat-splitting && python repeat_splitting.py   # ~30 s
cd code && python make_data_figures.py         # regenerates data/*.csv and figures/*.pdf
cd monograph && pdflatex main && pdflatex main && pdflatex main
```

Scripts locate their imports by relative path (`../10/pgl3_building.py`, `../bio3-repeat-splitting/repeat_splitting.py`), so keep the folder layout. Each paper folder's `README_ancillary.md` keys the script's checks to the paper's verification table and states the runtime. Every check is labelled *exact* (integer/symbolic), *sampled* (with a p-value) or *numeric* (with tolerance).

## Reading order

The monograph's introduction and concordance give the logical order; the papers were written in numerical order. Part I (Papers I, II, VI, III) is the foundation; Part II (IV, IX, XIV, V, VII, XIII) the towers; Part III (VIII, XVI) the operator algebras; Part IV (X, XI, XII, XV) rank two; Part V the applications. The closing chapter lists every result as proved / verified / open / impossible, with the corrections later papers made to earlier ones.

## Licence

Text of the papers and the monograph: CC BY 4.0. Code: MIT. See `LICENSE` and `LICENSE-TEXT`.

## Citing

Repository: https://github.com/Ruqing1963/anccft — Zenodo DOI: *(inserted at the first release)*.
See `CITATION.cff`. Cite the monograph as
*R. Chen, Algorithmic Non-commutative Class Field Theory: Hecke operators, sandpile groups and operator algebras on trees and buildings, with applications, 2026, Zenodo.*
Inside the book the papers cite one another as chapters; none of them has been published elsewhere, and the bibliography lists only the external literature.
