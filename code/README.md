# code/ — where every script lives and what it proves

Scripts stay beside their papers so that the relative imports the papers describe keep
working. This table maps chapters of the monograph (Volume I) to scripts.

| chapter | paper | script(s) | location | runtime |
|---|---|---|---|---|
| 2 | I | `hecke_snf.py`, `extra2.py`, `proofcheck.py`, `entropy.py`, `run_entropy.py` | `papers/04/ancillary/` | seconds–minutes |
| 3 | II | `tower.py`, `iwahori_tower.py` | `papers/04/ancillary/` | minutes |
| 4 | VI | `hecke_kirchberg.py` | `papers/06/` | seconds |
| 5 | III | `linedigraph_proof.py`, `trivialzero.py` | `papers/03/`, `papers/04/ancillary/` | seconds |
| 6 | IV | `limit_module.py` | `papers/04/` | seconds |
| 7 | IX | `exceptional_zero.py` | `papers/09/` | minutes |
| 8 | XIV | `cyclic_index.py` | `papers/14/` | ~1 min |
| 9 | XIX | `cyclic_obstruction.py` | `papers/19/` | seconds |
| 10 | V | `iw_stage1.py`, `iw_stage2.py`, `iw_stage3.py` | `papers/05/`, `papers/05/ancillary/` | minutes |
| 11 | VII | `iwahori_char.py` | `papers/07/` | seconds |
| 12 | XIII | `cohn_localization.py` | `papers/13/` | < 1 min |
| 13 | VIII | `tower_kirchberg.py` | `papers/08/` | 1–2 min |
| 14 | XVII | `composite_tail_norm.py` | `papers/17/` | seconds |
| 15 | XVI | `head_correspondence.py` | `papers/16/` | seconds |
| 16 | X | `pgl3_building.py` (also the shared exact Smith-form library) | `papers/10/` | 10–20 min |
| 17 | XI | `pgl3_dynamical.py` | `papers/11/` | minutes |
| 18 | XII | `pgl3_tower.py` | `papers/12/` | ~10 min |
| 19 | XV | `pgl3_parahoric.py` | `papers/15/` | several min |
| 20 | XVIII | `ad_unitarity.py` | `papers/18/` | seconds |
| 21 | XX | `ad_invariant_subspace.py` | `papers/20/` | seconds |
| 22 | Phys 1 | `bloch_mass.py` | `applications/phys-bloch-mass/` | seconds |
| 23 | Holo 1 | `boundary_blind.py` | `applications/holo-boundary-blind/` | seconds |
| 24 | Net 1 | `tower_lambda.py` | `applications/net-tower-lambda/` | ~7 min (27 towers) |
| 25 | Ctrl 1 | `hidden_eigs.py` | `applications/ctrl-hidden-eigenvalues/` | seconds |

Shared dependency: `papers/10/pgl3_building.py` provides exact integer Smith normal
forms (Bareiss determinant + p-adic elimination, torsion certified by
|Tor|·disc(ker) = pdet) and is imported by the scripts of XI, XII, XIII, XV, XVI, XVII,
XVIII, XX and by the applied notes through the relative path `../10/` or `../../papers/10/`
— see the first lines of each script.

`make_data_figures.py` tabulates the verified numbers into `data/*.csv` and draws
`figures/*.pdf`; it computes nothing new, every number is copied from the archived
`*_output.txt` beside the script that produced it.

`build_monograph.py` regenerates `monograph/chapters/*.tex`, `macros.tex`,
`bibliography.tex` and `frontmatter/concordance.tex` from the paper sources. It finds the
sources in either layout — the author's working folders (`Paper N`, `Phys 1`, …) or this
repository (`papers/NN`, `applications/…`) — so from the repository root

```
python code/build_monograph.py build monograph
cd monograph && pdflatex main && pdflatex main && pdflatex main
```

reproduces the book. The hand-written front matter (`preface.tex`, `introduction.tex`,
`state.tex`) is never touched by the builder. Titles containing mathematics are wrapped in
`\texorpdfstring` for the PDF bookmarks; the papers' self-citations (`\cite{II}`,
`\cite{ANCFT9}`, …) become chapter references and never enter the merged bibliography.

Archived outputs (`*_output.txt`) sit beside the scripts and are what the papers' tables
were typed from.
