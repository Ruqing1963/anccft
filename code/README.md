# code/ — where every script lives and what it proves

Scripts stay beside their papers so that the relative imports the papers describe keep
working. This table maps chapters of the monograph to scripts.

| chapter | paper | script(s) | location | runtime |
|---|---|---|---|---|
| 2 | I | `hecke_snf.py`, `extra2.py`, `proofcheck.py`, `entropy.py`, `run_entropy.py` | `papers/04/ancillary/` | seconds–minutes |
| 3 | II | `tower.py`, `iwahori_tower.py` | `papers/04/ancillary/` | minutes |
| 4 | VI | `hecke_kirchberg.py` | `papers/06/` | seconds |
| 5 | III | `linedigraph_proof.py`, `trivialzero.py` | `papers/03/`, `papers/04/ancillary/` | seconds |
| 6 | IV | `limit_module.py` | `papers/04/` | seconds |
| 7 | IX | `exceptional_zero.py` | `papers/09/` | minutes |
| 8 | XIV | `cyclic_index.py` | `papers/14/` | ~1 min |
| 9 | V | `iw_stage1.py`, `iw_stage2.py`, `iw_stage3.py` | `papers/05/`, `papers/05/ancillary/` | minutes |
| 10 | VII | `iwahori_char.py` | `papers/07/` | seconds |
| 11 | XIII | `cohn_localization.py` | `papers/13/` | < 1 min |
| 12 | VIII | `tower_kirchberg.py` | `papers/08/` | 1–2 min |
| 13 | XVI | `head_correspondence.py` | `papers/16/` | seconds |
| 14 | X | `pgl3_building.py` (also the shared exact Smith-form library) | `papers/10/` | 10–20 min |
| 15 | XI | `pgl3_dynamical.py` | `papers/11/` | minutes |
| 16 | XII | `pgl3_tower.py` | `papers/12/` | ~10 min |
| 17 | XV | `pgl3_parahoric.py` | `papers/15/` | several min |
| 18 | Phys 1 | `bloch_mass.py` | `applications/phys-bloch-mass/` | seconds |
| 19 | Holo 1 | `boundary_blind.py` | `applications/holo-boundary-blind/` | seconds |
| 20 | Net 1 | `tower_lambda.py` | `applications/net-tower-lambda/` | seconds |
| 21 | Ctrl 1 | `hidden_eigs.py` | `applications/ctrl-hidden-eigenvalues/` | seconds |
| 22 | Bio 1 | `dbg_sandpile.py` | `applications/bio1-dbg-sandpile/` | ~2 min |
| 23 | Bio 2 | `rotor_assembly.py` | `applications/bio2-rotor-assembly/` | ~5 s |
| 24 | Bio 3 | `repeat_splitting.py` | `applications/bio3-repeat-splitting/` | ~30 s |
| 25 | Bio 4 | `multicopy.py` | `applications/bio4-multicopy/` | ~40 s |

Shared dependency: `papers/10/pgl3_building.py` provides exact integer Smith normal
forms (Bareiss determinant + p-adic elimination, torsion certified by
|Tor|·disc(ker) = pdet) and is imported by the scripts of XI, XII, XIII, XV, XVI and
Bio 1–4 through the relative path `../10/` or `../../papers/10/` — see the first lines of
each script. Bio 4 imports Bio 3.

`make_data_figures.py` tabulates the verified numbers into `data/*.csv` and draws
`figures/*.pdf`; it computes nothing new. `build_monograph.py` regenerates
`monograph/chapters/*.tex`, `macros.tex`, `bibliography.tex` and `frontmatter/concordance.tex`
from the paper sources (run it from the folder that contains the `Paper N` / `Bio N`
source folders, or adapt `ROOT`).

Archived outputs (`*_output.txt`) sit beside the scripts and are what the papers' tables
were typed from.
