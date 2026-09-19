# Ancillary verification scripts for the ANCFT series (Papers I-II)

Every numerical or symbolic claim in the two papers is reproducible from these
scripts (Python 3, sympy, numpy, scipy; no network access required).

| Script            | Verifies                                                            | Paper, location            |
|-------------------|---------------------------------------------------------------------|----------------------------|
| hecke_snf.py      | SNF of Delta_p, I-T_p, I-A^nb; Bass identity (symbolic)             | I, Table 1, Prop. 3.5      |
| extra2.py         | 11-graph sweep: Z^r + Z/(r-1) boundary pattern; Jac(Y) structures   | I, Table 1, Thm 3.7        |
| proofcheck.py     | coker(I-A^nb) = coker(B) reduction used in the proof of Thm 3.7     | I, Thm 3.7                 |
| entropy.py        | Kesten-McKay integral vs Lyons closed form; graph statistics        | I, Thm 3.14 (baseline)     |
| run_entropy.py    | LPS X^{5,l}, PG(2,5), random regular: torsion density vs h(q), B(q) | I, Thm 3.14, Rem. 3.15     |
| tower.py          | Four abelian Z/l^k voltage towers: exact Iwasawa laws               | II, Table 2, Sec. 5        |
| iwahori_tower.py  | Iwahori path tower: Levine recursion; (mu,lambda,nu)=(6,-1,-4)      | II, Rem. 5.2               |

Additional one-off checks (trace formula Tr(C^2) = -n(q-1) on all Table-1 graphs;
Newton-polygon test mu_meas = mu(g), lambda_meas = lambda(g)-1; the Kasparov-product
counterexample; the M-vs-C spectral comparison on Petersen) were run as inline
sessions and are reproducible from the formulas stated in the papers; the load-bearing
general facts are proved in the text and do not depend on any script.

Bibliography status (final): every entry across Papers I-IV is verified against
its primary source or Crossref record, with DOIs/verification notes recorded in the
.tex files; Washington (Paper IV) verified via the Springer book page (2nd ed.,
GTM 83, 1997; the Chapter-13 pointer is a recollection-based reading aid). The only
remaining item is non-bibliographic: the Boyle-Handelman / Kim-Ormes-Roush name
spellings in Paper II's contextual (non-load-bearing) remark, to be confirmed at
source. No theorem depends on any external citation except as marked in the papers.

New in this bundle: limit_module.py -- the verification battery for Paper IV
(block-circulant identity, pushforward=projection, deck compatibility, torsion-norm
surjectivity with exact kernel orders), seven levels across two towers, all exact.

Paper V battery: iw_stage1.py (incidence calculus, tail-norm and mirror
intertwinings, non-normality, level-0 defect, tau eta^* = T_p), iw_stage2.py
(degeneracy complex: kernel = Z*1, cokernel = Z torsion-free), iw_stage3.py
(explicit tail norm on 2-primary sandpile torsion, transition 2->1, kernel 2^11).
Run all: iw_stage_all.py. All checks exact on the K4 tower, levels <= 3.
