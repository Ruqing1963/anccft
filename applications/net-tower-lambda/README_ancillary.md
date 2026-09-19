# Ancillary bundle for "The sign of lambda: an algebraic certificate for how a graph tower is generated" (`tower_lambda.tex`)

`tower_lambda.py` builds and measures 27 towers and runs 7 checks, all passing, in about
7 minutes (the cost is dominated by exact Bareiss determinants on the largest derived
graphs, up to a few hundred vertices). Requires `sympy`, reached through
`pgl3_building.py` in the `Paper 10` folder of the ANCFT series; the path is resolved
relative to this file, so keep the two folders as siblings. Run `python tower_lambda.py`.
The archived output is `tower_lambda_output.txt`.

Six base graphs (K_4, K_{3,3}, Q_3, Petersen, C_6 ring, C_4 x K_2 torus) and four
generating mechanisms.

Checks, keyed as in the paper's battery table (Table 4):

| key | content | arithmetic |
|-----|---------|-----------|
| (A) | abelian Z/2^k voltage towers: the law ord_l = mu*l^k + lambda*k + nu holds and lambda > 0, on 15 towers (6 bases x 3 voltage assignments, minus the disconnected covers). lambda values 1, 3, 5, 9. This VERIFIES published theory (McGown-Vallieres, Ann. Math. Quebec 48 (2024) 1-19; Gonet, Alg. Comb. 5 (2022) 827) — it is not a claim | exact |
| (L) | non-backtracking (line-digraph) towers: lambda = -1 on all 5 bases, and mu equals n(q+1)/q exactly in every case | exact |
| (E) | the -1 is an Euler characteristic: the one-step increment is n_k(q-1) - chi with chi = 1, verified level by level on K_4 (increments 11, 23, 47, 95) | exact |
| (K) | Prop. 3: the measured sequences equal the telescoped Knuth/Levine recursion term by term, on every base. **This is the honesty check**: it shows lambda = -1 is a consequence of Knuth, JCT 3 (1967) 309 / Levine, JCTA 118 (2011) 350, not a new computation | exact |
| (R) | rank-two (Z/2)^k towers: the ONE-variable law fails on 3/3, as a Z_l^2-extension requires (two-variable theory: Cuoco-Monsky; DuBose-Vallieres; Kleine-Mueller) | exact |
| (D) | dihedral D_(2^k) towers, i.e. non-abelian but still NORMAL: lambda = 3 > 0 on all 4. This is the finding that reframes the paper — the certificate detects normality, not commutativity | exact |
| (C) | the certificate assembled over the whole battery | exact |

Conventions and two traps.

- The fit uses the **last three levels**, not the first three. This is deliberate: the
  defects of the law are finitely supported [ANCFT IX], so early levels can be transient.
  Fitting from the head makes Petersen/"all ones" look like a violation when it is not.
  `show()` reports how many trailing levels the fit actually covers.
- A voltage assignment can produce a **disconnected** derived graph (then the spanning-tree
  count is 0 and ord is undefined). `connected()` screens for this; such levels are printed
  as "-" and excluded. Ignoring this silently produces nonsense.
- For mechanisms (a)-(c) the invariant is the spanning-tree count of an undirected graph;
  for (d), a digraph tower, it is the spanning-arborescence count. Different Laplacians;
  `spanning_trees` and `arborescences` respectively.

Attribution. Check (A) verifies established theory. Check (K) exists specifically to make
the classical origin of lambda = -1 explicit rather than letting a reader discover it. The
new content is (L) read as an Iwasawa law, (E), (D), and the certificate.
