# Ancillary bundle for "Effective mass from spanning trees, and the absence of Bloch phases on rank-two hyperbolic lattices" (`bloch_mass.tex`)

`bloch_mass.py` reproduces every computational claim of the paper: 90 checks, 90 passing,
about 13 s. Requires `sympy`, `numpy` and `mpmath`. Part II additionally imports
`pgl3_building.py` from the `Paper 10` folder of the ANCFT series (complex builder + exact
p-adic Smith normal form); the path is resolved relative to this file, so keep the two
folders as siblings. Run `python bloch_mass.py`. The archived output is
`bloch_mass_output.txt`.

Unit cells used in Part I: K4, K_{3,3}, Q_3 (cube), K5, Petersen, Heawood. Flux patterns:
a single link (d = 1), two links in independent directions (d = 2), and on K4/K5 the
multi-chord and triangle patterns of Papers IX and XIV for cross-validation.
Complexes used in Part II: Y_21, Y_24, Y_42, abelian covers of the three-site type quotient
of a triangle presentation over PG(2,2) compatible with a Singer cycle, deck groups Z/7,
(Z/2)^3, Z/14.

Checks, keyed as in the paper's battery table (Table IV):

| key | content | arithmetic |
|-----|---------|-----------|
| (S) | single-link flux: H = 1 - R_eff(e), and kappa_0 (1 - R_eff(e)) = tau(Y_0 - e) by an independent deletion-contraction spanning-tree count; Eq. (10), Thm. 2 | exact |
| (M) | Floquet determinant: Hess_k det Delta_k at k = 0 equals 2 kappa_0 H; Theta interpolated exactly as a Laurent polynomial over a rational grid, then differentiated; Eq. (9) | exact |
| (M) | the band itself: lambda_0(k) = k^T H k / n + O(\|k\|^4) at \|k\| = 1e-6 and 5e-7, 40-digit mpmath, one Richardson step; observed relative deviations 2e-27 .. 2e-26 | numeric |
| (Q) | quantization: kappa_0 * H is an integer matrix, and det(kappa_0 H) in Z; Cor. 3 | exact |
| (E) | Einstein relation: sigma^2 = H/m recomputed independently by a Green-Kubo sum over both orientations of every edge, not by substitution into the theorem; Eq. (13) | exact |
| (G) | the perturbation hypothesis: spectral gap lambda_1(0) > 0 and lambda_0(0) = 0 | numeric |
| (X) | cross-validation: H reproduces [IX, Table 1] values 1/2, 1, 3/2 on K4 chords (1,0,0), (1,1,1), (1,2,0) and [XIV, Table 1] value 7/5 on the K5 triangle | exact |
| (C) | loop flux: a cyclically oriented triangle on K5 gives H = 3 = the cycle length exactly, so m* = n/(2 t l); Cor. 4 | exact |
| (B) | rank two: the PG(2,2) triangle presentation and Singer cycle re-verified; b_1 = 0 by exact ranks and H_1(Y;Z) by exact p-adic Smith form on Y_21, Y_24, Y_42; Thm. 5, Table III | exact |

Conventions. Edges of K_n are listed as (a,b), a < b, oriented a -> b; the boundary is
d e = t(e) - o(e); the Bloch Laplacian is Delta_k = (q+1) I - A(k) with
A(k)_{uv} = sum_{e: u -> v} e^{i k . alpha(e)} + sum_{e: v -> u} e^{-i k . alpha(e)},
so Delta_k is Hermitian and Delta_{-k} = Delta_k^t. The single flux link is the first chord
of a BFS spanning tree, so the specific R_eff values in Table I depend on that choice of
link; the quantization statement (Q) does not.

Provenance discipline. Every floating-point row is labelled "numeric" in the script output;
every other row is exact integer or rational arithmetic. The two numeric rows are both
redundant: (M/numeric) confirms that the object differentiated in (M/exact) is indeed the
bottom band, and (G) confirms the non-degeneracy hypothesis of the perturbation argument.
No claim of the paper rests on a floating-point computation alone.
