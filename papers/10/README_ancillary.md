# Ancillary bundle for Paper X (`anccft10.tex`)

`pgl3_building.py` builds the test complexes and reproduces every computational claim of the
paper. Exact integer arithmetic throughout (Python integers and an in-house Smith
diagonalization); numpy is used only for floating-point spectra, and every such check is
labelled "numeric" in the output.

Run: `python pgl3_building.py` (10 to 20 minutes; the 294 x 294 Smith forms of Y_42 dominate).
Requires `sympy` and `numpy`. The archived output is `pgl3_building_output.txt`.

Complexes:

* thin A~_2 tori (q = 1): R^2 / mZ^2 with the equilateral tiling, m = 3, 6 (9 and 36 vertices);
* thick A~_2 complexes of order q = 2: abelian voltage covers of the 3-vertex type quotient of
  the triangle-presentation group over PG(2,2) with the Singer cycle lambda(x) = x+1, with
  deck groups Z/7, (Z/2)^3, Z/14, Z/2 x Z/2 x Z/7 and Z/2 x Z/2 x Z/14 (21, 24, 42, 84 and
  168 vertices). The triangle presentation is hard-coded (`ORBITS`) and its axioms are
  re-verified at start-up. Only Y_168 satisfies the flag condition (iv).

Smith forms are computed p-adically (Bareiss determinant, trial-division factoring,
elimination over Z/p^k with valuation-minimal pivots); torsion of singular symmetric
matrices is certified complete by |Tor| * disc(ker) = pdet.

Checks, keyed as in the paper's battery table:

| key | content |
|-----|---------|
| (P) | presentation axioms; simplicial; Heawood (resp. hexagon) links; regularity; flag condition Tr A_1^3 = 3F |
| (D) | d_1 d_2 = 0; Betti numbers; Euler characteristic |
| (H) | Smith forms of the Hodge Laplacians Delta_0, Delta_1, Delta_2; the DKM critical group Z_1 / im(d_2 d_2^t) |
| (L) | Hecke-Laplacian Delta^H: two-sided kernel, equal cofactors, |Jac^H| = (1/n) prod |P_j(1)| |
| (C) | explicit unimodular reduction U (I - C^(3)) V = P(1) + I + I |
| (T) | traces of the cubic companion, directly and by the Newton recursion; Tr C^3 = n(q-1)^2(q+1); net traces |
| (R) | spectrum of C^(3): Perron gap, Ramanujan test, threshold beyond which net traces are provably positive |
| (S) | signed shadows of Delta^H and Delta_1: unimodular identity, non-negativity, strong connectivity, double loops |

Conventions: vertices carry types in Z/3; every edge is oriented from type t to type t+1;
a triangle (a,b,c) has types (t,t+1,t+2) and boundary (a,b)+(b,c)+(c,a); A_1(x,y) = 1 iff
(x,y) is an edge, A_2 = A_1^t.
