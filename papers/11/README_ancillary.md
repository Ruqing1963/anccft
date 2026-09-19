# Ancillary bundle for Paper XI (`anccft11.tex`)

`pgl3_dynamical.py` reproduces every computational claim of the paper. It imports the complexes
and the exact linear algebra of Paper X (`../Paper 10/pgl3_building.py` must be present).
Exact integer arithmetic for all identities; numpy for floating spectra (labelled "numeric");
scipy's HiGHS for the linear-programming relaxation, whose infeasibility is then certified
exactly by a rational Farkas certificate where one was found.

Run: `python pgl3_dynamical.py` (a few minutes; the 3528 x 3528 chamber-flow spectrum dominates).
Requires `numpy`, `sympy`, `scipy`. The archived output is `pgl3_dynamical_output.txt`.

Checks, keyed as in the paper's battery table:

| key | content |
|-----|---------|
| (G) | the grading deg(i,v) = i - type(v) makes C^(3) a period-3 cyclic matrix; return maps P_d = C^3 restricted; spec P_0 = one cube per omega-orbit; det P_0 = q^(3n) |
| (Y) | sign data of P_0: negative entries, diagonal, Perron vectors (right positive, left not), no positive power, cycles in the negative digraph |
| (N) | LP relaxation of the nilpotent-cofactor conditions Y >= -P_0, Y P_0 <= 0, Tr Y = 0: infeasible; exact Farkas certificates for Y_21 (denominator 2) and Y_24 (denominator 128) |
| (L) | the cyclic lift of a shift equivalence (Lemma 2.2), checked on the trivial equivalence: the four equations hold with lag 3 |
| (E) | geodesic edge flow L_E on Y_168: row sums q^2, irreducible; 498 of 504 pencil roots in its spectrum, missing exactly the orbits of 1 and q; remainder of modulus q^(1/2); 1 is not an eigenvalue |
| (I) | the intertwiner Psi = [S^t T^t M^t] satisfies L_E Psi = Psi C' exactly; rank Psi = 3n - 6; C' = D C D^(-1) with D, C D^(-1) integral (lag-one shift equivalence) |
| (F) | chamber flow L_F on oriented chambers: irreducible, spectrum moduli 1, q^(1/4), q^(1/2), q |

Conventions as in Paper X: vertices typed in Z/3, edges oriented type t -> t+1, A_1 the
type-raising adjacency, A_2 = A_1^t; L_E((x,y),(y,z)) = 1 iff z is a type-raising neighbour
of y not adjacent to x; M(w,e) = 1 iff w is the third vertex of a chamber on the edge e.
