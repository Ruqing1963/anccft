# Ancillary bundle for Paper XV (`anccft15.tex`)

`pgl3_parahoric.py` reproduces every computational claim of the paper on the A~_2 complexes of Paper X:
Y_21, Y_24, Y_42, Y_168 (q = 2, abelian covers of the three-vertex complex X_3) and the thin torus T_6
(q = 1). It imports `../Paper 10/pgl3_building.py` (complex construction, companion matrix). Requires
`numpy`, `sympy`. Runtime: several minutes (dominated by the exact characteristic polynomials on Y_42 and the
modular determinants on Y_168). The archived output is `pgl3_parahoric_output.txt`.

Conventions. Typed edges e = (x -> y); L_E((x,y),(y,z)) = 1 iff (y,z) is an edge and {x,y,z} is not a chamber
(the chamber definition T^t S - K of Paper XII; on non-flag complexes it differs from the adjacency definition of
Paper XI, which is also computed for comparison). S, T, M are the tail, head and third-vertex incidences,
Psi = [S^t T^t M^t], W = (im Psi)^perp, R = L_E restricted to W, K = chamber continuation, P(u) the cubic pencil.
Pointed chambers are pairs (edge, third vertex); L_B is the Kang--Li operator (x,y;z) -> (y,z;w), w != x, which
is the chamber flow L_F of Paper XI; Phi = [sigma^* rho^* tau^*] pulls back edge functions along the three edges
of a pointed chamber; B = [[0,0,-I],[qI,K,K^t],[0,-I,0]] is the chamber companion.

| key | content | arithmetic |
|-----|---------|-----------|
| (D) | L_E L_E^t = qI + (q^2-q)T^tT, L_E^t L_E = qI + (q^2-q)S^tS; L_E Psi = Psi C', L_E^t Psi = Psi C'', K Psi = Psi Chat, K^t Psi = Psi Chat' | exact integers |
| (K) | rank Psi modulo two 31-bit primes plus the explicit rational kernel (certifies rank over Q); K(u) = det(I - uC' on ker Psi) | exact |
| (U) | R R^t = R^t R = qI on W; all eigenvalues of R of modulus sqrt(q) | numeric (theorem in the paper) |
| (E) | det(I - uL_E) K(u) = det P(u) G(u) in Z[u] via the type-cyclic reduction; deg G = dim W; functional equation of G; roots of G on the critical circle (Y_21, Y_24, Y_42, T_6 exact; Y_168 numeric) | exact / numeric |
| (B) | L_B Phi = Phi B exactly; det((1+qu^3)I + uK + u^2K^t) = (1-u^3)^n det(I-uL_E) det(I-u^2L_E)/det P(u) at 3\|E\|+1 points modulo two primes (certifies the identity in F_p[u]); det(I+uL_B) = det(I+uB) for q = 2; ratio (1-u^3)^(chi-n) on T_6 | modular |
| (V) | nontrivial pencil roots of modulus q (vertex Ramanujan property) | numeric |
