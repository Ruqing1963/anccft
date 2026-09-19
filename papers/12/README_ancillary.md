# Ancillary bundle for Paper XII (`anccft12.tex`)

`pgl3_tower.py` reproduces every computational claim of the paper. It imports the complexes
and exact linear algebra of Paper X (`../Paper 10/pgl3_building.py`) and the edge/chamber
flows of Paper XI (`../Paper 11/pgl3_dynamical.py`); both must be present. Exact integer
arithmetic for all identities and Smith forms; numpy for floating spectra ("numeric").

Run: `python pgl3_tower.py` (about ten minutes). Requires `numpy`, `sympy`, `scipy`.
The archived output is `pgl3_tower_output.txt`.

Checks, keyed as in the paper's battery table:

| key | content |
|-----|---------|
| (B) | b_1 = 0 and H_1(Y;Z) exact for Y_21, Y_24, Y_42 (prime set certified by the gcd of maximal minors of the coordinate matrix); b_1 = 0 for Y_84, Y_168 |
| (A) | covering pushforwards along X_3 <- Y_21 <- Y_84 <- Y_168, Y_21 <- Y_42 <- Y_168, Y_24 <- Y_168 intertwine A_1, A_2, Delta^H; Jac^H at five levels; kernel orders = ratios; character mechanism (product of new |P_j(1)|) and its Hodge analogue |
| (D) | two-dimensional incidence calculus on Y_168: S T^t, T S^t, S S^t, T T^t, S M^t, T M^t, M M^t, d_1 = T - S, Delta_0, d_0+d_1+d_2 = boundary, sum d_i d_i^t, d_i^t d_j, K, L_E = T^t S - K, A_1 A_2 = A_2 A_1 |
| (E) | geodesic edge digraph E_1 on Y_168: q^2-regular Eulerian, strongly connected; spectral accounting of kappa(E_1); tower constants mu_2 = 588, chi_2 = 2 |
| (F) | chamber digraph F_1: q-regular Eulerian; constants mu_1 = 1764, chi = 1 |
| (L) | [III, Thm 2.1] with d = 4 on the complete digraph K_5 |

Chambers are listed canonically (first vertex of type 0) when the face maps d_i are formed.
