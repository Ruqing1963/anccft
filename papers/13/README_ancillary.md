# Ancillary bundle for Paper XIII (`anccft13.tex`)

`cohn_localization.py` reproduces every computational claim of the paper on the K4 Iwahori
tower (q = 2, levels k = 1, 2, 3) and on the abelian voltage tower (1,1,1), l = 2. It imports
the exact linear algebra of Paper X (`../Paper 10/pgl3_building.py` must be present). All
arithmetic is exact except the pseudo-determinants marked "numeric".

Run: `python cohn_localization.py` (well under a minute). Requires `numpy`, `sympy`.
The archived output is `cohn_localization_output.txt`.

Checks, keyed as in the paper's battery table:

| key | content |
|-----|---------|
| (I) | incidence calculus of the degeneracy maps; A_k A_k^t != A_k^t A_k; degree compatibility (the four maps preserve the augmentation-zero lattices) |
| (S) | Delta_k restricted to ker(deg) is injective, invertible over Q, with determinant not a unit in Z_2 |
| (E) | det Delta_k^0 = n_k kappa_k (pseudo-determinant); ratio q^(n_k(q-1)) between consecutive levels; ord_2(n_k kappa_k) = n_k - 3 |
| (B) | 0 -> Z/n_k -> coker(Delta_k^0) -> Tor_k -> 0: orders, group structures, non-splitting; Delta_k = tau_*(tau^* - eta^*) |
| (A3) | abelian tower: pdet(Delta_k) = n kappa_0 prod_{chi != 1} det D(chi) exactly, k <= 3 |

Conventions follow Papers V and VII: vertices of Y_k are non-backtracking k-paths of K4,
tau = terminal truncation, eta = initial truncation, pushforwards sum over fibres, pullbacks
are transposes; ker(deg) has the basis e_i - e_0.
