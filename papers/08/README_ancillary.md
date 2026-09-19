# Ancillary bundle for Paper VIII (`anccft8.tex`)

`tower_kirchberg.py` reproduces every computational claim of the paper on the Iwahori
path tower over K4 (q = 2), levels k = 1, 2, 3 (n_k = 12, 24, 48; shadow matrices of
size 24, 48, 96). All arithmetic is exact (Python integers; sympy only for factoring,
a rational linear solve in check (P), and nothing else).

Run: `python tower_kirchberg.py` (about one to two minutes). Requires `sympy`.
The archived output is `tower_kirchberg_output.txt`.

Checks, keyed as in Table 2 of the paper:

| key | content |
|-----|---------|
| (S) | both unimodular factorizations of Theorem 1.3 with c = q-1 |
| (K) | SNF(I - B_k^t): one zero; torsion = Jac(Y_k); 2-adic orders 7, 18, 41; odd part Z/3 |
| (W) | the specification's c = q at a q-regular level: K_1 = 0, finite K_0 unrelated to Jac |
| (C) | tau^sharp is a surjective graph morphism, bijective on in-arcs, not on out-arcs |
| (U) | up-transfer diag(tau^*, tau^*): intertwining, unital, fixes the K_1 generator |
| (D) | down-transfer diag(tau_*, tau_*): intertwines I - B (Bowen-Franks side), not I - B^t; index q on ker(I - B) |
| (N) | induced map on torsion: surjective, kernel orders 2^11, 2^23, identical to tau_* on Jac |
| (Q) | Theorem 2.8: rank_F2 A_{k+1} = n_k; 2-rank of Jac(Y_{k+1}) = n_k - 1; tau_* kills Jac[2]; equality of orders |
| (P) | Pontryagin adjointness of tau_* and tau^* under the Q/Z pairing of Lemma 1.5 |
| (G) | Kirchberg witnesses at each level |
| (B) | plain tower: SNF(I - A_k^t) = Z^3 + Z/2 at k = 1, 2, 3 |
| (L) | degree scaling of tau^* (by q) and tau_* (by 1) |

Conventions follow Paper V (`iwahori_char.py`): vertices of Y_k are non-backtracking
k-paths of K4, arcs w[:-1] -> w[1:] for (k+1)-paths w, tau = terminal truncation.
