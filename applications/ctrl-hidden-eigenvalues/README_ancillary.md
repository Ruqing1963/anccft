# Ancillary bundle for "Power-sum lower bounds on hidden eigenvalues in positive realizations" (`hidden_eigs.tex`)

`hidden_eigs.py` reproduces every computational claim: 7 checks, 7 passing, about 23 s.
Requires `sympy` and `mpmath` only. Unlike the other bundles in this project it does NOT
import anything from the `Paper N` folders; it is fully self-contained. Run
`python hidden_eigs.py`. The archived output is `hidden_eigs_output.txt`.

Checks, keyed as in the paper's battery table (Table 2):

| key | content | arithmetic |
|-----|---------|-----------|
| (T) | Theorem 1 on 120 random nonnegative INTEGER matrices: Lambda is taken to be the root multiset of a rational factor of the characteristic polynomial containing the Perron root, power sums are computed exactly by Newton's identities from that factor, and the bound is compared with the true hidden count. **The spectral radius is replaced by a rational LOWER bound, which makes the tested inequality strictly stronger than the theorem.** 0 violations; attained in 0 of 120, which is the expected weak generic behaviour and is reported as such in the paper. | exact + numeric rho |
| (K) | at k = 1 the bound coincides with Benvenuti's zeta [ELA 36 (2020), Thm 4.7] on pole sets {1} u mult*{-a}, a in {1/2, 9/10}, mult in {2,3,5} | exact |
| (I) | the k = 2 bound strictly exceeds Benvenuti's Thm 4.7 lower bound on M_m(theta) in 15/15 instances (3 thetas x 5 sizes); the table printed by the script is Table 1 of the paper, including the rows where the new lower bound meets Benvenuti's own Thm 4.4 upper bound and N is thereby determined | exact |
| (S) | Theorem 5: on M_m = {1} u m{i,-i} the bound gives N >= 4m for m = 1..32, and the direct sum of m copies of the 4-cycle permutation matrix (verified nonnegative, charpoly x^4 - 1) attains it | exact |
| (E) | Remark 7: for the poles {1, +-0.9i} of [Benvenuti-Farina tutorial, Ex. 8] the bound gives s >= 1, attained by the explicit rational nonnegative circulant with first row (0, 19/20, 0, 1/20), charpoly (x-1)(x+1)(x^2 + 81/100) | exact |
| (C) | Proposition 3: the bound never exceeds n on 300 random pole sets inside the unit disc, so the method cannot certify N > 2n; worst ratio bound/n observed 0.572 | exact |
| (Z) | appending zeros leaves every power sum unchanged, so no bound of this shape can exist in the Boyle-Handelman / Kim-Ormes-Roush category | exact |

Conventions. p_k(Lambda) is the k-th power sum of the pole multiset; rho is the spectral
radius; s = N - n is the number of hidden eigenvalues. Power sums are computed from monic
polynomial coefficients by Newton's identities (`power_sums`), never from floating-point
roots. Benvenuti's Thm 4.4 and Thm 4.7 bounds are EVALUATED by us from the published
statements (`benvenuti_lower`), not reproved; the paper says so.

Caution for anyone extending this. `sympy.nsimplify` must not be applied to the polynomial
coefficients: on rationals with large denominators it returns radical expressions and
silently corrupts the power sums. Use `sympy.expand` / `sympy.sympify`. This bit us once.
