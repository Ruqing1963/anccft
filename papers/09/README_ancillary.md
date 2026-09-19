# Ancillary bundle for Paper IX (`anccft9.tex`)

`exceptional_zero.py` reproduces every computational claim of the paper on the base
graph K4 (q = 2, n = 4, r = 3, kappa(Y_0) = 16): the four abelian voltage towers of
[II, Table 2] and the Iwahori path tower at levels k <= 3. Exact arithmetic throughout
(sympy for symbolic determinants, resultants and series; a Bareiss determinant over
Python integers for the reduced Laplacians of derived graphs up to 128 vertices).

Run: `python exceptional_zero.py` (a few minutes). Requires `sympy`.
The archived output is `exceptional_zero_output.txt`.

Checks, keyed as in Table 1 of the paper:

| key | content |
|-----|---------|
| (Z) | Theta(0) = Theta'(0) = 0 and a_2 != 0: exceptional zero of order exactly two |
| (H) | a_2 = -kappa_0 * harmonic energy of the voltage cochain; L(Y_.) = a_2/kappa_0 |
| (H') | single-chord voltage: -a_2 = number of spanning trees of K4 avoiding the chord |
| (K) | critical partition function Z(1;t): double pole with leading coefficient n(q - 1/q)/L |
| (M) | Newton data (mu, lambda_h) of h = g/T^2; exceptional-zero regime iff ord_l a_2 = mu |
| (E) | defects eps_j for j <= 5 via resultants; finite support |
| (G) | growth law: direct determinant vs character product vs mu l^k + (lambda_h+1)k + nu_0 |
| (L) | nu_0 + ord_l L = (ord_l a_2 - mu) + sum eps_j; regime => nu_0 = -ord_l L |
| (B) | twisted Bass identity as a polynomial identity in (u, t) |
| (F) | functional equation; P(1/q, t) = q^-n Theta(t-1); derivatives at the critical point |
| (S) | rank D(1) = n-1, ord_T Theta = 2 |
| (I) | Iwahori tower: level-independent zeta determinant; residue constancy 9/8; kappa(Y_1) from level 0 |

Conventions: base edges of K4 ordered (0,1),(0,2),(0,3),(1,2),(1,3),(2,3) with voltages
(0,0,0,a12,a13,a23) as in Papers II, IV, VII; Iwahori levels as in Papers V, VII, VIII.
