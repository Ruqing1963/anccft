# Ancillary bundle for Paper XIV (`anccft14.tex`)

`cyclic_index.py` reproduces every computational claim of the paper on the base graphs K4 (q = 2,
voltage assignments (1,0,0), (1,1,1), (1,2,0) of [II, Table 2]; the (1,0,0) l = 3 tower has the
same voltage and the same L-invariant) and K5 (q = 3, a single-edge voltage and a triangle voltage).
Self-contained: requires `sympy` and `mpmath` only. Run `python cyclic_index.py` (about a minute).
The archived output is `cyclic_index_output.txt`.

Checks, keyed as in the paper's battery table:

| key | content | arithmetic |
|-----|---------|-----------|
| (C) | dim Z^1_lambda(A) = dim HC^1(A) = 0 for A = C(V) + C(E^+) (linear system on idempotents) | exact |
| (N) | the specification's pairing at t -> 0+: Tr_s([D_p, alpha] D_p) = 4 sum_e alpha(e), linear in alpha | exact |
| (H) | heat-kernel formulas: lim_{t->oo} <alpha, e^{-t d^t d} alpha> and the Duhamel integral both equal \|\|h_alpha\|\|^2; the full function of t is printed | exact (symbolic t) |
| (B) | bottom Bloch band: lambda_0(phi) = (\|\|h_alpha\|\|^2/n) phi^2 + O(phi^4); exact perturbation coefficient and 40-digit Richardson numerics | exact + numeric |
| (S) | spectral action: (n/2)(1/t) d^2/dphi^2 Tr e^{-t D(e^{i phi})} at phi = 0 for t = 5, 10, 20, 40 versus L | numeric (40 digits, h = 1e-7) |
| (V) | diffusion constant sigma^2 = \|\|h_alpha\|\|^2/m from the Poisson equation; exact Var(S_N) for N <= 24 | exact |
| (A) | period formula a_c^t G^{-1} a_c, gauge invariance, det G = kappa_0, kappa_0 L in Z | exact |
| (T) | [IX, Thm 1.2(ii)] a_2 = -kappa_0 \|\|h_alpha\|\|^2 on K5 | exact |

Conventions follow Paper IX: edges of K_n listed as (a,b), a < b, oriented a -> b; K4 voltages
(a12, a13, a23) sit on edges (1,2), (1,3), (2,3); the boundary is d e = t(e) - o(e); the twisted
Laplacian is D(t) = (q+1) I - A(t) with A(t)_{uv} = sum_{e: u -> v} t^{alpha(e)}.
