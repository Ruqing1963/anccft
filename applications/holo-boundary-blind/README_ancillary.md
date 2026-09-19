# Ancillary bundle for "What the boundary algebra of p-adic AdS/CFT does not see" (`boundary_blind.tex`)

`boundary_blind.py` reproduces every computational claim: 7 checks, 7 passing, about 0.3 s.
Requires `sympy`, reached through `pgl3_building.py` in the `Paper 10` folder of the ANCFT
series (exact integer Smith normal form); the path is resolved relative to this file, so
keep the two folders as siblings. Run `python boundary_blind.py`. The archived output is
`boundary_blind_output.txt`, and its first table is Table 1 of the paper verbatim.

Checks, keyed as in the paper's battery table (Table 3):

| key | content | arithmetic |
|-----|---------|-----------|
| (K) | K_0(C(dX) x Gamma) = Z^r + Z/(r-1) and K_1 = Z^r, computed as coker/ker of I - (A^nb)^T for the non-backtracking matrix on 2m oriented edges, on 11 standard quotients including q = 3 and q = 4. This is Robertson, Houston J. Math. 31 (2005), Thm 1 -- we verify, we do not claim it | exact |
| (B) | r - 1 = (q-1)n/2, so the boundary K-theory is a function of (q,n) alone | exact |
| (J) | witness families at (q,n) = (2,6), (2,8), (2,10): same boundary K-theory, different Jac(Y) | exact |
| (P) | **Theorem 3**: the two recorded 8-vertex cubic graphs have the same boundary K-theory (Z^5 + Z/4), the same von Neumann type (both non-bipartite, hence III_(1/q) by Robertson Thm 2), and coprime bulk orders 256 = 2^8 and 363 = 3*11^2. Hard-coded as `W_A`, `W_B` so the claim needs no search | exact |
| (S) | the spread inside one boundary class: 60 random cubic quotients at each n = 8,10,12,14 | sampled |
| (R) | the Busemann cocycle is q^Z-valued. **Normalisation check only** -- that the essential range is attained, which is what Krieger's theorem needs, is quoted from Ramagge-Robertson and Robertson, not verified | exact |

Conventions and one trap.

- `robertson_type(G)` returns III_(1/q^2) if the quotient graph is bipartite (equivalently
  Gamma is contained in PSL_2(F)) and III_(1/q) otherwise. This is Robertson, Thm 2 /
  Cor 1(2). **Do not state III_(1/q) unconditionally**: it is wrong for the bipartite case,
  which includes K_{3,3}, the cube, the Heawood graph and every even prism in Table 1.
- That value concerns the ORBIT relation of Gamma on the boundary. It is a different
  relation from the TAIL relation on the Cantor limit of a tower of covers, whose ratio set
  really is all of q^Z and which gives III_(1/q) unconditionally [ANCFT III, Thm 2.2].
  The two are easy to conflate and the paper flags this in Remark 1.
- The bulk channel is coker((q+1)I - A) = Z + Jac(Y) on n vertices; |Jac(Y)| is the number
  of spanning trees, and equals (1/n) prod_f (p+1-a_p(f)) when Gamma is the unit group of
  an Eichler order [ANCFT I].

Attribution. Checks (K) and (R) verify published theorems (Robertson 2005; Ramagge-
Robertson 1997) and are included as consistency checks on our conventions, not as claims.
The new content is check (P) and, through it, Theorem 3 of the paper.
