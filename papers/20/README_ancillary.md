# Ancillary files for ANCCFT XX, *The invariant subspace in rank d*

Ruqing Chen, 20 September 2026.

| file | what it is |
|---|---|
| `anccft20.tex` / `.pdf` | the paper (6 pp.) |
| `ad_invariant_subspace.py` | the verification battery, 5 checks |
| `ad_invariant_subspace_output.txt` | the transcript of one full run |

## Running it

```
python ad_invariant_subspace.py
```

About one second. Dependencies: `numpy`, `sympy`. It imports `pgl3_building.py` from
`../Paper 10` for the three thick Ã₂ complexes `Y_21`, `Y_24`, `Y_42` and their coloured
adjacency operators `A_1`, `A_2`.

## The five checks

| key | claim | cases | mode |
|---|---|---|---|
| F | in PG(d,q), the number of complete flags with `F_{i+1} = V` and `F_1` inside a hyperplane `H` takes exactly **two** values, according as `V ⊆ H` or not, and the jump is independent of `V` and `H` | 14 triples (d,q,i) | exact |
| K | the closed forms `f(m) = prod [j]_q`, `phi(m,m-1) = [m-1]_q f(m)/[m]_q` and `K_i = f(d-i) f(i+1) q^i / [i+1]_q` reproduce the enumeration, and `K_i` is an integer | 14 triples | exact |
| I | the three rank-2 intertwining identities and their three mirrors | 3 complexes | exact |
| P | `im Psi_2` absorbs `L_E` and `L_E^t`; `L_E Psi = Psi C` exactly; `C` obeys the cubic Hecke pencil `x^3 - A_1 x^2 + q A_2 x - q^3` | 3 complexes | exact |
| U | every eigenvalue of `L_E` restricted to `W_2 = (im Psi_2)^perp` has modulus `q^{1/2}` | 3 complexes | numeric |

Check (F) is **exhaustive**: it enumerates every complete flag of PG(d,q) and every pair
`(V,H)`, for the five pairs `(d,q)` with at most 20 000 flags — `(2,2)`, `(2,3)`, `(3,2)`,
`(3,3)`, `(4,2)`. There is no sampling anywhere in the battery. Check (U) is the only
floating-point item; the eigenvalue moduli agree with `sqrt(q)` to about 1e-12.

## What is and is not verified

The theorems are proved for **every** d ≥ 2. The *global* identities are verified only at
d = 2, because no Ã_d complex with d ≥ 3 small enough to build was available; the proof of
the intertwining uses only the link lemma of Paper XVIII and the flag count of check (F),
and both hold in every rank. The corank of `Psi_d` is 6 on all three complexes tested at
d = 2 and is **not** proved to be 6 in general — that is the paper's main remaining open
item, since `dim W_d = |E| - rank Psi_d`.

## Relation to the literature

The determinant identity at d = 2, with the cubic Hecke pencil
`det(I - A_1 u + q A_2 u^2 - q^3 u^3)`, is Kang–Li, *Zeta functions of complexes arising from
PGL(3)*, Adv. Math. **256** (2014), 46–103. What this note adds is the **operator-level**
statement `L_E Psi = Psi C` (not in the literature at any rank, as far as a targeted search
found), its reduction to one flag count, and the consequent unconditional statement for every
d. For general rank the current state of the art is the alternating-product Ihara identity of
Kang–Yu, arXiv:2607.21262 (2026), which carries no Riemann-hypothesis statement.

## Cross-references

* the identity `L_E L_E^t = q^{d-1} I + q^{d-1}(q-1) T^t T` and the link lemma: Paper XVIII;
* the thick Ã₂ complexes and the cubic Hecke pencil in companion form: Paper X;
* the conditional Corollary 2.7 that this note discharges: Paper XVIII, Cor. 2.7.
