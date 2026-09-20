# Ancillary bundle for Paper XVII (`anccft17.tex`)

`composite_tail_norm.py` reproduces every computational claim: 6 checks, 6 passing,
**about 5 s**. Requires `sympy`, and imports `pgl3_building.py` from the sibling
`Paper 10` folder (Bareiss determinant, p-adic Smith form); keep the folders as siblings.
Run `python composite_tail_norm.py`. Archived output:
`composite_tail_norm_output.txt`. No network access, no data files.

## What the paper proves, and what the script is for

The paper closes [VIII, Rem. 5.1(b)] — the tail norm and the splitting of
`0 → X_∞^∨ → K_0(O_∞) → Z[1/q] → 0` at composite `q`. The proof is short and rests on one
new lemma plus three facts already published for all `q`:

* **Lemma 1.1.** `A_{L(G)} = Hᵀ T` for the head and tail incidence matrices, whose rows are
  non-empty and disjointly supported; so `rank A_{L(G)} = n` over **every** field.
* known for all `q`: `τ_* η* = A_k`, `η* τ_* = A_{k+1}`, `τ_* τ* = η_* η* = qI` [V, Prop. 1.2];
  `|ker τ_*| = q^{r_k}` with `r_k = n_k(q−1)−1` [V, Thm. 3.1(ii)]; `ker τ_* ⊆ Jac[q]`
  [VIII, Prop. 2.7].

Counting then forces `q^{r_k} = |ker τ_*| ≤ |Jac[q]| ≤ q^{r_k}`, so **`ker τ_* = Jac(Y_{k+1})[q]
≅ (Z/q)^{r_k}` for every q**, and — because the second inequality is also an equality —
**every invariant factor of `Jac(Y_{k+1})` divisible by `p | q` is divisible by `p^{v_p(q)}`**.
The script exists mainly to test that last statement, which is falsifiable.

Checks, keyed as in the paper's battery table:

| key | content | arithmetic |
|-----|---------|-----------|
| (I) | the degeneracy identities `τ_*η* = A_k`, `η*τ_* = A_{k+1}`, `τ_*τ* = η_*η* = qI` hold verbatim at composite `q` (8 towers, level 1→2). This is the input to `ker τ_* ⊆ Jac[q]` | exact |
| (N) | the line-digraph recursion `κ(Y_{k+1}) = q^{r_k} κ(Y_k)` by exact Bareiss determinants, 10 transitions — this is `|ker τ_*| = q^{r_k}` | exact |
| (R) | Lemma 1.1: `rank A_{k+1} = n_k` over F_p for p = 2, 3, 5, 7 and 2147483647 — i.e. for p dividing q and p not dividing q alike. 8 towers × 5 primes | exact |
| (P) | Corollary 1.2: the p-rank of `Jac(Y_{k+1})` is `n_k(q−1)−1` for every `p | q`, 12 levels | exact |
| (E) | **the falsifiable one.** Every invariant factor of `Jac(Y_{k+1})` divisible by `p` is divisible by `p^{v_p(q)}`; equivalently `|Jac[q]| = q^{r_k}`. On the q=4 towers the 2-adic valuations are `2^84, 4, 5^4` (K₆ level 2), all ≥ 2 among 359 factors (K₆ level 3) and among 149 (K₅,₅ level 2). One factor of valuation 1 anywhere would refute Theorem 2.1 | exact |
| (T) | the published prime-q values are recovered: `Jac(Y_2) = (Z/2)^8+Z/4+Z/16+Z/48`, `Jac(Y_3) = (Z/2)^12+(Z/4)^8+Z/8+Z/32+Z/96`, `ord₂κ = 18, 41` | exact |

## The towers

| base | q | levels | sizes |
|---|---|---|---|
| K₄ | 2 | 3 | n₁=12, n₂=24, n₃=48 |
| Petersen | 2 | 3 | 30, 60, 120 |
| K₅ | 3 | 3 | 20, 60, 180 |
| K₄,₄ | 3 | 2 | 32, 96 |
| **K₆** | **4** | **3** | **30, 120, 480** |
| **K₅,₅** | **4** | **2** | **50, 200** |
| K₇ | 5 | 2 | 42, 210 |
| K₈ | 6 | 2 | 56, 336 |

`Y_0` must be `(q+1)`-regular, so the smallest base for `q = 4` is `K₆`. The two `q = 4`
towers (non-squarefree) and the `q = 6` tower (squarefree composite) are the cases the
prime-only theorem did not cover; `q = 5` is a prime control on a base that is not `K₄`,
and Petersen and `K₄,₄` check that nothing depends on the base being complete.

Two thresholds in the script: `DET_MAX = 210` vertices for the exact Bareiss determinant
of check (N), and `SNF_MAX = 520` for the p-adic valuations. Raising either is only a
matter of time; `pb.padic_valuations` is `O(n³)` with entries bounded by `p^60`, and the
480-vertex level-three `K₆` case takes about a second.

## Why prime powers are the interesting case

For squarefree `q`, `Jac[q] = ⊕_{p|q} Jac[p]` and the theorem assembles the prime picture
by CRT. For `q = 4` it says something with no prime analogue: `Jac(Y_{k+1})` has **no
element of order exactly 2**, so `Jac[2] ⊊ Jac[4] = ker τ_*` with index `2^{r_k}` — the
kernel of the tail norm strictly contains the p-torsion. That is what check (E) measures.
