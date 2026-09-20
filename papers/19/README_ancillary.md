# Ancillary bundle for Paper XIX (`anccft19.tex`)

`cyclic_obstruction.py` reproduces every number: 5 checks, 5 passing, **about 1 s**.
Requires `sympy` only; self-contained, no data files, no network.
Run `python cyclic_obstruction.py`. Archived output: `cyclic_obstruction_output.txt`.

## The three answers

Paper XIV proved `HC¹ = 0` for the finite-dimensional commutative triple and asked
(Rem. 9.1(a), = Ch. 27 Open 6) whether an infinite-dimensional non-commutative algebra
would carry an odd cyclic class pairing to the L-invariant. It does not, for three
reasons of increasing finality.

**1. The proposed algebra is the compacts.** A connected Z-cover has vertex set `V × Z`
with the deck group translating the second coordinate — a **free** action — so

```
C₀(Ỹ_α) ⋊ Z  ≅  c₀(V) ⊗ (c₀(Z) ⋊ Z)  ≅  ⊕_{v∈V} K(ℓ²(Z)),
```

with `K₀ = Z^n`, **`K₁ = 0`**, `HC^odd = 0`. There is no odd K-theory there at all. The
specification's own hint ("该代数具有非平凡的 K₁ 和 HC¹") is false.

**2. The natural class in the right algebra does not exist.** Bloch theory lives in the
*commutant*, `M_n(C(S¹))` by Floquet transform, where `K₁ = Z`. Its natural element is the
voltage pencil `D(t) = (q+1)I − A(t)` — but `D(1)` is the ordinary Laplacian, which is
singular, so `[D] ∈ K₁` is undefined. On a small circle about `t = 1` the pencil *is*
invertible, and the odd pairing there is the winding number of `det D`, which equals the
**order** of the exceptional zero — **2 for every voltage class** — while `κ₀·L` is its
leading **coefficient**, `½Θ″(0) = −κ₀‖h_α‖²`. Index theory reads orders, not
coefficients.

**3. Bilinearity (Theorem 3.1, the actual obstruction).** An index pairing is additive in
each slot. If `α ↦ τ_α ∈ HC¹(A)` and `α ↦ [u_α] ∈ K₁(A)` are both linear and `K₁(A)` has
**rank one**, then `⟨τ_α,[u_α]⟩ = ℓ(α)·μ(α)` is a product of two linear forms; for
`b₁ ≥ 2` the kernel of `ℓ` is nonzero, so the product vanishes at some `α` with `[α] ≠ 0`,
where `κ₀‖h_α‖² > 0`. Both candidate algebras have rank-one `K₁` — `M_n(C(S¹))` and the
shadow algebra (Paper VI, Thm 2.1) — so **both are excluded for every base graph with
b₁ ≥ 2**, in particular K₄ (b₁ = 3) and K₅ (b₁ = 6).

## Results

| | b₁ | κ₀ | ‖h_α‖² | κ₀‖h_α‖² | cleared pencil [T⁰,T¹,T²] | winding of det D |
|---|---|---|---|---|---|---|
| K₄, α=(1,0,0) | 3 | 16 | 1/2 | 8 | [0, 0, **−8**] | **2** |
| K₄, α=(1,1,1) | 3 | 16 | 1 | 16 | [0, 0, **−16**] | **2** |
| K₄, α=(1,2,0) | 3 | 16 | 3/2 | 24 | [0, 0, **−24**] | **2** |
| K₅, single edge | 6 | 125 | 3/5 | 75 | [0, 0, **−75**] | **2** |
| K₅, triangle | 6 | 125 | 7/5 | 175 | [0, 0, **−175**] | **2** |

The `‖h_α‖²` and κ₀ values reproduce Paper XIV; the T² coefficients reproduce Paper IX.

## Checks

| key | content | arithmetic |
|-----|---------|-----------|
| (A) | the deck action is free with exactly \|V\| orbits, so the crossed product is a sum of copies of the compacts | exact |
| (P) | the cleared voltage pencil has a double zero at t = 1 with `[T²] = −κ₀‖h_α‖²` | exact |
| (W) | the winding number of `det D` on `\|t−1\| = 10⁻³` is 2 for every voltage class — the order, not the coefficient | numeric |
| (I) | `κ₀‖h_α‖² = a_cᵀ adj(G) a_c ∈ Z`: integrality is **Cramer's rule**, not an index theorem | exact |
| (Q) | `adj(G)` is positive definite of rank exactly b₁ (3 for K₄, 6 for K₅), hence not a product of two linear forms — the hypothesis that makes Theorem 3.1 bite | exact |

## What is *not* proved

That no odd cyclic class in any algebra computes L. Two escapes are named in Remark 3.3:
find an algebra whose `K₁` has rank ≥ b₁ together with a natural map from `H¹(Y;Z)` (the
pairing may then be a genuine quadratic form), or abandon linearity — in which case
`α ↦ τ_α` is part of the input rather than a cohomological construction, and the exercise
is circular. The first is a real question; we do not settle it.

Also corrected: there is nothing for an index to "absorb". `κ₀‖h_α‖² = a_cᵀ adj(G) a_c`
is an integer because the adjugate of an integer matrix is an integer matrix.
