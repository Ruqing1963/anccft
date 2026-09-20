# Ancillary bundle for Paper XVIII (`anccft18.tex`)

`ad_unitarity.py` reproduces every computational claim: 4 checks, 4 passing, **about 3 s**.
Requires `numpy`, and imports `pgl3_building.py` from the sibling `Paper 10` folder to
build the thick Ã₂ complexes; keep the folders as siblings. Run `python ad_unitarity.py`.
Archived output: `ad_unitarity_output.txt`. No data files, no network.

## The result in one paragraph

Paper XV proved `L_E L_Eᵀ = qI + (q²−q) TᵀT` for Ã₂ complexes and left rank d ≥ 3 open
(Rem. 6.1(d)). It generalises, with exponent (d−1)/2, and the proof is one line of
projective geometry. In an Ã_d complex the Z/(d+1) typing forces the tails of the edges
into a vertex y to be the **hyperplanes** of the link PG(d,q) and the successors to be its
**points**, with `{x,y,z}` a simplex iff the point z lies in the hyperplane x. Two
distinct hyperplanes of PG(d,q) always meet in codimension two — for *every* d — so
inclusion–exclusion gives

```
diagonal    = θ_d − θ_{d−1}                = q^d
off-diagonal= θ_d − 2θ_{d−1} + θ_{d−2}     = q^{d−1}(q−1)      (θ_j = (q^{j+1}−1)/(q−1))
```

independently of the pair, hence

```
L_E L_Eᵀ = q^{d−1} I + q^{d−1}(q−1) TᵀT ,   L_Eᵀ L_E = q^{d−1} I + q^{d−1}(q−1) SᵀS ,
T L_E = q^d S ,                              S L_Eᵀ = q^d T .
```

At d = 2 these are Paper XV's Theorem 2.1. Unconditionally: **exactly |E| − n singular
values of L_E equal q^{(d−1)/2}**, and L_E is a bijective q^{(d−1)/2}-isometry from ker S
onto ker T. On any W ⊆ ker S ∩ ker T invariant under L_E and L_Eᵀ, `R = L_E|_W` satisfies
`RRᵀ = RᵀR = q^{d−1}I`, so the parahoric zeros lie on `|u| = q^{−(d−1)/2}` with no
Ramanujan hypothesis.

## Checks

| key | content | arithmetic |
|-----|---------|-----------|
| (H) | the block identity `BBᵀ = BᵀB = q^{d−1}I + q^{d−1}(q−1)J` for the point/hyperplane non-incidence matrix of PG(d,q), **exhaustively** for d = 2, 3, 4, 5 and q = 2, 3, 5, 7, 11 with θ_d ≤ 1300 — 14 pairs, up to PG(5,3) on 364 points and PG(4,5) on 781. This is the entire content of the theorem, since `L_E L_Eᵀ` is block diagonal by head and each block is `BBᵀ` | exact |
| (C) | **why the typing is not a convention**: run the same construction on lines against planes in PG(3,q) and the off-diagonal takes *two* values — 9/10 (q=2), 32/33 (q=3), 144/145 (q=5) — because two lines meet or are skew. With intermediate dimensions the correction would be a sum over the Grassmann scheme's classes rather than a multiple of TᵀT, and the remainder would carry several moduli | exact |
| (A) | the global identities on the thick Ã₂ complexes Y₂₁, Y₂₄, Y₄₂ of Paper X, with d = 2: `L_E L_Eᵀ`, `L_Eᵀ L_E`, `T L_E = q^d S`, `S L_Eᵀ = q^d T`, and that L_E carries ker S isometrically onto ker T. The first is Paper XV Thm 2.1, so the conventions agree | exact |
| (S) | 126, 144 and 252 singular values equal √2 out of 147, 168 and 294 — exactly \|E\| − n in each case, as rank TᵀT = n forces | numeric |

## What is not done

* **The invariant subspace.** Corollary 2.7 of the paper is conditional on a
  W ⊆ ker S ∩ ker T invariant under L_E and L_Eᵀ. For d = 2 that is Paper XV Thm 2.1, got
  from `L_E Ψ = Ψ C'` with `Ψ = [Sᵀ Tᵀ Mᵀ]`; in rank d the pencil has degree d+1 and the
  analogue of M — the "third vertex" of a triangle, canonical only when a chamber has
  three vertices — has no evident substitute. What survives unconditionally is
  Corollary 2.6: L_E is already q^{(d−1)/2} times an isometry **from ker S onto ker T**,
  and only the passage to a single invariant space is missing.
* **The Ã₃ worked example** (the specification's third task). No thick Ã₃ complex is at
  hand. The rank-two test complexes came from an exhaustive search for triangle
  presentations over PG(2,2) (Paper X, Constr. 1.4); there is no comparably small search
  in rank three, and known constructions of finite Ã_d quotients go through arithmetic
  lattices in PGL_{d+1} over a local field (Lubotzky–Samuels–Vishne). Separate problem.
* **The exact Ramanujan certificate**, the other half of Paper XV Rem. 6.1(d), is
  untouched.

## Note on the specification

It proposed building the identity from "the intersection counts of k-dimensional with
(d−k)-dimensional subspaces of PG(d,q)". No such formula is needed and none would serve —
see check (C). It also asked for `W_d = ⋂_{i=0}^d ker S_i`; only two operators occur, and
`ker S ∩ ker T` already carries `RRᵀ = q^{d−1}I`. Its predicted constant `q^{d−1}` and
circle `|u| = q^{−(d−1)/2}` are **correct**.
