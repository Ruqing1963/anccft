# -*- coding: utf-8 -*-
r"""
ad_unitarity.py -- verification battery for

  "Algorithmic non-commutative class field theory, XVIII: the unitary remainder in rank d"

THE QUESTION (Paper XV, Rem. 6.1(d)).  For A~_2 complexes of order q the geodesic edge
flow satisfies  L_E L_E^t = qI + (q^2-q) T^t T  and  L_E^t L_E = qI + (q^2-q) S^t S, so on
a subspace inside ker T \cap ker S the remainder R = L_E| is sqrt(q) times an orthogonal
matrix and the parahoric zeros lie on |u| = q^{-1/2} unconditionally.  Does this survive
in rank d >= 3?

THE ANSWER.  Yes, with the exponent (d-1)/2, and the proof is one line of projective
geometry.  In an A~_d complex an edge x -> y has type(y) = type(x)+1; in the link of y,
which is the flag complex of PG(d,q), a neighbour of relative type i is a subspace of
vector dimension i.  So the TAILS x are hyperplanes and the SUCCESSORS z are points, and
{x,y,z} is a simplex exactly when the point z lies in the hyperplane x.  Therefore
L_E L_E^t is block diagonal by head, and its block at y is B B^t with B[x][z] = [z not in x].
Two distinct hyperplanes of PG(d,q) always meet in a subspace of codimension two, so
inclusion-exclusion gives, with theta_j = (q^{j+1}-1)/(q-1),

    diagonal   = theta_d - theta_{d-1}                     = q^d,
    off-diag   = theta_d - 2 theta_{d-1} + theta_{d-2}     = q^{d-1}(q-1),

both independent of the pair and of d.  Hence  B B^t = q^{d-1} I + q^{d-1}(q-1) J  on every
link, i.e.

    L_E L_E^t = q^{d-1} I + q^{d-1}(q-1) T^t T,    L_E^t L_E = q^{d-1} I + q^{d-1}(q-1) S^t S,

which is Paper XV's identity at d = 2 and gives R R^t = q^{d-1} I on ker T \cap ker S.

WHY THE TYPING IS FORCED.  Had the tails and successors been subspaces of intermediate
dimension, the count would not be constant: two lines of PG(3,q) meet or are skew, and the
corresponding off-diagonal takes two values differing by one (check (C)).  It is the
Z/(d+1) typing of an A~_d complex that puts hyperplanes against points, and only for that
pairing is the intersection dimension forced.

 (H)  the block identity B B^t = q^{d-1} I + q^{d-1}(q-1) J and its transpose, exhaustively
      in PG(d,q) for d = 2..5                                                   exact
 (C)  the same construction at intermediate dimensions is NOT uniform            exact
 (A)  on the thick A~_2 complexes of Paper X the global identity holds with the
      predicted constants, reproducing Paper XV, Thm 2.1                         exact
 (S)  the singular-value consequence: EXACTLY |E| - n singular values of L_E equal
      q^{(d-1)/2}, since rank T^t T = n                                       numeric

Run: python ad_unitarity.py     (a few seconds)
Imports pgl3_building.py from the sibling Paper 10 folder.
"""

import sys, os, time, itertools, collections
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "Paper 10"))
import pgl3_building as pb

ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


def theta(j, q):
    return (q ** (j + 1) - 1) // (q - 1)


# ------------------------------------------------------------------ projective space

def points(d, q):
    """normalised representatives of the points of PG(d,q), q prime"""
    out = []
    for v in itertools.product(range(q), repeat=d + 1):
        if not any(v):
            continue
        lead = next(i for i, x in enumerate(v) if x)
        if v[lead] != 1:
            continue
        out.append(v)
    return out


def noninc(d, q):
    """B[x][z] = 1 iff the point z does NOT lie in the hyperplane x (self-dual labelling)"""
    P = np.array(points(d, q), dtype=np.int64)
    return ((P @ P.T) % q != 0).astype(np.float64)


def lines_pg3(q):
    """the lines of PG(3,q), each as a frozenset of its points"""
    P = points(3, q)
    idx = {p: i for i, p in enumerate(P)}
    A = np.array(P, dtype=np.int64)
    out = set()
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            span = set()
            for a in range(q):
                for b in range(q):
                    if a == 0 and b == 0:
                        continue
                    v = (a * A[i] + b * A[j]) % q
                    lead = next(t for t, x in enumerate(v) if x)
                    inv = pow(int(v[lead]), q - 2, q)
                    span.add(tuple((inv * v) % q))
            out.add(frozenset(idx[s] for s in span))
    return sorted(out, key=sorted)


# --------------------------------------------------------------------- A~_2 complexes

def edge_flow(X):
    """L_E, S (tail incidence) and T (head incidence) of a complex with .edges and .tris"""
    edges = list(X.edges)
    ei = {e: i for i, e in enumerate(edges)}
    tris = {frozenset(t) for t in X.tris}
    outs = collections.defaultdict(list)
    for i, (a, b) in enumerate(edges):
        outs[a].append(i)
    E, n = len(edges), X.n
    L = np.zeros((E, E), dtype=np.float64)
    for i, (a, b) in enumerate(edges):
        for j in outs[b]:
            c = edges[j][1]
            if frozenset((a, b, c)) not in tris:
                L[i, j] = 1.0
    S = np.zeros((n, E))
    T = np.zeros((n, E))
    for i, (a, b) in enumerate(edges):
        S[a, i] = 1.0
        T[b, i] = 1.0
    return L, S, T


# ----------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    print("=== the unitary remainder in rank d ===")

    # --------------------------------------------------------------- (H) the identity
    print("\n=== (H) B B^t = q^{d-1} I + q^{d-1}(q-1) J in PG(d,q) ===")
    print("      d  q   points     diagonal   predicted    off-diagonal   predicted")
    ok, cases = True, 0
    for d in (2, 3, 4, 5):
        for q in (2, 3, 5, 7, 11):
            n = theta(d, q)
            if n > 1300:
                continue
            B = noninc(d, q)
            off_mask = ~np.eye(n, dtype=bool)
            for Mx in (B @ B.T, B.T @ B):
                dg = set(np.rint(np.diag(Mx)).astype(int).tolist())
                of = set(np.rint(Mx[off_mask]).astype(int).tolist())
                good = dg == {q ** d} and of == {q ** (d - 1) * (q - 1)}
                ok &= good
            cases += 1
            print(f"      {d}  {q:2d} {n:8d} {str(sorted(dg)):>12s} {q**d:>11d} "
                  f"{str(sorted(of)):>15s} {q**(d-1)*(q-1):>11d}")
    row("H", "in PG(d,q) the point/hyperplane non-incidence matrix satisfies "
             "B B^t = B^t B = q^{d-1} I + q^{d-1}(q-1) J, for every d and q tested -- the whole "
             "content of the rank-d unitarity identity", f"{cases} pairs (d,q)", "exact", ok)

    # ------------------------------------------------------- (C) intermediate dimensions
    print("\n=== (C) the same construction at intermediate dimensions is not uniform ===")
    ok2 = True
    for q in (2, 3, 5):
        P = points(3, q)
        n = len(P)
        A = np.array(P, dtype=np.int64)
        inc = ((A @ A.T) % q == 0)
        L = lines_pg3(q)
        BL = np.array([[0.0 if all(inc[pl, p] for p in line) else 1.0 for pl in range(n)]
                       for line in L])
        M = BL @ BL.T
        off = sorted(set(np.rint(M[~np.eye(len(L), dtype=bool)]).astype(int).tolist()))
        ok2 &= (len(off) > 1)
        print(f"      q={q}: PG(3,q) has {len(L)} lines and {n} planes; the off-diagonal of "
              f"B B^t takes the values {off}")
    row("C", "with lines against planes the off-diagonal takes two values, differing by one "
             "(the lines meet or are skew): only the hyperplane/point pairing forced by the "
             "Z/(d+1) typing gives a constant", "3 values of q", "exact", ok2)

    # ------------------------------------------------------------- (A) real A~_2 complexes
    print("\n=== (A) the global identity on the thick A~_2 complexes of Paper X ===")
    ok3 = True
    sv_rows = []
    for moduli, name in (((1, 1, 7), "Y_21"), ((2, 2, 2), "Y_24"), ((1, 2, 7), "Y_42")):
        try:
            X = pb.thick_cover(moduli)
        except Exception as exc:
            print(f"      {name}: could not build ({exc})")
            ok3 = False
            continue
        q, n = X.q, X.n
        L, S, T = edge_flow(X)
        E = L.shape[0]
        d = 2
        c0, c1 = q ** (d - 1), q ** (d - 1) * (q - 1)
        good1 = np.allclose(L @ L.T, c0 * np.eye(E) + c1 * (T.T @ T))
        good2 = np.allclose(L.T @ L, c0 * np.eye(E) + c1 * (S.T @ S))
        # the two mixed identities: T L = q^d S and S L^t = q^d T
        good3 = np.allclose(T @ L, (q ** d) * S) and np.allclose(S @ L.T, (q ** d) * T)
        # L carries ker S into ker T, bijectively and as a q^{(d-1)/2}-isometry
        ns = np.linalg.svd(S)[2][np.linalg.matrix_rank(S):].T          # basis of ker S
        img = L @ ns
        good4 = (np.allclose(T @ img, 0)
                 and np.allclose(img.T @ img, c0 * np.eye(ns.shape[1])))
        ok3 &= good1 and good2 and good3 and good4
        sv = np.linalg.svd(L, compute_uv=False)
        nhit = int(np.sum(np.abs(sv - q ** ((d - 1) / 2)) < 1e-8))
        sv_rows.append((name, n, E, nhit, E - n))
        print(f"      {name}: n={n:4d} vertices, |E|={E:5d} edges, q={q};  "
              f"L L^t {good1}, L^t L {good2}, T L = q^d S and S L^t = q^d T {good3}, "
              f"L: ker S -> ker T an isometry {good4}")
    row("A", "on Y_21, Y_24 and Y_42 the rank-d identities with d = 2 hold, the first being "
             "exactly Paper XV's Theorem 2.1, so the conventions agree", "3 complexes",
        "exact", ok3)

    # ---------------------------------------------------------- (S) singular values
    print("\n=== (S) singular values of L_E equal to q^{(d-1)/2} ===")
    ok4 = True
    for name, n, E, nhit, exact in sv_rows:
        good = (nhit == exact)
        ok4 &= good
        print(f"      {name}: {nhit} of {E} singular values equal sqrt(q); "
              f"rank T^t T = n forces exactly |E|-n = {exact}   {'ok' if good else 'FAIL'}")
    row("S", "since T^t T has rank n, the identity forces EXACTLY |E| - n singular values of "
             "L_E to equal q^{(d-1)/2}, with no hypothesis on the complex",
        f"{len(sv_rows)} complexes", "numeric", ok4)

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
