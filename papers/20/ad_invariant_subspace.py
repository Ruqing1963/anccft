# -*- coding: utf-8 -*-
r"""
ad_invariant_subspace.py -- verification battery for

  "Algorithmic non-commutative class field theory, XX: the invariant subspace in rank d"

THE QUESTION (Paper XVIII, Rem. 4.1, Open (a)).  Paper XVIII proves that on an A~_d complex
of order q the geodesic edge flow satisfies L_E L_E^t = q^{d-1} I + q^{d-1}(q-1) T^t T and
its dual, so that on ANY subspace W inside ker S \cap ker T invariant under L_E and L_E^t the
restriction is q^{(d-1)/2} times an orthogonal matrix.  It leaves the existence of such a W
open, recording only a candidate intertwiner
      Psi_d = [ M_0^t  M_1^t  ...  M_d^t ],     M_i(v,e) = #{chambers C containing e
                                                  with C_{type(head e)+i} = v},
whose first and last blocks are T and S up to the constant number of chambers through an
edge.  Does im Psi_d absorb L_E and L_E^t?

THE ANSWER.  Yes, and the whole content is one flag count in PG(d,q).

THE MECHANISM.  Fix an edge e = (x,y).  In the link of y, which is the flag complex of
PG(d,q), the tail x is a HYPERPLANE H and the geodesic successors z are the points NOT on H
(Paper XVIII, Lem. 2.2).  Evaluating (L_E M_i^t)(e,v) sums M_i(v,(y,z)) over those z, i.e.
counts pairs (z, C) with C a chamber through {y,z} whose type-(type(y)+1+i) vertex is v and
z off H.  Writing V for the (i+1)-dimensional subspace of the link attached to v, the sum
over ALL z gives every flag with F_{i+1} = V, and the correction subtracts those with
F_1 \subseteq H.  So everything reduces to

    N(V,H) := #{complete flags F of PG(d,q) with F_{i+1} = V and F_1 \subseteq H},

and the point is that N takes exactly TWO values, according as V \subseteq H or not, with a
difference independent of V and of H:

    V \subseteq H :   N = f(i+1) f(d-i),
    V \not\subseteq H : N = phi(i+1,i) f(d-i),        f(m) = prod_{j<=m} [j]_q ,
    difference    K_i = f(d-i) f(i+1) q^i / [i+1]_q .

Since "V \subseteq H" is exactly the condition that {x,y,v} be a simplex, which is what
M_{i+1} records, this says

    L_E M_i^t = K_i ( T^t A_{i+1} - M_{i+1}^t / c_{i+1} ),     0 <= i <= d-1,
    L_E M_d^t = q^d M_0^t                                     (Paper XVIII, Prop. 2.4),

so im Psi_d is L_E-stable; the arrow-reversing argument gives L_E^t.  Hence
W_d := (im Psi_d)^perp is contained in ker S \cap ker T, is invariant under both, and
Paper XVIII's Corollary 2.7 holds with no hypothesis at all.

NORMALISATION.  Writing N_i for the 0/1 indicator "{x,y,v} lies in a chamber and v has
relative type i", so that M_i = c_i N_i with c_i = f(i) f(d-i) the number of chambers
through such a triple, the constant collapses to a bare power:  K_i / c_i = q^i.  The
recursion is then

    L_E N_i^t = q^i ( N_0^t A_{i+1} - N_{i+1}^t ),    N_0 = T, N_d = S, N_{d+1} = 0,

whose companion matrix has characteristic polynomial  sum_j (-1)^j q^{j(j-1)/2} A_j
x^{d+1-j}, the local Hecke polynomial of PGL_{d+1}: the cubic pencil of Kang-Li at d = 2,
and at d = 3 the quartic  x^4 - A_1 x^3 + q A_2 x^2 - q^3 A_3 x + q^6, obtained with no
global A~_3 complex.

 (F)  the flag count N(V,H) takes exactly two values in PG(d,q), with difference K_i,
      exhaustively over all (V,H) for d = 2,3,4 and q = 2,3                        exact
 (K)  the closed forms f(m), phi(m,m-1) and K_i, checked against the enumeration    exact
 (A)  K_i / c_i = q^i, symbolically in q, for 2 <= d <= 7 and every i              exact
 (H)  det(xI - C) is the local Hecke polynomial of PGL_{d+1}, for d = 2,3,4,5       exact
 (I)  on the thick A~_2 complexes of Paper X the three global identities hold:
      L_E T^t = T^t A_1 - M^t, L_E M^t = q(T^t A_2 - S^t), L_E S^t = q^2 T^t,
      together with their three mirrors                                            exact
 (P)  im Psi_2 absorbs L_E AND L_E^t, and L_E Psi = Psi C holds exactly             exact
 (U)  every eigenvalue of L_E restricted to W_2 = (im Psi_2)^perp has modulus
      q^{1/2}, unconditionally                                                    numeric
 (N)  N_i^t 1_{t+i} = b_i N_0^t 1_t with b_i = [d choose i]_q, forcing
      corank Psi_d >= d(d+1); at d = 2 the corank is exactly 6 and these span       exact

Run: python ad_invariant_subspace.py      (under ten seconds)
Imports pgl3_building.py from the sibling Paper 10 folder.
"""

import sys, os, time, itertools, collections
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "Paper 10"))
import pgl3_building as pb

ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# ------------------------------------------------------------------ q-combinatorics

def gauss(j, q):
    """[j]_q = (q^j - 1)/(q - 1), the number of points of a j-dimensional space"""
    return (q ** j - 1) // (q - 1)


def fact(m, q):
    """f(m) = number of complete flags of an m-dimensional F_q-space"""
    out = 1
    for j in range(1, m + 1):
        out *= gauss(j, q)
    return out


def K(d, i, q):
    """the predicted jump K_i = f(d-i) f(i+1) q^i / [i+1]_q"""
    num = fact(d - i, q) * fact(i + 1, q) * q ** i
    assert num % gauss(i + 1, q) == 0, "K_i is not an integer"
    return num // gauss(i + 1, q)


# ------------------------------------------------------------------ linear algebra over F_q

def subspaces(m, q):
    """all subspaces of F_q^m, as frozensets of vectors, grouped by dimension"""
    vecs = [tuple(v) for v in itertools.product(range(q), repeat=m)]
    zero = tuple([0] * m)
    bydim = collections.defaultdict(set)
    bydim[0].add(frozenset([zero]))
    cur = {frozenset([zero])}
    for dim in range(1, m + 1):
        nxt = set()
        for U in cur:
            for v in vecs:
                if v in U:
                    continue
                span = set()
                for u in U:
                    for a in range(q):
                        span.add(tuple((ui + a * vi) % q for ui, vi in zip(u, v)))
                nxt.add(frozenset(span))
        bydim[dim] = nxt
        cur = nxt
    return bydim


def flags(m, q, bydim):
    """all complete flags F_1 subset ... subset F_{m-1} of F_q^m, as tuples of subspaces"""
    if m == 0:
        return [()]
    out = [(U,) for U in bydim[1]]
    for dim in range(2, m):
        out = [F + (U,) for F in out for U in bydim[dim] if F[-1] <= U]
    return out


# --------------------------------------------------------------------- A~_2 complexes

def pieces(X):
    """L_E, S, T and M for a thick A~_2 complex"""
    edges = list(X.edges)
    ei = {e: i for i, e in enumerate(edges)}
    tris = {frozenset(t) for t in X.tris}
    outs = collections.defaultdict(list)
    for i, (a, b) in enumerate(edges):
        outs[a].append(i)
    E, n = len(edges), X.n
    L = np.zeros((E, E))
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
    M = np.zeros((n, E))
    for (a, b, c) in X.tris:            # types t, t+1, t+2: all three typed edges
        M[c, ei[(a, b)]] += 1.0
        M[a, ei[(b, c)]] += 1.0
        M[b, ei[(c, a)]] += 1.0
    return L, S, T, M, n, E


# ------------------------------------------------------------------------------- checks

def check_FK():
    print("\n=== (F) the flag count N(V,H) in PG(d,q) takes exactly two values ===")
    print("      d  q  i   N (V in H)  predicted   N (V not in H)  predicted    jump  K_i")
    okF = okK = True
    cases = 0
    for d in (2, 3, 4):
        for q in (2, 3):
            m = d + 1
            if fact(m, q) > 20000:
                continue
            bydim = subspaces(m, q)
            FL = flags(m, q, bydim)
            assert len(FL) == fact(m, q), f"flag count {len(FL)} != {fact(m, q)}"
            H_all = sorted(bydim[m - 1], key=sorted)
            for i in range(0, d):
                byV = collections.defaultdict(list)
                for F in FL:
                    byV[F[i]].append(F[0])          # F[i] is F_{i+1}; F[0] is F_1
                inn, out = set(), set()
                for V, firsts in byV.items():
                    for H in H_all:
                        N = sum(1 for p in firsts if p <= H)
                        (inn if V <= H else out).add(N)
                pin = fact(i + 1, q) * fact(d - i, q)
                pout = pin - K(d, i, q)
                good = (inn == {pin}) and (out == {pout})
                okF &= good
                okK &= (pin - pout == K(d, i, q))
                cases += 1
                print(f"      {d}  {q}  {i}  {str(sorted(inn)):>11s} {pin:10d}  "
                      f"{str(sorted(out)):>14s} {pout:10d}  {pin-pout:7d} {K(d,i,q):4d}"
                      + ("" if good else "   MISMATCH"))
    row("F", "in PG(d,q) the number of complete flags with F_{i+1} = V and F_1 inside a "
             "hyperplane H depends only on whether V lies in H, and the jump between the two "
             "values is independent of V and H -- the whole content of the rank-d "
             "intertwining", f"{cases} triples (d,q,i)", "exact", okF)
    row("K", "the closed forms f(m) = prod [j]_q, phi(m,m-1) = [m-1]_q f(m)/[m]_q and "
             "K_i = f(d-i) f(i+1) q^i / [i+1]_q reproduce the enumeration exactly, and K_i "
             "is an integer for every (d,i,q)", f"{cases} triples", "exact", okK)
    return okF and okK


def check_HA():
    """(A) kappa_i = K_i/c_i = q^i ; (H) det(xI-C) is the PGL_{d+1} Hecke polynomial"""
    q = sp.Symbol("q", positive=True)
    gq = lambda j: sum(q ** t for t in range(j))
    fq = lambda m: sp.prod([gq(j) for j in range(1, m + 1)])
    print("\n=== (A) the constant collapses: K_i / c_i = q^i with c_i = f(i) f(d-i) ===")
    okA, cases = True, 0
    for d in range(2, 8):
        for i in range(d + 1):
            kap = sp.simplify(fq(d - i) * fq(i + 1) * q ** i / gq(i + 1)
                              / (fq(i) * fq(d - i)))
            okA &= (sp.simplify(kap - q ** i) == 0)
            cases += 1
    print(f"      K_i / c_i = q^i symbolically in q, for 2 <= d <= 7 and every i "
          f"({cases} pairs): {okA}")
    row("A", "the coefficient in the intertwining is the bare power q^i, because a chamber "
             "through {x,y,v} is a flag through a chain of dimensions i < d, so "
             "c_i = f(i) f(d-i)", f"{cases} pairs (d,i)", "exact", okA)

    print("\n=== (H) det(xI - C) is the local Hecke polynomial of PGL_{d+1} ===")
    x = sp.Symbol("x")
    okH = True
    for d in (2, 3, 4, 5):
        A = [sp.Integer(1)] + [sp.Symbol(f"A{j}") for j in range(1, d + 1)] + [sp.Integer(1)]
        m = d + 1
        C = sp.zeros(m, m)
        for i in range(m):
            C[0, i] = q ** i * A[i + 1]
            if i + 1 < m:
                C[i + 1, i] = -q ** i
        got = sp.expand(sp.det(x * sp.eye(m) - C))
        want = sp.expand(sum((-1) ** j * q ** (j * (j - 1) // 2) * A[j] * x ** (m - j)
                             for j in range(m + 1)))
        good = sp.simplify(got - want) == 0
        okH &= good
        print(f"      d={d}: {sp.collect(got, x)}   {'ok' if good else 'MISMATCH'}")
    row("H", "the companion matrix read off the intertwining has characteristic polynomial "
             "sum_j (-1)^j q^{j(j-1)/2} A_j x^{d+1-j}: the cubic Hecke pencil at d = 2 and "
             "the quartic x^4 - A_1x^3 + qA_2x^2 - q^3A_3x + q^6 of PGL_4 at d = 3",
        "4 values of d", "exact", okH)
    return okA and okH


okNs = []


def check_IPU():
    print("\n=== (I) the three global identities on the thick A~_2 complexes of Paper X ===")
    okI = okP = okU = True
    for moduli, name in (((1, 1, 7), "Y_21"), ((2, 2, 2), "Y_24"), ((1, 2, 7), "Y_42")):
        X = pb.thick_cover(moduli)
        q = X.q
        L, S, T, M, n, E = pieces(X)
        A1 = np.array(X.A1, dtype=float)
        A2 = np.array(X.A2, dtype=float)
        i1 = np.allclose(L @ S.T, q * q * T.T)
        i2 = np.allclose(L @ T.T, T.T @ A1 - M.T)
        i3 = np.allclose(L @ M.T, q * (T.T @ A2 - S.T))
        j1 = np.allclose(L.T @ T.T, q * q * S.T)
        j2 = np.allclose(L.T @ S.T, S.T @ A2 - M.T)
        j3 = np.allclose(L.T @ M.T, q * (S.T @ A1 - T.T))
        okI &= i1 and i2 and i3 and j1 and j2 and j3
        print(f"      {name}: n={n:3d}, |E|={E:4d}, q={q};  "
              f"L_E S^t = q^2 T^t {i1};  L_E T^t = T^t A_1 - M^t {i2};  "
              f"L_E M^t = q(T^t A_2 - S^t) {i3}")
        print(f"            the mirror, with S and T and A_1 and A_2 exchanged: "
              f"{j1}, {j2}, {j3}")
        # (P) two-sided absorption and the companion matrix, in the order [N_0 N_1 N_2]
        P = np.hstack([T.T, M.T, S.T])
        r = np.linalg.matrix_rank(P, tol=1e-8)
        U = np.linalg.svd(P)[0][:, :r]
        esc = lambda Op: np.linalg.norm((Op @ P) - U @ (U.T @ (Op @ P))) \
            / max(np.linalg.norm(Op @ P), 1e-12)
        e1, e2 = esc(L), esc(L.T)
        Z = np.zeros((n, n))
        I = np.eye(n)
        Cn = np.block([[A1, q * A2, q * q * I],
                       [-I, Z, Z],
                       [Z, -q * I, Z]])
        C = Cn
        good = np.allclose(L @ P, P @ Cn)
        C2, C3 = Cn @ Cn, Cn @ Cn @ Cn
        A1b = np.block([[A1, Z, Z], [Z, A1, Z], [Z, Z, A1]])
        A2b = np.block([[A2, Z, Z], [Z, A2, Z], [Z, Z, A2]])
        Ib = np.eye(3 * n)
        pencil = np.allclose(C3, A1b @ C2 - q * (A2b @ Cn) + q ** 3 * Ib)
        okP &= good and e1 < 1e-8 and e2 < 1e-8 and pencil
        print(f"            rank Psi_2 = {r} of {3*n};  escape under L_E {e1:.1e}, "
              f"under L_E^t {e2:.1e};  L_E Psi = Psi C {good};  "
              f"C^3 = A_1 C^2 - q A_2 C + q^3 {pencil}")
        # (U) unitarity on W
        Wb = np.linalg.svd(P)[0][:, r:]
        assert np.allclose(S @ Wb, 0) and np.allclose(T @ Wb, 0), "W not inside ker S cap ker T"
        R = Wb.T @ L @ Wb
        stab = np.linalg.norm(L @ Wb - Wb @ R) / np.linalg.norm(L @ Wb)
        ev = np.abs(np.linalg.eigvals(R))
        unit = np.allclose(ev, np.sqrt(q))
        okU &= unit and stab < 1e-8
        print(f"            dim W = {Wb.shape[1]};  L_E W = W R with residual {stab:.1e};  "
              f"|eig R| in [{ev.min():.6f}, {ev.max():.6f}], sqrt(q) = {np.sqrt(q):.6f}  "
              f"{'ok' if unit else 'FAIL'}")
        # (N) the corank and the type-indicator relations of Theorem 6
        typ = list(X.typ)
        ind = lambda s: np.array([1.0 if typ[v] == s else 0.0 for v in range(n)])
        # b_i = #{v of relative type i forming a simplex with e} = Gaussian binomial [d,i]_q
        b = [1.0, float(q + 1), 1.0]
        assert all(np.allclose(N.sum(0), bi) for N, bi in zip((T, M, S), b)), \
            "b_i is not the number of vertices of relative type i on an edge"
        rel = []
        for t in range(3):
            for a in ((1, -1, 0), (0, 1, -1)):
                w = np.zeros(3 * n)
                for i, coef in enumerate(a):
                    if coef:
                        w[i * n:(i + 1) * n] += coef * ind((t + i) % 3) / b[i]
                rel.append(w)
        Rm = np.array(rel).T
        relres = np.linalg.norm(P @ Rm)
        ker = np.linalg.svd(P)[2][r:].T
        spanned = np.linalg.matrix_rank(np.hstack([ker, Rm]), tol=1e-8)
        okN = (relres < 1e-8) and (3 * n - r == 2 * 3) and (spanned == 3 * n - r)
        okNs.append(okN)
        print(f"            corank Psi_2 = {3*n-r} = d(d+1);  the relations "
              f"N_i^t 1_{{t+i}} = b_i N_0^t 1_t, b = {[int(x) for x in b]}, give "
              f"||Psi w|| = {relres:.1e} and span the kernel: {spanned == 3*n-r}")
    row("I", "the three rank-2 intertwining identities and their three mirrors hold on Y_21, "
             "Y_24 and Y_42; each is the flag count of check (F) at i = 0, 1, 2, the mirrors "
             "being the same count read through the duality of PG(d,q)",
        "3 complexes, 6 identities each", "exact", okI)
    row("P", "im Psi_2 absorbs L_E and L_E^t, the companion matrix assembled from the three "
             "identities satisfies L_E Psi = Psi C exactly, and C obeys the cubic Hecke "
             "pencil of Paper X", "3 complexes", "exact", okP)
    row("U", "W_2 = (im Psi_2)^perp lies in ker S cap ker T, is L_E-invariant, and every "
             "eigenvalue of L_E there has modulus q^{1/2} -- Paper XVIII's Cor. 2.7 with its "
             "hypothesis discharged", "3 complexes", "numeric", okU)
    row("N", "the relations N_i^t 1_{t+i} = b_i N_0^t 1_t force corank Psi_d >= d(d+1); at "
             "d = 2 the corank is exactly 6 = d(d+1) and these relations span the kernel, so "
             "dim W_2 = |E| - 3n + 6", "3 complexes", "exact", all(okNs) and len(okNs) == 3)
    return okI and okP and okU


def main():
    t0 = time.time()
    print("=== the invariant subspace in rank d ===")
    check_FK()
    check_HA()
    check_IPU()
    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
