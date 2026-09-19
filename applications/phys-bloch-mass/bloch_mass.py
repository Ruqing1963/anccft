# -*- coding: utf-8 -*-
"""
bloch_mass.py -- verification battery for

    "Effective mass from spanning trees, and the absence of Bloch phases
     on rank-two hyperbolic lattices"

PART I  (rank one: Z^d-periodic lattices built as abelian covers of a finite graph)

 (M)  effective-mass tensor = Albanese metric.  With Delta_k the Peierls-twisted
      Laplacian of the Z^d cover and H_ab = alpha_a^T Pi alpha_b the Gram matrix of the
      harmonic projections,
          Hess_k det Delta_k |_{k=0} = 2 kappa_0 H          [exact, Laurent interpolation]
          lambda_0(k) = (1/n) k^T H k + O(|k|^4)            [numeric, Richardson]
 (S)  single-link flux: ||h_e||^2 = 1 - R_eff(e) and kappa_0 (1 - R_eff(e)) = tau(Y - e),
      the number of spanning trees avoiding the link.       [exact]
 (Q)  quantization: kappa_0 * H is an INTEGER matrix, so every entry of the inverse
      effective-mass tensor is an integer multiple of 2t/(n kappa_0).  [exact]
 (E)  Einstein relation: sigma^2 = H/m by an independent Green-Kubo edge sum. [exact]
 (G)  parabolicity window: the gap lambda_1(0) and the size of the quartic term. [numeric]

PART II (rank two: thick A~_2 complexes -- the 2-dimensional hyperbolic lattices)

 (B)  b_1(Y) = 0 by exact ranks, and H_1(Y;Z) finite: no Z^d cover exists, hence no
      Brillouin zone; the only threadable fluxes form the finite group H_1(Y;Z). [exact]

Exact arithmetic is sympy Rational/Integer throughout; every floating-point row is
labelled "numeric".  Part II imports the complex builder of Paper X.
"""

import sys, os, time, itertools
import sympy as sp
import numpy as np
import mpmath as mp

mp.mp.dps = 40

PAPER10 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Paper 10")

EXACT, NUM = "exact", "numeric"
ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# ----------------------------------------------------------------------------- graphs

def graph(name, n, edges):
    """edges: list of (u,v), u != v, each undirected edge once.  Returns a dict."""
    deg = [0] * n
    for (u, v) in edges:
        deg[u] += 1
        deg[v] += 1
    assert len(set(deg)) == 1, f"{name} is not regular: {deg}"
    return dict(name=name, n=n, edges=list(edges), m=len(edges), q=deg[0] - 1)


def K(n):
    return graph(f"K_{n}", n, [(i, j) for i in range(n) for j in range(i + 1, n)])


def K33():
    return graph("K_{3,3}", 6, [(i, 3 + j) for i in range(3) for j in range(3)])


def cube():
    vs = [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    idx = {v: k for k, v in enumerate(vs)}
    es = set()
    for v in vs:
        for i in range(3):
            w = list(v)
            w[i] ^= 1
            es.add(tuple(sorted((idx[v], idx[tuple(w)]))))
    return graph("Q_3", 8, sorted(es))


def petersen():
    out = [(i, (i + 1) % 5) for i in range(5)]
    inn = [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    spk = [(i, 5 + i) for i in range(5)]
    return graph("Petersen", 10, out + inn + spk)


def heawood():
    """Incidence graph of PG(2,2): 7 points + 7 lines, 21 edges, 3-regular, girth 6.
    This is exactly the vertex link of a thick A~_2 building of order 2 (Part II)."""
    LINES = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2)]
    es = [(p, 7 + li) for li, L in enumerate(LINES) for p in L]
    return graph("Heawood", 14, es)


# ------------------------------------------------------------------- exact linear algebra

def incidence(G):
    """B : n x m over Z, B[t(e)][e] = +1, B[o(e)][e] = -1 for the stored orientation."""
    B = sp.zeros(G["n"], G["m"])
    for k, (u, v) in enumerate(G["edges"]):
        B[v, k] += 1
        B[u, k] -= 1
    return B


def laplacian_pinv(L, n):
    """Moore-Penrose inverse of a connected graph Laplacian: (L + J/n)^{-1} - J/n."""
    J = sp.ones(n, n) / n
    return (L + J).inv() - J


def spanning_trees(G):
    """Kirchhoff: determinant of the reduced Laplacian (delete vertex 0)."""
    n = G["n"]
    L = sp.zeros(n, n)
    for (u, v) in G["edges"]:
        L[u, u] += 1
        L[v, v] += 1
        L[u, v] -= 1
        L[v, u] -= 1
    return sp.Integer(sp.Matrix(L[1:, 1:]).det())


def delete_edge(G, k):
    es = [e for i, e in enumerate(G["edges"]) if i != k]
    Gd = dict(G)
    Gd["edges"] = es
    Gd["m"] = len(es)
    return Gd


def spanning_tree_and_chords(G):
    """A spanning tree by BFS; returns (tree_edge_indices, chord_indices)."""
    adj = {i: [] for i in range(G["n"])}
    for k, (u, v) in enumerate(G["edges"]):
        adj[u].append((v, k))
        adj[v].append((u, k))
    seen, tree, order = {0}, set(), [0]
    while order:
        u = order.pop(0)
        for (v, k) in adj[u]:
            if v not in seen:
                seen.add(v)
                tree.add(k)
                order.append(v)
    chords = [k for k in range(G["m"]) if k not in tree]
    return sorted(tree), chords


# ------------------------------------------------- the Laurent determinant and its Hessian

def laurent_coeffs(G, alphas):
    """Exact Laurent coefficients of  Theta(t) = det Delta(t),
       Delta(t)[u][v] = -sum_{e: u->v} t^{alpha(e)},  diagonal q+1,
    where alphas is a list of d integer vectors (one per Z-direction), each of length m.

    Support: exponent j_a ranges in [-M_a, M_a], M_a = sum_e |alpha_a(e)|.  We evaluate the
    determinant exactly at a product grid of rational points and solve the (small)
    Vandermonde system.  Returns {multi-index j : coefficient}.
    """
    n, m, q = G["n"], G["m"], G["q"]
    d = len(alphas)
    Ms = [sum(abs(int(a[k])) for k in range(m)) for a in alphas]
    axes = [list(range(1, 2 * M + 2)) for M in Ms]          # distinct positive integers
    idxs = [list(range(-M, M + 1)) for M in Ms]
    grid = list(itertools.product(*axes))
    basis = list(itertools.product(*idxs))

    def theta_at(tvals):
        D = sp.zeros(n, n)
        for i in range(n):
            D[i, i] = q + 1
        for k, (u, v) in enumerate(G["edges"]):
            e = sp.Integer(1)
            for a in range(d):
                e *= sp.Rational(tvals[a]) ** int(alphas[a][k])
            D[u, v] -= e
            D[v, u] -= 1 / e
        return sp.together(D.det())

    A = sp.zeros(len(grid), len(basis))
    b = sp.zeros(len(grid), 1)
    for r, gpt in enumerate(grid):
        for c, j in enumerate(basis):
            val = sp.Integer(1)
            for a in range(d):
                val *= sp.Rational(gpt[a]) ** int(j[a])
            A[r, c] = val
        b[r] = theta_at(gpt)
    sol = A.solve(b)
    return {basis[c]: sp.nsimplify(sol[c]) for c in range(len(basis))}


def hessian_from_coeffs(coeffs, d):
    """Hess_k Theta(e^{ik})|_{k=0} = - sum_j c_j j_a j_b   (all first derivatives vanish)."""
    Hs = sp.zeros(d, d)
    for j, c in coeffs.items():
        for a in range(d):
            for bb in range(d):
                Hs[a, bb] -= c * sp.Integer(j[a]) * sp.Integer(j[bb])
    return sp.simplify(Hs)


# ------------------------------------------------------------------- numeric Bloch bands

def lambda0_numeric(G, alphas, kvec):
    """Float64 bands; used only for the spectral gap lambda_1(0)."""
    n, q, d = G["n"], G["q"], len(alphas)
    D = np.zeros((n, n), dtype=complex)
    for i in range(n):
        D[i, i] = q + 1
    for k, (u, v) in enumerate(G["edges"]):
        ph = np.exp(1j * sum(float(kvec[a]) * int(alphas[a][k]) for a in range(d)))
        D[u, v] -= ph
        D[v, u] -= np.conj(ph)
    w = np.linalg.eigvalsh(D)
    return float(w[0]), float(w[1])


def lambda0_mp(G, alphas, kvec):
    """Smallest Bloch band at high precision (mpmath, dps set by the caller)."""
    n, q, d = G["n"], G["q"], len(alphas)
    D = mp.zeros(n, n)
    for i in range(n):
        D[i, i] = mp.mpf(q + 1)
    for k, (u, v) in enumerate(G["edges"]):
        arg = sum(kvec[a] * int(alphas[a][k]) for a in range(d))
        ph = mp.exp(mp.mpc(0, 1) * arg)
        D[u, v] -= ph
        D[v, u] -= mp.conj(ph)
    E = mp.eighe(D, eigvals_only=True)
    return min(E, key=lambda x: mp.mpf(x))


def richardson_curvature(G, alphas, u, h=None):
    """lambda_0(h u)/h^2 -> u^T H u / n.  High precision + one Richardson step:
    truncation O(h^4) relative, roundoff far below it, so the residual is ~1e-24."""
    if h is None:
        h = mp.mpf(10) ** (-6)
    r1 = lambda0_mp(G, alphas, [h * x for x in u]) / h ** 2
    r2 = lambda0_mp(G, alphas, [h / 2 * x for x in u]) / (h / 2) ** 2
    return (4 * r2 - r1) / 3


# --------------------------------------------------------------------------- Part I driver

def analyse(G, alphas, label):
    n, m, q, d = G["n"], G["m"], G["q"], len(alphas)
    B = incidence(G)
    L = B * B.T
    Lp = laplacian_pinv(L, n)
    Pi = sp.eye(m) - B.T * Lp * B
    kappa0 = spanning_trees(G)

    Amat = sp.Matrix([[sp.Integer(int(alphas[a][k])) for a in range(d)] for k in range(m)])
    H = sp.simplify(Amat.T * Pi * Amat)

    print(f"\n--- {G['name']}  (n={n}, m={m}, q={q}, kappa_0={kappa0})   voltage: {label}")
    print(f"    H = Albanese/Gram matrix of harmonic projections = {sp.nsimplify(H).tolist()}")
    print(f"    kappa_0 * H = {(kappa0 * H).tolist()}")

    # (Q) quantization
    kH = sp.simplify(kappa0 * H)
    ok_q = all(sp.Integer(kH[i, j]) == kH[i, j] for i in range(d) for j in range(d))
    row("Q", f"kappa_0*H integral (det = {sp.simplify((kappa0*H).det())})",
        f"{G['name']}/{label}", EXACT, ok_q)

    # (M-exact) Hessian of the Floquet determinant
    coeffs = laurent_coeffs(G, alphas)
    Hs = hessian_from_coeffs(coeffs, d)
    ok_m = sp.simplify(Hs - 2 * kappa0 * H) == sp.zeros(d, d)
    row("M", f"Hess_k det Delta_k = 2*kappa_0*H  (= {sp.nsimplify(Hs).tolist()})",
        f"{G['name']}/{label}", EXACT, ok_m)

    # (M-num) the band itself
    errs = []
    dirs = [[1 if b == a else 0 for b in range(d)] for a in range(d)]
    if d > 1:
        dirs.append([1] * d)
    for u in dirs:
        got = richardson_curvature(G, alphas, [mp.mpf(x) for x in u])
        want = mp.mpf(sp.Rational(sum(H[a, b] * u[a] * u[b]
                                      for a in range(d) for b in range(d)), n).p) / \
               mp.mpf(sp.Rational(sum(H[a, b] * u[a] * u[b]
                                      for a in range(d) for b in range(d)), n).q)
        errs.append(abs(got - want) / abs(want))
    worst = max(errs)
    ok_b = worst < mp.mpf(10) ** (-15)
    row("M", f"lambda_0(k) = k^T H k / n + O(|k|^4), max rel err {mp.nstr(worst, 3)}",
        f"{G['name']}/{label}", NUM, ok_b)

    # (E) Einstein relation via an independent Green-Kubo edge sum
    Hc = sp.Matrix(Pi * Amat)                      # harmonic representatives, m x d
    sig = sp.zeros(d, d)
    for a in range(d):
        for b in range(d):
            s = sp.Integer(0)
            for k in range(m):                     # both orientations of every edge
                s += 2 * Hc[k, a] * Hc[k, b]
            sig[a, b] = sp.nsimplify(s / (n * (q + 1)))
    ok_e = sp.simplify(sig - H / m) == sp.zeros(d, d)
    row("E", f"sigma^2 = H/m (Green-Kubo) = {sp.nsimplify(sig).tolist()}",
        f"{G['name']}/{label}", EXACT, ok_e)

    # (G) parabolicity window
    l0, l1 = lambda0_numeric(G, alphas, [0.0] * d)
    row("G", f"spectral gap lambda_1(0) = {l1:.6f} > 0, lambda_0(0) = {l0:.2e}",
        f"{G['name']}/{label}", NUM, l1 > 1e-9 and abs(l0) < 1e-9)

    return H, kappa0


def single_link(G):
    """(S) single-link flux: ||h_e||^2 = 1 - R_eff(e), kappa_0(1-R_eff) = tau(Y - e)."""
    n, m = G["n"], G["m"]
    B = incidence(G)
    L = B * B.T
    Lp = laplacian_pinv(L, n)
    Pi = sp.eye(m) - B.T * Lp * B
    kappa0 = spanning_trees(G)
    _, chords = spanning_tree_and_chords(G)
    k0 = chords[0]
    Reff = sp.nsimplify(1 - Pi[k0, k0])
    tau_avoid = spanning_trees(delete_edge(G, k0))
    ok = (sp.nsimplify(kappa0 * (1 - Reff)) == tau_avoid) and (sp.nsimplify(Pi[k0, k0]) == 1 - Reff)
    mstar = sp.nsimplify(sp.Rational(n, 2) / (1 - Reff))          # in units of 1/t
    row("S", f"R_eff={Reff}, ||h||^2={1-Reff}, kappa_0(1-R_eff)={sp.nsimplify(kappa0*(1-Reff))}"
             f" = tau(Y-e)={tau_avoid}, m* = {mstar}/t",
        G["name"], EXACT, ok)
    return k0, Reff, tau_avoid, kappa0, mstar


# -------------------------------------------------------------------------- Part II driver

def rank_two():
    print("\n=== PART II: rank two (thick A~_2 complexes) ===")
    sys.path.insert(0, PAPER10)
    import pgl3_building as pb

    ok_pres = pb.check_presentation()
    row("B", "triangle presentation over PG(2,2) with Singer cycle re-verified",
        "X_3", EXACT, ok_pres)

    results = []
    for moduli, want_n in [((1, 1, 7), 21), ((2, 2, 2), 24), ((1, 1, 14), 42)]:
        t0 = time.time()
        X = pb.thick_cover(moduli)
        E, F = len(X.edges), len(X.tris)
        r1, r2 = pb.exact_rank(X.d1), pb.exact_rank(X.d2)
        b1 = E - r1 - r2
        # torsion of H_1 = ker d1 / im d2, in coordinates on ker d1
        dd, _, Qi, _ = pb.diagonalize(X.d1, False, True)
        rank = sum(1 for x in dd if x != 0)
        Lc = pb.mul([Qi[i] for i in range(rank, E)], X.d2)
        inv, free, _ = pb.coker_structure(Lc, primes=pb.SMALL_PRIMES)
        grp = pb.group_str(inv, free)
        dt = time.time() - t0
        print(f"    {X.name}: n={X.n}, |E|={E}, |F|={F}, b_1={b1}, H_1(Y;Z)={grp}  ({dt:.1f}s)")
        results.append((X.name, X.n, b1, grp, free))
        row("B", f"b_1 = {b1} and H_1(Y;Z) = {grp} finite -> no Z-cover, no Brillouin zone",
            X.name, EXACT, b1 == 0 and free == 0)
    return results


# --------------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    print("=== PART I: Z^d-periodic lattices from abelian covers of a finite graph ===")

    graphs = [K(4), K33(), cube(), K(5), petersen(), heawood()]

    print("\n--- (S) single-link flux: effective mass from spanning trees ---")
    sl = {}
    for G in graphs:
        sl[G["name"]] = single_link(G)

    # d = 1, single chord, on every graph
    for G in graphs:
        k0 = sl[G["name"]][0]
        a = [0] * G["m"]
        a[k0] = 1
        analyse(G, [a], f"single link e_{k0}")

    # cross-validation against [ANCFT IX, Table 1] and [XIV, Table 1]
    print("\n--- cross-validation against the rank-one theory [IX, XIV] ---")
    G = K(4)
    _, chords = spanning_tree_and_chords(G)
    for lab in [(1, 0, 0), (1, 1, 1), (1, 2, 0)]:
        a = [0] * G["m"]
        for i, c in enumerate(chords):
            a[c] = lab[i]
        H, _ = analyse(G, [a], f"chords {lab}")
        want = {(1, 0, 0): sp.Rational(1, 2), (1, 1, 1): sp.Integer(1),
                (1, 2, 0): sp.Rational(3, 2)}[lab]
        row("X", f"||h||^2 = {sp.nsimplify(H[0,0])} matches [IX, Tab.1] value {want}",
            f"K_4/{lab}", EXACT, sp.simplify(H[0, 0] - want) == 0)

    G = K(5)
    tri = [k for k, (u, v) in enumerate(G["edges"]) if {u, v} <= {0, 1, 2}]
    # (a) the [XIV] convention: +1 on each triangle edge in the stored orientation
    a = [0] * G["m"]
    for k in tri:
        a[k] = 1
    H, _ = analyse(G, [a], "triangle 0-1-2, stored orientation")
    row("X", f"||h||^2 = {sp.nsimplify(H[0,0])} matches [XIV, Tab.1] value 7/5",
        "K_5/triangle", EXACT, sp.simplify(H[0, 0] - sp.Rational(7, 5)) == 0)
    # (b) the cyclically oriented triangle: alpha lies in the cycle space, so h = alpha
    a = [0] * G["m"]
    sgn = {(0, 1): 1, (1, 2): 1, (0, 2): -1}
    for k in tri:
        a[k] = sgn[G["edges"][k]]
    H, _ = analyse(G, [a], "triangle 0-1-2, cyclically oriented")
    row("C", f"flux pattern a cycle of length 3 => ||h||^2 = 3 exactly, m* = n/(2*3)",
        "K_5/3-cycle", EXACT, sp.simplify(H[0, 0] - 3) == 0)

    # d = 2: the effective-mass TENSOR (new beyond the rank-one theory)
    print("\n--- (M/Q) d = 2: the effective-mass tensor and its quantized determinant ---")
    for G in [K(4), cube(), K(5), heawood()]:
        _, chords = spanning_tree_and_chords(G)
        a1 = [0] * G["m"]
        a1[chords[0]] = 1
        a2 = [0] * G["m"]
        a2[chords[1]] = 1
        analyse(G, [a1, a2], f"two links e_{chords[0]}, e_{chords[1]}")

    rank_two()

    print("\n=== battery summary ===")
    bad = [r for r in ROWS if not r[4]]
    for key in ("S", "M", "Q", "E", "G", "X", "C", "B"):
        sub = [r for r in ROWS if r[0] == key]
        print(f"  ({key})  {sum(1 for r in sub if r[4])}/{len(sub)} passed")
    print(f"  TOTAL {sum(1 for r in ROWS if r[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(r[0], r[2]) for r in bad]}"))
    print(f"  elapsed {time.time()-t0:.1f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
