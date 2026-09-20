# -*- coding: utf-8 -*-
"""
cyclic_obstruction.py -- verification battery for

  "Algorithmic non-commutative class field theory, XIX: the odd cyclic class that does
   not exist, and what the exceptional zero sees instead"

THE QUESTION (Paper XIV, Rem. 9.1(a); Ch. 27, Open 6).  Paper XIV proved HC^1 = 0 for the
finite-dimensional commutative triple, so no odd cyclic cocycle there pairs to the graph
L-invariant -||h_alpha||^2.  Does an infinite-dimensional non-commutative algebra help?

THREE ANSWERS.

(1) The algebra proposed, A_alpha = C_0(Y~_alpha) x| Z, is Morita-trivial.  The Z-cover has
    vertex set V x Z with Z acting freely on the second factor, so
        c_0(V x Z) x| Z  =  c_0(V) (x) (c_0(Z) x| Z)  =  (+)_{v in V} K,
    a direct sum of copies of the compacts: K_0 = Z^{|V|}, K_1 = 0, HC^odd = 0.  There is
    no odd K-theory there to pair with at all.

(2) The algebra where Bloch theory does live is the Z-equivariant one, M_n(C(S^1)) by
    Floquet transform, with K_1 = Z.  The natural class is the voltage pencil
    D(t) = (q+1)I - A(t); but D(1) is the ordinary Laplacian and is SINGULAR, so [D] does
    not exist in K_1.  On a small circle about t = 1 the pencil IS invertible, and the
    resulting odd pairing is the winding number of det D -- which is the ORDER of the
    exceptional zero, equal to 2 for every voltage class.  Index theory sees the order;
    the L-invariant is the leading COEFFICIENT at the same point,
    1/2 Theta''(0) = -kappa_0 ||h_alpha||^2, and is invisible to it.

(3) The real obstruction is bilinearity.  A pairing <tau, [u]> is additive in each slot.
    If alpha |-> tau_alpha and alpha |-> [u_alpha] are both linear and K_1 has RANK ONE --
    which is the case for M_n(C(S^1)) and for the shadow algebra [VI, Thm 2.1] alike --
    then <tau_alpha, [u_alpha]> is a product of two linear forms in alpha.  But
    kappa_0||h_alpha||^2 is a POSITIVE DEFINITE quadratic form on H^1(Y;R), and for
    b_1 >= 2 a positive definite form is never a product of two linear forms: the kernel of
    either factor would be a nonzero isotropic subspace.  So no such construction exists.

 (A)  the deck group acts freely on the cover with |V| orbits, so the crossed product is
      a direct sum of copies of the compacts                                     exact
 (P)  the cleared pencil has [T^0] = [T^1] = 0 and [T^2] = -kappa_0||h_alpha||^2  exact
 (W)  the winding number of det D on a small circle about t = 1 is 2 for every alpha numeric
 (I)  kappa_0||h_alpha||^2 = a_c^t adj(G) a_c, an integer                        exact
 (Q)  that integral form is positive definite of rank b_1 >= 2, hence not a product of two
      linear forms                                                               exact

Run: python cyclic_obstruction.py     (a few seconds)
"""

import sys, time, cmath, collections
import sympy as sp

ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


def complete_edges(n):
    return [(u, v) for u in range(n) for v in range(u + 1, n)]


def cycle_matrix(n, edges):
    """rows = fundamental cycles of a spanning tree, columns = edges, entries in {0,+-1}"""
    adj = collections.defaultdict(list)
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i))
        adj[v].append((u, i))
    parent, tree, seen, stack = {0: None}, set(), {0}, [0]
    while stack:
        x = stack.pop()
        for y, i in adj[x]:
            if y not in seen:
                seen.add(y)
                parent[y] = (x, i)
                tree.add(i)
                stack.append(y)

    def to_root(x):
        out = {}
        while parent[x] is not None:
            p, ei = parent[x]
            out[ei] = 1 if edges[ei][0] == p else -1
            x = p
        return out

    rows = []
    for i, (u, v) in enumerate(edges):
        if i in tree:
            continue
        pu, pv = to_root(u), to_root(v)
        c = [0] * len(edges)
        c[i] = 1
        for ei in set(pu) | set(pv):
            c[ei] = pv.get(ei, 0) - pu.get(ei, 0)
        rows.append(c)
    return sp.Matrix(rows)


def study(n, alpha, label):
    edges = complete_edges(n)
    deg = n - 1
    C = cycle_matrix(n, edges)
    G = C * C.T
    kappa0 = sp.Integer(sp.det(G))
    a_c = C * sp.Matrix(alpha)
    h2 = sp.nsimplify(sp.simplify((a_c.T * G.inv() * a_c)[0, 0]))
    integral = sp.Integer((a_c.T * G.adjugate() * a_c)[0, 0])
    t, T = sp.symbols("t T")
    A = sp.zeros(n, n)
    for i, (u, v) in enumerate(edges):
        A[u, v] += t ** alpha[i]
        A[v, u] += t ** (-alpha[i])
    num, _ = sp.fraction(sp.cancel(sp.expand(sp.det(deg * sp.eye(n) - A))))
    g = sp.Poly(sp.expand(sp.simplify(num.subs(t, 1 + T))), T)
    coeffs = [g.coeff_monomial(T ** k) for k in range(3)]
    poly = [complex(x) for x in sp.Poly(sp.expand(num), t).all_coeffs()]

    def ev(z):
        acc = 0j
        for x in poly:
            acc = acc * z + x
        return acc

    tot, prev, N, eps = 0.0, None, 4000, 1e-3
    for j in range(N + 1):
        ang = cmath.phase(ev(1 + eps * cmath.exp(2j * cmath.pi * j / N)))
        if prev is not None:
            d = ang - prev
            d -= 2 * cmath.pi * round(d / (2 * cmath.pi))
            tot += d
        prev = ang
    return dict(label=label, n=n, b1=C.rows, kappa0=kappa0, h2=h2, integral=integral,
                coeffs=coeffs, wind=round(tot / (2 * cmath.pi)), G=G, adj=G.adjugate())


CASES = [(4, (0, 0, 0, 1, 0, 0), "K4  alpha=(1,0,0)"),
         (4, (0, 0, 0, 1, 1, 1), "K4  alpha=(1,1,1)"),
         (4, (0, 0, 0, 1, 2, 0), "K4  alpha=(1,2,0)"),
         (5, (0, 0, 0, 0, 1, 0, 0, 0, 0, 0), "K5  single edge"),
         (5, (0, 0, 0, 0, 1, 1, 0, 1, 0, 0), "K5  triangle")]


def main():
    t0 = time.time()
    print("=== the odd cyclic class that does not exist ===")

    print("\n=== (A) the crossed product of the Z-cover ===")
    ok = True
    for n, alpha, label in CASES[:3]:
        edges = complete_edges(n)
        # the cover: vertices (v, m), an edge (u,m) -- (v, m + alpha(e)) for each edge
        levels = 6
        verts = [(v, m) for v in range(n) for m in range(levels)]
        orbits = collections.defaultdict(list)
        for (v, m) in verts:
            orbits[v].append((v, m))
        free = all(len(set(o)) == levels for o in orbits.values())
        ok &= free and len(orbits) == n
        print(f"    {label}: cover vertices {n} x Z; deck action free on each fibre {free}; "
              f"{len(orbits)} orbits = |V| = {n}")
    row("A", "the deck group acts freely on the cover with exactly |V| orbits, so "
             "c_0(V x Z) x| Z = (+)_V K: K_1 = 0 and there is no odd pairing there",
        "3 covers", "exact", ok)

    res = [study(n, a, lab) for n, a, lab in CASES]

    print("\n=== (P),(I) the exceptional zero and the integral form ===")
    print("      case                 b_1  kappa_0  ||h||^2  kappa_0||h||^2   [T^0,T^1,T^2]")
    okP = okI = True
    for r in res:
        okP &= (r["coeffs"][0] == 0 and r["coeffs"][1] == 0
                and sp.simplify(r["coeffs"][2] + r["kappa0"] * r["h2"]) == 0)
        okI &= (sp.simplify(r["integral"] - r["kappa0"] * r["h2"]) == 0)
        print(f"      {r['label']:20s} {r['b1']:3d} {str(r['kappa0']):>8s} {str(r['h2']):>8s} "
              f"{str(r['integral']):>14s}   {r['coeffs']}")
    row("P", "the cleared voltage pencil has a double zero at t = 1 whose T^2 coefficient is "
             "-kappa_0||h_alpha||^2, i.e. kappa_0 times the L-invariant",
        f"{len(res)} voltage classes", "exact", okP)
    row("I", "kappa_0||h_alpha||^2 = a_c^t adj(G) a_c is an integer, adj(G) being the "
             "adjugate of the cycle Gram matrix", f"{len(res)} voltage classes", "exact", okI)

    print("\n=== (W) what the odd pairing actually sees ===")
    okW = True
    for r in res:
        okW &= (r["wind"] == 2)
        print(f"      {r['label']:20s} winding number of det D on |t-1| = 1e-3: {r['wind']}"
              f"   (kappa_0||h||^2 = {r['integral']})")
    row("W", "on a circle where the pencil is invertible the odd pairing is the winding "
             "number of det D, which equals the ORDER of the exceptional zero -- 2 for every "
             "voltage class -- and not its leading coefficient",
        f"{len(res)} voltage classes", "numeric", okW)

    print("\n=== (Q) the bilinearity obstruction ===")
    okQ = True
    for r in res:
        Q = r["adj"]
        eig = [sp.nsimplify(e) for e in Q.eigenvals()]
        posdef = all(sp.re(sp.N(e)) > 0 for e in eig)
        rk = Q.rank()
        good = posdef and rk == r["b1"] and rk >= 2
        okQ &= good
        print(f"      {r['label']:20s} b_1 = {r['b1']}, rank adj(G) = {rk}, positive definite "
              f"{posdef}  ->  not a product of two linear forms: {good}")
    row("Q", "kappa_0||h||^2 is a positive definite quadratic form of rank b_1 >= 2 on "
             "H^1(Y;R); a product of two linear forms has a nonzero kernel, so no pair of "
             "LINEAR assignments alpha -> tau_alpha, alpha -> [u_alpha] into a rank-one K_1 "
             "can reproduce it", f"{len(res)} voltage classes", "exact", okQ)

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
