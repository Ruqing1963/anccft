# -*- coding: utf-8 -*-
"""
boundary_blind.py -- verification battery for

    "The boundary algebra of p-adic AdS/CFT does not see the bulk arithmetic"

Setting.  X is the (q+1)-regular Bruhat-Tits tree of GL_2 over a non-archimedean local
field with residue field of size q; Gamma is a torsion-free cocompact lattice, free of
rank r; Y = Gamma\\X is a finite (q+1)-regular graph on n vertices, m = (q+1)n/2 edges,
r = m - n + 1.  Two algebras are attached to the same data:

  BOUNDARY   C(dX) x Gamma  =  O_{A^nb},  A^nb the non-backtracking matrix on 2m oriented
             edges.  K_0 = coker(I - (A^nb)^T), K_1 = ker(I - (A^nb)^T).
  BULK       the Hecke-Laplacian channel, Delta = (q+1)I - A on n vertices.
             K_0 = coker(Delta) = Z + Jac(Y), |Jac(Y)| = (1/n) prod_f (p+1-a_p(f)).

 (K)  K_0(boundary) = Z^r + Z/(r-1) and K_1 = Z^r, on every test graph.       exact
 (B)  r - 1 = (q-1)n/2 depends only on (q,n): the boundary K-theory is CONSTANT
      on the set of (q+1)-regular graphs with n vertices.                     exact
 (J)  the bulk critical group Jac(Y) is NOT constant on that set: explicit witness
      pairs with identical boundary K-theory and different bulk arithmetic.   exact
 (S)  the spread: sampling cubic graphs on n vertices, how many distinct |Jac| occur
      inside a single boundary K-theory class.                                exact / sampled
 (R)  the Radon-Nikodym cocycle of the harmonic measure on cylinder sets takes values
      in q^Z -- the finite input to Krieger's type III_{1/q} classification.  exact

Exact integer arithmetic throughout. Requires sympy. Run: python boundary_blind.py
"""

import sys, os, time, random, itertools, collections
from fractions import Fraction
import sympy as sp

PAPER10 = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "..", "papers", "10")
sys.path.insert(0, PAPER10)
import pgl3_building as pb

EXACT, SAMPLED = "exact", "sampled"
ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# --------------------------------------------------------------------------- graphs

def graph(name, n, edges):
    deg = [0] * n
    for (u, v) in edges:
        deg[u] += 1
        deg[v] += 1
    assert len(set(deg)) == 1, f"{name} not regular: {deg}"
    return dict(name=name, n=n, edges=sorted(edges), m=len(edges), q=deg[0] - 1)


def bipartite(G):
    adj = collections.defaultdict(list)
    for (u, v) in G["edges"]:
        adj[u].append(v)
        adj[v].append(u)
    col, stack = {0: 0}, [0]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in col:
                col[y] = 1 - col[x]
                stack.append(y)
            elif col[y] == col[x]:
                return False
    return True


def robertson_type(G):
    """Robertson, Houston J. Math. 31 (2005), Thm 2 / Cor 1(2): L^inf(dX) x Gamma is the
    hyperfinite factor of type III_lambda with lambda = 1/q^2 if the quotient graph is
    bipartite (equivalently Gamma is contained in PSL_2(F)), and 1/q otherwise."""
    return "III_(1/q^2)" if bipartite(G) else "III_(1/q)"


def K(n):
    return graph(f"K_{n}", n, [(i, j) for i in range(n) for j in range(i + 1, n)])


def K33():
    return graph("K_{3,3}", 6, [(i, 3 + j) for i in range(3) for j in range(3)])


def prism(k):
    """C_k x K_2, cubic, n = 2k."""
    e = [(i, (i + 1) % k) for i in range(k)]
    e += [(k + i, k + (i + 1) % k) for i in range(k)]
    e += [(i, k + i) for i in range(k)]
    return graph(f"prism C_{k}xK_2", 2 * k, e)


def mobius(k):
    """Moebius ladder M_{2k}: cycle C_{2k} plus the k diameters. Cubic, n = 2k."""
    n = 2 * k
    e = [(i, (i + 1) % n) for i in range(n)]
    e += [(i, i + k) for i in range(k)]
    return graph(f"Moebius M_{n}", n, e)


def petersen():
    out = [(i, (i + 1) % 5) for i in range(5)]
    inn = [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    spk = [(i, 5 + i) for i in range(5)]
    return graph("Petersen", 10, out + inn + spk)


def heawood():
    LINES = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2)]
    return graph("Heawood", 14, [(p, 7 + li) for li, L in enumerate(LINES) for p in L])


def random_cubic(n, rng, tries=400):
    """Configuration model with rejection; returns a simple connected cubic graph."""
    for _ in range(tries):
        stubs = [v for v in range(n) for _ in range(3)]
        rng.shuffle(stubs)
        es = set()
        ok = True
        for i in range(0, len(stubs), 2):
            u, v = stubs[i], stubs[i + 1]
            if u == v or (min(u, v), max(u, v)) in es:
                ok = False
                break
            es.add((min(u, v), max(u, v)))
        if not ok or len(es) != 3 * n // 2:
            continue
        # connectivity
        adj = collections.defaultdict(list)
        for (u, v) in es:
            adj[u].append(v)
            adj[v].append(u)
        seen, stack = {0}, [0]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        if len(seen) != n:
            continue
        return graph(f"random cubic n={n}", n, sorted(es))
    return None


# --------------------------------------------------------------- exact linear algebra

def snf_invariants(M):
    """Invariant factors and free rank of coker(M) for an integer matrix M.
    Uses the exact integer Smith normal form of Paper X of the ANCFT series."""
    if not M or not M[0]:
        return [], 0
    d = pb.diagonalize([list(map(int, r)) for r in M], False, False)[0]
    ncols = len(M[0])
    nz = sum(1 for x in d if x != 0)
    free = ncols - nz
    tor = [abs(int(x)) for x in d if x != 0 and abs(x) != 1]
    return tor, free


def bulk_jacobian(G):
    """coker((q+1)I - A) = Z + Jac(Y); returns (invariant factors of Jac, order)."""
    n, q = G["n"], G["q"]
    L = [[0] * n for _ in range(n)]
    for i in range(n):
        L[i][i] = q + 1
    for (u, v) in G["edges"]:
        L[u][v] -= 1
        L[v][u] -= 1
    tor, free = snf_invariants(L)
    order = 1
    for x in tor:
        order *= x
    return tor, free, order


def oriented_edges(G):
    oe = []
    for (u, v) in G["edges"]:
        oe.append((u, v))
        oe.append((v, u))
    return oe


def nonbacktracking(G):
    oe = oriented_edges(G)
    N = len(oe)
    A = [[0] * N for _ in range(N)]
    for i, (a, b) in enumerate(oe):
        for j, (c, d) in enumerate(oe):
            if b == c and not (c == b and d == a):
                A[i][j] = 1
    return A, N


def boundary_ktheory(G):
    """K_0 = coker(I - (A^nb)^T), K_1 = ker(I - (A^nb)^T)."""
    A, N = nonbacktracking(G)
    M = [[(1 if i == j else 0) - A[j][i] for j in range(N)] for i in range(N)]
    tor, free = snf_invariants(M)
    order = 1
    for x in tor:
        order *= x
    return tor, free, order


# ------------------------------------------------------------------------------ (K)(B)

NAMED = [K(4), K33(), prism(3), prism(4), mobius(4), mobius(5), prism(5),
         petersen(), heawood(), K(5), K(6)]

# Two cubic graphs on 8 vertices found by search and recorded here, so that the coprimality
# claim of check (P) needs no search to re-verify.  Both have boundary K_0 = Z^5 + Z/4.
W_A = graph("witness A (n=8)", 8,
            [(0, 3), (0, 6), (0, 7), (1, 2), (1, 5), (1, 6),
             (2, 4), (2, 5), (3, 4), (3, 7), (4, 5), (6, 7)])
W_B = graph("witness B (n=8)", 8,
            [(0, 2), (0, 3), (0, 7), (1, 4), (1, 6), (1, 7),
             (2, 4), (2, 5), (3, 4), (3, 5), (5, 6), (6, 7)])


def check_boundary_formula():
    print("\n=== (K) boundary K-theory = Z^r + Z/(r-1), K_1 = Z^r ===")
    print("    graph             q   n   m    r   K_0(boundary)         K_1    vN type"
          "        |Jac(Y)|")
    ok_all = True
    for G in NAMED:
        n, m, q = G["n"], G["m"], G["q"]
        r = m - n + 1
        tor, free, order = boundary_ktheory(G)
        want_tor = [r - 1] if r - 1 > 1 else []
        ok = (free == r and tor == want_tor)
        ok_all &= ok
        k0 = f"Z^{free}" + ("" if not tor else " + " + " + ".join(f"Z/{t}" for t in tor))
        _, _, jorder = bulk_jacobian(G)
        print(f"    {G['name']:<17s} {q:3d} {n:3d} {m:3d} {r:4d}   {k0:<20s}  Z^{free:<4d} "
              f"{robertson_type(G):<13s} {jorder:8d}" + ("" if ok else "   <-- MISMATCH"))
    row("K", "K_0(C(dX) x Gamma) = Z^r + Z/(r-1) and K_1 = Z^r on every test graph",
        f"{len(NAMED)} graphs", EXACT, ok_all)

    # (B) r-1 depends only on (q,n)
    ok_b = all((G["m"] - G["n"] + 1) - 1 == (G["q"] - 1) * G["n"] // 2 for G in NAMED)
    row("B", "r - 1 = (q-1)n/2, so the boundary K-theory is a function of (q,n) alone",
        f"{len(NAMED)} graphs", EXACT, ok_b)
    return ok_all and ok_b


# ---------------------------------------------------------------------------- (J)

def check_witnesses():
    print("\n=== (J) witness pairs: identical boundary, different bulk ===")
    groups = collections.defaultdict(list)
    for G in NAMED:
        groups[(G["q"], G["n"])].append(G)
    ok_any = False
    for (q, n), gs in sorted(groups.items()):
        if len(gs) < 2:
            continue
        r = gs[0]["m"] - gs[0]["n"] + 1
        btor, bfree, _ = boundary_ktheory(gs[0])
        same_boundary = all(boundary_ktheory(G)[:2] == (btor, bfree) for G in gs)
        print(f"\n    q={q}, n={n}:  boundary K_0 = Z^{bfree} + Z/{r-1} for ALL of them")
        jacs = []
        for G in gs:
            tor, free, order = bulk_jacobian(G)
            j = " + ".join(f"Z/{t}" for t in tor) if tor else "0"
            jacs.append(order)
            print(f"      {G['name']:<17s}  Jac(Y) = {j:<28s} |Jac| = {order}")
        differ = len(set(jacs)) > 1
        ok_any |= (same_boundary and differ)
        if same_boundary and differ:
            print(f"      -> same boundary K-theory, |Jac| ranges over {sorted(set(jacs))}"
                  f"  (ratio {max(jacs)/min(jacs):.3f})")
    row("J", "explicit graphs with identical boundary K-theory and different Jac(Y)",
        "q=2 at n=6,8,10", EXACT, ok_any)
    return ok_any


# ---------------------------------------------------------------------------- (P)

def prime_support(k):
    ps, d = set(), 2
    while d * d <= k:
        while k % d == 0:
            ps.add(d)
            k //= d
        d += 1
    if k > 1:
        ps.add(k)
    return ps


def check_prime_supports():
    """The strongest form: inside one boundary class the bulk critical groups need not
    even share a prime. |Jac(Y)| = (1/n) prod_f (p+1-a_p(f)) is a product of Frobenius
    traces, so this says the boundary K-theory does not determine which primes divide
    that product."""
    print("\n=== (P) prime supports of Jac(Y) inside one boundary class ===")
    groups = collections.defaultdict(list)
    for G in NAMED:
        groups[(G["q"], G["n"])].append(G)
    found_disjoint = False
    for (q, n), gs in sorted(groups.items()):
        if len(gs) < 2:
            continue
        sups = []
        for G in gs:
            _, _, order = bulk_jacobian(G)
            sups.append((G["name"], order, prime_support(order)))
        print(f"    q={q}, n={n}  (one boundary class):")
        for name, order, sup in sups:
            print(f"      {name:<17s} |Jac| = {order:<6d} primes {sorted(sup)}")
        for (n1, o1, s1), (n2, o2, s2) in itertools.combinations(sups, 2):
            if not (s1 & s2):
                print(f"      -> {n1} and {n2}: prime supports {sorted(s1)} and "
                      f"{sorted(s2)} are DISJOINT")
                found_disjoint = True
    if not found_disjoint:
        print("    among the named graphs no two have disjoint prime supports.")

    # the recorded witnesses: same K-theory, same von Neumann type, coprime bulk orders.
    # This is strictly stronger than Robertson's own counterexample (theta vs dumbbell),
    # which the von Neumann type DOES separate (III_(1/4) against III_(1/2)).
    import math
    print("\n    recorded witnesses, both cubic on n = 8:")
    data = []
    for G in (W_A, W_B):
        btor, bfree, _ = boundary_ktheory(G)
        jtor, _, order = bulk_jacobian(G)
        data.append((G["name"], bfree, btor, robertson_type(G), order))
        j = " + ".join(f"Z/{t}" for t in jtor) if jtor else "0"
        fac = '*'.join(f"{p}^{e}" for p, e in sorted(sp.factorint(order).items()))
        print(f"      {G['name']:<18s} K_0 = Z^{bfree} + Z/{btor[0]},  type "
              f"{robertson_type(G):<12s} Jac = {j:<22s} |Jac| = {order} = {fac}")
    same_boundary = (data[0][1], data[0][2]) == (data[1][1], data[1][2])
    same_type = data[0][3] == data[1][3]
    g = math.gcd(data[0][4], data[1][4])
    print(f"      same K-theory: {same_boundary};  same vN type: {same_type};  "
          f"gcd of bulk orders = {g}" + ("  -> COPRIME" if g == 1 else ""))
    print("      => neither the boundary K-theory nor the von Neumann type determines the")
    print("         bulk arithmetic, not even which primes divide it.")
    ok = same_boundary and same_type and g == 1
    row("P", f"two cubic graphs with the same boundary K-theory AND the same von Neumann "
             f"type but COPRIME bulk orders ({data[0][4]}, {data[1][4]})",
        "recorded witnesses, n=8", EXACT, ok)
    return ok


# ---------------------------------------------------------------------------- (S)

def check_spread(samples=60):
    print("\n=== (S) how much bulk arithmetic hides inside one boundary class? ===")
    rng = random.Random(2026)
    ok = True
    print("      n    r-1   boundary K_0        distinct |Jac| found   min .. max")
    for n in (8, 10, 12, 14):
        seen = {}
        for _ in range(samples):
            G = random_cubic(n, rng)
            if G is None:
                continue
            _, _, order = bulk_jacobian(G)
            seen.setdefault(order, G)
        r = 3 * n // 2 - n + 1
        vals = sorted(seen)
        print(f"     {n:3d}  {r-1:5d}   Z^{r} + Z/{r-1:<8d}  {len(vals):14d}      "
              f"{min(vals)} .. {max(vals)}")
        if len(vals) < 2:
            ok = False
    row("S", "many distinct bulk critical groups share one boundary K-theory class",
        f"{samples} random cubic graphs per n, n = 8,10,12,14", SAMPLED, ok)
    return ok


# ---------------------------------------------------------------------------- (R)

def check_cocycle(q=2, depth=8):
    """Harmonic measure on the boundary of the (q+1)-regular tree: a cylinder set given by
    a non-backtracking ray of length L has measure c*q^{-L}.  Moving the basepoint by one
    step multiplies the measure of a cylinder by q^{+-1}.  Hence the Radon-Nikodym cocycle
    of the boundary action takes values in q^Z, which is the input to Krieger's theorem."""
    print("\n=== (R) the Radon-Nikodym cocycle takes values in q^Z ===")
    # measures of cylinders at depth L from a basepoint: nu(C) = 1/((q+1) q^{L-1})
    def nu(L):
        return Fraction(1, (q + 1) * q ** (L - 1)) if L >= 1 else Fraction(1, 1)
    ratios = set()
    for L in range(1, depth + 1):
        for L2 in range(1, depth + 1):
            rr = nu(L) / nu(L2)
            ratios.add(rr)
    # every ratio must be an integer power of q
    def is_power(x, q):
        num, den = x.numerator, x.denominator
        for e in range(-depth - 2, depth + 3):
            if Fraction(q) ** e == x:
                return True
        return False
    ok = all(is_power(x, q) for x in ratios)
    exps = sorted({next(e for e in range(-depth - 2, depth + 3) if Fraction(q) ** e == x)
                   for x in ratios})
    print(f"    q = {q}: cylinder measures nu(L) = 1/((q+1)q^(L-1)); all pairwise ratios are")
    print(f"    integer powers of q, exponents observed: {exps[0]} .. {exps[-1]}")
    print(f"    => the Busemann cocycle is q^Z-valued.  TWO RELATIONS MUST BE DISTINGUISHED:")
    print(f"       - the TAIL relation on the Cantor limit of the tower realises all of q^Z,")
    print(f"         giving type III_(1/q)  [ANCFT III, Thm 2.2];")
    print(f"       - the ORBIT relation of a lattice Gamma realises only q^(2Z) when Gamma")
    print(f"         is type-preserving (quotient bipartite), giving III_(1/q^2), and q^Z")
    print(f"         otherwise  [Robertson, Houston J. Math. 31 (2005), Cor. 1(2)].")
    print(f"       NOTE: this row is a consistency check on the normalisation only. That the")
    print(f"       essential range is attained is a quoted theorem, not a finite computation.")
    row("R", f"cylinder-measure ratios lie in q^Z (normalisation check only; both ratio-set "
             f"theorems are quoted, q={q}, depth {depth})", f"{len(ratios)} ratios", EXACT, ok)
    return ok


# --------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    check_boundary_formula()
    check_witnesses()
    check_prime_supports()
    check_spread()
    check_cocycle(2)
    check_cocycle(3)

    print("\n=== battery summary ===")
    bad = [r for r in ROWS if not r[4]]
    for key in ("K", "B", "J", "P", "S", "R"):
        sub = [r for r in ROWS if r[0] == key]
        if sub:
            print(f"  ({key})  {sum(1 for r in sub if r[4])}/{len(sub)} passed")
    print(f"  TOTAL {sum(1 for r in ROWS if r[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(r[0], r[2]) for r in bad]}"))
    print(f"  elapsed {time.time()-t0:.1f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
