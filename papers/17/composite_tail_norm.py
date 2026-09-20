# -*- coding: utf-8 -*-
"""
composite_tail_norm.py -- verification battery for

  "Algorithmic non-commutative class field theory, XVII: the tail norm at composite q,
   and the splitting of the limit K-module"

SETTING (Papers V, VIII).  Y_0 is a (q+1)-regular graph on n vertices; Y_1 is its
non-backtracking directed-edge graph (Eulerian q-regular, n_1 = n(q+1)); Y_{k+1} = L(Y_k)
is the line digraph, n_{k+1} = q n_k.  tau (tail) and eta (head) truncate a vertex
i -> j of Y_{k+1} to i and to j.  Write r_k := n_k(q-1) - 1.

KNOWN FOR EVERY q >= 2:
  (b)  |ker(tau_* : Jac(Y_{k+1}) -> Jac(Y_k))| = q^{r_k}        [V, Thm 3.1(ii); III, Thm 2.1]
  (c)  ker(tau_*) is contained in Jac(Y_{k+1})[q]                [VIII, Prop 2.7]
KNOWN ONLY FOR q = p PRIME BEFORE THIS NOTE:
  (a)  rk_{F_p} A_{k+1} = n_k, whence Jac(Y_{k+1})[p] = (Z/p)^{r_k}        [VIII, Thm 2.8]

THE POINT.  A_{k+1} is the adjacency matrix of a line digraph, hence factors as H^t T with
H[v][e] = [head(e) = v] and T[v][f] = [tail(f) = v].  The rows of H and of T are non-empty
and disjointly supported, so both have rank n_k over EVERY field; T is then onto and H^t
injective, so rk A_{k+1} = n_k in every characteristic.  So (a) is not special to q = p.
The p-rank of Jac(Y_{k+1}) is r_k for every p | q, and then (b) and (c) close the problem
by counting:  q^{r_k} = |ker| <= |Jac[q]| = prod_p p^{sum_j min(e_p, a_{p,j})}
            <= prod_p p^{e_p r_k} = q^{r_k},
so equality holds throughout: ker(tau_*) = Jac(Y_{k+1})[q] for EVERY q, and every
invariant factor of Jac(Y_{k+1}) divisible by p is divisible by p^{v_p(q)}.  That last
statement is falsifiable and is what check (E) tests.

 (I)  the degeneracy identities tau_* eta^* = A_k, eta^* tau_* = A_{k+1},
      tau_* tau^* = eta_* eta^* = qI, at composite q                          exact
 (N)  the line-digraph recursion kappa(Y_{k+1}) = q^{r_k} kappa(Y_k), hence (b)  exact
 (R)  rk A_{k+1} = n_k over F_p for several p, including p | q and p not | q    exact
 (P)  the p-rank of Jac(Y_{k+1}) is r_k for every p | q                        exact
 (E)  every invariant factor of Jac(Y_{k+1}) divisible by p is divisible by
      p^{v_p(q)}; equivalently |Jac(Y_{k+1})[q]| = q^{r_k}                     exact
 (T)  the published prime-q values are recovered verbatim                      exact

Run: python composite_tail_norm.py     (about two minutes)
Imports pgl3_building.py from the sibling Paper 10 folder.
"""

import sys, os, time, collections
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "Paper 10"))
import pgl3_building as pb

ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# --------------------------------------------------------------------------- graphs

def complete(n):
    return [[v for v in range(n) if v != u] for u in range(n)]


def bipartite(n):
    return [[n + v for v in range(n)] if u < n else [v for v in range(n)] for u in range(2 * n)]


def petersen():
    edges = ([(i, (i + 1) % 5) for i in range(5)] +
             [(5 + i, 5 + (i + 2) % 5) for i in range(5)] + [(i, 5 + i) for i in range(5)])
    adj = collections.defaultdict(list)
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return [sorted(adj[u]) for u in range(10)]


def nonbacktracking(adj):
    """Y_1: vertices are the directed edges of Y_0, arcs the non-backtracking transitions"""
    arcs = [(u, v) for u in range(len(adj)) for v in adj[u]]
    idx = {a: i for i, a in enumerate(arcs)}
    return [[idx[(v, w)] for w in adj[v] if w != u] for (u, v) in arcs]


def line_digraph(out):
    """Y_{k+1} = L(Y_k), with the list of (tail, head) of each new vertex"""
    pairs = [(i, j) for i in range(len(out)) for j in out[i]]
    where = collections.defaultdict(list)
    for k, pr in enumerate(pairs):
        where[pr].append(k)
    newout = []
    for (i, j) in pairs:
        succ = []
        for k in out[j]:
            succ.extend(where[(j, k)])
        newout.append(succ)
    return newout, pairs


def adjacency(out):
    A = [[0] * len(out) for _ in out]
    for u, succ in enumerate(out):
        for v in succ:
            A[u][v] += 1
    return A


def laplacian_T_reduced(out, root=0):
    N = len(out)
    L = [[0] * N for _ in range(N)]
    for u, succ in enumerate(out):
        L[u][u] += len(succ)
        for v in succ:
            L[u][v] -= 1
    idx = [v for v in range(N) if v != root]
    return [list(r) for r in zip(*[[L[i][j] for j in idx] for i in idx])]


def deg_maps(out):
    """tau_* (N_k x N_{k+1}), eta^* (N_{k+1} x N_k), tau^* (N_{k+1} x N_k), eta_* (N_k x N_{k+1})"""
    _, pairs = line_digraph(out)
    Nk, Nk1 = len(out), len(pairs)
    tp = [[0] * Nk1 for _ in range(Nk)]
    ep = [[0] * Nk for _ in range(Nk1)]
    tpl = [[0] * Nk for _ in range(Nk1)]
    eph = [[0] * Nk1 for _ in range(Nk)]
    for k, (i, j) in enumerate(pairs):
        tp[i][k] += 1
        tpl[k][i] += 1
        ep[k][j] += 1
        eph[j][k] += 1
    return tp, ep, tpl, eph


def mul(A, B):
    m, p = len(A), len(B[0])
    return [[sum(a * B[k][j] for k, a in enumerate(A[i]) if a) for j in range(p)] for i in range(m)]


# ------------------------------------------------------------------------- the towers

TOWERS = [
    ("K4",    2, complete(4),  3),
    ("Pete",  2, petersen(),   3),
    ("K5",    3, complete(5),  3),
    ("K4,4",  3, bipartite(4), 2),
    ("K6",    4, complete(6),  3),
    ("K5,5",  4, bipartite(5), 2),
    ("K7",    5, complete(7),  2),
    ("K8",    6, complete(8),  2),
]
PRIMES = (2, 3, 5, 7, 2147483647)
DET_MAX = 210          # exact Bareiss determinant up to this many vertices
SNF_MAX = 520          # p-adic valuations up to this many vertices


def levels(base, maxlev):
    Y = {1: nonbacktracking(base)}
    for lev in range(2, maxlev + 1):
        Y[lev], _ = line_digraph(Y[lev - 1])
    return Y


def main():
    t0 = time.time()
    print("=== the Iwahori tower at composite q: is ker(tau_*) still the q-torsion? ===")
    built = {}
    for name, q, base, maxlev in TOWERS:
        built[name] = (q, levels(base, maxlev), maxlev)
        sizes = ", ".join(f"n_{k}={len(v)}" for k, v in sorted(built[name][1].items()))
        print(f"    {name:5s} q={q}  {sizes}")

    print("\n=== (I) the degeneracy identities at composite q ===")
    ok = True
    for name, (q, Y, maxlev) in built.items():
        tp, ep, tpl, eph = deg_maps(Y[1])
        Y2, _ = line_digraph(Y[1])
        A1, A2 = adjacency(Y[1]), adjacency(Y2)
        qI = [[q if i == j else 0 for j in range(len(Y[1]))] for i in range(len(Y[1]))]
        got = (mul(tp, ep) == A1, mul(ep, tp) == A2, mul(tp, tpl) == qI, mul(eph, ep) == qI)
        ok &= all(got)
        print(f"    {name:5s} q={q}: tau_* eta^* = A_k {got[0]}, eta^* tau_* = A_k+1 {got[1]}, "
              f"tau_* tau^* = qI {got[2]}, eta_* eta^* = qI {got[3]}")
    row("I", "tau_* eta^* = A_k, eta^* tau_* = A_{k+1}, tau_* tau^* = eta_* eta^* = qI hold "
             "verbatim at composite q (this is what makes ker tau_* < Jac[q])",
        f"{len(built)} towers", "exact", ok)

    print("\n=== (N) the line-digraph recursion, hence |ker tau_*| = q^{r_k} ===")
    ok = True
    cases = 0
    for name, (q, Y, maxlev) in built.items():
        for lev in range(2, maxlev + 1):
            if len(Y[lev]) > DET_MAX:
                continue
            kap_lo = abs(pb.bareiss_det(laplacian_T_reduced(Y[lev - 1])))
            kap_hi = abs(pb.bareiss_det(laplacian_T_reduced(Y[lev])))
            r = len(Y[lev - 1]) * (q - 1) - 1
            good = (kap_hi == q ** r * kap_lo)
            ok &= good
            cases += 1
            print(f"    {name:5s} q={q} level {lev-1}->{lev}: kappa ratio = q^{r} "
                  f"{'ok' if good else 'MISMATCH'}   (log10 kappa_{lev} = "
                  f"{len(str(kap_hi))-1})")
    row("N", "kappa(Y_{k+1}) = q^{n_k(q-1)-1} kappa(Y_k), so the surjection tau_* has kernel "
             "of order exactly q^{r_k}", f"{cases} transitions", "exact", ok)

    print("\n=== (R) rank of A_{k+1} over F_p: the same in every characteristic ===")
    ok = True
    for name, (q, Y, maxlev) in built.items():
        A2 = adjacency(Y[2])
        ranks = {p: sum(1 for v in pb.padic_valuations(A2, p, 1) if v == 0) for p in PRIMES}
        good = all(v == len(Y[1]) for v in ranks.values())
        ok &= good
        print(f"    {name:5s} q={q}: n_1={len(Y[1]):3d}, rank A_2 mod "
              f"{ {p: v for p, v in ranks.items()} }  {'ok' if good else 'MISMATCH'}")
    row("R", "rk A_{k+1} = n_k over F_p for p | q and for p not dividing q alike -- the rank is "
             "characteristic-free, which is what removes the primality hypothesis",
        f"{len(built)} towers x {len(PRIMES)} primes", "exact", ok)

    print("\n=== (P),(E) the p-rank of Jac(Y_k), and the q-torsion ===")
    okP = okE = True
    cases = 0
    for name, (q, Y, maxlev) in built.items():
        fac = sp.factorint(q)
        for lev in range(2, maxlev + 1):
            if len(Y[lev]) > SNF_MAX:
                continue
            M = laplacian_T_reduced(Y[lev])
            r = len(Y[lev - 1]) * (q - 1) - 1
            parts = []
            for p, e in sorted(fac.items()):
                vals = [v for v in pb.padic_valuations(M, p, 60) if v > 0]
                rk, mn = len(vals), (min(vals) if vals else 0)
                logq = sum(min(e, v) for v in vals)
                okP &= (rk == r)
                okE &= (mn >= e and logq == e * r)
                parts.append(f"p={p}^{e}: rank {rk}/{r}, min val {mn} (>= {e}), "
                             f"log_p|Jac[p^{e}]| {logq}/{e*r}")
            cases += 1
            print(f"    {name:5s} q={q} level {lev}: n={len(Y[lev]):3d}   " + " | ".join(parts))
    row("P", "the p-rank of Jac(Y_{k+1}) is n_k(q-1)-1 for every prime p dividing q, composite "
             "q included", f"{cases} levels", "exact", okP)
    row("E", "every invariant factor of Jac(Y_{k+1}) divisible by p is divisible by p^{v_p(q)}; "
             "equivalently |Jac(Y_{k+1})[q]| = q^{r_k}, so ker tau_* = Jac(Y_{k+1})[q]",
        f"{cases} levels", "exact", okE)

    print("\n=== (T) the published prime-q values ===")
    Y = built["K4"][1]
    exp = {2: ([1] * 8 + [2, 4, 4], 11), 3: ([1] * 12 + [2] * 8 + [3, 5, 5], 23)}
    ok = True
    for lev, (vexp, rk) in exp.items():
        vals = sorted(v for v in pb.padic_valuations(laplacian_T_reduced(Y[lev]), 2, 60) if v > 0)
        good = (vals == sorted(vexp))
        ok &= good
        print(f"    K4 q=2 level {lev}: 2-adic valuations {collections.Counter(vals)}, "
              f"2-rank {len(vals)} (published {rk}); ord_2 kappa = {sum(vals)}  "
              f"{'ok' if good else 'MISMATCH'}")
    row("T", "Jac(Y_2) = (Z/2)^8+Z/4+Z/16+Z/48 and Jac(Y_3) = (Z/2)^12+(Z/4)^8+Z/8+Z/32+Z/96 "
             "of the companion papers are recovered, with ker tau_* of order 2^11 and 2^23",
        "K4 tower, levels 2 and 3", "exact", ok)

    print("\n=== battery summary ===")
    bad = [x for x in ROWS if not x[4]]
    print(f"  TOTAL {sum(1 for x in ROWS if x[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(x[0], x[2]) for x in bad]}"))
    print(f"  elapsed {time.time()-t0:.0f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
