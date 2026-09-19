#!/usr/bin/env python3
"""Verification of the Hecke -> Smith-normal-form K-theory claims (Sec. 3).

For a finite connected (q+1)-regular graph Y with vertex adjacency (Hecke) matrix T:
  Delta = (q+1)I - T                (Hecke-Laplacian; = I - ([E_p] - q))
  A_NB  = non-backtracking edge matrix (2m x 2m)
We compute SNF of Delta, of I - T, and of I - A_NB^t, and compare with
  kappa(Y) = (1/n) prod_{lambda != q+1} (q+1 - lambda)   (Kirchhoff)
  det(I-T) = prod (1 - lambda)
  Bass-Ihara: det(I - u A_NB) = (1-u^2)^{m-n} det((1+q u^2)I - u T)
"""
import itertools
from fractions import Fraction
import numpy as np
import sympy as sp


def snf_invariants(M):
    """Smith normal form invariant factors of an integer sympy Matrix."""
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(sp.Matrix(M))
    d = [abs(S[i, i]) for i in range(min(S.rows, S.cols))]
    return d


def cokernel(M):
    d = snf_invariants(M)
    n = sp.Matrix(M).cols
    tors = [x for x in d if x not in (0, 1)]
    free = n - len([x for x in d if x != 0])
    return free, tors


def nonbacktracking(edges, n):
    """edges: list of (u,v) oriented edges (both directions present)."""
    m = len(edges)
    B = sp.zeros(m, m)
    for i, (u, v) in enumerate(edges):
        for j, (x, y) in enumerate(edges):
            if v == x and not (y == u and x == v):
                if not (x == v and y == u):
                    B[i, j] = 1
    return B


def oriented_edges(adj):
    e = []
    n = len(adj)
    for u in range(n):
        for v in range(n):
            for _ in range(adj[u][v]):
                e.append((u, v))
    return e


def report(name, adj, q):
    n = len(adj)
    T = sp.Matrix(adj)
    I = sp.eye(n)
    lam = sorted([sp.nsimplify(x) for x in sp.Matrix(adj).eigenvals(multiple=True)],
                 key=lambda z: -float(z))
    Delta = (q + 1) * I - T
    edges = oriented_edges(adj)
    m2 = len(edges)          # 2m oriented edges
    ANB = nonbacktracking(edges, n)

    kappa = sp.prod([(q + 1) - l for l in lam[1:]]) / n
    freeD, torsD = cokernel(Delta)
    freeT, torsT = cokernel(I - T)
    freeB, torsB = cokernel((sp.eye(m2) - ANB).T)

    # Bass-Ihara check
    u = sp.symbols('u')
    lhs = sp.factor(sp.expand((sp.eye(m2) - u * ANB).det()))
    rhs = sp.factor(sp.expand(((1 - u**2) ** (m2 // 2 - n)) *
                              (((1 + q * u**2) * sp.eye(n) - u * T).det())))
    bass_ok = sp.simplify(lhs - rhs) == 0

    print(f"--- {name}:  n={n}, q+1={q+1}, oriented edges={m2}")
    print(f"    spectrum T           : {lam}")
    print(f"    kappa (spanning trees): {sp.nsimplify(kappa)}")
    print(f"    SNF(Delta)  coker = Z^{freeD} (+) {torsD}   |tors|={sp.prod(torsD) if torsD else 1}")
    print(f"    SNF(I-T)    coker = Z^{freeT} (+) {torsT}   det(I-T)={ (I-T).det() }")
    print(f"    SNF(I-A_NB^t) coker = Z^{freeB} (+) {torsB}")
    print(f"    Bass-Ihara identity verified: {bass_ok}")
    print()
    return dict(kappa=kappa, torsD=torsD, torsB=torsB, freeB=freeB)


def complete(n):
    return [[1 if i != j else 0 for j in range(n)] for i in range(n)]


def petersen():
    # standard Kneser graph K(5,2)
    verts = list(itertools.combinations(range(5), 2))
    n = len(verts)
    adj = [[0] * n for _ in range(n)]
    for i, a in enumerate(verts):
        for j, b in enumerate(verts):
            if i != j and not set(a) & set(b):
                adj[i][j] = 1
    return adj


def cube():
    n = 8
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for b in range(3):
            j = i ^ (1 << b)
            adj[i][j] = 1
    return adj


def k33():
    n = 6
    adj = [[0] * n for _ in range(n)]
    for i in range(3):
        for j in range(3, 6):
            adj[i][j] = adj[j][i] = 1
    return adj


def heawood():
    # incidence graph of the Fano plane, 3-regular, 14 vertices
    pts = list(range(7))
    lines = [(0, 1, 3), (1, 2, 4), (2, 3, 5), (3, 4, 6), (4, 5, 0), (5, 6, 1), (6, 0, 2)]
    n = 14
    adj = [[0] * n for _ in range(n)]
    for li, L in enumerate(lines):
        for p in L:
            adj[p][7 + li] = adj[7 + li][p] = 1
    return adj


if __name__ == "__main__":
    report("K_4  (q=2)", complete(4), 2)
    report("K_5  (q=3)", complete(5), 3)
    report("K_{3,3} (q=2)", k33(), 2)
    report("Q_3 cube (q=2)", cube(), 2)
    report("Petersen (q=2)", petersen(), 2)
    report("Heawood (q=2)", heawood(), 2)

    # arithmetic illustration: Brandt matrix D=11, p=3
    B = sp.Matrix([[2, 2], [3, 1]])
    print("Brandt B(3), D=11 :", B.eigenvals(), " 4 - a_3 =", 4 - (-1),
          " #E_{11a}(F_3) = 5")
