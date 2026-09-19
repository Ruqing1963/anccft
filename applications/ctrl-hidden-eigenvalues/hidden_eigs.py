# -*- coding: utf-8 -*-
"""
hidden_eigs.py -- verification battery for

    "Power-sum lower bounds on hidden eigenvalues in positive realizations"

Setting.  A transfer function of order n with pole multiset Lambda (|Lambda| = n,
spectral radius rho attained at a positive real dominant pole) has a positive
realization of size N iff there is a nonnegative N x N matrix A with
Lambda contained in spec(A) and rho(A) = rho.  The s = N - n eigenvalues of A outside
Lambda are the HIDDEN eigenvalues; characterising them is an open problem
[Benvenuti-Farina, IEEE TAC 49(5):651-664, 2004, Sec. IX].

Theorem 1 (this paper).  s >= ceil( max_{k>=1} [ -p_k(Lambda) ]_+ / rho^k ),
where p_k(Lambda) = sum of the k-th powers of the poles.

 (T)  validity: on nonnegative integer matrices with Lambda a rational factor of the
      characteristic polynomial, the bound never exceeds the true s.     exact / numeric rho
 (K)  k = 1 reproduces Benvenuti's zeta [ELA 36:367-384 (2020), Thm 4.7]. exact
 (I)  k = 2 strictly improves it on the family M_m(theta):  their bound is
      vacuous (N >= n) while ours rises to their own upper bound 4m.     exact
 (S)  sharpness: on M_m = {1} u m{i,-i} the bound gives N >= 4m and the direct sum of
      m copies of the 4-cycle permutation attains it, so N_min = 4m.     exact
 (E)  Benvenuti-Farina tutorial Example 8: poles {1, +-0.9i}; the bound gives s >= 1 and
      an explicit rational nonnegative 4x4 circulant attains it.         exact
 (C)  the ceiling: the bound never exceeds n, so this method cannot prove more than
      N <= 2n -- it is complementary to digraph bounds, not competitive.  exact
 (Z)  why no such bound exists in the Boyle-Handelman/Kim-Ormes-Roush category:
      there one appends only ZEROS, which change no power sum.            exact

Exact rational arithmetic throughout; the only numeric quantity is the spectral radius in
check (T), and it is labelled. Requires sympy and mpmath. Run: python hidden_eigs.py
"""

import sys, time, itertools
from fractions import Fraction
import sympy as sp
import mpmath as mp

mp.mp.dps = 40

EXACT, NUM = "exact", "exact + numeric rho"
ROWS = []


def row(key, desc, cases, mode, ok):
    ROWS.append((key, desc, cases, mode, bool(ok)))
    print(("  [ok]   " if ok else "  [FAIL] ") + f"({key}) {desc}   [{cases}] [{mode}]")
    return bool(ok)


# ------------------------------------------------------------------ power sums (exact)

def power_sums(coeffs, kmax):
    """Newton's identities. coeffs = [1, c1, ..., cd] of a monic degree-d polynomial
    f(x) = x^d + c1 x^(d-1) + ... + cd.  Returns [p_1, ..., p_kmax] as exact rationals."""
    d = len(coeffs) - 1
    c = [sp.sympify(x) for x in coeffs]
    p = []
    for k in range(1, kmax + 1):
        if k <= d:
            s = k * c[k] + sum(c[i] * p[k - i - 1] for i in range(1, k))
        else:
            s = sum(c[i] * p[k - i - 1] for i in range(1, d + 1))
        p.append(sp.expand(-s))
    return p


def bound_from_coeffs(coeffs, rho, kmax=40):
    """max_k [-p_k]_+ / rho^k  and the k attaining it."""
    ps = power_sums(coeffs, kmax)
    best, bk = sp.Integer(0), None
    for k in range(1, kmax + 1):
        v = sp.Max(-ps[k - 1], 0) / rho ** k
        if sp.simplify(v - best) > 0:
            best, bk = sp.simplify(v), k
    return best, bk, ps


def poly_from_roots(roots):
    x = sp.Symbol('x')
    f = sp.expand(sp.prod([x - r for r in roots]))
    p = sp.Poly(f, x)
    return [sp.expand(cc) for cc in p.all_coeffs()]


# ---------------------------------------------------------------------------- (T)

def check_validity(trials=120):
    print("\n=== (T) validity of Theorem 1 on nonnegative integer matrices ===")
    import random
    random.seed(2026)
    x = sp.Symbol('x')
    viol = 0
    tested = 0
    tight = 0
    while tested < trials:
        N = random.randint(3, 7)
        B = sp.Matrix(N, N, lambda i, j: random.choice([0, 0, 1, 1, 2, 3]))
        if all(v == 0 for v in B):
            continue
        cp = sp.Poly(B.charpoly(x).as_expr(), x)
        facs = sp.factor_list(cp.as_expr())[1]
        if len(facs) < 2:
            continue
        # the spectral radius, numerically, and the factor carrying it
        ev = mp.polyroots([mp.mpf(str(c)) for c in cp.all_coeffs()],
                          maxsteps=200, extraprec=200)
        rho_num = max(abs(e) for e in ev)
        if rho_num < 1e-9:
            continue
        perron_fac = None
        for f, mult in facs:
            fp = sp.Poly(f, x)
            if fp.degree() == 0:
                continue
            rts = mp.polyroots([mp.mpf(str(c)) for c in fp.all_coeffs()],
                               maxsteps=200, extraprec=200)
            if any(abs(abs(r) - rho_num) < 1e-20 for r in rts):
                perron_fac = (f, mult)
                break
        if perron_fac is None:
            continue
        # Lambda = roots of the Perron factor, times a random subset of the other factors
        chosen = [perron_fac[0]]
        for f, mult in facs:
            if f == perron_fac[0] or sp.Poly(f, x).degree() == 0:
                continue
            if random.random() < 0.5:
                chosen.append(f)
        g = sp.Poly(sp.expand(sp.prod(chosen)), x)
        g = g.monic()
        d = g.degree()
        if d >= N:
            continue
        s_true = N - d
        rho = sp.nsimplify(mp.nstr(rho_num, 30), rational=False)
        # use a rational LOWER bound for rho: that makes the bound LARGER, i.e. harder
        rho_lo = sp.Rational(int(mp.floor(rho_num * 10**12)), 10**12)
        if rho_lo <= 0:
            continue
        b, k, _ = bound_from_coeffs(g.all_coeffs(), rho_lo)
        tested += 1
        if b > s_true:
            viol += 1
            print(f"   VIOLATION: N={N} deg={d} s={s_true} bound={sp.nsimplify(b)} k={k}")
        if sp.ceiling(b) == s_true and s_true > 0:
            tight += 1
    row("T", f"bound <= true hidden count in {tested} instances ({viol} violations, "
             f"{tight} attained exactly)", f"{tested} integer matrices", NUM, viol == 0)
    return viol == 0


# ---------------------------------------------------------------------------- (K)(I)

def benvenuti_lower(m, theta):
    """Benvenuti, ELA 36 (2020), Thm 4.7, for M_m(theta) = {1} u m{i*theta, -i*theta}:
    N >= N_D + sum_{|s|<1} m(s) + zeta, with N_D = 1 (only 1 is dominant, r = 1),
    sum_{|s|<1} m(s) = 2m, and zeta = 0 because m_{R1} + sum s*m(s) = 1 + 0 >= 0."""
    N_D = 1
    tail = 2 * m
    mR1 = 1
    s_sum = sp.Integer(0)          # i*theta*m - i*theta*m = 0
    zeta = sp.ceiling(-mR1 - s_sum) if (mR1 + s_sum) < 0 else sp.Integer(0)
    return N_D + tail + zeta


def ours_lower(m, theta):
    """Theorem 1 for M_m(theta). rho = 1."""
    roots = [sp.Integer(1)] + [sp.I * theta, -sp.I * theta] * m
    coeffs = poly_from_roots(roots)
    b, k, ps = bound_from_coeffs(coeffs, sp.Integer(1), kmax=8)
    n = len(roots)
    return n + int(sp.ceiling(b)), b, k, ps


def check_k1_and_improvement():
    print("\n=== (K) k = 1 reproduces Benvenuti's zeta; (I) k = 2 strictly improves ===")
    # (K): a family where zeta actually bites -- put mass on a negative real pole
    ok_k1 = True
    for a in (sp.Rational(1, 2), sp.Rational(9, 10)):
        for mult in (2, 3, 5):
            roots = [sp.Integer(1)] + [-a] * mult
            coeffs = poly_from_roots(roots)
            ps = power_sums(coeffs, 1)
            p1 = ps[0]                                    # 1 - mult*a
            ours = sp.Max(-p1, 0)                         # rho = 1
            benv = sp.Max(sp.ceiling(-1 + mult * a), 0)   # -m_R1 - sum s*m(s) = -1 + mult*a
            same = sp.ceiling(ours) == benv
            ok_k1 &= bool(same)
    row("K", "at k = 1 the bound equals Benvenuti's zeta [ELA 2020, Thm 4.7] on every "
             "tested pole set", "a in {1/2, 9/10}, multiplicity 2,3,5", EXACT, ok_k1)

    print("\n    family M_m(theta) = {1} u m copies of {i*theta, -i*theta},  n = 2m+1")
    print("      theta   m    n   Benvenuti Thm 4.7   ours (k=2)   Benvenuti Thm 4.4 (upper)")
    strictly_better = 0
    total = 0
    for theta in (sp.Rational(9, 10), sp.Rational(99, 100), sp.Rational(999, 1000)):
        for m in (1, 2, 4, 8, 16):
            n = 2 * m + 1
            lo_b = benvenuti_lower(m, theta)
            lo_o, b, k, ps = ours_lower(m, theta)
            up = 4 * m                      # Thm 4.4 with kappa(i*theta) = 4
            total += 1
            if lo_o > lo_b:
                strictly_better += 1
            print(f"      {str(theta):>9s} {m:3d} {n:4d}   {int(lo_b):11d}       "
                  f"{lo_o:6d}       {up:10d}")
    row("I", f"our k=2 bound strictly exceeds Benvenuti's Thm 4.7 lower bound in "
             f"{strictly_better}/{total} cases", "3 thetas x 5 sizes", EXACT,
        strictly_better == total)
    return ok_k1 and strictly_better == total


# ---------------------------------------------------------------------------- (S)

def check_sharpness():
    print("\n=== (S) sharpness: M_m = {1} u m{i,-i}, the bound is attained ===")
    print("      m    n   p_2   bound s>=   N >=    (+)C_4 gives N    equal?")
    ok = True
    for m in (1, 2, 3, 4, 6, 8, 16, 32):
        roots = [sp.Integer(1)] + [sp.I, -sp.I] * m
        n = len(roots)
        coeffs = poly_from_roots(roots)
        b, k, ps = bound_from_coeffs(coeffs, sp.Integer(1), kmax=6)
        s_lo = int(sp.ceiling(b))
        N_lo = n + s_lo
        N_constr = 4 * m
        eq = (N_lo == N_constr)
        ok &= eq
        print(f"      {m:3d} {n:4d} {int(ps[1]):5d}   {s_lo:8d}  {N_lo:5d}    "
              f"{N_constr:12d}      {'YES' if eq else 'no'}")
    # the construction really is nonnegative and really has those eigenvalues
    x = sp.Symbol('x')
    P = sp.Matrix(4, 4, lambda i, j: 1 if j == (i + 1) % 4 else 0)
    cp = sp.factor(P.charpoly(x).as_expr())
    nonneg = all(v >= 0 for v in P)
    row("S", f"bound attained for every m tested; the 4-cycle has charpoly {cp} and "
             f"nonnegative entries", "m = 1..32", EXACT, ok and nonneg)
    return ok and nonneg


# ---------------------------------------------------------------------------- (E)

def check_example8():
    print("\n=== (E) Benvenuti-Farina tutorial, Example 8: poles {1, +-0.9i} ===")
    th = sp.Rational(9, 10)
    roots = [sp.Integer(1), sp.I * th, -sp.I * th]
    coeffs = poly_from_roots(roots)
    b, k, ps = bound_from_coeffs(coeffs, sp.Integer(1), kmax=8)
    print(f"    p_2 = {ps[1]} < 0,  bound s >= {b} = {float(b):.3f}  ->  ceil = "
          f"{int(sp.ceiling(b))} at k = {k}")
    # an explicit rational nonnegative realizer: circulant with first row (0, 19/20, 0, 1/20)
    c0, c1, c2, c3 = 0, sp.Rational(19, 20), 0, sp.Rational(1, 20)
    A = sp.Matrix(4, 4, lambda i, j: [c0, c1, c2, c3][(j - i) % 4])
    x = sp.Symbol('x')
    cp = sp.factor(sp.simplify(A.charpoly(x).as_expr()))
    spec_ok = sp.simplify(A.charpoly(x).as_expr() - (x - 1) * (x + 1) * (x**2 + th**2)) == 0
    nonneg = all(v >= 0 for v in A)
    print(f"    realizer: circulant first row (0, 19/20, 0, 1/20), charpoly = {cp}")
    print(f"    spectrum {{1, -1, +-0.9i}} contains the poles, so s = 1 -- the bound is tight")
    row("E", "bound = 1 and an explicit rational nonnegative 4x4 attains it",
        "Example 8", EXACT, int(sp.ceiling(b)) == 1 and spec_ok and nonneg)
    return spec_ok and nonneg


# ---------------------------------------------------------------------------- (C)(Z)

def check_ceiling_and_zeros():
    print("\n=== (C) the ceiling of the method, and (Z) why it is empty for KOR ===")
    import random
    random.seed(11)
    ok_c = True
    worst = sp.Integer(0)
    for _ in range(300):
        d = random.randint(2, 8)
        roots = [sp.Integer(1)]
        while len(roots) < d:
            a = sp.Rational(random.randint(-9, 9), 10)
            bb = sp.Rational(random.randint(0, 9), 10)
            if a**2 + bb**2 > 1:
                continue
            if bb == 0:
                roots.append(a)
            elif len(roots) + 2 <= d:
                roots += [a + sp.I * bb, a - sp.I * bb]
            else:
                roots.append(a)
        coeffs = poly_from_roots(roots)
        b, k, _ = bound_from_coeffs(coeffs, sp.Integer(1), kmax=25)
        n = len(roots)
        if b > n:
            ok_c = False
        worst = sp.Max(worst, sp.nsimplify(b / n))
    row("C", f"bound <= n always, so the method cannot prove more than N <= 2n "
             f"(worst ratio bound/n observed: {float(worst):.3f})", "300 random pole sets",
        EXACT, ok_c)

    # (Z) appending zeros changes no power sum
    roots = [sp.Integer(1), sp.I, -sp.I]
    c1 = poly_from_roots(roots)
    c2 = poly_from_roots(roots + [sp.Integer(0)] * 7)
    p1 = power_sums(c1, 10)
    p2 = power_sums(c2, 10)
    same = all(sp.simplify(a - b2) == 0 for a, b2 in zip(p1, p2))
    row("Z", "appending any number of zeros leaves every power sum unchanged: in the "
             "Boyle-Handelman/KOR category no trace bound on the dilation can exist",
        "7 zeros appended", EXACT, same)
    return ok_c and same


# --------------------------------------------------------------------------- main

def main():
    t0 = time.time()
    check_validity()
    check_k1_and_improvement()
    check_sharpness()
    check_example8()
    check_ceiling_and_zeros()

    print("\n=== battery summary ===")
    bad = [r for r in ROWS if not r[4]]
    for key in ("T", "K", "I", "S", "E", "C", "Z"):
        sub = [r for r in ROWS if r[0] == key]
        if sub:
            print(f"  ({key})  {sum(1 for r in sub if r[4])}/{len(sub)} passed")
    print(f"  TOTAL {sum(1 for r in ROWS if r[4])}/{len(ROWS)} passed"
          + ("" if not bad else f"   FAILURES: {[(r[0], r[2]) for r in bad]}"))
    print(f"  elapsed {time.time()-t0:.1f}s")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
