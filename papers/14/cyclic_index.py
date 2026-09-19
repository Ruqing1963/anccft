"""
Paper XIV battery (cyclic_index.py): the graph L-invariant as a second variation of the spectral action,
as the effective mass of the bottom Bloch band, as the diffusion constant of the Z-cover, and as the Albanese
length of the voltage class.  Base graphs K4 (q = 2) with the four voltage assignments of [II, Table 2] and
K5 (q = 3) with two voltage assignments.

 (C)  no-go, cyclic side:  for A = C(V) + C(E^+) (commutative, semisimple) the space of cyclic 1-cocycles
      Z^1_lambda(A) is zero, so HC^1(A) = 0 and there is no class Ch^1(D_p) in HC^1(A) to pair with
 (N)  no-go, the specification's formula:  Tr_s([D_p, alpha] D_p e^{-t D_p^2}) -> Tr_s([D_p, alpha] D_p) = 4 sum_e alpha(e)
      as t -> 0+ : linear in alpha, not the quadratic quantity -||h_alpha||^2
 (H)  heat-kernel formulas (exact, symbolic t):  ||h_alpha||^2 = lim_{t->oo} <alpha, e^{-t d^t d} alpha>
                                                 = ||alpha||^2 - int_0^oo <d alpha, e^{-s Delta_0} d alpha> ds
 (B)  Bloch band: the bottom eigenvalue of the twisted Laplacian D(e^{i phi}) is
      lambda_0(phi) = (||h_alpha||^2 / n) phi^2 + O(phi^4);  exact second-order perturbation coefficient and
      high-precision numerics (Richardson)
 (S)  spectral action: (n/2) (1/t) d^2/dphi^2 Tr e^{-t D(e^{i phi})} |_{phi=0} -> L(Y_.) = -||h_alpha||^2 as t -> oo
 (V)  diffusion: the asymptotic variance of the voltage functional S_N = sum alpha(e_i) along simple random walk
      is sigma^2 = ||h_alpha||^2 / m = -L / m  (exact Poisson-equation formula), and Var(S_N) - N sigma^2 converges
      (exact Var(S_N) for N <= 24 from the truncated moment generating matrix)
 (A)  Albanese length: ||h_alpha||^2 = a_c^t G^{-1} a_c with a_c the periods of alpha on a cycle basis and G the Gram
      matrix of the basis; gauge invariance under alpha -> alpha + d^t phi
 (T)  [IX, Thm 1.2(ii)] on K5:  a_2 = Theta''(0)/2 = -kappa(Y_0) ||h_alpha||^2
"""
import sympy as sp, mpmath as mp
from fractions import Fraction
import itertools, random

t, s, phi = sp.symbols('t s phi', real=True)

def complete_graph(nv):
    return nv, [(a, b) for a in range(nv) for b in range(a+1, nv)]

class Base:
    def __init__(self, nv, BASE, alpha, name):
        self.n, self.E, self.alpha, self.name = nv, BASE, list(alpha), name
        self.m = len(BASE); self.q = (2*self.m)//nv - 1
        n, m = nv, self.m
        DEL = sp.zeros(n, m)
        for j, (a, b) in enumerate(BASE): DEL[b, j] += 1; DEL[a, j] -= 1
        self.DEL = DEL; self.LAP = DEL*DEL.T
        J = sp.ones(n, n); self.LAPP = (self.LAP + J/n).inv() - J/n
        self.PI = sp.eye(m) - DEL.T*self.LAPP*DEL
        a = sp.Matrix(alpha); self.a = a
        self.hE = sp.Rational((a.T*self.PI*a)[0, 0]); self.L = -self.hE
        self.kappa0 = abs(self.LAP[1:, 1:].det())
    def pencil(self, var):
        A = sp.zeros(self.n, self.n)
        for (a, b), v in zip(self.E, self.alpha):
            A[a, b] += var**v; A[b, a] += var**(-v)
        return (self.q+1)*sp.eye(self.n) - A
    def twisted(self, ph):          # Hermitian twisted Laplacian D(e^{i ph}) as mpmath matrix
        n = self.n; M = mp.matrix(n, n)
        for i in range(n): M[i, i] = self.q + 1
        for (a, b), v in zip(self.E, self.alpha):
            z = mp.expj(ph*v); M[a, b] -= z; M[b, a] -= mp.conj(z)
        return M

def cyclic_1_cocycles(N):
    """dimension of Z^1_lambda for the commutative semisimple algebra C^N (basis of minimal idempotents)."""
    idx = lambda i, j: i*N + j
    rows = []
    for i in range(N):
        for j in range(N):
            r = [0]*(N*N); r[idx(i, j)] += 1; r[idx(j, i)] += 1; rows.append(r)          # cyclicity
    for i, j, k in itertools.product(range(N), repeat=3):                                   # Hochschild
        r = [0]*(N*N)
        if i == j: r[idx(i, k)] += 1
        if j == k: r[idx(i, j)] -= 1
        if k == i: r[idx(i, j)] += 1
        if any(r): rows.append(r)
    Mx = sp.Matrix(rows)
    return N*N - Mx.rank()

def brief_formula(B):
    """Tr_s([D, alpha] D) with alpha acting as multiplication on l^2(E^+):  the t -> 0+ limit of the
    specification's pairing."""
    n, m = B.n, B.m; DEL = B.DEL
    D = sp.zeros(n+m, n+m); D[:n, n:] = DEL; D[n:, :n] = DEL.T
    Al = sp.zeros(n+m, n+m); Al[n:, n:] = sp.diag(*B.alpha)
    gam = sp.diag(*([1]*n + [-1]*m))
    X = (D*Al - Al*D)*D
    return (gam*X).trace()

def heat_formulas(B):
    E1 = (-t*B.DEL.T*B.DEL).exp()
    qf = sp.simplify((B.a.T*E1*B.a)[0, 0])
    lim = sp.limit(qf, t, sp.oo)
    da = B.DEL*B.a
    E0 = (-s*B.LAP).exp()
    integrand = sp.simplify((da.T*E0*da)[0, 0])
    integ = sp.integrate(integrand, (s, 0, sp.oo))
    return lim, sp.nsimplify(B.a.dot(B.a) - integ), qf

def bloch_exact(B):
    """second-order Rayleigh-Schroedinger coefficient of the ground eigenvalue of D(e^{i phi})."""
    n = B.n; one = sp.ones(n, 1)
    D1 = sp.zeros(n, n); D2 = sp.zeros(n, n)
    for (a, b), v in zip(B.E, B.alpha):
        D1[a, b] -= sp.I*v; D1[b, a] -= -sp.I*v            # -A_1,  A_1 = d/dphi A(e^{i phi})
        D2[a, b] -= -sp.Rational(v*v, 2); D2[b, a] -= -sp.Rational(v*v, 2)   # -A_2
    first = sp.simplify((one.T*D1*one)[0, 0])/n
    c = sp.simplify(((one.T*D2*one)[0, 0] - ((D1*one).H*B.LAPP*(D1*one))[0, 0]))/n
    return first, sp.nsimplify(c)

def bloch_numeric(B, h=mp.mpf('1e-3')):
    lam = lambda ph: min(mp.eighe(B.twisted(ph))[0])
    f = lambda ph: lam(ph)/ph**2
    return (4*f(h/2) - f(h))/3            # Richardson: lambda_0 even in phi, error O(h^4)

def spectral_action(B, tt, h=mp.mpf('1e-7')):
    g = lambda ph: sum(mp.e**(-tt*e) for e in mp.eighe(B.twisted(ph))[0])
    d2 = (g(h) - 2*g(0) + g(-h))/h**2
    return mp.mpf(B.n)/2/tt*d2

def diffusion(B):
    n, q = B.n, B.q; P = sp.Rational(1, q+1)
    f = {}
    for (a, b), v in zip(B.E, B.alpha): f[(a, b)] = v; f[(b, a)] = -v
    F = sp.Matrix([sum(P*f[(x, y)] for y in range(n) if (x, y) in f) for x in range(n)])
    g = (q+1)*B.LAPP*F                     # (I - P) g = F
    sig2 = sum(sp.Rational(1, n)*P*(f[(x, y)] + g[y] - g[x])**2 for (x, y) in f)
    # exact Var(S_N) from the truncated moment generating matrix over Q[s]/(s^3)
    def tmul(X, Y):
        N = len(X); Z = [[(Fraction(0),)*3 for _ in range(N)] for _ in range(N)]
        for i in range(N):
            for j in range(N):
                c0 = c1 = c2 = Fraction(0)
                for k in range(N):
                    x, y = X[i][k], Y[k][j]
                    c0 += x[0]*y[0]; c1 += x[0]*y[1] + x[1]*y[0]; c2 += x[0]*y[2] + x[1]*y[1] + x[2]*y[0]
                Z[i][j] = (c0, c1, c2)
        return Z
    Mx = [[(Fraction(0),)*3 for _ in range(n)] for _ in range(n)]
    for (x, y), v in f.items():
        Mx[x][y] = (Fraction(1, q+1), Fraction(v, q+1), Fraction(v*v, 2*(q+1)))
    var = {}; Pw = Mx; N = 1
    while N <= 24:
        tot = [sum((Pw[i][j][k] for i in range(n) for j in range(n)), Fraction(0)) for k in range(3)]
        mean = tot[1]/n; var[N] = 2*tot[2]/n - mean**2
        Pw = tmul(Pw, Mx); N += 1
    return sp.nsimplify(sig2), var

def albanese(B):
    n, m = B.n, B.m
    # spanning tree by BFS, fundamental cycles
    adj = {v: [] for v in range(n)}
    for j, (a, b) in enumerate(B.E): adj[a].append((b, j, 1)); adj[b].append((a, j, -1))
    parent = {0: None}; order = [0]; tree = set()
    for v in order:
        for (w, j, sgn) in adj[v]:
            if w not in parent: parent[w] = (v, j, sgn); order.append(w); tree.add(j)
    def path_to_root(v):
        c = [0]*m
        while parent[v] is not None:
            u, j, sgn = parent[v]; c[j] += sgn; v = u          # edge traversed from u to v
        return c
    cycles = []
    for j, (a, b) in enumerate(B.E):
        if j in tree: continue
        c = [0]*m; c[j] = 1
        pa, pb = path_to_root(a), path_to_root(b)
        for i in range(m): c[i] += pa[i] - pb[i]               # root -> a, a -> b, b -> root (p_* = root-to-vertex chains)
        cycles.append(c)
    C = sp.Matrix(cycles).T                                    # m x r
    assert (B.DEL*C).is_zero_matrix
    G = C.T*C; ac = C.T*B.a
    period_norm = sp.Rational((ac.T*G.inv()*ac)[0, 0])
    random.seed(1); ph = sp.Matrix([random.randint(-3, 3) for _ in range(n)])
    a2 = B.a + B.DEL.T*ph; gauge = sp.Rational((a2.T*B.PI*a2)[0, 0])
    return period_norm, gauge, C.shape[1], G.det()

if __name__ == "__main__":
    mp.mp.dps = 40
    n4, E4 = complete_graph(4); n5, E5 = complete_graph(5)
    bases = [Base(n4, E4, [0, 0, 0, 1, 0, 0], "K4 (1,0,0)"),
             Base(n4, E4, [0, 0, 0, 1, 1, 1], "K4 (1,1,1)"),
             Base(n4, E4, [0, 0, 0, 1, 2, 0], "K4 (1,2,0)"),
             Base(n5, E5, [1] + [0]*9, "K5 single edge (0,1)"),
             Base(n5, E5, [0, 0, 0, 0, 1, 1, 0, 1, 0, 0], "K5 triangle {1,2,3}")]
    print("===== (C) cyclic 1-cocycles of A = C(V) + C(E^+) =====")
    for nv, ne in ((4, 6), (5, 10)):
        print(f"  A = C^{nv+ne} (K{nv}): dim Z^1_lambda(A) = dim HC^1(A) = {cyclic_1_cocycles(nv+ne)}")
    print()
    for B in bases:
        print(f"===== {B.name}: n={B.n}, m={B.m}, q={B.q}, kappa_0={B.kappa0}; alpha={B.alpha} =====")
        print(f"  ||h_alpha||^2 = {B.hE},  L(Y_.) = {B.L}")
        print(f"  (N) Tr_s([D,alpha] D) = {brief_formula(B)} = 4 sum alpha(e) = {4*sum(B.alpha)}: {brief_formula(B) == 4*sum(B.alpha)}")
        lim, integral_form, qf = heat_formulas(B)
        print(f"  (H) lim_{{t->oo}} <alpha, e^(-t d^t d) alpha> = {lim};  ||alpha||^2 - int_0^oo <d alpha, e^(-s Delta_0) d alpha> ds = {integral_form};  both = ||h_alpha||^2: {lim == B.hE and integral_form == B.hE}")
        print(f"      <alpha, e^(-t d^t d) alpha> = {qf}")
        first, c = bloch_exact(B)
        cn = bloch_numeric(B)
        print(f"  (B) lambda_0(phi) = {first}*phi + ({c})*phi^2 + O(phi^3);  ||h||^2/n = {B.hE/B.n}: {c == B.hE/B.n};  numeric lambda_0/phi^2 (Richardson) = {mp.nstr(cn, 15)}")
        for tt in (5, 10, 20, 40):
            val = spectral_action(B, mp.mpf(tt))
            print(f"  (S) t={tt:>3}: (n/2)(1/t) d^2/dphi^2 Tr e^(-t D(e^(i phi)))|_0 = {mp.nstr(val, 15)}   (L = {B.L}, error {mp.nstr(val - mp.mpf(B.L.p)/B.L.q, 3)})")
        sig2, var = diffusion(B)
        print(f"  (V) sigma^2 (Poisson equation) = {sig2};  ||h||^2/m = {B.hE/B.m}: {sig2 == B.hE/B.m}")
        print("      Var(S_N) - N sigma^2 for N = 1,2,4,8,16,24: " + ", ".join(f"{sp.nsimplify(var[N]) - N*sig2}" for N in (1, 2, 4, 8, 16, 24)))
        pn, gauge, r, detG = albanese(B)
        print(f"  (A) rank H_1 = {r};  period formula a_c^t G^-1 a_c = {pn}: {pn == B.hE};  gauge-shifted energy = {gauge}: {gauge == B.hE};  det G = {detG} = kappa_0: {detG == B.kappa0};  kappa_0 L in Z: {(B.kappa0*B.L).is_integer}")
        if B.n == 5:
            D = B.pencil(t); f = sp.together(sp.simplify(D.det())); num, den = sp.fraction(f)
            T = sp.symbols('T'); g = sp.Poly(sp.expand(num.subs(t, 1+T)), T); dpoly = sp.Poly(sp.expand(den.subs(t, 1+T)), T)
            a2 = sp.Rational(g.coeff_monomial(T**2), dpoly.coeff_monomial(1))
            print(f"  (T) Theta(T) = det D(1+T): a_0 = {g.coeff_monomial(1)}, a_1 = {g.coeff_monomial(T)}, a_2 = {a2};  -kappa_0 ||h||^2 = {-B.kappa0*B.hE}: {a2 == -B.kappa0*B.hE}")
        print()
    print("done.")
