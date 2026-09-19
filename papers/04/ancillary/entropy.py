"""Tree entropy vs the Ramanujan torsion bound of Corollary 3.10.

|Jac(Y)| = kappa(Y) = (1/n) prod_{i>=2} (d - lambda_i)   [Kirchhoff]
so   (1/n) log|Jac| = (1/n)[ sum_{i>=2} log(d-lambda_i) - log n ].
Compare with
  h_d  = int log(d-x) dmu_KM(x)      (Lyons tree entropy, limit for BS-convergence to T_d)
  B_d  = log(d + 2 sqrt(d-1)) = 2 log(1+sqrt q),  q=d-1   (Cor 3.10 upper bound)
"""
import numpy as np, math, random
from scipy.integrate import quad

# ---------- Kesten-McKay ----------
def km_density(x, d):
    R = 2*math.sqrt(d-1)
    if abs(x) >= R: return 0.0
    return d*math.sqrt(4*(d-1)-x*x) / (2*math.pi*(d*d - x*x))

def km_cdf(x, d):
    R = 2*math.sqrt(d-1)
    if x <= -R: return 0.0
    if x >=  R: return 1.0
    return quad(km_density, -R, x, args=(d,), limit=200)[0]

def h_tree_integral(d):
    R = 2*math.sqrt(d-1)
    return quad(lambda x: math.log(d-x)*km_density(x,d), -R, R, limit=400)[0]

def h_tree_closed(d):                      # Lyons:  (d-1)^{d-1} / (d^2-2d)^{d/2-1}
    return (d-1)*math.log(d-1) - (d/2-1)*math.log(d*d-2*d)

def bound(d):                              # Corollary 3.10
    return math.log(d + 2*math.sqrt(d-1))

# ---------- graph -> statistics ----------
def stats(A, d, name):
    n = A.shape[0]
    lam = np.linalg.eigvalsh(A.astype(float))[::-1]      # descending
    triv = [lam[0]]                                       # the Perron eigenvalue d
    bip = abs(lam[-1] + d) < 1e-8
    bulk = lam[1:-1] if bip else lam[1:]                  # drop +-d for the KM comparison
    lam2 = max(abs(lam[1]), abs(lam[-1]) if not bip else abs(lam[-2]))
    ram = lam2 <= 2*math.sqrt(d-1) + 1e-8
    logkappa = (np.log(d - lam[1:]).sum() - math.log(n)) / n
    xs = np.sort(bulk); m = len(xs)
    emp = (np.arange(1, m+1))/m
    ks = max(abs(km_cdf(x, d) - e) for x, e in zip(xs[::max(1, m//400)], emp[::max(1, m//400)]))
    return dict(name=name, n=n, d=d, lam2=lam2, ram=ram, ks=ks, ent=logkappa, bip=bip)

# ---------- LPS Ramanujan graphs ----------
def four_squares(p):
    sols = []
    r = int(math.isqrt(p))
    for a in range(1, r+1, 2):
        for b in range(-r, r+1, 2):
            for c in range(-r, r+1, 2):
                t = p - a*a - b*b - c*c
                if t < 0: continue
                s = int(math.isqrt(t))
                if s*s == t and s % 2 == 0:
                    for dd in ({s, -s} if s else {0}):
                        sols.append((a, b, c, dd))
    return sols

def norm_mat(M, l):
    M = tuple(x % l for x in M)
    for x in M:
        if x: 
            inv = pow(x, l-2, l)
            return tuple((y*inv) % l for y in M)
    return M

def lps(p, l):
    i = next(k for k in range(1, l) if (k*k) % l == l-1)
    gens = [norm_mat((a+b*i, c+d*i, -c+d*i, a-b*i), l) for (a,b,c,d) in four_squares(p)]
    assert len(gens) == p+1, (len(gens), p+1)
    idm = norm_mat((1,0,0,1), l)
    idx = {idm: 0}; order = [idm]; qu = [idm]
    def mul(X, Y):
        return norm_mat((X[0]*Y[0]+X[1]*Y[2], X[0]*Y[1]+X[1]*Y[3],
                         X[2]*Y[0]+X[3]*Y[2], X[2]*Y[1]+X[3]*Y[3]), l)
    while qu:
        X = qu.pop()
        for g in gens:
            Z = mul(g, X)
            if Z not in idx:
                idx[Z] = len(order); order.append(Z); qu.append(Z)
    n = len(order)
    A = np.zeros((n, n), dtype=np.int16)
    for X in order:
        for g in gens:
            A[idx[X], idx[mul(g, X)]] += 1
    return A

# ---------- random regular graphs ----------
def random_regular(n, d, seed=0):
    rng = random.Random(seed)
    for _ in range(400):
        stubs = [v for v in range(n) for _ in range(d)]
        rng.shuffle(stubs)
        E = set(); bad = False
        for k in range(0, len(stubs), 2):
            u, v = stubs[k], stubs[k+1]
            if u == v or (min(u,v), max(u,v)) in E: bad = True; break
            E.add((min(u,v), max(u,v)))
        if not bad:
            A = np.zeros((n,n), dtype=np.int16)
            for u,v in E: A[u,v] = A[v,u] = 1
            return A
    return None

def incidence_pg(qq):                     # bipartite incidence graph of PG(2,q), q prime
    pts = []
    for x in range(qq):
        for y in range(qq): pts.append((x,y,1))
    for x in range(qq): pts.append((x,1,0))
    pts.append((1,0,0))
    N = len(pts); A = np.zeros((2*N,2*N), dtype=np.int16)
    for a,P in enumerate(pts):
        for b,L in enumerate(pts):
            if sum(P[k]*L[k] for k in range(3)) % qq == 0:
                A[a, N+b] = A[N+b, a] = 1
    return A
