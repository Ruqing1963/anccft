"""
Paper IX battery: exceptional zeros, the graph L-invariant, and the twisted Bass determinant
(exceptional_zero.py).  Base graph K4 (q = 2, n = 4, m = 6, r = 3, kappa(Y_0) = 16).

Abelian voltage towers of [II, Table 2]:  (1,0,0) l=2;  (1,1,1) l=2;  (1,2,0) l=2;  (1,0,0) l=3.
 (Z)  Theta(T) = det D(1+T):  Theta(0) = Theta'(0) = 0,  a_2 = Theta''(0)/2 != 0  (double zero, order exactly 2)
 (H)  a_2 = - kappa(Y_0) * ||h_alpha||^2,  h_alpha = harmonic (cycle-space) projection of the voltage cochain
      = projection onto ker of the Hecke--Dirac operator of [II] on edges;  L(Y_.) := a_2/kappa(Y_0) = -||h_alpha||^2
 (H') single-chord voltage: -a_2 = kappa_0 (1 - R_eff(chord)) = number of spanning trees of Y_0 avoiding the chord
 (K)  critical partition function Z(1;t) = q(1-q^{-2}) tr D(t)^{-1}: double pole at t = 1 with leading
      coefficient n (q - q^{-1}) / L(Y_.)
 (M)  Newton data of h = g/T^2: mu, lambda_h;  regime lambda_h = 0  <=>  ord_l a_2 = mu
 (E)  defects eps_j = ord_l prod_{prim zeta} h(zeta-1) - (mu*phi(l^j) + lambda_h), j = 1..5: finite support
 (G)  growth law: ord_l kappa(Y_k) computed directly (reduced derived Laplacian determinant) and from the
      character product; = mu l^k + (lambda_h+1) k + [ord_l kappa_0 - mu + sum eps_j] for k >= k_0;
      matches [II, Table 2]
 (L)  nu_0 + ord_l L(Y_.) = (ord_l a_2 - mu) + sum_j eps_j;  in the lambda_h = 0 regime nu_0 = -ord_l L(Y_.)
 (B)  twisted Bass identity  det(I - u B Lambda) = (1-u^2)^{r-1} det(I - u A(t) + q u^2)  as a polynomial identity
 (F)  functional equation P(1/(qu),t) = (qu^2)^{-n} P(u,t);  P(1/q,t) = q^{-n} Theta(t-1);
      d_u P(1/q,1) = -n(q-1) q^{1-n} kappa_0;  (1/2) d_t^2 P(1/q,1) = q^{-n} a_2
 (S)  rank D(1) = n-1 while ord_T Theta = 2: the augmentation-prime length is 2 on a single invariant factor

Iwahori path tower (K4), levels k = 1,2,3:
 (I)  det(zI - A_{k+1}) = z^{n_{k+1}-n_k} det(zI - A_k);  the zeta determinant det(I - uA_k) is level-independent;
      u = 1/q is a simple zero with derivative -q^{2-n_k} n_k kappa(Y_k);  n_k kappa(Y_k) q^{-n_k} is a tower constant;
      kappa(Y_1) = (1-q^{-2})^{r-1} (q-1) n kappa_0 q^{n_1-1-n} / n_1  (level-0 formula)
"""
import sympy as sp
from fractions import Fraction

t, T, u, z = sp.symbols('t T u z')
Q_ = 2; n = 4
BASE = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]          # edges of K4, oriented u -> v
Tp = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]

def voltages(volt):
    a12,a13,a23 = volt
    return [0,0,0,a12,a13,a23]

# ---------------------------------------------------------------- exact integer determinant (Bareiss)
def bareiss(M):
    A = [list(map(int,r)) for r in M]; N = len(A); sign = 1; prev = 1
    for k in range(N-1):
        if A[k][k] == 0:
            sw = next((i for i in range(k+1,N) if A[i][k] != 0), None)
            if sw is None: return 0
            A[k],A[sw] = A[sw],A[k]; sign = -sign
        for i in range(k+1,N):
            for j in range(k+1,N):
                A[i][j] = (A[i][j]*A[k][k] - A[i][k]*A[k][j]) // prev
        prev = A[k][k]
    return sign*A[N-1][N-1]

def vl(x, l):
    x = abs(int(x))
    if x == 0: return None
    v = 0
    while x % l == 0: x //= l; v += 1
    return v

# ---------------------------------------------------------------- pencil, Theta, harmonic energy
def pencil(volt):
    A = sp.zeros(n,n)
    for (a,b),v in zip(BASE, voltages(volt)):
        A[a,b] += t**v; A[b,a] += t**(-v)
    return (Q_+1)*sp.eye(n) - A

def theta_data(volt):
    D = pencil(volt); f = sp.together(sp.simplify(D.det())); num, den = sp.fraction(f)
    num = sp.Poly(sp.expand(num), t); den = sp.Poly(sp.expand(den), t)
    g = sp.Poly(sp.expand(num.as_expr().subs(t, 1+T)), T)
    coeffs = [int(g.coeff_monomial(T**i)) for i in range(g.degree()+1)]
    return D, f, num, den, g, coeffs

# boundary  d: Z^E -> Z^V,  d(e) = e_v - e_u ;  Laplacian = d d^t ;  cycle-space projection  I - d^t Lap^+ d
DEL = sp.zeros(n, len(BASE))
for j,(a,b) in enumerate(BASE): DEL[b,j] += 1; DEL[a,j] -= 1
LAP = DEL*DEL.T
assert LAP == (Q_+1)*sp.eye(n) - sp.Matrix(Tp)
J = sp.ones(n,n)
LAPPLUS = (LAP + J/n).inv() - J/n
PCYC = sp.eye(len(BASE)) - DEL.T*LAPPLUS*DEL
KAPPA0 = abs(bareiss(LAP[1:,1:].tolist()))

def harmonic_energy(volt):
    a = sp.Matrix(voltages(volt))
    return sp.Rational((a.T*PCYC*a)[0,0])

# ---------------------------------------------------------------- derived graphs
def derived_laplacian(volt, N):
    size = n*N; L = [[0]*size for _ in range(size)]
    idx = lambda v,g: v*N + (g % N)
    for (a,b),v in zip(BASE, voltages(volt)):
        for g in range(N):
            i,j = idx(a,g), idx(b,g+v)
            L[i][j] -= 1; L[j][i] -= 1; L[i][i] += 1; L[j][j] += 1
    return L

def kappa_direct(volt, N):
    L = derived_laplacian(volt, N)
    return abs(bareiss([r[1:] for r in L[1:]]))

# ---------------------------------------------------------------- twisted Bass identity data
def edge_data(volt):
    vol = voltages(volt)
    E = [(a,b,v) for (a,b),v in zip(BASE,vol)] + [(b,a,-v) for (a,b),v in zip(BASE,vol)]   # 2m oriented
    m2 = len(E)
    S = sp.zeros(n,m2); Tm = sp.zeros(n,m2); Jr = sp.zeros(m2,m2); Lam = sp.zeros(m2,m2)
    for i,(a,b,v) in enumerate(E):
        S[a,i] = 1; Tm[b,i] = 1; Lam[i,i] = t**v
        Jr[i,(i+len(BASE)) % m2] = 1
    B = Tm.T*S - Jr
    return S, Tm, Jr, Lam, B

# ================================================================= abelian towers
TOWERS = [((1,0,0),2,(3,1,1)), ((1,1,1),2,(0,3,4)), ((1,2,0),2,(2,3,3)), ((1,0,0),3,(0,1,0))]
KMAX = {2:5, 3:3}
print(f"K4: n={n}, q={Q_}, r={len(BASE)-n+1}, kappa(Y_0)={KAPPA0}\n")
for volt, l, (mu_m, lam_m, nu_m) in TOWERS:
    D, f, num, den, g, c = theta_data(volt)
    a2 = c[2]
    print(f"===== tower {volt}, l={l} =====")
    print(f"  (Z) g(T) = (1+T)^{den.degree()} Theta(T) = {g.as_expr()}")
    print(f"      Theta(0)={c[0]}, Theta'(0)={c[1]}, a_2=Theta''(0)/2={a2}  -> double zero of order exactly 2: {c[0]==c[1]==0 and a2!=0}")
    hE = harmonic_energy(volt); Linv = sp.Rational(a2, KAPPA0)
    print(f"  (H) ||h_alpha||^2 = {hE};  -kappa_0*||h||^2 = {-KAPPA0*hE};  a_2 = -kappa_0 ||h_alpha||^2: {a2 == -KAPPA0*hE};  L(Y_.) = a_2/kappa_0 = {Linv} = -||h||^2: {Linv == -hE}")
    if volt == (1,0,0):
        Lm = LAP.copy(); Lm[1,2] += 1; Lm[2,1] += 1; Lm[1,1] -= 1; Lm[2,2] -= 1        # K4 minus the chord (1,2)
        trees_avoiding = abs(bareiss(Lm[1:,1:].tolist()))
        Reff = sp.Rational((DEL.T*LAPPLUS*DEL)[3,3])
        print(f"  (H') single chord: R_eff = {Reff}, kappa_0(1-R_eff) = {KAPPA0*(1-Reff)}, spanning trees of K4 avoiding the chord = {trees_avoiding};  = -a_2: {trees_avoiding == -a2}")
    Zc = Q_*(1 - sp.Rational(1,Q_**2))*(D.inv().trace())
    lead = sp.limit(sp.simplify(Zc*(t-1)**2), t, 1)
    print(f"  (K) leading coefficient of the double pole of Z(1;t) at t=1: {lead};  n(q-1/q)/L = {sp.Rational(n)*(Q_-sp.Rational(1,Q_))/Linv}: {lead == n*(Q_-sp.Rational(1,Q_))/Linv}")
    h = [c[i] for i in range(2, len(c))]
    ords = [vl(x,l) for x in h]
    mu = min(o for o in ords if o is not None); lam_h = next(i for i,o in enumerate(ords) if o == mu)
    print(f"  (M) h coefficients {h};  mu={mu}, lambda_h={lam_h};  regime lambda_h=0: {lam_h==0};  ord_l a_2 = {vl(a2,l)} (= mu iff regime)")
    hpoly = sp.Poly(sum(sp.Integer(h[i])*T**i for i in range(len(h))), T)
    eps = {}; Nj = {}
    for j in range(1, 6):
        Phi = sp.Poly(sp.expand(sum((1+T)**(i*l**(j-1)) for i in range(l))), T)     # Phi_{l^j}(1+T), monic
        R = int(sp.resultant(Phi, hpoly))
        Nj[j] = vl(R, l); eps[j] = Nj[j] - (mu*(l-1)*l**(j-1) + lam_h)
    print(f"  (E) eps_j = {eps}   (finite support: {all(eps[j]==0 for j in (4,5))})")
    ord_k0 = vl(KAPPA0, l) or 0
    nu0 = ord_k0 - mu + sum(eps.values())
    print(f"  (G) predicted law for k >= k_0: ord_l kappa(Y_k) = {mu}*{l}^k + {lam_h+1} k + {nu0};  [II, Table 2]: ({mu_m},{lam_m},{nu_m});  "
          f"match: {(mu,lam_h+1,nu0)==(mu_m,lam_m,nu_m)}")
    rows = []
    for k in range(0, KMAX[l]+1):
        direct = vl(kappa_direct(volt, l**k), l) or 0
        formula = ord_k0 - k + sum(2 + Nj[j] for j in range(1, k+1))
        law = mu*l**k + (lam_h+1)*k + nu0
        rows.append((k, direct, formula, law))
    print("      k : direct, character-formula, law  ->", rows)
    print(f"      direct == formula at all k: {all(r[1]==r[2] for r in rows)};  law exact for k>=2: {all(r[1]==r[3] for r in rows if r[0]>=2)};  "
          f"law exact for k>=0: {all(r[1]==r[3] for r in rows)}")
    ordL = (vl(a2,l) or 0) - ord_k0
    print(f"  (L) ord_l L(Y_.) = {ordL};  nu_0 + ord_l L = {nu0 + ordL} = (ord_l a_2 - mu) + sum eps = {(vl(a2,l)-mu) + sum(eps.values())}: "
          f"{nu0+ordL == (vl(a2,l)-mu)+sum(eps.values())};  regime => nu_0 = -ord_l L: {(not lam_h==0) or nu0 == -ordL}")
    # (B) twisted Bass identity
    S, Tm, Jr, Lam, B = edge_data(volt)
    m2 = B.shape[0]; r = len(BASE) - n + 1
    lhs = sp.expand((sp.eye(m2) - u*B*Lam).det())
    A_t = (Q_+1)*sp.eye(n) - D
    P = sp.expand((sp.eye(n) - u*A_t + Q_*u**2*sp.eye(n)).det())
    rhs = sp.expand((1-u**2)**(r-1)*P)
    print(f"  (B) det(I - uB Lambda) = (1-u^2)^(r-1) det(I - uA(t) + q u^2): {sp.simplify(lhs - rhs) == 0}")
    # (F) functional equation and derivatives at the critical point
    fe = sp.simplify(P.subs(u, 1/(Q_*u)) - (Q_*u**2)**(-n)*P)
    Pc = sp.simplify(P.subs(u, sp.Rational(1,Q_)))
    crit = sp.simplify(Pc - Q_**(-n)*f)
    dP_u = sp.simplify(sp.diff(P,u).subs({u:sp.Rational(1,Q_), t:1}))
    d2P_t = sp.simplify(sp.diff(P,t,2).subs({u:sp.Rational(1,Q_), t:1})/2)
    print(f"  (F) P(1/(qu),t) = (qu^2)^(-n) P(u,t): {fe==0};  P(1/q,t) = q^(-n) Theta(t-1): {crit==0};  "
          f"d_u P(1/q,1) = {dP_u} = -n(q-1)q^(1-n)kappa_0: {dP_u == -n*(Q_-1)*sp.Rational(1,Q_**(n-1))*KAPPA0};  "
          f"(1/2) d_t^2 P(1/q,1) = {d2P_t} = q^(-n) a_2: {d2P_t == sp.Rational(a2, Q_**n)}")
    print(f"  (S) rank D(1) = {D.subs(t,1).rank()} = n-1;  ord_T Theta = 2\n")

# ================================================================= Iwahori tower
def nb_paths(k):
    if k == 0: return [(v,) for v in range(4)]
    ps = [(a,b) for a in range(4) for b in range(4) if Tp[a][b]]
    for _ in range(k-1):
        ps = [p+(w,) for p in ps for w in range(4) if Tp[p[-1]][w] and w != p[-2]]
    return ps
def level(k):
    V = nb_paths(k); idx = {p:i for i,p in enumerate(V)}; N = len(V)
    A = [[0]*N for _ in range(N)]
    for w in nb_paths(k+1): A[idx[w[:-1]]][idx[w[1:]]] += 1
    return A

print("===== Iwahori path tower over K4 =====")
r = len(BASE) - n + 1
chi = {}; kap = {}; nk = {}
for k in (1,2,3):
    A = level(k); N = len(A); nk[k] = N
    M = sp.Matrix(A)
    chi[k] = sp.Poly(M.charpoly(z).as_expr(), z)
    Lk = [[(Q_ if i==j else 0) - A[i][j] for j in range(N)] for i in range(N)]
    kap[k] = abs(bareiss([row[1:] for row in Lk[1:]]))
    print(f"  k={k}: n_k={N}, kappa(Y_k) = {kap[k]} = 2^{vl(kap[k],2)}*{kap[k]>>vl(kap[k],2)}")
for k in (1,2):
    ok = sp.expand(chi[k+1].as_expr() - z**(nk[k+1]-nk[k])*chi[k].as_expr()) == 0
    print(f"  (I) det(zI - A_{k+1}) = z^(n_{k+1}-n_k) det(zI - A_{k}): {ok}")
const = {k: sp.Rational(nk[k]*kap[k], Q_**nk[k]) for k in (1,2,3)}
print(f"  (I) n_k kappa_k q^(-n_k) = {const}  constant: {len(set(const.values()))==1}")
for k in (1,2,3):
    zeta_inv = sp.expand(u**nk[k]*chi[k].as_expr().subs(z, 1/u))       # det(I - uA_k)
    val = zeta_inv.subs(u, sp.Rational(1,Q_)); der = sp.diff(zeta_inv,u).subs(u, sp.Rational(1,Q_))
    pred = -sp.Rational(nk[k]*kap[k], Q_**(nk[k]-2))
    print(f"  (I) level {k}: det(I-uA_k)|_(u=1/q) = {val};  d/du = {der} = -q^(2-n_k) n_k kappa_k: {der==pred}")
n1 = nk[1]
kappa1_pred = sp.Rational((Q_**2-1)**(r-1)*(Q_-1)*n*KAPPA0*Q_**(n1-1-n), Q_**(2*(r-1))*n1)
print(f"  (I) kappa(Y_1) from level 0: (1-q^-2)^(r-1)(q-1) n kappa_0 q^(n_1-1-n)/n_1 = {kappa1_pred} = {kap[1]}: {kappa1_pred==kap[1]}")
print("done.")
