"""
Paper XIII battery: the degeneracy algebra of the Iwahori tower, the augmentation-zero lattice, and the
level evaluations of the K_1 class  (cohn_localization.py).  K4 tower, q = 2, levels k = 1,2,3.

 (I)  incidence calculus [V, Prop 1.2]: tau_* eta^* = A_k, eta^* tau_* = A_{k+1}, tau_* tau^* = eta_* eta^* = qI;
      non-commutativity A_k A_k^t != A_k^t A_k; all four degeneracy maps preserve the augmentation-zero lattices
      (deg tau_* = deg, deg eta_* = deg, deg tau^* = q deg, deg eta^* = q deg)
 (S)  Sigma-inverting: Delta_k^0 := Delta_k restricted to ker(deg) is injective with finite cokernel; det Delta_k^0 is
      not a unit in Z_p (so the localization is genuinely needed) but Delta_k^0 is invertible over Q
 (E)  level evaluation: |det Delta_k^0| = n_k kappa_k = pseudo-determinant (product of nonzero eigenvalues), and
      det Delta_{k+1}^0 / det Delta_k^0 = q^{n_k(q-1)} exactly: pure arc growth, no Euler term
 (B)  boundary class at level k: 0 -> Z/n_k -> coker(Delta_k^0) -> Tor_k -> 0 exact; orders; whether it splits
 (A3) abelian tower (1,1,1), l = 2: pdet(Delta_k) = n kappa_0 prod_{chi != 1} det D(chi) exactly (resultants)
"""
import sys, math
sys.path.insert(0, r"C:\Users\LAPPIE\Desktop\KMS\算法化非交换类域论\Paper 10")
import numpy as np, sympy as sp
import pgl3_building as pb

q = 2
Tp = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]

def nb_paths(k):
    if k == 0: return [(v,) for v in range(4)]
    ps = [(u, v) for u in range(4) for v in range(4) if Tp[u][v]]
    for _ in range(k-1):
        ps = [p+(w,) for p in ps for w in range(4) if Tp[p[-1]][w] and w != p[-2]]
    return ps
def level(k):
    V = nb_paths(k); idx = {p: i for i, p in enumerate(V)}; n = len(V)
    A = [[0]*n for _ in range(n)]
    for w in nb_paths(k+1): A[idx[w[:-1]]][idx[w[1:]]] += 1
    return V, idx, A
def degeneracy(k):
    """tau_*, eta_* : Z^{n_{k+1}} -> Z^{n_k}  (tail / head truncation), pullbacks are transposes."""
    Vk, ik, _ = level(k); Vk1, ik1, _ = level(k+1)
    tau = [[0]*len(Vk1) for _ in range(len(Vk))]; eta = [[0]*len(Vk1) for _ in range(len(Vk))]
    for w in Vk1: tau[ik[w[:-1]]][ik1[w]] += 1; eta[ik[w[1:]]][ik1[w]] += 1
    return tau, eta

def M(x): return np.array(x, dtype=object)
def laplacian(A): n = len(A); return [[(q if i == j else 0) - A[i][j] for j in range(n)] for i in range(n)]
def reduced(Mx): n = len(Mx); return [[Mx[i][j] for j in range(1, n)] for i in range(1, n)]
def aug_basis(n):
    """columns e_i - e_0 (i = 1..n-1): a saturated basis of ker(deg)."""
    B = [[0]*(n-1) for _ in range(n)]
    for i in range(1, n): B[i][i-1] = 1; B[0][i-1] = -1
    return B
def restrict_to_aug(Mx):
    """matrix of Mx on ker(deg) in the basis e_i - e_0: since Mx maps Z^n into ker(deg) when its column sums
    vanish, the coordinates are just the rows 1..n-1 of Mx B."""
    n = len(Mx); B = aug_basis(n); MB = pb.mul(Mx, B)
    assert all(sum(MB[i][j] for i in range(n)) == 0 for j in range(n-1))
    return [MB[i] for i in range(1, n)]

if __name__ == "__main__":
    lev = {k: level(k) for k in (1, 2, 3, 4)}
    print("===== (I) incidence calculus and the augmentation-zero lattices =====")
    for k in (1, 2, 3):
        _, _, Ak = lev[k]; _, _, Ak1 = lev[k+1]; tau, eta = degeneracy(k)
        tauM, etaM, AkM, Ak1M = M(tau), M(eta), M(Ak), M(Ak1); nk, nk1 = len(Ak), len(Ak1)
        I = np.eye(nk, dtype=object)
        ok = ((tauM @ etaM.T == AkM).all(), (etaM.T @ tauM == Ak1M).all(), (tauM @ tauM.T == q*I).all(), (etaM @ etaM.T == q*I).all())
        noncomm = not (AkM @ AkM.T == AkM.T @ AkM).all()
        one1 = np.ones((nk1, 1), dtype=object); onek = np.ones((nk, 1), dtype=object)
        degs = ((onek.T @ tauM == one1.T).all(), (onek.T @ etaM == one1.T).all(), (one1.T @ tauM.T == q*onek.T).all(), (one1.T @ etaM.T == q*onek.T).all())
        print(f"  k={k}: n_k={nk}: tau eta^* = A_k {ok[0]}, eta^* tau = A_(k+1) {ok[1]}, tau tau^* = qI {ok[2]}, eta eta^* = qI {ok[3]}; "
              f"A_k A_k^t != A_k^t A_k: {noncomm}; deg-compatibility (deg tau_* = deg, deg eta_* = deg, deg tau^* = q deg, deg eta^* = q deg): {all(degs)}")
    print()
    print("===== (S),(E),(B) the restricted Laplacian on ker(deg) =====")
    dets = {}
    for k in (1, 2, 3):
        _, _, Ak = lev[k]; nk = len(Ak); D = laplacian(Ak); D0 = restrict_to_aug(D)
        det0 = pb.bareiss_det(D0); dets[k] = abs(det0)
        kappa = abs(pb.bareiss_det(reduced(D)))
        ev = np.linalg.eigvals(np.array(D, dtype=float)); pdet = float(np.prod([z for z in ev if abs(z) > 1e-8]).real)
        invZ, _, remZ = pb.coker_structure(D0)
        invT, _, remT = pb.coker_structure(reduced(D))
        v2 = pb.sp.factorint(dets[k]).get(2, 0)
        print(f"  k={k}: n_k={nk}, kappa_k={kappa} (2-adic order {sp.factorint(kappa).get(2,0)}); det Delta_k^0 = {det0} = ±n_k kappa_k: {dets[k]==nk*kappa}; "
              f"pdet (numeric) = {pdet:.1f}; 2-adic order of det = {v2} (not a unit in Z_2: {v2>0}; invertible over Q: {det0!=0})")
        # splitting test: is coker(Delta_k^0) isomorphic to Tor_k + Z/n_k ?  compare invariant factors
        split_inv, _ = pb.invariant_factors(list(invT) + [nk], 0)
        # factorization Delta_k = tau_* (tau^* - eta^*)
        tau, eta = degeneracy(k); tauM, etaM = M(tau), M(eta)
        fact = (tauM @ (tauM.T - etaM.T) == M(D)).all()
        print(f"        coker(Delta_k^0) = {pb.group_str(invZ, 0)};  Tor_k = Jac(Y_k) = {pb.group_str(invT, 0)};  |coker|/|Tor| = {dets[k]//kappa} = n_k: {dets[k]//kappa == nk}; "
              f"extension 0 -> Z/n_k -> coker -> Tor_k -> 0 splits: {sorted(invZ) == sorted(split_inv)};  Delta_k = tau_*(tau^* - eta^*): {fact}")
    for k in (1, 2):
        nk = len(lev[k][2])
        print(f"  ratio det Delta_(k+1)^0 / det Delta_k^0 = {sp.Rational(dets[k+1], dets[k])} = q^(n_k(q-1)) = 2^{nk*(q-1)}: {dets[k+1] == dets[k]*q**(nk*(q-1))}")
    print(f"  law: ord_2(n_k kappa_k) = {[sp.factorint(dets[k]).get(2,0) for k in (1,2,3)]} vs n_k + c: n_k = {[len(lev[k][2]) for k in (1,2,3)]} -> constant c = {[sp.factorint(dets[k]).get(2,0) - len(lev[k][2]) for k in (1,2,3)]}")
    print()
    print("===== (A3) abelian voltage tower (1,1,1), l = 2: pdet(Delta_k) = n kappa_0 prod_{chi != 1} det D(chi) =====")
    t = sp.symbols('t')
    E = [((0,1),0),((0,2),0),((0,3),0),((1,2),1),((1,3),1),((2,3),1)]
    Dt = sp.zeros(4, 4)
    for (u, v), a in E: Dt[u, v] += t**a; Dt[v, u] += t**(-a)
    Dt = 3*sp.eye(4) - Dt
    f = sp.together(sp.simplify(Dt.det())); num, den = sp.fraction(f); num = sp.Poly(sp.expand(num), t); den = sp.Poly(sp.expand(den), t)
    def derived(N):
        L = [[0]*(4*N) for _ in range(4*N)]; idx = lambda v, g: v*N + (g % N)
        for (u, v), a in E:
            for g in range(N):
                i, j = idx(u, g), idx(v, g+a); L[i][j] -= 1; L[j][i] -= 1; L[i][i] += 1; L[j][j] += 1
        return L
    kappa0 = abs(pb.bareiss_det(reduced(derived(1))))
    for k in (1, 2, 3):
        N = 2**k; Lk = derived(N); nk = 4*N
        kap = abs(pb.bareiss_det(reduced(Lk))); pdet_exact = nk*kap
        # prod over nontrivial characters: prod_{zeta^N = 1, zeta != 1} det D(zeta) = prod_{j=1}^{k} Res(Phi_{2^j}, num)/Res(Phi_{2^j}, den)
        prod = sp.Integer(1)
        for j in range(1, k+1):
            Phi = sp.Poly(sum(t**(i*2**(j-1)) for i in range(2)), t)
            prod *= sp.Rational(sp.resultant(Phi, num), sp.resultant(Phi, den))
        rhs = 4*kappa0*prod
        print(f"  k={k}: n_k={nk}, kappa_k={kap}; pdet = n_k kappa_k = {pdet_exact}; n kappa_0 prod det D(chi) = {rhs}; equal up to sign: {abs(rhs) == pdet_exact}")
    print("done.")
