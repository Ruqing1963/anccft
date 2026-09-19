"""
Paper XVI battery (head_correspondence.py): the out-covering correspondence of the head map on the shadow
graphs of the Iwahori tower over K4 (q = 2), levels k = 1,2,3 (shadow graphs on 2 n_k = 24, 48, 96 vertices).

Conventions of Paper VIII: (CK1) s_e^* s_e = p_{r(e)}, (CK2) p_v = sum_{s(e)=v} s_e s_e^*;  K_0(C^*(E)) = coker(1 - A_E^t),
K_1 = ker(1 - A_E^t).  Shadow graph Y^sharp: vertices v, v'; arcs of Y; v -> v'; c = q-1 labelled return arcs v' -> v;
two loops at v'.  B = [[A, I],[cI, 2I]].  U = [[I,-cI],[0,I]], V = [[I,0],[-I,I]]:  U (I-B^t) V = Delta^t (+) (-I).

 (O)  eta^sharp : Y_{k+1}^sharp -> Y_k^sharp is a graph morphism, bijective on OUT-arcs at every vertex (out-covering),
      q-to-1 on the in-arcs of Y-vertices (not an in-covering);  tau^sharp is the mirror (in-covering, [VIII, Prop 2.3])
 (R)  reversal symmetry of the path tower: the path-reversal permutation rho_k conjugates A_k to A_k^t and satisfies
      rho_k tau_* = eta_* rho_{k+1}: everything proved for (tau, Delta) in [VIII] mirrors to (eta, Delta^t)
 (K)  the K_0-map of the head correspondence, H_k = diag(eta_*, eta_*):  H_k (I - B_{k+1}^t) = (I - B_k^t) H_k (it descends,
      unlike diag(tau_*,tau_*), [VIII, Thm 2.5(iv)]);  H_k 1 = q 1  ([P_k] = q [1]);  H_k (1,-1) = q (1,-1)  (K_1 map = q)
 (H)  Hecke composites:  H_k T^k = diag(A_k^t, A_k^t),  T^k H_k = diag(A_{k+1}^t, A_{k+1}^t)  (T^k = diag(tau^*,tau^*) =
      K_0 of the shadow transfer);  both act as multiplication by q on K_0 = coker(I - B^t): explicit integer preimage
      diag(A^t - q, A^t - q) = (I - B^t) V N,  N = [[-I, cI],[0, Delta^t]]
 (B)  Bowen--Franks side (opposite algebras): D_k = diag(tau_*,tau_*) descends on coker(I - B); composites
      D_k E^k = diag(A_k,A_k), E^k D_k = diag(A_{k+1},A_{k+1}) with E^k = diag(eta^*,eta^*); both act as q on coker(I - B):
      diag(A - q, A - q) = (I - B) U^t N',  N' = [[-I, I],[0, Delta]]
 (N)  the in-module transitions eta_* : Tor coker Delta_{k+1}^t -> Tor coker Delta_k^t : surjective, kernel = the 2-torsion,
      of order 2^{n_k(q-1)-1} (mirror of [VIII, Thm 2.8]); group structures
 (X)  UCT bookkeeping: Ext(K_0(O_{k+1}), K_1(O_k)) = Ext(Z (+) Jac_{k+1}^vee, Z) = Jac(Y_{k+1}); orders; the KK-lift
      ambiguity of [VIII, Rem 2.6]; Ext(K_0(O_k), K_1(O_k)) = Jac(Y_k) hosts the Eisenstein class [H_k] - q[1]
"""
import sys
sys.path.insert(0, r"C:\Users\LAPPIE\Desktop\KMS\算法化非交换类域论\Paper 10")
import numpy as np, sympy as sp
import pgl3_building as pb

q = 2; c = q - 1
Tp = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]

def nb_paths(k):
    if k == 0: return [(v,) for v in range(4)]
    ps = [(u, v) for u in range(4) for v in range(4) if Tp[u][v]]
    for _ in range(k-1):
        ps = [p+(w,) for p in ps for w in range(4) if Tp[p[-1]][w] and w != p[-2]]
    return ps

class Level:
    def __init__(self, k):
        self.k = k; self.V = nb_paths(k); self.idx = {p: i for i, p in enumerate(self.V)}; n = self.n = len(self.V)
        self.arcs = nb_paths(k+1)                       # arc p : p[:-1] -> p[1:]
        A = np.zeros((n, n), dtype=np.int64)
        for w in self.arcs: A[self.idx[w[:-1]], self.idx[w[1:]]] += 1
        self.A = A; self.Delta = q*np.eye(n, dtype=np.int64) - A
        I = np.eye(n, dtype=np.int64)
        self.B = np.block([[A, I], [c*I, 2*I]])
        self.U = np.block([[I, -c*I], [0*I, I]]); self.Vm = np.block([[I, 0*I], [-I, I]])
        # shadow graph as explicit arc list: vertices ('y', p) and ('s', p)
        self.svert = [('y', p) for p in self.V] + [('s', p) for p in self.V]
        self.sidx = {v: i for i, v in enumerate(self.svert)}
        arcs = []
        for w in self.arcs: arcs.append((('Y', w), ('y', w[:-1]), ('y', w[1:])))
        for p in self.V:
            arcs.append((('up', p), ('y', p), ('s', p)))
            for i in range(c): arcs.append((('ret', p, i), ('s', p), ('y', p)))
            for j in range(2): arcs.append((('loop', p, j), ('s', p), ('s', p)))
        self.sarcs = arcs                                 # (label, source, range)
        Bcheck = np.zeros((2*n, 2*n), dtype=np.int64)
        for (lab, s, r) in arcs: Bcheck[self.sidx[s], self.sidx[r]] += 1
        assert (Bcheck == self.B).all()

def maps(Lk, Lk1):
    """tau_*, eta_* : Z^{n_{k+1}} -> Z^{n_k} ; and the shadow maps eta^sharp, tau^sharp on vertices and arcs."""
    tau = np.zeros((Lk.n, Lk1.n), dtype=np.int64); eta = tau.copy()
    for w in Lk1.V: tau[Lk.idx[w[:-1]], Lk1.idx[w]] += 1; eta[Lk.idx[w[1:]], Lk1.idx[w]] += 1
    def sharp(which):
        cut = (lambda p: p[1:]) if which == 'eta' else (lambda p: p[:-1])
        vmap = {v: (v[0], cut(v[1])) for v in Lk1.svert}
        def amap(a):
            lab = a[0]
            if lab[0] == 'Y': return ('Y', cut(lab[1]))
            if lab[0] == 'up': return ('up', cut(lab[1]))
            return (lab[0], cut(lab[1]), lab[2])
        return vmap, amap
    return tau, eta, sharp('eta'), sharp('tau')

def covering_report(Lk, Lk1, vmap, amap):
    """Returns (is_morphism, out_bijective, in_bijective, max in-fibre multiplicity over Y-vertices)."""
    arcsF = {a[0]: a for a in Lk.sarcs}
    morphism = all(vmap[s] == arcsF[amap(a)][1] and vmap[r] == arcsF[amap(a)][2] for a in Lk1.sarcs for (_, s, r) in [a])
    outE = {}; inE = {}
    for a in Lk1.sarcs: outE.setdefault(a[1], []).append(a); inE.setdefault(a[2], []).append(a)
    outF = {}; inF = {}
    for a in Lk.sarcs: outF.setdefault(a[1], []).append(a[0]); inF.setdefault(a[2], []).append(a[0])
    out_ok = in_ok = True; maxfib = 1
    for w in Lk1.svert:
        # out-arcs: bijective onto out-arcs at the image  <=>  image multiset = image vertex's out-arc set
        if sorted(amap(a) for a in outE[w]) != sorted(outF[vmap[w]]): out_ok = False
        imgs_in = [amap(a) for a in inE[w]]
        if len(imgs_in) != len(set(imgs_in)) or sorted(set(imgs_in)) != sorted(inF[vmap[w]]): in_ok = False
        if w[0] == 'y':
            fib = {}
            for x in imgs_in: fib[x] = fib.get(x, 0) + 1
            maxfib = max(maxfib, max(fib.values()))
    return morphism, out_ok, in_ok, maxfib

def group(M):
    """Tor(coker M) for a singular Laplacian M: delete row 0 and column 0 (reduced Laplacian), whose
    cokernel is the sandpile group Jac and whose determinant is +/- kappa."""
    n = len(M); red = [[int(M[i][j]) for j in range(1, n)] for i in range(1, n)]
    inv, free, _ = pb.coker_structure(red)
    return inv, free

def torsion_order(inv):
    o = 1
    for d in inv: o *= d
    return o

if __name__ == "__main__":
    L = {k: Level(k) for k in (1, 2, 3, 4)}
    for k in (1, 2, 3):
        Lk, Lk1 = L[k], L[k+1]; n, n1 = Lk.n, Lk1.n
        tau, eta, (ve, ae), (vt, at) = maps(Lk, Lk1)
        print(f"===== transition k+1 -> k with k={k}: n_k={n}, n_(k+1)={n1}, shadow sizes {2*n}, {2*n1} =====")
        m, out_ok, in_ok, mf = covering_report(Lk, Lk1, ve, ae)
        print(f"  (O) eta^sharp: morphism {m}; bijective on OUT-arcs (out-covering) {out_ok}; bijective on in-arcs {in_ok} (max Y in-fibre {mf} = q: {mf == q})")
        m2, out_ok2, in_ok2, mf2 = covering_report(Lk, Lk1, vt, at)
        print(f"      tau^sharp: morphism {m2}; bijective on out-arcs {out_ok2}; bijective on IN-arcs (in-covering, [VIII, Prop 2.3]) {in_ok2}")
        # (R) reversal symmetry
        rho_k = np.zeros((n, n), dtype=np.int64); rho_k1 = np.zeros((n1, n1), dtype=np.int64)
        for p in Lk.V: rho_k[Lk.idx[p[::-1]], Lk.idx[p]] = 1
        for p in Lk1.V: rho_k1[Lk1.idx[p[::-1]], Lk1.idx[p]] = 1
        print(f"  (R) reversal: rho A_k rho^-1 = A_k^t: {(rho_k @ Lk.A @ rho_k.T == Lk.A.T).all()};  rho_k tau_* = eta_* rho_(k+1): {(rho_k @ tau == eta @ rho_k1).all()}")
        # (K) head correspondence K_0 map
        I2n = np.eye(2*n, dtype=np.int64); I2n1 = np.eye(2*n1, dtype=np.int64)
        H = np.block([[eta, 0*eta], [0*eta, eta]]); Tk = np.block([[tau.T, 0*tau.T], [0*tau.T, tau.T]])
        D = np.block([[tau, 0*tau], [0*tau, tau]]); Ek = np.block([[eta.T, 0*eta.T], [0*eta.T, eta.T]])
        one1 = np.ones(2*n1, dtype=np.int64); onek = np.ones(2*n, dtype=np.int64)
        pm1 = np.concatenate([np.ones(n1, dtype=np.int64), -np.ones(n1, dtype=np.int64)]); pmk = np.concatenate([np.ones(n, dtype=np.int64), -np.ones(n, dtype=np.int64)])
        print(f"  (K) H_k (I-B_(k+1)^t) = (I-B_k^t) H_k: {(H @ (I2n1 - Lk1.B.T) == (I2n - Lk.B.T) @ H).all()};"
              f"  diag(tau_*,tau_*) descends on the K_0 side: {(D @ (I2n1 - Lk1.B.T) == (I2n - Lk.B.T) @ D).all()} (false, [VIII, Thm 2.5(iv)]);"
              f"  H_k 1 = q 1: {(H @ one1 == q*onek).all()};  H_k (1,-1) = q (1,-1): {(H @ pm1 == q*pmk).all()}")
        # (H) Hecke composites and multiplication by q on K_0
        HT = H @ Tk; TH = Tk @ H
        okH = (HT == np.block([[Lk.A.T, 0*Lk.A], [0*Lk.A, Lk.A.T]])).all() and (TH == np.block([[Lk1.A.T, 0*Lk1.A], [0*Lk1.A, Lk1.A.T]])).all()
        In = np.eye(n, dtype=np.int64); Nmat = np.block([[-In, c*In], [0*In, Lk.Delta.T]])
        M = np.block([[Lk.A.T - q*In, 0*In], [0*In, Lk.A.T - q*In]])
        okq = ((I2n - Lk.B.T) @ Lk.Vm @ Nmat == M).all()
        print(f"  (H) H_k T^k = diag(A_k^t,A_k^t) and T^k H_k = diag(A_(k+1)^t,A_(k+1)^t): {okH};  diag(A_k^t - q, A_k^t - q) = (I-B_k^t) V N (so both Hecke composites act as q on K_0): {okq}")
        # (B) Bowen--Franks side
        okD = (D @ (I2n1 - Lk1.B) == (I2n - Lk.B) @ D).all()
        okDE = (D @ Ek == np.block([[Lk.A, 0*Lk.A], [0*Lk.A, Lk.A]])).all() and (Ek @ D == np.block([[Lk1.A, 0*Lk1.A], [0*Lk1.A, Lk1.A]])).all()
        Np = np.block([[-In, In], [0*In, Lk.Delta]]); Mp = np.block([[Lk.A - q*In, 0*In], [0*In, Lk.A - q*In]])
        okqB = ((I2n - Lk.B) @ Lk.U.T @ Np == Mp).all()
        okE = (Ek @ (I2n - Lk.B) == (I2n1 - Lk1.B) @ Ek).all()
        print(f"  (B) opposite algebras: D_k = diag(tau_*,tau_*) descends on coker(I-B): {okD}; E^k = diag(eta^*,eta^*) descends: {okE};"
              f" D_k E^k = diag(A_k,A_k), E^k D_k = diag(A_(k+1),A_(k+1)): {okDE};  diag(A_k - q, A_k - q) = (I-B_k) U^t N': {okqB}")
        # (N) in-module transitions
        invT1, _ = group(Lk1.Delta.T); invT0, _ = group(Lk.Delta.T)
        o1, o0 = torsion_order(invT1), torsion_order(invT0)
        # eta_* on torsion: kernel via 2-torsion: eta_* kills Tor[2]?  test: for x with 2x in im Delta^t ... use the mirror
        # exact statement: eta_* rho_{k+1} = rho_k tau_* and rho conjugates Delta_{k+1} to Delta_{k+1}^t, so ker(eta_* on Tor coker
        # Delta^t) = rho(ker tau_* on Tor coker Delta) = rho(Jac[2]) of order 2^{n_k(q-1)-1}
        r2 = pb.exact_rank([[int(x) % 2 for x in r] for r in Lk1.A])
        print(f"  (N) Tor coker Delta_(k+1)^t = {pb.group_str(invT1, 0)} (order {o1}), Tor coker Delta_k^t = {pb.group_str(invT0, 0)} (order {o0});"
              f" ratio {o1//o0} = 2^{{n_k(q-1)-1}} = 2^{n*(q-1)-1}: {o1 == o0 * 2**(n*(q-1)-1)};  rank_F2 A_(k+1) = {r2} = n_k: {r2 == n}"
              f" (kernel of eta_* on torsion = the 2-torsion, by reversal from [VIII, Thm 2.8])")
        # (X) UCT bookkeeping
        print(f"  (X) Ext(K_0(O_(k+1)), K_1(O_k)) = Ext(Z + Jac_(k+1)^vee, Z) = Jac(Y_(k+1)) of order {o1}: KK-lift ambiguity of [VIII, Rem 2.6];"
              f"  Ext(K_0(O_k), K_1(O_k)) = Jac(Y_k) of order {o0}: home of the Eisenstein class [H_k] - q[1]")
        print()
    print("done.")
