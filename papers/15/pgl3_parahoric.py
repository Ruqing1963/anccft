"""
Paper XV battery (pgl3_parahoric.py): the parahoric factorization of the geodesic edge flow of an A~_2 complex,
the unitarity of the parahoric remainder, and the chamber companion.  Complexes of Paper X: Y_21, Y_24, Y_42,
Y_168 (q = 2, abelian covers of X_3) and the thin torus T_6 (q = 1).

Notation.  E = typed edges (x -> y), L_E((x,y),(y,z)) = 1 iff (y,z) in E and {x,y,z} is NOT a chamber (this is
T^t S - K of [XII, Thm 3.1]; it agrees with the adjacency definition of [XI, Def 3.1] exactly when Y is a flag
complex, and differs on Y_21, Y_24, Y_42).  S, T, M = tail / head / third-vertex incidences, Psi = [S^t T^t M^t],
W = (im Psi)^perp = ker S cap ker T cap ker M,  K = chamber continuation,  P(u) = I - A_1 u + q A_2 u^2 - q^3 u^3.

 (D)  exact integer identities:  L_E L_E^t = qI + (q^2-q) T^t T,  L_E^t L_E = qI + (q^2-q) S^t S,
      L_E Psi = Psi C',  L_E^t Psi = Psi C'',  K Psi = Psi Chat,  K^t Psi = Psi Chat'
 (K)  rank Psi (modular, certified by the explicit 6-dimensional type-constant kernel when it is 3n-6); dim W;
      K(u) = det(I - u C'|ker Psi) exactly
 (U)  unitarity of the remainder R = L_E|_W:  R R^t = R^t R = q I (numeric on an orthonormal basis of W),
      all eigenvalues of R of modulus sqrt(q)
 (E)  exact factorization det(I - u L_E) K(u) = det P(u) G(u) in Z[u] (small complexes: exact characteristic
      polynomials through the type-cyclic reduction; Y_168: numeric spectra), deg G = dim W, functional equation
      G(u) = (-u)^d det(R) G(1/(qu))
 (B)  chamber companion: pointed chambers, L_B ((x,y;z) -> (y,z;w), w != x), Phi = [sigma^* rho^* tau^*],
      L_B Phi = Phi B exactly with B = [[0,0,-I],[qI,K,K^t],[0,-I,0]];  the identity
      det((1+qu^3)I + uK + u^2 K^t) (1-u^3)^{-n} det P(u) = det(I - u L_E) det(I - u^2 L_E)  checked modulo
      two 31-bit primes at 3|E|+1 points (small complexes; a polynomial identity of degree <= 3|E| is thereby
      certified mod p) and at random points on Y_168;  det(I + u L_B) = det(I + uB) for q = 2 (points mod p);
      det(I + u L_B)/det(I + uB) = (1-u^3)^{chi-n} on T_6
 (V)  vertex Ramanujan check: all nontrivial roots of det P(u) have |u| = 1/q (numeric)
"""
import sys, time
sys.path.insert(0, r"C:\Users\LAPPIE\Desktop\KMS\算法化非交换类域论\Paper 10")
import numpy as np, sympy as sp
import pgl3_building as pb

PRIMES = [2147483629, 2147483587]
u = sp.symbols('u')

# ------------------------------------------------------------------ modular linear algebra (numpy int64, p < 2^31)
def modred(M, p):
    return np.array(M, dtype=np.int64) % p
def moddet(M, p):
    A = modred(M, p).copy(); n = A.shape[0]; det = 1
    for i in range(n):
        piv = np.nonzero(A[i:, i])[0]
        if len(piv) == 0: return 0
        k = i + piv[0]
        if k != i: A[[i, k]] = A[[k, i]]; det = (-det) % p
        det = (det * int(A[i, i])) % p
        inv = pow(int(A[i, i]), p - 2, p)
        A[i] = (A[i] * inv) % p
        if i + 1 < n:
            col = A[i+1:, i].copy()
            nz = np.nonzero(col)[0]
            if len(nz): A[i+1:][nz] = (A[i+1:][nz] - np.outer(col[nz], A[i])) % p
    return det
def modrank(M, p):
    A = modred(M, p).copy(); m, n = A.shape; r = 0
    for c in range(n):
        if r == m: break
        piv = np.nonzero(A[r:, c])[0]
        if len(piv) == 0: continue
        k = r + piv[0]
        if k != r: A[[r, k]] = A[[k, r]]
        inv = pow(int(A[r, c]), p - 2, p); A[r] = (A[r] * inv) % p
        col = A[:, c].copy(); col[r] = 0; nz = np.nonzero(col)[0]
        if len(nz): A[nz] = (A[nz] - np.outer(col[nz], A[r])) % p
        r += 1
    return r

# ------------------------------------------------------------------ the operators
class Data:
    def __init__(self, X):
        self.X = X; n, q = X.n, X.q; E = X.edges; eidx = X.eidx; nE = len(E)
        self.n, self.q, self.nE, self.nF = n, q, nE, len(X.tris)
        chambers = {frozenset(t) for t in X.tris}
        out = {}
        for (a, b) in E: out.setdefault(a, []).append(b)
        LE = np.zeros((nE, nE), dtype=np.int64)
        for k, (x, y) in enumerate(E):
            for z in out[y]:
                if frozenset((x, y, z)) not in chambers: LE[k, eidx[(y, z)]] = 1
        S = np.zeros((n, nE), dtype=np.int64); T = S.copy(); M = S.copy()
        for k, (a, b) in enumerate(E): S[a, k] = 1; T[b, k] = 1
        for (a, b, c) in X.tris:
            M[c, eidx[(a, b)]] += 1; M[a, eidx[(b, c)]] += 1; M[b, eidx[(c, a)]] += 1
        A1 = np.array(X.A1, dtype=np.int64); A2 = A1.T; I = np.eye(n, dtype=np.int64); Z = np.zeros((n, n), dtype=np.int64)
        self.LE, self.S, self.T, self.M, self.A1, self.A2 = LE, S, T, M, A1, A2
        self.K = T.T @ S - LE
        self.Psi = np.hstack([S.T, T.T, M.T])
        self.Cp = np.block([[Z, Z, -q*I], [q*q*I, A1, q*A2], [Z, -I, Z]])
        self.Cpp = np.block([[A2, q*q*I, q*A1], [Z, Z, -q*I], [-I, Z, Z]])
        self.Chat = np.block([[Z, Z, q*I], [(q+1)*I, Z, A2], [Z, I, Z]])
        self.Chatp = np.block([[Z, (q+1)*I, A1], [Z, Z, q*I], [I, Z, Z]])
        # pointed chambers
        pointed = []
        for (a, b, c) in X.tris: pointed += [((a, b), c), ((b, c), a), ((c, a), b)]
        pidx = {p: k for k, p in enumerate(pointed)}; third = {}
        for (e, z) in pointed: third.setdefault(e, []).append(z)
        LB = np.zeros((len(pointed), len(pointed)), dtype=np.int64)
        sig = np.zeros((len(pointed), nE), dtype=np.int64); rho = sig.copy(); tau = sig.copy()
        for k, ((x, y), z) in enumerate(pointed):
            for w in third[(y, z)]:
                if w != x: LB[k, pidx[((y, z), w)]] = 1
            sig[k, eidx[(x, y)]] = 1; rho[k, eidx[(y, z)]] = 1; tau[k, eidx[(z, x)]] = 1
        self.pointed, self.LB = pointed, LB
        self.Phi = np.hstack([sig, rho, tau])
        IE = np.eye(nE, dtype=np.int64); ZE = np.zeros((nE, nE), dtype=np.int64)
        self.B = np.block([[ZE, ZE, -IE], [q*IE, self.K, self.K.T], [ZE, -IE, ZE]])
        self.chi = n - nE + len(X.tris)
    def type_constant_kernel(self):
        """the six vectors (f,-f,0), ((q+1)f,0,-f), f = omega^{k tau} 1, as rational vectors: real/imag parts."""
        n = self.n; typ = self.X.typ; q = self.q
        ind = lambda i: [1 if typ[v] == i else 0 for v in range(n)]          # indicator of type i
        vecs = []
        for i in range(3):
            # (f, g, h) type-constant with f_j + g_{j+1} + (q+1) h_{j+2} = 0:  g = delta_i  =>  f = -delta_{i-1}
            vecs.append(sp.Matrix([-x for x in ind((i-1) % 3)] + ind(i) + [0]*n))
            # h = delta_i  =>  f = -(q+1) delta_{i-2}
            vecs.append(sp.Matrix([-(q+1)*x for x in ind((i-2) % 3)] + [0]*n + ind(i)))
        return sp.Matrix.hstack(*vecs)

def charpoly_det_I_minus(Mint):
    """det(I - v M) as a sympy Poly in v for an integer matrix M (sympy Berkowitz)."""
    v = sp.symbols('v')
    cp = sp.Matrix(Mint.tolist()).charpoly()          # det(x I - M)
    coeffs = cp.all_coeffs()                          # x^d + c1 x^{d-1} + ...
    d = len(coeffs) - 1
    return sp.Poly(sum(int(coeffs[k]) * v**k for k in range(d+1)), v), v

def cyclic_reduced_det(LE, typ_of_edge):
    """det(I - u L_E) = det(I - u^3 L_2 L_1 L_0) using the type grading of edges."""
    blocks = {i: [k for k, t in enumerate(typ_of_edge) if t == i] for i in range(3)}
    L = {i: LE[np.ix_(blocks[i], blocks[(i+1) % 3])] for i in range(3)}
    Pi = L[0] @ L[1] @ L[2]                            # E_0 -> E_1 -> E_2 -> E_0 (row-vector convention)
    return Pi

def poly_from_v(polyv, v, power):
    return sp.Poly(polyv.as_expr().subs(v, u**power), u)

def run(X, exact=True, points_full=True):
    t0 = time.time(); D = Data(X); n, q, nE = D.n, D.q, D.nE
    print(f"===== {X.name}: n={n}, |E|={nE}, |F|={D.nF}, q={q}, chi={D.chi}, |P|={len(D.pointed)} =====")
    LE, S, T, M, Psi, K = D.LE, D.S, D.T, D.M, D.Psi, D.K
    IE = np.eye(nE, dtype=np.int64)
    print("  (D) L_E = T^tS - K (chamber definition), row sums", sorted(set(LE.sum(1))), "|",
          "L_E L_E^t = qI+(q^2-q)T^tT:", (LE @ LE.T == q*IE + (q*q-q)*(T.T @ T)).all(),
          "L_E^t L_E = qI+(q^2-q)S^tS:", (LE.T @ LE == q*IE + (q*q-q)*(S.T @ S)).all())
    print("      L_E Psi = Psi C':", (LE @ Psi == Psi @ D.Cp).all(), " L_E^t Psi = Psi C'':", (LE.T @ Psi == Psi @ D.Cpp).all(),
          " K Psi = Psi Chat:", (K @ Psi == Psi @ D.Chat).all(), " K^t Psi = Psi Chat':", (K.T @ Psi == Psi @ D.Chatp).all())
    # adjacency-based L_E of [XI] for comparison
    adj = {(a, b) for (a, b) in X.edges} | {(b, a) for (a, b) in X.edges}
    out = {}
    for (a, b) in X.edges: out.setdefault(a, []).append(b)
    LEadj = np.zeros((nE, nE), dtype=np.int64)
    for k, (x, y) in enumerate(X.edges):
        for z in out[y]:
            if (x, z) not in adj: LEadj[k, X.eidx[(y, z)]] = 1
    print("      [XI, Def 3.1] adjacency version equals chamber version:", (LEadj == LE).all(), "(row sums", sorted(set(LEadj.sum(1))), ")")
    # (K) rank Psi
    rk = [modrank(Psi, p) for p in PRIMES]
    Kv = D.type_constant_kernel()
    Cp_s = sp.Matrix(D.Cp.tolist()); PsiS = sp.Matrix(Psi.tolist())
    kill = (PsiS * Kv).is_zero_matrix
    if rk[0] == 3*n - 6 and Kv.shape[1] == 6 and kill:
        # C' on the kernel: solve Kv X = C' Kv
        Xc = Kv.solve_least_squares(Cp_s * Kv); assert (Kv*Xc - Cp_s*Kv).is_zero_matrix
        Kpoly = sp.Poly((sp.eye(6) - u*Xc).det(), u)
        print(f"  (K) rank Psi = {rk[0]} = 3n-6 (mod {PRIMES[0]}, {PRIMES[1]}: {rk}); explicit 6-dim type-constant kernel Psi K = 0: {kill};"
              f" hence rank_Q Psi = 3n-6 exactly.  K(u) = det(I - u C'|ker) = {sp.factor(Kpoly.as_expr())}")
        dimW = nE - (3*n - 6); ND = True
    else:
        ns = PsiS.nullspace(); Kv = sp.Matrix.hstack(*ns)
        Xc = Kv.solve_least_squares(Cp_s * Kv); assert (Kv*Xc - Cp_s*Kv).is_zero_matrix
        Kpoly = sp.Poly((sp.eye(len(ns)) - u*Xc).det(), u)
        print(f"  (K) rank Psi = {rk[0]} (exact nullspace dim {len(ns)} != 6: nondegeneracy (ND) FAILS); K(u) = {sp.factor(Kpoly.as_expr())}")
        dimW = nE - rk[0]; ND = False
    print(f"      dim W = |E| - rank Psi = {dimW}  (4n+6 = {4*n+6})")
    # (U) unitarity on W (numeric)
    Uo, s, Vt = np.linalg.svd(Psi.astype(float)); r = (s > 1e-8).sum(); Nb = Uo[:, r:]
    R = Nb.T @ LE.astype(float) @ Nb
    evR = np.linalg.eigvals(R)
    print(f"  (U) numeric rank Psi = {r}; on W: R R^t = qI: {np.allclose(R @ R.T, q*np.eye(dimW), atol=1e-9)}, R^t R = qI: {np.allclose(R.T @ R, q*np.eye(dimW), atol=1e-9)};"
          f" L_E W in W: {np.allclose(Psi.T @ (LE @ Nb), 0)}; L_E^t W in W: {np.allclose(Psi.T @ (LE.T @ Nb), 0)};"
          f" |spec R| in [{np.abs(evR).min():.12f}, {np.abs(evR).max():.12f}], sqrt(q) = {np.sqrt(q):.12f}")
    # (V) vertex Ramanujan
    evA = np.linalg.eigvals(D.A1.astype(float)); N = q*q+q+1
    nontriv = [a for a in evA if abs(abs(a) - N) > 1e-6]
    mods = [abs(z) for a in nontriv for z in np.roots([1, -a, q*np.conj(a), -q**3])]
    print(f"  (V) nontrivial pencil roots (1/u): {len(mods)} values, modulus in [{min(mods):.10f}, {max(mods):.10f}] (Ramanujan iff all = q = {q})")
    # (E) exact polynomials
    typ_edge = [X.typ[a] for (a, b) in X.edges]
    if exact:
        Pi = cyclic_reduced_det(LE, typ_edge)
        polyv, v = charpoly_det_I_minus(Pi); detLE = poly_from_v(polyv, v, 3)            # det(I - u L_E)
        # det P(u) = det(I - u C) = det(I - u^3 P_0), P_0 return map of [XI]
        C = np.array(pb.companion(X), dtype=np.int64)
        pieces = {d: [i*n + w for i in range(3) for w in range(n) if (i - X.typ[w]) % 3 == d] for d in range(3)}
        Xd = {d: C[np.ix_(pieces[(d+1) % 3], pieces[d])] for d in range(3)}
        P0 = Xd[2] @ Xd[1] @ Xd[0]
        pv, v2 = charpoly_det_I_minus(P0); detP = poly_from_v(pv, v2, 3)
        # sanity: det(I - uC) at u=1 vs det of full companion
        lhs = detLE * Kpoly; Gq, Gr = sp.div(lhs, detP)
        print(f"  (E) exact: deg det(I-uL_E) = {detLE.degree()}, deg det P = {detP.degree()}, deg K = {Kpoly.degree()};"
              f" det P | det(I-uL_E) K: {Gr.is_zero};  G := quotient, deg G = {Gq.degree()} = dim W: {Gq.degree() == dimW}; G in Z[u]: {all(c.is_integer for c in Gq.all_coeffs())}")
        Gc = Gq.all_coeffs()
        print(f"      G(u) leading coeff {Gc[0]}, constant {Gc[-1]}, first coefficients {[int(c) for c in Gq.all_coeffs()[::-1][:6]]} ...")
        d = Gq.degree(); detR = sp.Integer(Gc[0]) * (-1)**d       # G(u) = det(I - uR): coefficient of u^d is (-1)^d det R
        fe = sp.expand(Gq.as_expr() - (-u)**d * detR * Gq.as_expr().subs(u, 1/(q*u)))
        print(f"      det R = {detR} (q^(d/2) = {sp.sqrt(sp.Integer(q)**d)});  functional equation G(u) = (-u)^d det(R) G(1/(qu)): {sp.simplify(fe) == 0}")
        # roots of G = 1/eigenvalues of R (the coefficient vector of G is far too ill-conditioned for np.roots)
        print(f"      roots of G = 1/spec R: |u| in [{1/np.abs(evR).max():.12f}, {1/np.abs(evR).min():.12f}], q^(-1/2) = {q**-0.5:.12f};"
              f" G(u) at u = 1/q^(1/2) e^(i t) for t = 0.7: |G| = {abs(sum(float(c)*(np.exp(0.7j)/np.sqrt(q))**k for k, c in enumerate(Gc[::-1]))):.4e} (nonzero: roots are not all at one point)")
        self_detLE = detLE
    else:
        detLE = detP = None
        print("  (E) exact characteristic polynomials skipped (size); spectra numeric: |spec R| above; det(I-uL_E) has the pencil roots minus the six trivial ones:")
        evLE = np.linalg.eigvals(LE.astype(float))
        pr = np.array([z for a in evA for z in np.roots([1, -a, q*np.conj(a), -q**3])])
        rem = list(evLE); miss = []
        for z in pr:
            k = min(range(len(rem)), key=lambda k: abs(rem[k]-z))
            if abs(rem[k]-z) < 1e-5: rem.pop(k)
            else: miss.append(z)
        print(f"      pencil roots missing from spec L_E: {np.round(sorted(miss, key=lambda z: (abs(z), np.angle(z))), 6)}; remainder count {len(rem)} = dim W: {len(rem) == dimW};"
              f" remainder moduli in [{min(abs(z) for z in rem):.10f}, {max(abs(z) for z in rem):.10f}]")
    # (B) chamber companion
    LB, Phi, B = D.LB, D.Phi, D.B
    print(f"  (B) L_B Phi = Phi B: {(LB @ Phi == Phi @ B).all()};  row sums L_B {sorted(set(LB.sum(1)))};  rank Phi (mod p) = {modrank(Phi, PRIMES[0])} of {Phi.shape}")
    def rhs_mod(uu, p):
        # (1-u^3)^n det(I - u L_E) det(I - u^2 L_E) * inv(det P(u))   mod p
        In = np.eye(n, dtype=np.int64)
        dP = moddet(In - uu*D.A1 + q*uu*uu*D.A2 - q**3*uu**3*In, p)
        d1 = moddet(IE - uu*LE, p); d2 = moddet(IE - uu*uu*LE, p)
        return (pow((1 - uu**3) % p, n, p) * d1 * d2 * pow(dP, p-2, p)) % p
    def lhsB_mod(uu, p):
        return moddet((1 + q*uu**3)*IE + uu*K + uu*uu*K.T, p)
    rng = np.random.default_rng(15)
    if points_full:
        pts = list(range(1, 3*nE + 2))
    else:
        pts = [int(x) for x in rng.integers(2, 10**6, size=6)]
    okB = True
    for p in PRIMES:
        for uu in pts:
            if lhsB_mod(uu, p) != rhs_mod(uu, p): okB = False; break
    print(f"      det((1+qu^3)I+uK+u^2K^t) = (1-u^3)^n det(I-uL_E) det(I-u^2L_E)/det P(u): {okB}  "
          f"({'all ' + str(len(pts)) + ' points 1..3|E|+1' if points_full else str(len(pts)) + ' random points'} mod two primes{'; certifies the polynomial identity mod p' if points_full else ''})")
    if len(D.pointed) <= 1000:
        ptsB = [int(x) for x in rng.integers(2, 10**6, size=8)]
        IP = np.eye(len(D.pointed), dtype=np.int64)
        if q == 2:
            ok2 = all(moddet(IP + uu*LB, p) == lhsB_mod(uu, p) for p in PRIMES[:1] for uu in ptsB)
            print(f"      q = 2: det(I + u L_B) = det(I + uB) = det((1+qu^3)I+uK+u^2K^t): {ok2} ({len(ptsB)} random points mod p)")
        else:
            p = PRIMES[0]
            ok1 = all((moddet(IP + uu*LB, p) * pow(lhsB_mod(uu, p), p-2, p)) % p == pow((1 - uu**3) % p, (D.chi - n) % (p-1), p) for uu in ptsB)
            print(f"      det(I + u L_B)/det(I + uB) = (1-u^3)^(chi-n) with chi-n = {D.chi-n}: {ok1} ({len(ptsB)} random points mod p)")
    print(f"      Kang--Li form: det(I+uL_B) = (1-u^3)^chi det(I-uL_E) det(I-u^2 L_E^t)/det P(u) is the previous two lines combined.")
    print(f"  [{time.time()-t0:.1f}s]\n")

if __name__ == "__main__":
    run(pb.thick_cover((1, 1, 7)))
    run(pb.thick_cover((2, 2, 2)))
    run(pb.thin_torus(6))
    run(pb.thick_cover((1, 1, 14)))
    run(pb.thick_cover((2, 2, 14)), exact=False, points_full=False)
    print("done.")
