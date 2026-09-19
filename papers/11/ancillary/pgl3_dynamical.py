import os
"""
Paper XI battery: shift equivalence over Z for the PGL_3 companion, the graded return map, positive
presentations, and the geodesic edge flow  (pgl3_shift.py).  Reuses the complexes of Paper X.

 (G)  grading: deg(i,v) = i - type(v) mod 3 makes C^(3) a degree-one (period-3 cyclic) matrix; pieces of size n;
      return maps P_d = C^3|_piece d; spectrum of P_0 = cubes of one representative per omega-orbit; det P_0 = q^{3n}
 (Y)  sign data of P_0: negative entries, diagonal, right/left Perron vectors, eventual positivity, cycles in the
      negative digraph (obstruction to a strictly triangular nilpotent cofactor)
 (N)  nilpotent-cofactor criterion (Prop. 3.3): LP relaxation  Y >= -P_0,  Y P_0 <= 0,  Tr Y = 0  (necessary
      linear conditions); feasibility reported; if infeasible no nilpotent cofactor exists
 (E)  geodesic edge flow L_E on typed edges ((x,y)->(y,z), z not adjacent to x): row sums q^2, irreducibility,
      spectrum vs the pencil roots (which trivial orbits are missing), remainder moduli, log|det(I-L_E)|,
      eigenvalue 1 absent (K_1(O_{L_E}) = 0, K_0 finite)
 (F)  chamber flow L_F on oriented chambers: spectrum moduli
 (L)  cyclic lift of a shift equivalence (Lemma 2.2) checked on a synthetic example: given P_0 and B_0 = P_0
      (trivial SE) the lifted R,S satisfy the four equations with lag 3
"""
import sys, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "10"))
import numpy as np, sympy as sp
import pgl3_building as pb

W3 = np.exp(2j*np.pi/3)

def graded_pieces(X):
    n = X.n; typ = X.typ
    return {d: [i*n + v for i in range(3) for v in range(n) if (i - typ[v]) % 3 == d] for d in range(3)}

def graded_data(X):
    C = np.array(pb.companion(X), dtype=object); pieces = graded_pieces(X)
    deg1 = all(not (C[np.ix_(pieces[(d+k) % 3], pieces[d])] != 0).any() for d in range(3) for k in (0, 2))
    Xd = {d: C[np.ix_(pieces[(d+1) % 3], pieces[d])] for d in range(3)}          # X_d : piece d -> piece d+1
    P = {d: Xd[(d+2) % 3] @ Xd[(d+1) % 3] @ Xd[d] for d in range(3)}             # return maps
    C3 = C @ C @ C
    P_check = all((C3[np.ix_(pieces[d], pieces[d])] == P[d]).all() for d in range(3))
    return C, pieces, Xd, P, deg1 and P_check

def perron(Pf):
    w, V = np.linalg.eig(Pf); i = np.argmax(w.real); v = V[:, i].real; v = v/np.sign(v[np.argmax(np.abs(v))])
    return w[i].real, v

def edge_flow(X):
    E = X.edges; eidx = X.eidx
    adj = set()
    for (a, b) in E: adj.add((a, b)); adj.add((b, a))
    out = {}
    for (a, b) in E: out.setdefault(a, []).append(b)
    L = np.zeros((len(E), len(E)), dtype=np.int64)
    for k, (x, y) in enumerate(E):
        for z in out[y]:
            if (x, z) not in adj: L[k, eidx[(y, z)]] = 1
    return L

def chamber_flow(X):
    orient = []
    for (a, b, c) in X.tris: orient += [(a, b, c), (b, c, a), (c, a, b)]
    oidx = {o: k for k, o in enumerate(orient)}; on_edge = {}
    for (a, b, c) in orient: on_edge.setdefault((a, b), []).append(c)
    L = np.zeros((len(orient), len(orient)), dtype=np.int64)
    for k, (a, b, c) in enumerate(orient):
        for d in on_edge[(b, c)]:
            if d != a: L[k, oidx[(b, c, d)]] = 1
    return L

def strongly_connected(M):
    N = len(M); adj = [np.nonzero(M[i])[0].tolist() for i in range(N)]; radj = [np.nonzero(M[:, j])[0].tolist() for j in range(N)]
    def reach(a):
        seen = {0}; st = [0]
        while st:
            u = st.pop()
            for w in a[u]:
                if w not in seen: seen.add(w); st.append(w)
        return len(seen) == N
    return reach(adj) and reach(radj)

def match_spectra(target, pool, tol=1e-4):
    rem = list(pool); found = []; missing = []
    for z in target:
        if rem:
            k = min(range(len(rem)), key=lambda k: abs(rem[k]-z))
            if abs(rem[k]-z) < tol: found.append(z); rem.pop(k); continue
        missing.append(z)
    return found, missing, rem

def lp_nilpotent_relaxation(P0):
    """Feasibility of  Y >= -P_0 (entrywise),  Y P_0 <= 0 (entrywise),  Tr Y = 0  over the rationals (LP)."""
    from scipy.optimize import linprog
    n = P0.shape[0]; P = np.array(P0, dtype=float)
    N = n*n
    # variables y_{ij} flattened row-major; (YP)_{ik} = sum_j y_{ij} P_{jk}
    A_ub = []; b_ub = []
    for i in range(n):
        for k in range(n):
            row = np.zeros(N); row[i*n:(i+1)*n] = P[:, k]; A_ub.append(row); b_ub.append(0.0)
    A_eq = np.zeros((1, N)); A_eq[0, [i*n+i for i in range(n)]] = 1.0
    bounds = [(-P[i, j], None) for i in range(n) for j in range(n)]
    res = linprog(c=np.zeros(N), A_ub=np.array(A_ub), b_ub=np.array(b_ub), A_eq=A_eq, b_eq=[0.0], bounds=bounds, method="highs")
    return res.status, res.message

def cyclic_lift_check(Xd, P0):
    """Lemma 2.2 with the trivial SE  B_0 = P_0, R_0 = S_0' = I, l = 0:  R = diag(R_0, X_0 R_0, X_1 X_0 R_0),
    S = diag(S_2 X_1 X_0, S_2 X_1, S_2) with S_2 = S_0' X_2, B = cyc(I, I, B_0). Check CR = RB, SC = BS, RS = C^3, SR = B^3."""
    n = P0.shape[0]; I = np.eye(n, dtype=object); Z = np.zeros((n, n), dtype=object)
    X0, X1, X2 = Xd[0], Xd[1], Xd[2]
    R0 = I.copy(); S0p = I.copy(); B0 = P0
    R1 = X0 @ R0; R2 = X1 @ X0 @ R0
    S2 = S0p @ X2; S1 = S2 @ X1; S0 = S2 @ X1 @ X0
    def blockdiag(*Ms):
        N = sum(M.shape[0] for M in Ms); out = np.zeros((N, N), dtype=object); o = 0
        for M in Ms: out[o:o+M.shape[0], o:o+M.shape[1]] = M; o += M.shape[0]
        return out
    def cyc(Y0, Y1, Y2):   # Y_d : piece d -> piece d+1 ; matrix maps piece d block to piece d+1 block
        out = np.zeros((3*n, 3*n), dtype=object)
        out[n:2*n, 0:n] = Y0; out[2*n:3*n, n:2*n] = Y1; out[0:n, 2*n:3*n] = Y2
        return out
    Cg = cyc(X0, X1, X2); B = cyc(I, I, B0); R = blockdiag(R0, R1, R2); S = blockdiag(S0, S1, S2)
    e1 = (Cg @ R == R @ B).all(); e2 = (S @ Cg == B @ S).all()
    C3 = Cg @ Cg @ Cg; B3 = B @ B @ B
    e3 = (R @ S == C3).all(); e4 = (S @ R == B3).all()
    return e1, e2, e3, e4

def incidence(X):
    """S (start), T (end), M (third chamber vertex) : n x E integer matrices."""
    n, E = X.n, len(X.edges); S = np.zeros((n, E), dtype=object); Tm = np.zeros((n, E), dtype=object); M = np.zeros((n, E), dtype=object)
    for k, (a, b) in enumerate(X.edges): S[a, k] = 1; Tm[b, k] = 1
    for (a, b, c) in X.tris:
        M[c, X.eidx[(a, b)]] += 1; M[a, X.eidx[(b, c)]] += 1; M[b, X.eidx[(c, a)]] += 1
    return S, Tm, M

def intertwiner_check(X, LE):
    """Psi = [S^t | T^t | M^t] : Z^{3n} -> Z^E intertwines C' = [[0,0,-qI],[q^2 I, A_1, qA_2],[0,-I,0]] with L_E:
    L_E Psi = Psi C'.  Also C' = D C D^{-1} with D = [[0,0,qI],[I,0,0],[0,-I,0]], R = D, S = C D^{-1} integral (lag-1 SE)."""
    n, q = X.n, X.q; S, Tm, M = incidence(X); A1 = np.array(X.A1, dtype=object); A2 = np.array(X.A2, dtype=object)
    I = np.eye(n, dtype=object); Z = np.zeros((n, n), dtype=object)
    Psi = np.hstack([S.T, Tm.T, M.T])
    Cp = np.block([[Z, Z, -q*I], [q*q*I, A1, q*A2], [Z, -I, Z]])
    LEo = np.array(LE, dtype=object)
    inter = (LEo @ Psi == Psi @ Cp).all()
    C = np.array(pb.companion(X), dtype=object)
    D = np.block([[Z, Z, q*I], [I, Z, Z], [Z, -I, Z]])
    Dinv_q = np.block([[Z, q*I, Z], [Z, Z, -q*I], [I, Z, Z]])          # = q * D^{-1}
    Sm_q = C @ Dinv_q                                                   # q * (C D^{-1})
    S_int = (Sm_q % q == 0).all(); Sm = Sm_q // q
    e = ((Cp @ D == D @ C).all(), (Sm @ Cp == C @ Sm).all(), (D @ Sm == Cp).all(), (Sm @ D == C).all())
    rankPsi = np.linalg.matrix_rank(np.array(Psi, dtype=float))
    colsums = sorted(set(int(x) for x in M.sum(axis=0)))
    return inter, S_int, e, rankPsi, colsums

def farkas_certificate(P0):
    """Exact rational Farkas certificate for infeasibility of {Y >= -P0, Y P0 <= 0, Tr Y = 0}.
    With Y = -P0 + Z, Z >= 0: constraints  Z P0 <= P0^2 (entrywise),  Tr Z = Tr P0.
    Dual: find lambda >= 0 (n x n), mu with  lambda P0^t + mu I >= 0 entrywise  and  <lambda, P0^2> + mu Tr P0 < 0."""
    from scipy.optimize import linprog
    from fractions import Fraction
    n = P0.shape[0]; P = np.array(P0, dtype=float); P2 = P @ P; trP = float(np.trace(P))
    # variables: lambda (n*n, >=0), mu (free); objective: minimize <lambda,P2> + mu*trP ; normalization sum lambda + t = 1 via bound
    # constraint: (lambda P^t)_{ij} + mu*[i==j] >= 0  for all i,j   <=>  -(lambda P^t)_{ij} - mu*[i==j] <= 0
    N = n*n
    A_ub = []; b_ub = []
    for i in range(n):
        for j in range(n):
            row = np.zeros(N+1)
            for k in range(n): row[i*n+k] = -P[j, k]        # (lambda P^t)_{ij} = sum_k lambda_{ik} P_{jk}
            if i == j: row[N] = -1.0
            A_ub.append(row); b_ub.append(0.0)
    c = np.concatenate([P2.flatten(), [trP]])
    bounds = [(0, 1)]*N + [(-1e3, 1e3)]
    res = linprog(c=c, A_ub=np.array(A_ub), b_ub=np.array(b_ub), bounds=bounds, method="highs")
    if res.status != 0 or res.fun >= -1e-9: return None
    lam = res.x[:N]; mu = res.x[N]
    # rationalize and verify exactly
    for den in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 256, 512, 1024):
        L = [[Fraction(round(lam[i*n+j]*den), den) for j in range(n)] for i in range(n)]
        m = Fraction(round(mu*den), den)
        if any(x < 0 for r in L for x in r): continue
        P0i = [[int(x) for x in r] for r in P0.tolist()]
        ok = True
        for i in range(n):
            for j in range(n):
                val = sum(L[i][k]*P0i[j][k] for k in range(n)) + (m if i == j else 0)
                if val < 0: ok = False; break
            if not ok: break
        if not ok: continue
        P2i = (P0 @ P0).tolist()
        obj = sum(L[i][j]*int(P2i[i][j]) for i in range(n) for j in range(n)) + m*int(np.trace(P0))
        if obj < 0: return den, obj
    return "float-only"

if __name__ == "__main__":
    for mod in [(1, 1, 7), (2, 2, 2), (1, 1, 14), (2, 2, 14)]:
        X = pb.thick_cover(mod); n, q = X.n, X.q
        C, pieces, Xd, P, ok = graded_data(X)
        P0 = P[0]; P0f = np.array(P0, dtype=float)
        print(f"===== {X.name}: n={n} =====")
        print(f"  (G) C is degree one for deg = i - type: {ok}; pieces of size {[len(pieces[d]) for d in range(3)]}; P_d = C^3|piece d")
        evC = np.linalg.eigvals(np.array(C, dtype=float)); evP = np.linalg.eigvals(P0f)
        cubes = sorted(np.abs(evC**3)); cubes_rep = cubes[::3]
        print(f"      spec P_0 = cubes of C-spectrum, one per orbit (moduli): {np.allclose(sorted(np.abs(evP)), cubes_rep)}; det P_0 = 2^{sp.factorint(abs(pb.bareiss_det(P0.tolist()))).get(2,0)} (q^(3n) = 2^{3*n})")
        neg = [(i, j) for i in range(n) for j in range(n) if P0[i, j] < 0]
        rho, v = perron(P0f); rhoL, u = perron(P0f.T)
        Pk = P0.copy(); K0 = None
        for k in range(1, 31):
            if (Pk > 0).all(): K0 = k; break
            Pk = Pk @ P0
        adjn = {i: [j for (a, j) in neg if a == i] for i in range(n)}; color = {}
        def dfs(s):
            color[s] = 1
            for t in adjn[s]:
                if color.get(t) == 1: return True
                if t not in color and dfs(t): return True
            color[s] = 2; return False
        cyc_neg = any(dfs(s) for s in range(n) if s not in color)
        print(f"  (Y) P_0: entries in [{int(P0.min())}, {int(P0.max())}], {len(neg)} negative of {n*n}, diagonal values {sorted(set(int(P0[i,i]) for i in range(n)))}")
        print(f"      Perron root {rho:.6f} = q^6; right Perron vector > 0: {(v > 0).all()}; left Perron vector > 0: {(u > 0).all()}; "
              f"smallest k <= 30 with P_0^k > 0: {K0}; negative digraph has a cycle: {cyc_neg}")
        if n <= 50:
            st, msg = lp_nilpotent_relaxation(P0)
            cert = farkas_certificate(P0)
            print(f"  (N) LP relaxation of the nilpotent-cofactor conditions: status {st} ({'feasible' if st == 0 else 'infeasible' if st == 2 else msg}); "
                  f"exact rational Farkas certificate: {cert if cert is None or cert == 'float-only' else f'denominator {cert[0]}, dual value {cert[1]}'}")
        e = cyclic_lift_check(Xd, P0)
        print(f"  (L) cyclic lift of the trivial shift equivalence: CR=RB {e[0]}, SC=BS {e[1]}, RS=C^3 {e[2]}, SR=B^3 {e[3]}")
        if n >= 84:
            LE = edge_flow(X)
            rs = sorted(set(LE.sum(axis=1).tolist()))
            evE = np.linalg.eigvals(LE.astype(float)); nz = [z for z in evE if abs(z) > 1e-6]
            found, missing, rem = match_spectra(list(evC), nz)
            print(f"  (E) L_E: {LE.shape[0]} typed edges, row sums {rs} (q^2 = {q*q}), strongly connected: {strongly_connected(LE)}, rank-deficiency {LE.shape[0]-len(nz)}")
            print(f"      pencil roots in spec(L_E): {len(found)} of {len(evC)}; missing (rounded): {sorted(set(np.round(missing, 3).tolist()), key=lambda z: (abs(z), np.angle(z)))}")
            print(f"      remainder: {len(rem)} eigenvalues, moduli {sorted(set(np.round(np.abs(rem), 4).tolist()))} (q^(1/2) = {math.sqrt(q):.4f}); "
                  f"remainder arguments (deg, rounded): {sorted(set(np.round(np.degrees(np.angle(rem))).astype(int).tolist()))[:12]}...")
            one = min(abs(z-1) for z in evE); logdet = float(np.sum(np.log(np.abs(1-evE))))
            print(f"      eigenvalue 1 of L_E: distance {one:.4f} (so ker(I-L_E) = 0, K_1(O_LE) = 0); log|det(I-L_E)| = {logdet:.4f}  (|K_0(O_LE)| = exp of that)")
            # traces of L_E against traces of C^(3) (the specification's check): they differ by the parahoric part
            Lk = LE.copy(); trL = []
            for k in range(1, 11):
                trL.append(int(np.trace(Lk))); Lk = Lk @ LE
            trC, _, _ = pb.traces(X, kmax=10, kdirect=1)
            print(f"  (T) Tr L_E^k, k=1..10: {trL}")
            print(f"      Tr C^k,  k=1..10: {[trC[k] for k in range(1,11)]};  equal: {trL == [trC[k] for k in range(1,11)]} "
                  f"(they cannot be: L_E has the extra modulus-sqrt(q) spectrum)")
            inter, S_int, e4, rankPsi, colsums = intertwiner_check(X, LE)
            print(f"  (I) L_E Psi = Psi C' exactly: {inter}; column sums of M (chambers per edge) {colsums}; rank Psi = {rankPsi} = 3n-6 = {3*n-6}: {rankPsi == 3*n-6}")
            print(f"      C' = D C D^-1 with integral D and S = C D^-1 integral: {S_int}; lag-1 shift equivalence C'D = DC, SC' = CS, DS = C', SD = C: {e4}")
            LF = chamber_flow(X); evF = np.linalg.eigvals(LF.astype(float)); nzF = [z for z in evF if abs(z) > 1e-6]
            foundF, _, remF = match_spectra(list(evC), nzF)
            print(f"  (F) L_F: {LF.shape[0]} oriented chambers, row sums {sorted(set(LF.sum(axis=1).tolist()))}, strongly connected: {strongly_connected(LF)}; "
                  f"pencil roots in spec: {len(foundF)}; moduli of spectrum: {sorted(set(np.round(np.abs(nzF), 4).tolist()))} (q^(1/4) = {q**0.25:.4f})")
        print()
    print("done.")
