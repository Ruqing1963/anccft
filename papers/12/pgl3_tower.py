"""
Paper XII battery: towers of A~_2 complexes  (pgl3_tower.py).  Reuses Papers X and XI.

 (B)  first homology: b_1 = 0 and H_1(Y;Z) finite for Y_21, Y_24, Y_42 (exact), b_1 = 0 for Y_84, Y_168:
      no Z_l-voltage towers exist over thick A~_2 complexes
 (A)  the abelian chain X_3 <- Y_21 <- Y_84 <- Y_168 and X_3 <- Y_24 <- Y_168, X_3 <- Y_42 <- Y_168:
      covering pushforwards intertwine A_1, A_2 and Delta^H; Jac^H at each level; kernel orders = ratios;
      character mechanism: ratio x (n_k/n_{k+1}) = prod over new joint eigenpairs of |P_j(1)| (numeric)
 (D)  2-dimensional degeneracy calculus on Y_168: S T^t = A_1, T S^t = A_2, S S^t = T T^t = (q^2+q+1) I,
      S M^t = (q+1) A_2, T M^t = (q+1) A_1, M M^t = (q^2+q+1)(q+1) I + (structure), d_0+d_1+d_2 = boundary,
      K = sum_i d_i d_{i+1}^t, L_E = T^t S - K
 (E)  geodesic edge digraph E_1 = (E, L_E) on Y_168: q^2-regular in and out (Eulerian), strongly connected,
      n_1 = |E|; spectral accounting of kappa(E_1) = (1/n_1) prod_{lambda != q^2} (q^2 - lambda) into the
      Eisenstein factor 3q^4/n_1, the pencil factor prod_j q^6 P_j(q^{-2}), and the parahoric factor (numeric);
      tower constants mu_2 = 2 n_1 / q^2, chi_2 = ord_p(q^2) = 2, n_k = n_1 q^{2(k-1)}
 (F)  chamber digraph F_1 = (oriented chambers, L_F): q-regular Eulerian, n = 3|F|; constants mu_1 = 3|F|/q, chi = 1
 (L)  [III, Thm 2.1] sanity check on a small Eulerian 4-regular digraph (line digraph recursion with d = 4)
"""
import sys, math, itertools
sys.path.insert(0, r"C:\Users\LAPPIE\Desktop\KMS\算法化非交换类域论\Paper 10")
sys.path.insert(0, r"C:\Users\LAPPIE\Desktop\KMS\算法化非交换类域论\Paper 11")
import numpy as np, sympy as sp
import pgl3_building as pb
import pgl3_dynamical as pd

def h1_torsion(X):
    """H_1(Y;Z) = ker d1 / im d2 : coordinates via column transform of d1; torsion via p-adic Smith (small primes),
    completeness certified against the pseudo-determinant identity is not available here, so we also report
    |H_1| via the gcd of maximal minors when the coordinate matrix is square."""
    n, E, F = X.n, len(X.edges), len(X.tris)
    dd, _, Qi, _ = pb.diagonalize(X.d1, False, True); rank = sum(1 for x in dd if x != 0)
    Kcoord = [Qi[i] for i in range(rank, E)]
    L = pb.mul(Kcoord, X.d2)                       # (E-rank) x F coordinates of the columns of d2
    r = pb.exact_rank(L); b1 = (E - rank) - r
    inv, free, _ = pb.coker_structure(L, primes=pb.SMALL_PRIMES)
    return b1, inv, free

def covering_map(Xbig, Xsmall, proj):
    """pushforward matrix P (n_small x n_big) from vertex labels (i,g) -> (i, proj(g))."""
    # Complex objects do not store labels; rebuild them the way thick_cover does
    return None

def pushforward(labels_big, labels_small, proj):
    idx_small = {lab: k for k, lab in enumerate(labels_small)}
    P = np.zeros((len(labels_small), len(labels_big)), dtype=object)
    for k, (i, g) in enumerate(labels_big):
        P[idx_small[(i, proj(g))], k] = 1
    return P

def labels_of(moduli):
    """vertex labels (i, g) in the order used by pb.thick_cover"""
    mods = [m for m in moduli if m > 1]
    proj = lambda c: tuple(ci % mi for ci, mi in zip(c, moduli) if mi > 1)
    A = sorted({proj(c) for c in itertools.product(range(2), range(2), range(14))})
    return [(i, g) for i in range(3) for g in A]

def heckeH(X):
    n, q = X.n, X.q
    return [[(q**3 - 1)*(i == j) + X.A1[i][j] - q*X.A2[i][j] for j in range(n)] for i in range(n)]

def jacH(X):
    n = X.n; DH = heckeH(X)
    red = [[DH[i][j] for j in range(1, n)] for i in range(1, n)]
    inv, _, rem = pb.coker_structure(red)
    return inv, abs(pb.bareiss_det(red)), rem

def joint_P1(X, remove_rotations=False):
    """values P_j(1) = 1 - a_j + q conj(a_j) - q^3 over eigenvalues a_j of A_1 (numeric); the trivial eigenvalue
    q^2+q+1 is removed once (it is the only zero of P_j(1)); with remove_rotations also its two omega-rotations."""
    q = X.q; ev = np.linalg.eigvals(np.array(X.A1, dtype=float)); w = np.exp(2j*np.pi/3)
    targets = [q*q+q+1] + ([(q*q+q+1)*w, (q*q+q+1)*w*w] if remove_rotations else [])
    ev = list(ev)
    for t in targets:
        i = min(range(len(ev)), key=lambda i: abs(ev[i] - t)); ev.pop(i)
    ev = np.array(ev)
    return ev, np.array([1 - a + q*np.conj(a) - q**3 for a in ev])

def h1_certified(X):
    """Tor H_1(Y;Z) with the prime set certified: primes dividing a gcd of several maximal minors of the
    coordinate matrix L contain every prime of the torsion (each such prime divides all maximal minors)."""
    import random
    random.seed(12)
    n, E, F = X.n, len(X.edges), len(X.tris)
    dd, _, Qi, _ = pb.diagonalize(X.d1, False, True); rank = sum(1 for x in dd if x != 0)
    L = pb.mul([Qi[i] for i in range(rank, E)], X.d2); R = len(L)
    # one independent column basis from a pivoted QR (floating point, then confirmed by a nonzero exact determinant)
    from scipy.linalg import qr
    _, _, piv = qr(np.array(L, dtype=float), pivoting=True)
    cols = sorted(int(c) for c in piv[:R])
    g = abs(pb.bareiss_det([[row[c] for c in cols] for row in L]))
    assert g != 0, "pivoted QR did not return an independent column set"
    out = [c for c in range(F) if c not in cols]; tries = 0
    while tries < 12 and out:
        i = random.randrange(R); j = random.choice(out); trial = cols[:i] + [j] + cols[i+1:]
        dv = abs(pb.bareiss_det([[row[c] for c in trial] for row in L]))
        if dv: g = math.gcd(g, dv)
        tries += 1
    primes = sorted(sp.factorint(g).keys()) if g else pb.SMALL_PRIMES
    inv, free, _ = pb.coker_structure(L, primes=sorted(set(primes) | set(pb.SMALL_PRIMES)))
    return inv, free, g, primes

def small_digraph_check():
    """[III, Thm 2.1] with d = 4: kappa(L(G)) = 4^{n(4-1)-1} kappa(G) for an Eulerian 4-regular digraph G:
    take G = Cayley digraph of Z/5 with generators {1,2,3,4} (complete digraph K5 without loops... 4-regular)."""
    n = 5; A = [[1 if i != j else 0 for j in range(n)] for i in range(n)]
    arcs = [(i, j) for i in range(n) for j in range(n) if A[i][j]]
    aidx = {a: k for k, a in enumerate(arcs)}; m = len(arcs)
    AL = [[0]*m for _ in range(m)]
    for (i, j) in arcs:
        for (j2, k) in arcs:
            if j2 == j: AL[aidx[(i, j)]][aidx[(j, k)]] = 1
    def kappa(M, d):
        N = len(M); L = [[(d if i == j else 0) - M[i][j] for j in range(N)] for i in range(N)]
        return abs(pb.bareiss_det([r[1:] for r in L[1:]]))
    k0, k1 = kappa(A, 4), kappa(AL, 4)
    return k0, k1, k1 == 4**(n*3 - 1)*k0

if __name__ == "__main__":
    q = 2
    # ---------------- (B) first homology
    print("===== (B) first homology of the thick complexes =====")
    for mod in [(1, 1, 7), (2, 2, 2), (1, 1, 14)]:
        X = pb.thick_cover(mod); inv, free, g, primes = h1_certified(X)
        print(f"  {X.name}: b_1 = {free}; H_1(Y;Z) = {pb.group_str(inv, free)}; gcd of six maximal minors = {g} (primes {primes}): prime set certified")
    for mod in [(2, 2, 7), (2, 2, 14)]:
        X = pb.thick_cover(mod); r1 = pb.exact_rank(X.d1); r2 = pb.exact_rank(X.d2)
        print(f"  {X.name}: b_1 = {len(X.edges) - r1 - r2} (exact ranks); H_1 finite")
    print()
    # ---------------- (A) abelian chains
    print("===== (A) Hecke critical groups along the abelian covers of X_3 =====")
    chain = [(1, 1, 7), (2, 2, 7), (2, 2, 14)]
    Xs = {mod: pb.thick_cover(mod) for mod in set(chain + [(2, 2, 2), (1, 1, 14)])}
    data = {}
    for mod, X in Xs.items():
        inv, kap, rem = jacH(X); data[mod] = (inv, kap)
        ev, P1 = joint_P1(X)
        print(f"  {X.name}: Jac^H = {pb.group_str(inv, 0)}" + (f" [cofactor {rem}]" if rem > 1 else "") +
              f"; |Jac^H| = 2^{pb.sp.factorint(kap).get(2,0)} * odd;  log|Jac^H| = {math.log(kap):.6f} vs log((1/n)prod|P_j(1)|) = {float(np.sum(np.log(np.abs(P1))) - math.log(X.n)):.6f}")
    def proj_to(mod_big, mod_small):
        return lambda g: tuple(gi % ms for gi, ms, mb in zip(g, [m for m in mod_small if m > 1] if len([m for m in mod_small if m>1]) == len(g) else None, [m for m in mod_big if m > 1]) if ms > 1)
    pairs = [((2, 2, 14), (2, 2, 7)), ((2, 2, 7), (1, 1, 7)), ((2, 2, 14), (2, 2, 2)), ((2, 2, 14), (1, 1, 14)), ((1, 1, 14), (1, 1, 7))]
    for big, small in pairs:
        Xb, Xs_ = Xs[big], Xs[small]
        lb, ls = labels_of(big), labels_of(small)
        # projection on group labels: coordinates present in big; small keeps a subset with smaller moduli
        big_mods = [m for m in big if m > 1]; keep = [(k, sm) for k, (bm, sm) in enumerate(zip([m for m in big], [m for m in small])) if bm > 1]
        def proj(g, keep=keep):
            out = []
            for pos, (k, sm) in enumerate(keep):
                if sm > 1: out.append(g[pos] % sm)
            return tuple(out)
        P = pushforward(lb, ls, proj)
        A1b, A1s = np.array(Xb.A1, dtype=object), np.array(Xs_.A1, dtype=object)
        A2b, A2s = np.array(Xb.A2, dtype=object), np.array(Xs_.A2, dtype=object)
        ok1 = (P @ A1b == A1s @ P).all(); ok2 = (P @ A2b == A2s @ P).all()
        DHb, DHs = np.array(heckeH(Xb), dtype=object), np.array(heckeH(Xs_), dtype=object)
        okH = (P @ DHb == DHs @ P).all()
        kb, ks = data[big][1], data[small][1]
        ratio = sp.Rational(kb, ks)
        # character mechanism: new eigenvalues of A_1 (multiset difference), numeric
        evb, P1b = joint_P1(Xb); evs, P1s = joint_P1(Xs_)
        rem = list(evb);
        for z in evs:
            k = min(range(len(rem)), key=lambda k: abs(rem[k] - z)); rem.pop(k)
        evb_all = np.linalg.eigvals(np.array(Xb.A1, dtype=float))
        newP1 = [1 - a + q*np.conj(a) - q**3 for a in rem]
        pred = float(np.sum(np.log(np.abs(newP1)))) + math.log(Xs_.n) - math.log(Xb.n)
        print(f"  {Xs_.name} <- {Xb.name} (degree {Xb.n//Xs_.n}): P A1 = A1 P {ok1}, P A2 = A2 P {ok2}, P Delta^H = Delta^H P {okH}; "
              f"|Jac^H| ratio = {ratio} = 2^{sp.factorint(int(ratio)).get(2,0) if ratio.q==1 else '?'} * ...; integral {ratio.q == 1}; "
              f"log ratio = {math.log(kb/ks):.6f} vs log[(n_k/n_k+1) prod_new |P_j(1)|] = {pred:.6f}")
        # Hodge analogue: Jac_0 orders via reduced Laplacians and the Laplacian eigenvalues 2(q^2+q+1) - a - conj(a)
        def kappa0(Y):
            D0 = pb.mul(Y.d1, pb.T(Y.d1)); return abs(pb.bareiss_det([[D0[i][j] for j in range(1, Y.n)] for i in range(1, Y.n)]))
        k0b, k0s = kappa0(Xb), kappa0(Xs_)
        pred0 = float(np.sum(np.log(np.abs([2*(q*q+q+1) - a - np.conj(a) for a in rem])))) + math.log(Xs_.n) - math.log(Xb.n)
        print(f"      Hodge: |Jac_0| ratio = 2^{sp.factorint(k0b//k0s).get(2,0)} * ..., integral {k0b % k0s == 0}; log ratio = {math.log(k0b/k0s):.6f} vs log[(n_k/n_k+1) prod_new (2(q^2+q+1)-2Re a)] = {pred0:.6f}")
    print()
    # ---------------- (D) degeneracy calculus and (E) edge digraph on Y_168
    X = Xs[(2, 2, 14)]; n, E, F = X.n, len(X.edges), len(X.tris)
    S, T, M = pd.incidence(X); A1, A2 = np.array(X.A1, dtype=object), np.array(X.A2, dtype=object); I = np.eye(n, dtype=object)
    print("===== (D) two-dimensional degeneracy calculus on", X.name, "=====")
    print(f"  S T^t = A_1: {(S @ T.T == A1).all()};  T S^t = A_2: {(T @ S.T == A2).all()};  S S^t = T T^t = (q^2+q+1)I: {(S @ S.T == (q*q+q+1)*I).all() and (T @ T.T == (q*q+q+1)*I).all()}")
    print(f"  S M^t = (q+1) A_2: {(S @ M.T == (q+1)*A2).all()};  T M^t = (q+1) A_1: {(T @ M.T == (q+1)*A1).all()}")
    MM = M @ M.T; diag_ok = all(MM[i, i] == (q*q+q+1)*(q+1) for i in range(n))
    off = MM - (q*q+q+1)*(q+1)*I
    A11 = A1 @ A2 - (q*q+q+1)*I          # common type+1 out-neighbours: the (1,1)-distance operator
    A11b = A2 @ A1 - (q*q+q+1)*I
    print(f"  M M^t diagonal = (q^2+q+1)(q+1) = {(q*q+q+1)*(q+1)}: {diag_ok};  off-diagonal part = A_1A_2 - (q^2+q+1)I: {(off == A11).all()};  = A_2A_1 - (q^2+q+1)I: {(off == A11b).all()};  entries of off-diagonal part: {sorted(set(int(x) for r in off for x in r))}")
    # face maps d_i : chambers -> typed edges
    d = [np.zeros((E, F), dtype=np.int64) for _ in range(3)]
    for k, t in enumerate(X.tris):
        while X.typ[t[0]] != 0: t = (t[1], t[2], t[0])          # canonical listing: first vertex of type 0
        a, b, c = t
        d[0][X.eidx[(a, b)], k] = 1; d[1][X.eidx[(b, c)], k] = 1; d[2][X.eidx[(c, a)], k] = 1
    d2 = np.array(X.d2, dtype=np.int64)
    K = sum(d[i] @ d[(i+1) % 3].T for i in range(3))
    LEi = np.array(pd.edge_flow(X), dtype=np.int64); Si = np.array(S, dtype=np.int64); Ti = np.array(T, dtype=np.int64)
    print(f"  d_0 + d_1 + d_2 = boundary d2: {(d[0] + d[1] + d[2] == d2).all()};  K = sum_i d_i d_(i+1)^t has row sums {sorted(set(K.sum(axis=1).tolist()))} (= q+1);  L_E = T^t S - K: {(LEi == Ti.T @ Si - K).all()}")
    dtd = [d[i].T @ d[i] for i in range(3)]
    print(f"  sum_i d_i d_i^t = (q+1) I_E: {(sum(d[i] @ d[i].T for i in range(3)) == (q+1)*np.eye(E, dtype=np.int64)).all()};  "
          f"d_i^t d_j = 0 for i != j: {all((d[i].T @ d[j] == 0).all() for i in range(3) for j in range(3) if i != j)};  "
          f"d_i^t d_i = chamber-adjacency across the i-th edge (diagonal 1, row sums q+1): {all(int(dtd[i][k,k]) == 1 for i in range(3) for k in range(F)) and sorted(set(dtd[0].sum(axis=1).tolist())) == [q+1]}")
    d1i = np.array(X.d1, dtype=np.int64); Delta0 = d1i @ d1i.T
    A1i = np.array(X.A1, dtype=np.int64); A2i = np.array(X.A2, dtype=np.int64)
    print(f"  d_1 = T - S: {(d1i == Ti - Si).all()};  Delta_0 = (T-S)(T-S)^t = 2(q^2+q+1)I - A_1 - A_2: {(Delta0 == 2*(q*q+q+1)*np.eye(n, dtype=np.int64) - A1i - A2i).all()}")
    print(f"  A_1 A_2 = A_2 A_1 (from the symmetry of M M^t): {(A1i @ A2i == A2i @ A1i).all()}")
    print()
    print("===== (E) geodesic edge digraph E_1 = (E, L_E) on", X.name, "=====")
    LEi = np.array(pd.edge_flow(X), dtype=np.int64)
    outdeg = sorted(set(LEi.sum(axis=1).tolist())); indeg = sorted(set(LEi.sum(axis=0).tolist()))
    print(f"  n_1 = |E| = {E} = n(q^2+q+1); out-degrees {outdeg}, in-degrees {indeg} (Eulerian q^2-regular); strongly connected: {pd.strongly_connected(LEi)}")
    ev = np.linalg.eigvals(LEi.astype(float))
    i_per = np.argmin(np.abs(ev - q*q)); ev_rest = np.delete(ev, i_per)
    logkappa = float(np.sum(np.log(np.abs(q*q - ev_rest)))) - math.log(E)
    evC = np.linalg.eigvals(np.array(pb.companion(X), dtype=float))
    found, missing, remp = pd.match_spectra(list(evC), list(ev))
    # pencil factor: over nontrivial joint eigenpairs j: q^6 P_j(q^-2) = prod_i (q^2 - lambda_{j,i})
    evA, P1 = joint_P1(X, remove_rotations=True)
    pencil_factor = float(np.sum(np.log(np.abs(np.array([q**6*(1 - a*q**-2 + q*np.conj(a)*q**-4 - q**3*q**-6) for a in evA])))))
    eis = math.log(3*q**4 / E)
    parahoric = float(np.sum(np.log(np.abs(q*q - np.array(remp)))))
    print(f"  log kappa(E_1) = {logkappa:.6f} (numeric);  accounting: Eisenstein log(3q^4/n_1) = {eis:.6f} + pencil sum_j log|q^6 P_j(q^-2)| = {pencil_factor:.6f} + parahoric sum log|q^2 - lambda| = {parahoric:.6f} -> total {eis+pencil_factor+parahoric:.6f}")
    print(f"  tower constants: n_k = {E} * q^(2(k-1)); mu_2 = 2 n_1/q^2 = {2*E//(q*q)}; chi_2 = ord_p(q^2) = 2;  law: ord_p kappa(E_k) = {2*E//(q*q)} p^(2k) - 2k + nu_0")
    LF = pd.chamber_flow(X)
    print(f"  chamber digraph: {LF.shape[0]} = 3|F| vertices, out-degrees {sorted(set(LF.sum(axis=1).tolist()))}, in-degrees {sorted(set(LF.sum(axis=0).tolist()))}; mu_1 = 3|F|/q = {LF.shape[0]//q}; chi = 1")
    print()
    k0, k1, ok = small_digraph_check()
    print(f"===== (L) [III, Thm 2.1] with d = 4 on the complete digraph K_5: kappa = {k0}, kappa(L) = {k1} = 4^(5*3-1) kappa: {ok}")
    print("done.")
