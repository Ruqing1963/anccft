"""
Paper VIII battery: the tower-level Hecke--Kirchberg limit  (tower_kirchberg.py)

Iwahori path tower over K4 (q = 2), levels k = 1,2,3 (n_k = 12, 24, 48), level-k shadow
matrix  B_k = [[A_k, I],[cI, 2I]]  of size 2 n_k = 24, 48, 96.  Exact integer arithmetic.

 (S)  normalization: with c = q-1,  U (I-B_k^t) V = Delta_k^t (+) (-I)  and
      V^t (I-B_k) U^t = Delta_k (+) (-I),  U = [[I,-cI],[0,I]], V = [[I,0],[-I,I]].
 (K)  SNF(I-B_k^t): exactly one zero (K_1 = Z); torsion = Jac(Y_k); 2-adic orders 7, 18, 41;
      the prime-to-q part is the same at every level.
 (W)  the specification's c = q at a q-regular level: no zero invariant (K_1 = 0), finite K_0
      unrelated to Jac(Y_k)  -- the correction of Remark 1.3.
 (C)  tau^sharp : Y_{k+1}^sharp -> Y_k^sharp is bijective on IN-arcs at every vertex and is
      NOT bijective on out-arcs  -- the in-covering hypothesis of Theorem 2.2.
 (U)  up-transfer  T^k = diag(tau^*,tau^*):  T^k (I-B_k^t) = (I-B_{k+1}^t) T^k;  T^k 1 = 1;
      T^k (1,-1) = (1,-1)  (K_1 map = id).
 (D)  down-transfer T_k = diag(tau_*,tau_*):  T_k (I-B_{k+1}) = (I-B_k) T_k  (Bowen--Franks side);
      T_k (I-B_{k+1}^t) != (I-B_k^t) T_k  (it does NOT descend on the K_0 side);
      T_k (1,-c1) = q (1,-c1)  (Eisenstein index q on ker(I-B)).
 (N)  induced map coker(I-B_{k+1}) -> coker(I-B_k) on torsion: surjective, |ker| = 2^11, 2^23,
      and identical to tau_* : Tor coker Delta_{k+1} -> Tor coker Delta_k.
 (P)  Pontryagin adjointness  <tau_* x, y>_k = <x, tau^* y>_{k+1}  in Q/Z on random torsion classes.
 (Q)  kernel = p-torsion (Theorem 2.8, q = p = 2):  rank_{F_2} A_{k+1} = n_k;  the 2-rank of
      Jac(Y_{k+1}) is n_k(q-1)-1 = 11, 23;  tau_* kills Jac(Y_{k+1})[2];  hence ker tau_* = Jac(Y_{k+1})[2].
 (G)  Kirchberg witnesses at each level (strong connectivity, double loop at a shadow vertex).
 (B)  plain tower is blind: SNF(I-A_k^t) is the same at k = 1,2,3  (= Z^r (+) Z/(r-1), r = 3).
 (L)  limit data: deg(tau^* y) = q deg(y)  (free part of K_0 -> Z[1/q]);  K_1 constant.
"""
import itertools, random
from fractions import Fraction
import sympy as sp

random.seed(8)
Q_ = 2
Tp = [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]

# ---------------------------------------------------------------- tower
def nb_paths(k):
    if k == 0: return [(v,) for v in range(4)]
    ps = [(u,v) for u in range(4) for v in range(4) if Tp[u][v]]
    for _ in range(k-1):
        ps = [p+(w,) for p in ps for w in range(4) if Tp[p[-1]][w] and w != p[-2]]
    return ps

def level(k):
    V = nb_paths(k); idx = {p:i for i,p in enumerate(V)}; n = len(V)
    A = [[0]*n for _ in range(n)]
    for w in nb_paths(k+1): A[idx[w[:-1]]][idx[w[1:]]] += 1
    D = [[(Q_ if i==j else 0) - A[i][j] for j in range(n)] for i in range(n)]
    return V, idx, A, D

def tail(k):                       # tau_* : Z^{n_{k+1}} -> Z^{n_k}
    Vk, ik, _, _ = level(k); Vk1, ik1, _, _ = level(k+1)
    M = [[0]*len(Vk1) for _ in range(len(Vk))]
    for w in Vk1: M[ik[w[:-1]]][ik1[w]] += 1
    return M

# ---------------------------------------------------------------- integer linear algebra
def eye(n): return [[int(i==j) for j in range(n)] for i in range(n)]
def zeros(m,n): return [[0]*n for _ in range(m)]
def T(M): return [list(r) for r in zip(*M)]
def mul(A,B):
    Bt = T(B)
    return [[sum(a*b for a,b in zip(r,c)) for c in Bt] for r in A]
def add(A,B): return [[a+b for a,b in zip(r,s)] for r,s in zip(A,B)]
def sub(A,B): return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]
def scal(c,A): return [[c*a for a in r] for r in A]
def block(TL,TR,BL,BR): return [r+s for r,s in zip(TL,TR)] + [r+s for r,s in zip(BL,BR)]
def blockdiag(X,Y):
    m,n = len(X),len(X[0]); p,q = len(Y),len(Y[0])
    return block(X, zeros(m,q), zeros(p,n), Y)
def eq(A,B): return A == B
def vec(v): return [[x] for x in v]

def shadow(A, c):
    n = len(A); I = eye(n)
    return block(A, I, scal(c,I), scal(2,I))

def diagonalize(M, transforms=False):
    """P M Q = D diagonal (d_i not necessarily divisibility-ordered). Returns (d, P, Pinv)."""
    A = [list(r) for r in M]; m = len(A); n = len(A[0])
    P = eye(m) if transforms else None; Pinv = eye(m) if transforms else None
    def rowop(r,i,q):            # row_r -= q*row_i  ;  Pinv: col_i += q*col_r
        A[r] = [a-q*b for a,b in zip(A[r],A[i])]
        if transforms:
            P[r] = [a-q*b for a,b in zip(P[r],P[i])]
            for t in range(m): Pinv[t][i] += q*Pinv[t][r]
    def rowswap(r,i):
        A[r],A[i] = A[i],A[r]
        if transforms:
            P[r],P[i] = P[i],P[r]
            for t in range(m): Pinv[t][r],Pinv[t][i] = Pinv[t][i],Pinv[t][r]
    def colop(c,j,q):
        for t in range(m): A[t][c] -= q*A[t][j]
    def colswap(c,j):
        for t in range(m): A[t][c],A[t][j] = A[t][j],A[t][c]
    t = 0
    while t < min(m,n):
        piv = None
        for r in range(t,m):
            for c in range(t,n):
                if A[r][c] != 0 and (piv is None or abs(A[r][c]) < abs(A[piv[0]][piv[1]])): piv = (r,c)
        if piv is None: break
        rowswap(t,piv[0]); colswap(t,piv[1])
        while True:
            changed = False
            for r in range(t+1,m):
                if A[r][t] != 0:
                    q = A[r][t] // A[t][t]; rowop(r,t,q)
                    if A[r][t] != 0: rowswap(t,r); changed = True; break
            if changed: continue
            for c in range(t+1,n):
                if A[t][c] != 0:
                    q = A[t][c] // A[t][t]; colop(c,t,q)
                    if A[t][c] != 0: colswap(t,c); changed = True; break
            if not changed: break
        t += 1
    d = [abs(A[i][i]) for i in range(min(m,n))]
    return d, P, Pinv

def invariant_factors(d, ncols):
    """Smith invariants (>1) and free rank from a diagonal form."""
    primes = {}
    for x in d:
        if x in (0,1): continue
        f = sp.factorint(x)
        for p,e in f.items(): primes.setdefault(p,[]).append(e)
    k = max((len(v) for v in primes.values()), default=0)
    inv = [1]*k
    for p,es in primes.items():
        es = sorted(es, reverse=True)
        for i,e in enumerate(es): inv[k-1-i] *= p**e   # largest invariant gets largest power
    inv = sorted(x for x in inv if x > 1)
    free = sum(1 for x in d if x == 0) + (ncols - len(d))
    return inv, free

def group_str(inv, free):
    parts = ([f"Z^{free}"] if free > 1 else ["Z"] if free == 1 else [])
    from collections import Counter
    for v,mult in sorted(Counter(inv).items()):
        parts.append(f"(Z/{v})^{mult}" if mult > 1 else f"Z/{v}")
    return " + ".join(parts) if parts else "0"

def v2(x):
    x = abs(int(x)); v = 0
    while x and x % 2 == 0: x //= 2; v += 1
    return v

def torsion_order(d):
    o = 1
    for x in d:
        if x not in (0,1): o *= x
    return o

def coker_map_on_torsion(Mk, Mk1, F):
    """F : Z^{rows Mk1} -> Z^{rows Mk} with F Mk1 = Mk F'.  Returns (|Tor_k1|, |Tor_k|, |im|, |ker|) of the
    induced map on torsion subgroups of the cokernels."""
    dk, Pk, _ = diagonalize(Mk, True); dk1, _, Pk1inv = diagonalize(Mk1, True)
    G = mul(mul(Pk, F), Pk1inv)
    tk = [i for i,x in enumerate(dk) if x not in (0,1)]; tk1 = [j for j,x in enumerate(dk1) if x not in (0,1)]
    fk = [i for i,x in enumerate(dk) if x == 0]
    # torsion generators must land in torsion: free coordinates of their images vanish
    assert all(G[i][j] == 0 for i in fk for j in tk1), "torsion not preserved"
    Ft = [[G[i][j] % dk[i] for j in tk1] for i in tk]
    stacked = [Ft[r] + [dk[tk[r]] if s == r else 0 for s in range(len(tk))] for r in range(len(tk))]
    dd, _, _ = diagonalize(stacked)
    idx = 1
    for x in dd: idx *= x
    Tk = torsion_order(dk); Tk1 = torsion_order(dk1)
    im = Tk // idx
    return Tk1, Tk, im, Tk1 // im

def kills_p_torsion(Mk, Mk1, F, p):
    """Does the induced map coker Mk1 -> coker Mk annihilate the p-torsion subgroup of Tor coker Mk1?
    Returns (p-rank of Tor coker Mk1, killed?)."""
    dk, Pk, _ = diagonalize(Mk, True); dk1, _, Pk1inv = diagonalize(Mk1, True)
    G = mul(mul(Pk, F), Pk1inv)
    tk = [i for i,x in enumerate(dk) if x not in (0,1)]
    pj = [j for j,x in enumerate(dk1) if x not in (0,1) and x % p == 0]
    killed = True
    for j in pj:
        g = dk1[j] // p                       # generator of the p-torsion of the j-th cyclic factor
        for i in range(len(dk)):
            val = g * G[i][j]
            ok = (val == 0) if dk[i] == 0 else (val % dk[i] == 0)
            killed &= ok
    return len(pj), killed

def rank_mod_p(M, p):
    d,_,_ = diagonalize(M)
    return sum(1 for x in d if x % p != 0)

# ---------------------------------------------------------------- shadow graphs as arc lists
def shadow_graph(k, c):
    """vertices ('V',p) / ('S',p) for p in V_k ; arcs as (label, source, range)."""
    Vk = nb_paths(k); arcs = []
    for x in nb_paths(k+1):                       # arcs of Y_k : x[:-1] -> x[1:]
        arcs.append((('g',x), ('V',x[:-1]), ('V',x[1:])))
    for p in Vk:
        arcs.append((('in',p), ('V',p), ('S',p)))
        for i in range(c): arcs.append((('ret',p,i), ('S',p), ('V',p)))
        for j in range(2): arcs.append((('loop',p,j), ('S',p), ('S',p)))
    verts = [('V',p) for p in Vk] + [('S',p) for p in Vk]
    return verts, arcs

def tau_sharp_arc(a):
    lab = a[0]
    if lab[0] == 'g': return ('g', lab[1][:-1])
    if lab[0] == 'in': return ('in', lab[1][:-1])
    return (lab[0], lab[1][:-1], lab[2])
def tau_sharp_vert(v): return (v[0], v[1][:-1])

def covering_check(k, c):
    Vk, Ak = shadow_graph(k, c); Vk1, Ak1 = shadow_graph(k+1, c)
    labels_k = {a[0] for a in Ak}
    assert all(tau_sharp_arc(a) in labels_k for a in Ak1)
    src_k = {a[0]:a[1] for a in Ak}; rng_k = {a[0]:a[2] for a in Ak}
    # graph morphism
    morph = all(src_k[tau_sharp_arc(a)] == tau_sharp_vert(a[1]) and rng_k[tau_sharp_arc(a)] == tau_sharp_vert(a[2]) for a in Ak1)
    in_k = {}; out_k = {}
    for a in Ak: in_k.setdefault(a[2],[]).append(a[0]); out_k.setdefault(a[1],[]).append(a[0])
    in_k1 = {}; out_k1 = {}
    for a in Ak1: in_k1.setdefault(a[2],[]).append(a[0]); out_k1.setdefault(a[1],[]).append(a[0])
    in_bij = all(sorted(tau_sharp_arc((l,)) for l in in_k1[w]) == sorted(in_k[tau_sharp_vert(w)]) for w in Vk1)
    out_bij = all(sorted(tau_sharp_arc((l,)) for l in out_k1[w]) == sorted(out_k[tau_sharp_vert(w)]) for w in Vk1)
    surj = {tau_sharp_vert(w) for w in Vk1} == set(Vk)
    return morph, surj, in_bij, out_bij

def strongly_connected(B):
    n = len(B); reach = [[B[i][j] > 0 or i == j for j in range(n)] for i in range(n)]
    for k_ in range(n):
        for i in range(n):
            if reach[i][k_]:
                for j in range(n):
                    if reach[k_][j]: reach[i][j] = True
    return all(all(r) for r in reach)

# ---------------------------------------------------------------- pairing  Tor coker M  x  Tor coker M^t -> Q/Z
def pairing(M, x, y):
    """<[x],[y]> = y^t z mod 1, where M z = x over Q; requires deg-zero (torsion) classes."""
    Ms = sp.Matrix(M); sol, params = Ms.gauss_jordan_solve(sp.Matrix(x))
    z = sol.subs({p:0 for p in params})
    val = sum(sp.Rational(int(yi))*zi for yi,zi in zip(y, z))
    return sp.Rational(val) % 1

# ================================================================= run
c = Q_ - 1
print(f"K4 Iwahori tower, q = {Q_}, shadow return arcs c = q-1 = {c}\n")
LV = {k: level(k) for k in (1,2,3)}
ords2 = {}
plain = {}
for k in (1,2,3):
    V, idx, A, D = LV[k]; n = len(V); I = eye(n)
    B = shadow(A, c); M = sub(eye(2*n), T(B)); Mn = sub(eye(2*n), B)
    U = block(I, scal(-c,I), zeros(n,n), I); Vm = block(I, zeros(n,n), scal(-1,I), I)
    S1 = eq(mul(mul(U,M),Vm), blockdiag(T(D), scal(-1,I)))
    S2 = eq(mul(mul(T(Vm),Mn),T(U)), blockdiag(D, scal(-1,I)))
    d,_,_ = diagonalize(M); inv, free = invariant_factors(d, 2*n)
    dD,_,_ = diagonalize(D); invD, freeD = invariant_factors(dD, n)
    tor = torsion_order(d); ords2[k] = v2(tor); odd = tor >> v2(tor)
    # wrong normalization c = q
    Bw = shadow(A, Q_); dw,_,_ = diagonalize(sub(eye(2*n), T(Bw))); invw, freew = invariant_factors(dw, 2*n)
    # plain tower
    dp,_,_ = diagonalize(sub(I, T(A))); plain[k] = invariant_factors(dp, n)
    G = strongly_connected(B) and B[n][n] == 2
    print(f"level k={k}: n_k={n}, 2n_k={2*n}")
    print(f"  (S) U(I-B^t)V = Delta^t(+)(-I): {S1};  V^t(I-B)U^t = Delta(+)(-I): {S2}")
    print(f"  (K) K_0(O_k) = coker(I-B^t) = {group_str(inv,free)};  coker(Delta_k) = {group_str(invD,freeD)};  match: {inv==invD and free==freeD==1}")
    print(f"      |Tor_k| = 2^{v2(tor)} * {odd}   (2-adic order {v2(tor)}; odd part {odd})")
    print(f"  (W) c = q instead: coker(I-B^t) = {group_str(invw,freew)}  -> K_1 = 0, finite K_0, not Z(+)Jac")
    print(f"  (G) Kirchberg witnesses (strongly connected, double loop at shadow): {G}")
    print(f"  (B) plain C*(Y_k): coker(I-A_k^t) = {group_str(*plain[k])}")
print(f"\n2-adic orders (7,18,41 predicted by mu p^k - k + nu_0, mu=6, nu_0=-4): {[ords2[k] for k in (1,2,3)]}  "
      f"{[ords2[k]==6*2**k-k-4 for k in (1,2,3)]}")
print(f"plain tower blind: SNF(I-A_k^t) constant over k=1,2,3: {plain[1]==plain[2]==plain[3]}  "
      f"(= Z^r (+) Z/(r-1), r=3 [I, Thm 3.7])\n")

for k in (1,2):
    V, idx, A, D = LV[k]; V1, idx1, A1, D1 = LV[k+1]; n, n1 = len(V), len(V1)
    tau = tail(k); taut = T(tau)
    B = shadow(A,c); B1 = shadow(A1,c)
    M, M1 = sub(eye(2*n), T(B)), sub(eye(2*n1), T(B1))        # K_0 side (transposed)
    Mn, Mn1 = sub(eye(2*n), B), sub(eye(2*n1), B1)           # Bowen--Franks side
    Tup = blockdiag(taut, taut); Tdn = blockdiag(tau, tau)
    # (C)
    morph, surj, inb, outb = covering_check(k, c)
    print(f"transition {k}->{k+1}:")
    print(f"  (C) tau^sharp graph morphism: {morph}; surjective: {surj}; bijective on in-arcs: {inb}; bijective on out-arcs: {outb}")
    # (U)
    U1 = eq(mul(Tup, M), mul(M1, Tup))
    one = vec([1]*(2*n)); one1 = vec([1]*(2*n1))
    U2 = eq(mul(Tup, one), one1)
    kvec = vec([1]*n + [-1]*n); kvec1 = vec([1]*n1 + [-1]*n1)
    U3 = eq(mul(M, kvec), zeros(2*n,1)) and eq(mul(Tup, kvec), kvec1)
    print(f"  (U) T^k(I-B_k^t) = (I-B_(k+1)^t)T^k: {U1};  T^k 1 = 1: {U2};  ker(I-B^t)=Z(1,-1), T^k(1,-1)=(1,-1): {U3}")
    # (D)
    D1ok = eq(mul(Tdn, Mn1), mul(Mn, Tdn))
    D2bad = eq(mul(Tdn, M1), mul(M, Tdn))
    kv = vec([1]*n1 + [-c]*n1); kvk = vec([1]*n + [-c]*n)
    D3 = eq(mul(Mn1, kv), zeros(2*n1,1)) and eq(mul(Tdn, kv), scal(Q_, kvk))
    print(f"  (D) T_k(I-B_(k+1)) = (I-B_k)T_k: {D1ok};  same with transposes: {D2bad} (must be False);  T_k(1,-c1) = q(1,-c1): {D3}")
    # (N)
    Tk1, Tk, im, ker = coker_map_on_torsion(Mn, Mn1, Tdn)
    Tk1d, Tkd, imd, kerd = coker_map_on_torsion(D, D1, tau)
    print(f"  (N) shadow Bowen--Franks: |Tor_(k+1)|=2^{v2(Tk1)}*{Tk1>>v2(Tk1)}, |Tor_k|=2^{v2(Tk)}*{Tk>>v2(Tk)}, "
          f"surjective: {im==Tk}, |ker| = 2^{v2(ker)} (odd part {ker>>v2(ker)})")
    print(f"      sandpile tau_*   : surjective: {imd==Tkd}, |ker| = 2^{v2(kerd)};  identical to shadow transfer: {(im,ker)==(imd,kerd)};  "
          f"predicted 2^{n*(Q_-1)-1}: {ker == 2**(n*(Q_-1)-1)}")
    # (Q) kernel of the tail norm = p-torsion subgroup  (q = p = 2)
    rkA = rank_mod_p(A1, Q_)
    prank, killed = kills_p_torsion(D, D1, tau, Q_)
    print(f"  (Q) rank_F2(A_(k+1)) = {rkA} = n_k: {rkA==n};  2-rank of Jac(Y_(k+1)) = {prank} = n_k(q-1)-1: {prank==n*(Q_-1)-1};  "
          f"tau_* kills Jac(Y_(k+1))[2]: {killed};  |Jac[2]| = 2^{prank} = |ker tau_*|: {2**prank==kerd}")
    # (P) adjointness on random torsion classes
    ok = True
    for _ in range(4):
        a,b = random.sample(range(n1),2); cc,dd = random.sample(range(n),2)
        x = [0]*n1; x[a] += 1; x[b] -= 1
        y = [0]*n; y[cc] += 1; y[dd] -= 1
        lhs = pairing(D, [r[0] for r in mul(tau, vec(x))], y)
        rhs = pairing(D1, x, [r[0] for r in mul(taut, vec(y))])
        ok &= (lhs == rhs)
    print(f"  (P) <tau_* x, y>_k = <x, tau^* y>_(k+1) in Q/Z on 4 random torsion pairs: {ok}")
    # (L)
    degs = all(sum(col) == Q_ for col in T(taut)) and all(sum(col) == 1 for col in T(tau))
    print(f"  (L) deg(tau^* y) = q deg(y), deg(tau_* x) = deg(x): {degs}\n")
print("done.")
