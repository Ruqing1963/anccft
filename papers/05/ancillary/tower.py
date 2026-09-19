"""
Z/l^k voltage-cover tower over K_4: exact |Jac(Y_k)|, its l-part, and the
Iwasawa-type growth  ord_l|Jac(Y_k)| = mu*l^k + lambda*k + nu (k >> 0).

Method A (direct, small k): build the derived graph, reduced-Laplacian
determinant = kappa (Kirchhoff), and SNF for the group structure.
Method B (character factorization, any k):
   kappa(Y_k) = kappa(Y_0)/N * prod_{j=1}^{N-1} det D(zeta_N^j),   N = l^k,
computed exactly as  Res(t^N-1, f(t)) / det D(1)  with f = t^M det D(t) in Z[t].
Method B is verified against Method A before being trusted.
"""
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

t = sp.symbols('t')

# base graph K4, spanning tree (0,1),(0,2),(0,3) voltage 0; chords with voltages a
def K4_pencil(volt):   # volt = (a12, a13, a23) on chords (1,2),(1,3),(2,3)
    a12,a13,a23 = volt
    edges = [((0,1),0), ((0,2),0), ((0,3),0), ((1,2),a12), ((1,3),a13), ((2,3),a23)]
    n=4
    A = sp.zeros(n,n)
    for (u,v),a in edges:
        A[u,v] += t**a
        A[v,u] += t**(-a)
    D = 3*sp.eye(n) - A          # (q+1)I - A(t),  q=2
    return edges, D

def derived_graph_laplacian(edges, N, n=4):
    size = n*N
    L = sp.zeros(size,size)
    def idx(u,g): return u*N + (g % N)
    for (u,v),a in edges:
        for g in range(N):
            i,j = idx(u,g), idx(v,g+a)
            L[i,j] -= 1; L[j,i] -= 1
            L[i,i] += 1; L[j,j] += 1
    return L

def kappa_direct(L):
    M = L[1:,1:]
    return M.det()

def jac_structure(L):
    S = smith_normal_form(sp.Matrix(L))
    d = [abs(S[i,i]) for i in range(min(S.rows,S.cols))]
    return [x for x in d if x not in (0,1)]

def kappa_resultant(D, N):
    f_laurent = sp.simplify(D.det())
    # clear denominators: f(t) = t^M * det D(t)
    num, den = sp.fraction(sp.together(f_laurent))
    poly_num = sp.Poly(sp.expand(num), t)
    poly_den = sp.Poly(sp.expand(den), t)   # should be t^M
    F = poly_num
    detD1 = sp.simplify(f_laurent.subs(t,1))
    # prod over all N-th roots zeta of detD(zeta) = prod f_laurent(zeta)
    # = Res(t^N-1, num)/Res(t^N-1, den);  Res(t^N-1, t^M) = prod zeta^M = ((-1)^(N+1))^M
    R_num = sp.resultant(sp.Poly(t**N-1, t), poly_num)
    R_den = sp.resultant(sp.Poly(t**N-1, t), poly_den)
    total = sp.Rational(R_num, R_den)
    # kappa(Y_0) with these q: base Laplacian D(1)
    prod_nontriv = sp.simplify(total / detD1)   # this is prod_{j!=0} det D(zeta^j) ... but detD(1)=0!
    return total, detD1

# careful: det D(1) = 0 (base Laplacian singular). Use kappa formula instead:
# kappa(Y_k) = (kappa(Y_0)/N) * prod_{j=1..N-1} det D(zeta^j)
# and prod_{j=1..N-1} det D(zeta^j) = [ Res(t^N-1, f)/Res(t^N-1,t^M) ] / detD(1) fails; instead
# divide polynomials: (t^N-1)/(t-1) = 1+t+...+t^{N-1}; use Res(cyclN(t), f)/adjust:
def kappa_via_res(D, N, kappa0):
    f_laurent = sp.together(sp.simplify(D.det()))
    num, den = sp.fraction(f_laurent)
    pnum = sp.Poly(sp.expand(num), t); pden = sp.Poly(sp.expand(den), t)
    Phi = sp.Poly(sum(t**i for i in range(N)), t)     # (t^N-1)/(t-1)
    Rn = sp.resultant(Phi, pnum)
    Rd = sp.resultant(Phi, pden)
    prod_nontrivial = sp.Rational(Rn, Rd)             # prod_{j=1..N-1} detD(zeta^j)
    return sp.nsimplify(sp.Rational(kappa0, N) * prod_nontrivial)

def ordl(x, l):
    x = int(x); c=0
    while x % l == 0 and x != 0:
        x //= l; c += 1
    return c

def run(volt, l, kmax_direct=3, kmax_res=9):
    edges, D = K4_pencil(volt)
    L0 = derived_graph_laplacian(edges, 1)
    kappa0 = int(kappa_direct(L0))
    print(f"\n=== voltages {volt}, l={l};  base kappa(K4) = {kappa0}, Jac = {jac_structure(L0)} ===")
    # verify Method B against Method A
    for k in range(0, kmax_direct+1):
        N = l**k
        L = derived_graph_laplacian(edges, N)
        kd = int(kappa_direct(L))
        kr = int(kappa_via_res(D, N, kappa0)) if N>1 else kappa0
        js = jac_structure(L) if N <= 16//1 and 4*N<=32 else None
        print(f"k={k}: n={4*N:4d}  kappa_direct={kd}  kappa_resultant={kr}  match={kd==kr}"
              + (f"  Jac={js}" if js is not None else ""))
    # extend with Method B
    print("growth of ord_l kappa:")
    prev=None
    for k in range(0, kmax_res+1):
        N=l**k
        kr = kappa0 if N==1 else int(kappa_via_res(D, N, kappa0))
        o = ordl(kr, l)
        inc = "" if prev is None else f"  increment={o-prev}"
        print(f"  k={k:2d}: ord_{l}(kappa) = {o}{inc}")
        prev=o

run((1,0,0), 2)
run((1,1,1), 2)
run((1,2,0), 2)
run((1,0,0), 3, kmax_direct=2, kmax_res=6)
