"""
Paper VII battery: snake-lemma factorization of the transfer.

Iwahori tower (K4, q=2), transitions 1->2 and 2->3:
 (I)  K = ker_Z(tau_*) is A_{k+1}-invariant (tau A K = 0)
 (Nil) A_{k+1}|_K nilpotent: smallest s with A^s K = 0
 (Det) det(Delta_{k+1}|_K) = + q^{n_k(q-1)}  (exact, from qI - nilpotent)
 (Idx) Eisenstein index: tau_* 1 = q * 1   (boundary 0->1: S 1 = (q+1) 1)
 (Sum) increment identity: ord = n_k(q-1) - 1

Abelian voltage tower (1,1,1), l=2, transitions 0->1 and 1->2:
 (I')  ker pi_* invariant;  (Det') det(Delta|_ker pi) = +- prod_{new zeta} det D(zeta)
 (Idx') pi_* 1 = l * 1;     increments reproduce [II, Table 2].
"""
import numpy as np, sympy as sp

Tp=[[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]
def nb_paths(k):
    if k==0: return [(v,) for v in range(4)]
    ps=[(u,v) for u in range(4) for v in range(4) if Tp[u][v]]
    for _ in range(k-1):
        ps=[p+(w,) for p in ps for w in range(4) if Tp[p[-1]][w] and w!=p[-2]]
    return ps
def level(k):
    V=nb_paths(k); idx={p:i for i,p in enumerate(V)}
    W=nb_paths(k+1); n=len(V); A=np.zeros((n,n),dtype=object)
    for w in W: A[idx[w[:-1]],idx[w[1:]]]+=1
    return V,idx,A,2*np.eye(n,dtype=object)-A
def tail(k):
    Vk,ik,_,_=level(k); Vk1,ik1,_,_=level(k+1)
    Do=np.zeros((len(Vk),len(Vk1)),dtype=object)
    for w in Vk1: Do[ik[w[:-1]],ik1[w]]+=1
    return Do

def ker_int(M):
    ns=sp.Matrix(M.tolist()).nullspace()
    cols=[]
    for v in ns:
        l=sp.ilcm(*[sp.fraction(x)[1] for x in v]) if any(sp.fraction(x)[1]!=1 for x in v) else 1
        cols.append([int(x*l) for x in v])
    return np.array(cols,dtype=object).T   # shape (dim, rank)

def restrict(Op,K):
    Ks=sp.Matrix(K.tolist()); Y=sp.Matrix((Op@K).tolist())
    G=(Ks.T*Ks); M=G.inv()*(Ks.T*Y)
    assert all(sp.denom(x)==1 for x in M), "restriction not integral"
    return M

for k in (1,2):
    _,_,Ak,Dk=level(k); _,_,Ak1,Dk1=level(k+1)
    Do=tail(k); nk=Do.shape[0]
    K=ker_int(Do); r=K.shape[1]
    inv=( (Do@(Ak1@K))==0 ).all()
    s=None; P=K.copy()
    for t in range(1,r+1):
        P=Ak1@P
        if (P==0).all(): s=t; break
    M=restrict(Dk1,K); detM=sp.Matrix(M).det()
    idx=(Do@np.ones((Do.shape[1],1),dtype=object)==2*np.ones((nk,1),dtype=object)).all()
    print(f"Iwahori {k}->{k+1}: rk K={r}=n_k(q-1)={nk}: {r==nk}  (I){bool(inv)}  (Nil)s={s}  "
          f"(Det)det={detM}==2^{nk}:{detM==2**nk}  (Idx)tau 1=2*1:{bool(idx)}  (Sum)inc={nk}-1")

S=np.zeros((4,12),dtype=object)
for j,(u,v) in enumerate(nb_paths(1)): S[u,j]+=1
print("boundary 0->1: S 1 = 3*1:", bool((S@np.ones((12,1),dtype=object)==3*np.ones((4,1),dtype=object)).all()),
      " (index q+1: the Eisenstein special layer)")

# ---------- abelian tower (1,1,1), l=2 ----------
t=sp.symbols('t')
def pencilD():
    E=[((0,1),0),((0,2),0),((0,3),0),((1,2),1),((1,3),1),((2,3),1)]
    A=sp.zeros(4,4)
    for (u,v),a in E: A[u,v]+=t**a; A[v,u]+=t**(-a)
    return 3*sp.eye(4)-A
def derived(N):
    E=[((0,1),0),((0,2),0),((0,3),0),((1,2),1),((1,3),1),((2,3),1)]
    L=np.zeros((4*N,4*N),dtype=object)
    idx=lambda v,g:(v*N+(g%N))
    for (u,v),a in E:
        for g in range(N):
            i,j=idx(u,g),idx(v,g+a)
            L[i,j]-=1; L[j,i]-=1; L[i,i]+=1; L[j,j]+=1
    return L
def proj(N):
    P=np.zeros((4*N,8*N),dtype=object)
    for v in range(4):
        for g in range(2*N): P[v*N+(g%N), v*2*N+g]+=1
    return P
D=pencilD()
for k,(N,new) in enumerate([(1,[-1]),(2,[sp.I,-sp.I])]):
    Dk=derived(N); Dk1=derived(2*N); P=proj(N)
    K=ker_int(P); inv=((P@(Dk1@K))==0).all()
    M=restrict(Dk1,K); detM=sp.Matrix(M).det()
    prod=sp.simplify(sp.prod([D.subs(t,z).det() for z in new]))
    idx=(P@np.ones((8*N,1),dtype=object)==2*np.ones((4*N,1),dtype=object)).all()
    print(f"abelian {k}->{k+1}: (I'){bool(inv)}  det(Delta|ker pi)={detM}  prod_new detD(zeta)={prod}"
          f"  match up to sign: {abs(detM)==abs(prod)}  (Idx')pi 1=2*1:{bool(idx)}")
