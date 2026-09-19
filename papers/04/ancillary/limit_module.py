"""
Verification battery for the limit-module theorem (abelian voltage towers over K4):

  (B) block-circulant identity:  Delta_k = D(S_N)  in the ordering idx(v,g)=v*N+g
  (S) deck equivariance:         R_k Delta_k = Delta_k R_k,  R_k = I_n (x) shift
  (P) norm intertwining:         P Delta_{k+1} = Delta_k P,  P = I_n (x) (g -> g mod N)
  (C) deck-norm compatibility:   P R_{k+1} = R_k P
  (N) induced norm on l-primary torsion: surjective, kernel order = l^{ord_{k+1}-ord_k}
"""
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

def base_edges(volt):
    a12,a13,a23=volt
    return [((0,1),0),((0,2),0),((0,3),0),((1,2),a12),((1,3),a13),((2,3),a23)]

def shift(N):
    S=sp.zeros(N,N)
    for g in range(N): S[(g+1)%N,g]=1
    return S

def delta_level(volt,N,n=4):
    E=base_edges(volt); size=n*N
    A=sp.zeros(size,size)
    idx=lambda v,g: v*N+(g%N)
    for (u,v),a in E:
        for g in range(N):
            A[idx(u,g),idx(v,g+a)]+=1
            A[idx(v,g+a),idx(u,g)]+=1
    return 3*sp.eye(size)-A

def delta_block(volt,N,n=4):
    # assemble D(S_N) from the pencil coefficients
    E=base_edges(volt); S=shift(N)
    A=sp.zeros(n*N,n*N)
    def put(u,v,M):
        for i in range(N):
            for j in range(N):
                if M[i,j]:
                    A[u*N+i,v*N+j]+=M[i,j]
    for (u,v),a in E:
        Sa=S**(a%N) if a>=0 else (S.T)**((-a)%N)
        put(u,v,Sa.T)   # careful with convention; test both if mismatch
        put(v,u,Sa)
    return 3*sp.eye(n*N)-A

def coker_data(M):
    """Return (invariants>1 list, U) with U*M*V = SNF; coker = Z^m / diag via y -> U y."""
    M=sp.Matrix(M)
    # sympy smith_normal_form doesn't return transforms; do our own via hermite-ish:
    # use sympy's matrix normal form with transformation: implement simple algorithm
    A=M.copy(); m,nn=A.rows,A.cols
    U=sp.eye(m); V=sp.eye(nn)
    def improve():
        # standard SNF with transforms (small matrices only)
        i=0; j=0
        while i<m and j<nn:
            # find pivot = nonzero min abs
            piv=None; pv=None
            for r in range(i,m):
                for c in range(j,nn):
                    x=A[r,c]
                    if x!=0 and (pv is None or abs(x)<pv): piv=(r,c); pv=abs(x)
            if piv is None: break
            r,c=piv
            A.row_swap(i,r); U.row_swap(i,r)
            A.col_swap(j,c); V.col_swap(j,c)
            done=False
            while not done:
                done=True
                for r in range(i+1,m):
                    q=A[r,j]//A[i,j]
                    if q: A.row_op(r,lambda x,k:x-q*A[i,k]); U.row_op(r,lambda x,k:x-q*U[i,k])
                for r in range(i+1,m):
                    if A[r,j]!=0:
                        A.row_swap(i,r); U.row_swap(i,r); done=False; break
                if not done: continue
                for c in range(j+1,nn):
                    q=A[i,c]//A[i,j]
                    if q: A.col_op(c,lambda x,k:x-q*A[k,j]); V.col_op(c,lambda x,k:x-q*V[k,j])
                for c in range(j+1,nn):
                    if A[i,c]!=0: done=False; break
                if not done:
                    # restart pivot cleanup on this position
                    continue
            i+=1; j+=1
        return A,U,V
    A,U,V=improve()
    # ensure divisibility not needed for our order computations; collect diagonal
    d=[abs(A[i,i]) for i in range(min(m,nn))]
    return d,U

def vl(x,l):
    x=abs(int(x)); v=0
    if x==0: return None
    while x%l==0: x//=l; v+=1
    return v

def lpart_order(d,l):
    return sum(vl(x,l) for x in d if x not in (0,) and vl(x,l) is not None and x!=0)

def subgroup_order_in(prims, gens, l):
    """prims: exponents e_i of ⊕ Z/l^{e_i}; gens: list of vectors (ints).
       order of generated subgroup = |⊕| / |coker([gens | diag(l^e)])|."""
    m=len(prims)
    if m==0: return 1
    cols=[sp.Matrix(m,1,[g[i] for i in range(m)]) for g in gens]
    D=sp.diag(*[l**e for e in prims])
    Mat=sp.Matrix.hstack(*(cols+[D[:,i] for i in range(m)])) if cols else D
    S=smith_normal_form(Mat)
    dd=[abs(S[i,i]) for i in range(min(S.rows,S.cols))]
    cok=1
    for x in dd[:m]:
        cok*= (x if x!=0 else 0)
    total=1
    for e in prims: total*= l**e
    return total//cok

def run(volt,l,kmax):
    print(f"\n===== tower voltages {volt}, l={l} =====")
    deltas={}; orders={}
    for k in range(kmax+1):
        N=l**k
        Dk=delta_level(volt,N)
        deltas[k]=Dk
        # (B) block identity (try both edge conventions)
        Db=delta_block(volt,N)
        okB = (Dk==Db)
        if not okB:
            # swap convention
            E=base_edges(volt); S=shift(N); n=4
            A=sp.zeros(4*N,4*N)
            for (u,v),a in E:
                Sa=S**(a%N)
                for i in range(N):
                    for j in range(N):
                        if Sa[i,j]: A[u*N+j,v*N+i]+=1; A[v*N+i,u*N+j]+=1
            okB = (Dk == 3*sp.eye(4*N)-A)
        # (S) deck equivariance
        R=sp.zeros(4*N,4*N)
        for v in range(4):
            for g in range(N): R[v*N+(g+1)%N, v*N+g]=1
        okS = (R*Dk == Dk*R)
        d,U=coker_data(Dk)
        orders[k]=lpart_order(d,l)
        print(f" k={k}: n={4*N:3d}  (B) block-circulant: {okB}   (S) sigma-equivariance: {okS}   ord_{l}|Tor|={orders[k]}")
    for k in range(kmax):
        N=l**k; M=l**(k+1)
        P=sp.zeros(4*N,4*M)
        for v in range(4):
            for g in range(M): P[v*N+(g%N), v*M+g]=1
        okP = (P*deltas[k+1] == deltas[k]*P)
        Rk=sp.zeros(4*N,4*N); Rk1=sp.zeros(4*M,4*M)
        for v in range(4):
            for g in range(N): Rk[v*N+(g+1)%N, v*N+g]=1
            for g in range(M): Rk1[v*M+(g+1)%M, v*M+g]=1
        okC = (P*Rk1 == Rk*P)
        # (N) induced norm on l-primary parts
        dk,Uk = coker_data(deltas[k])
        dk1,Uk1 = coker_data(deltas[k+1])
        # generators of coker_{k+1}: classes of columns of Uk1^{-1}; induced map matrix:
        T = Uk*P*(Uk1.inv())
        idx_k  =[i for i,x in enumerate(dk)  if x!=0 and vl(x,l)]
        idx_k1 =[i for i,x in enumerate(dk1) if x!=0 and vl(x,l)]
        prims=[vl(dk[i],l) for i in idx_k]
        gens=[]
        for j in idx_k1:
            oddj = dk1[j]//(l**vl(dk1[j],l))
            vec=[]
            for pos,i in enumerate(idx_k):
                e=prims[pos]; di=dk[i]; oddi=di//(l**e)
                val=int(T[i,j])*int(oddj)
                inv=pow(int(oddi), -1, l**e)
                vec.append((val*inv) % (l**e))
            gens.append(vec)
        img = subgroup_order_in(prims,gens,l)
        full = 1
        for e in prims: full*= l**e
        surj = (img==full)
        ker_ord = (l**orders[k+1])//img if img else None
        print(f"  norm {k+1}->{k}: (P) intertwine: {okP}  (C) deck-compat: {okC}  "
              f"(N) surjective on l-part: {surj}  |ker|=l^{vl(ker_ord,l) if ker_ord else '?'} "
              f"(predicted l^{orders[k+1]-orders[k]})")

run((1,1,1),2,3)
run((1,0,0),3,2)
