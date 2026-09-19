"""
Paper VI verification: the shadow (Hecke double) graph  B(Y) = [[T, I],[qI, 2I]].

 (S)  Schur/unimodular identity:  U (I - B^t) V = Delta_p (+) (-I)   with explicit
      integer unimodular U, V  =>  coker(I-B^t) = Z (+) Jac(Y),  ker = Z.
 (K)  SNF(I - B^t) vs SNF(Delta_p): invariant factors match Table 1 of [I].
 (G)  Kirchberg witnesses: B strongly connected; shadow vertices carry two loops
      (two distinct cycles at one vertex + aperiodicity).
 (N)  Non-contradiction with [II]: spec(B) is NOT that of the companion C
      (traces differ); the only shared invariant is Bowen-Franks:
      det(I - uB) = prod_lambda [ (1-2u)(1-u lambda) - q u^2 ],  at u=1 equal to
      (-1)^n det(Delta_p)  -- the exact reconciliation identity.
 (F)  functoriality spot-check: automorphism equivariance on K4.
"""
import numpy as np, sympy as sp, itertools
from sympy.matrices.normalforms import smith_normal_form

def shadow(T,q):
    n=T.shape[0]
    B=np.zeros((2*n,2*n),dtype=np.int64)
    B[:n,:n]=T; B[:n,n:]=np.eye(n,dtype=np.int64)
    B[n:,:n]=q*np.eye(n,dtype=np.int64); B[n:,n:]=2*np.eye(n,dtype=np.int64)
    return B

def snf_inv(M):
    S=smith_normal_form(sp.Matrix(M.tolist()))
    d=[abs(S[i,i]) for i in range(min(S.rows,S.cols))]
    return [x for x in d if x not in (0,1)], sum(1 for x in d if x==0)+ (M.shape[1]-len(d))

def strongly_connected(B):
    n=B.shape[0]; R=(B>0).astype(np.int64); P=np.eye(n,dtype=np.int64)
    acc=np.zeros_like(R)
    for _ in range(n): P=P@R; acc=acc|(P>0)
    return bool(acc.all())

def check(name, T, q, jac_expected):
    T=np.array(T,dtype=np.int64); n=T.shape[0]
    B=shadow(T,q); M=np.eye(2*n,dtype=np.int64)-B.T
    tors,free=snf_inv(M)
    D=(q+1)*np.eye(n,dtype=np.int64)-T
    torsD,freeD=snf_inv(D)
    K = (tors==torsD==jac_expected) and free==freeD==1
    # explicit unimodular factorization:  M ~ Delta (+) (-I)
    # I-B^t = [[I-T, -qI],[-I, -I]];  U=[[I, -qI],[0, I]], V=[[I,0],[-I,I]]  (check)
    I=np.eye(n,dtype=np.int64); Z=np.zeros((n,n),dtype=np.int64)
    U=np.block([[I,-q*I],[Z,I]]); V=np.block([[I,Z],[-I,I]])
    S_=U@M@V
    Sok = (S_[:n,:n]==D).all() and (S_[:n,n:]==0).all() and (S_[n:,:n]==0).all() and (S_[n:,n:]==-I).all()
    G = strongly_connected(B) and B[n,n]==2
    trB=int(B.trace()); trC=0   # tr(C)=tr(T)=0 for simple graphs
    N1 = trB!=trC
    # reconciliation at u=1: det(I-B) == (-1)^n det(Delta)
    detI_B=int(round(np.linalg.det((np.eye(2*n)-B))))
    detD=int(round(np.linalg.det(D.astype(float))))
    N2 = detI_B == ((-1)**n)*detD
    print(f"{name:9s} n={n:2d} q={q}: (S){Sok} (K){K} Jac={tors} (G){G} (N)tr(B)={trB}!=0:{N1}, det(I-B)=(-1)^n detD:{N2}")
    return B

check("K4", [[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]], 2, [4,4])
check("K5", [[1 if i!=j else 0 for j in range(5)] for i in range(5)], 3, [5,5,5])
k33=[[0]*6 for _ in range(6)]
for i in range(3):
    for j in range(3,6): k33[i][j]=k33[j][i]=1
check("K33", k33, 2, [3,3,9])
V=list(itertools.combinations(range(5),2)); pet=[[0]*10 for _ in range(10)]
for i,a in enumerate(V):
    for j,b in enumerate(V):
        if i!=j and not set(a)&set(b): pet[i][j]=1
check("Petersen", pet, 2, [2,10,10,10])
# (F) equivariance on K4: transposition (0 1)
T=np.array([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]],dtype=np.int64); B=shadow(T,2)
pi=np.eye(4,dtype=np.int64)[[1,0,2,3]]; P=np.block([[pi,np.zeros((4,4),dtype=np.int64)],[np.zeros((4,4),dtype=np.int64),pi]])
print("(F) automorphism equivariance P B P^-1 = B:", bool((P@B@P.T==B).all()))
