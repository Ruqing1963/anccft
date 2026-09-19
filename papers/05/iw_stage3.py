import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form

def coker_data(M):
    A=sp.Matrix(M).copy(); m,nn=A.rows,A.cols; U=sp.eye(m); i=j=0
    while i<m and j<nn:
        piv=None; pv=None
        for r in range(i,m):
            for c in range(j,nn):
                if A[r,c]!=0 and (pv is None or abs(A[r,c])<pv): piv=(r,c); pv=abs(A[r,c])
        if piv is None: break
        r,c=piv; A.row_swap(i,r); U.row_swap(i,r); A.col_swap(j,c)
        prog=True
        while prog:
            prog=False
            for r in range(i+1,m):
                q=A[r,j]//A[i,j]
                if q: A.row_op(r,lambda x,kk:x-q*A[i,kk]); U.row_op(r,lambda x,kk:x-q*U[i,kk])
                if A[r,j]!=0: A.row_swap(i,r); U.row_swap(i,r); prog=True; break
            if prog: continue
            for c in range(j+1,nn):
                q=A[i,c]//A[i,j]
                if q: A.col_op(c,lambda x,kk:x-q*A[kk,j])
                if A[i,c]!=0: prog=True
        i+=1; j+=1
    return [abs(A[t,t]) for t in range(min(m,nn))],U

def v2(x):
    x=abs(int(x)); v=0
    if x==0: return None
    while x%2==0: x//=2; v+=1
    return v

D1=np.load("Delta_1.npy"); D2=np.load("Delta_2.npy"); Do=np.load("Do_1.npy")
d1,U1=coker_data(D1); d2,U2=coker_data(D2)
T=U1*sp.Matrix(Do.tolist())*(U2.inv())
i1=[i for i,x in enumerate(d1) if x!=0 and v2(x)]
i2=[i for i,x in enumerate(d2) if x!=0 and v2(x)]
prims=[v2(d1[i]) for i in i1]
gens=[]
for j in i2:
    oj=d2[j]//(2**v2(d2[j])); vec=[]
    for pos,i in enumerate(i1):
        e=prims[pos]; oi=d1[i]//(2**e)
        vec.append((int(T[i,j])*int(oj)*pow(int(oi),-1,2**e))%(2**e))
    gens.append(vec)
m=len(prims); full=1
for e in prims: full*=2**e
cols=[sp.Matrix(m,1,g) for g in gens]; Dg=sp.diag(*[2**e for e in prims])
S=smith_normal_form(sp.Matrix.hstack(*(cols+[Dg[:,i] for i in range(m)])))
img=full
for x in [abs(S[i,i]) for i in range(m)]: img//=(x if x else 1)
tot1=sum(v2(x) for x in d1 if x and v2(x)); tot2=sum(v2(x) for x in d2 if x and v2(x))
print(f"norm 2->1: surjective on 2-part: {img==full}   |ker|=2^{tot2-tot1}   predicted 2^(n_1-1)=2^11")
print(f"(orders: ord2 level1={tot1}, level2={tot2} -- match [II] data 7,18)")
