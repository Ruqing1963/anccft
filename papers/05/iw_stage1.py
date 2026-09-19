import numpy as np, json
Tp=np.array([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]],dtype=np.int64)
def nb_paths(k):
    if k==0: return [(v,) for v in range(4)]
    ps=[(u,v) for u in range(4) for v in range(4) if Tp[u][v]]
    for _ in range(k-1):
        ps=[p+(w,) for p in ps for w in range(4) if Tp[p[-1]][w] and w!=p[-2]]
    return ps
def level(k):
    V=nb_paths(k); idx={p:i for i,p in enumerate(V)}
    if k==0: return V,idx,Tp.copy(),3*np.eye(4,dtype=np.int64)-Tp
    W=nb_paths(k+1); n=len(V); A=np.zeros((n,n),dtype=np.int64)
    for w in W: A[idx[w[:-1]],idx[w[1:]]]+=1
    return V,idx,A,2*np.eye(n,dtype=np.int64)-A
def degen(k):
    Vk,ik,_,_=level(k); Vk1,ik1,_,_=level(k+1)
    Do=np.zeros((len(Vk),len(Vk1)),dtype=np.int64); Di=np.zeros_like(Do)
    for w in Vk1: Do[ik[w[:-1]],ik1[w]]=1; Di[ik[w[1:]],ik1[w]]=1
    return Do,Di
out={}
for k in (1,2):
    Vk,ik,Ak,Dk=level(k); Vk1,ik1,Ak1,Dk1=level(k+1)
    Do,Di=degen(k); n=len(Vk)
    F = (Do@Di.T==Ak).all() and (Di.T@Do==Ak1).all() and (Di@Di.T==2*np.eye(n,dtype=np.int64)).all() and (Do@Do.T==2*np.eye(n,dtype=np.int64)).all()
    Tt= (Do@Dk1==Dk@Do).all() and (Do@Ak1==Ak@Do).all()
    Din_in=(Di@(2*np.eye(len(Vk1),dtype=np.int64)-Ak1.T)==(2*np.eye(n,dtype=np.int64)-Ak.T)@Di).all()
    Hasym=(Di@Dk1==2*(Di-Do)).all()
    Cnn=not (Ak@Ak.T==Ak.T@Ak).all()
    out[k]=dict(n=n,F=bool(F),T=bool(Tt),Hmirror=bool(Din_in),Hasym=bool(Hasym),C=bool(Cnn))
    np.save(f"Do_{k}.npy",Do); np.save(f"Di_{k}.npy",Di); np.save(f"Delta_{k}.npy",Dk); np.save(f"Delta_{k+1}.npy",Dk1)
_,_,A3,_=level(3)
out[3]=dict(n=A3.shape[0],C=bool(not (A3@A3.T==A3.T@A3).all()))
# level-0 defect + T_p composite
V0,i0,A0,D0=level(0); E1=nb_paths(1)
S=np.zeros((4,12),dtype=np.int64); T=np.zeros((4,12),dtype=np.int64)
for j,(u,v) in enumerate(E1): S[u,j]=1; T[v,j]=1
_,_,A1,D1=level(1)
out['Z']=bool((S@D1-D0@S==-(S-T)).all()); out['Tp']=bool((S@T.T==Tp).all())
print(json.dumps(out,indent=1))
