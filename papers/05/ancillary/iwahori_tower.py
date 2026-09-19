"""
The Iwahori path tower over K4 (q=2): the canonical combinatorial model of
Gamma_0(p^k)\X.  Level-k vertices = non-backtracking paths of length k in Y0
(k=0: vertices; k=1: oriented edges with the Hashimoto adjacency A^nb; k>=2:
iterated line digraphs of level 1).  |V(Y_k)| = n(q+1)q^{k-1}, matching the
index [Gamma(1):Gamma_0(p^k)] = (p+1)p^{k-1} exactly.

For k>=1 these are q-out-regular Eulerian digraphs; the sandpile group is
coker of the reduced out-degree Laplacian qI - A_k (sink-independent for
Eulerian digraphs).  We compute exact orders and 2-adic valuations.
"""
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

def K4_hashimoto():
    Tp=[[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]]
    E=[(u,v) for u in range(4) for v in range(4) if Tp[u][v]]
    arcs={}
    for i,(u,v) in enumerate(E):
        arcs[i]=[j for j,(x,y) in enumerate(E) if x==v and not(y==u)]
    return arcs   # adjacency dict: vertex -> list of out-neighbors

def line_digraph(adj):
    arcs=[(u,w) for u in adj for w in adj[u]]
    idx={a:i for i,a in enumerate(arcs)}
    new={}
    for (u,w) in arcs:
        new[idx[(u,w)]]=[idx[(w,x)] for x in adj[w]]
    return new

def sandpile_order_and_val(adj, q, do_snf=False):
    n=len(adj)
    L=sp.zeros(n,n)
    for u in adj:
        L[u,u]=q
        for w in adj[u]:
            L[u,w]-=1
    M=L[1:,1:]
    det=int(M.det())
    v2=0; d=abs(det)
    while d%2==0 and d: d//=2; v2+=1
    snf=None
    if do_snf:
        S=smith_normal_form(sp.Matrix(M))
        snf=[abs(S[i,i]) for i in range(S.rows) if abs(S[i,i]) not in (0,1)]
    return n, abs(det), v2, snf

# level 0: undirected K4 sandpile, q+1=3
Tp=sp.Matrix([[0,1,1,1],[1,0,1,1],[1,1,0,1],[1,1,1,0]])
L0=3*sp.eye(4)-Tp
det0=int(L0[1:,1:].det()); v20=0; d=det0
while d%2==0 and d: d//=2; v20+=1
print(f"k=0: n=   4  |Jac|={det0}  v2={v20}")

adj=K4_hashimoto(); q=2
prev_v2=v20
rows=[(0,4,det0,v20)]
for k in range(1,7):
    do_snf = (k<=3)
    n,order,v2,snf = sandpile_order_and_val(adj,q,do_snf)
    inc=v2-prev_v2
    print(f"k={k}: n={n:4d}  |sandpile|={order}  v2={v2}  increment={inc}"
          + (f"  SNF-torsion={snf}" if snf else ""))
    rows.append((k,n,order,v2)); prev_v2=v2
    if k<6: adj=line_digraph(adj)
