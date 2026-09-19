import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

def coker(M):
    S=smith_normal_form(sp.Matrix(M)); d=[abs(S[i,i]) for i in range(min(S.rows,S.cols))]
    return sp.Matrix(M).rows-len([x for x in d if x!=0]), [x for x in d if x not in (0,1)]

def oriented(adj):
    return [(u,v) for u in range(len(adj)) for v in range(len(adj)) if adj[u][v]]

def anb(adj):
    E=oriented(adj); m=len(E); B=sp.zeros(m,m)
    for i,(u,v) in enumerate(E):
        for j,(x,y) in enumerate(E):
            if v==x and not(x==v and y==u): B[i,j]=1
    return B,E

def reduced_B(adj):
    """The (m+1) x n presentation matrix B = [ D^t ; -lambda ] of the proof."""
    n=len(adj); Ep=[(u,v) for u in range(n) for v in range(u+1,n) if adj[u][v]]
    m=len(Ep)
    B=sp.zeros(m+1,n)
    for i,(o,t) in enumerate(Ep):
        B[i,t]+=1; B[i,o]-=1                       # D^t
    dminus=[sum(1 for (o,t) in Ep if o==v) for v in range(n)]
    for v in range(n): B[m,v]=-(1-dminus[v])       # -lambda
    return B,m,n

def check(name,adj,q):
    A,E=anb(adj); f1,t1=coker((sp.eye(len(E))-A).T)
    B,m,n=reduced_B(adj); f2,t2=coker(B)
    r=m-n+1
    print(f"{name:22s} r={r:2d} | coker(I-A^nb): Z^{f1}+{t1} | coker(B): Z^{f2}+{t2} | "
          f"predicted Z^{r}+[{r-1}] | {'OK' if (f1,t1)==(f2,t2)==(r,[r-1] if r-1>1 else []) else 'MISMATCH'}")

def complete(n): return [[1 if i!=j else 0 for j in range(n)] for i in range(n)]
def prism(k):
    n=2*k; a=[[0]*n for _ in range(n)]
    for i in range(k):
        for (x,y) in [(i,(i+1)%k),(k+i,k+(i+1)%k),(i,k+i)]: a[x][y]=a[y][x]=1
    return a
def petersen():
    import itertools
    V=list(itertools.combinations(range(5),2)); n=len(V); a=[[0]*n for _ in range(n)]
    for i,x in enumerate(V):
        for j,y in enumerate(V):
            if i!=j and not set(x)&set(y): a[i][j]=1
    return a
def heawood():
    L=[(0,1,3),(1,2,4),(2,3,5),(3,4,6),(4,5,0),(5,6,1),(6,0,2)]
    a=[[0]*14 for _ in range(14)]
    for li,Ln in enumerate(L):
        for p in Ln: a[p][7+li]=a[7+li][p]=1
    return a

check("K_4 (q=2)",complete(4),2); check("K_5 (q=3)",complete(5),3)
check("K_6 (q=4)",complete(6),4); check("prism C_4xK_2",prism(4),2)
check("prism C_5xK_2",prism(5),2); check("Petersen",petersen(),2)
check("Heawood",heawood(),2)
