import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

def coker(M):
    S=smith_normal_form(sp.Matrix(M)); d=[abs(S[i,i]) for i in range(min(S.rows,S.cols))]
    tors=[x for x in d if x not in (0,1)]; free=sp.Matrix(M).cols-len([x for x in d if x!=0])
    return free,tors

def oriented(adj):
    e=[]
    for u in range(len(adj)):
        for v in range(len(adj)):
            if adj[u][v]: e.append((u,v))
    return e

def anb(adj):
    E=oriented(adj); m=len(E); B=sp.zeros(m,m)
    for i,(u,v) in enumerate(E):
        for j,(x,y) in enumerate(E):
            if v==x and not(x==v and y==u): B[i,j]=1
    return B,m

def rep(name,adj,q):
    n=len(adj); T=sp.Matrix(adj); I=sp.eye(n)
    lam=sorted(T.eigenvals(multiple=True),key=lambda z:-float(z))
    kappa=sp.nsimplify(sp.prod([(q+1)-l for l in lam[1:]])/n)
    fD,tD=coker((q+1)*I-T)
    B,m2=anb(adj); fB,tB=coker((sp.eye(m2)-B).T)
    r=m2//2-n+1
    print(f"{name:26s} n={n:3d} r={r:3d} | SNF(Delta): Z^{fD}+{tD} kappa={kappa} | "
          f"SNF(I-ANB^t): Z^{fB}+{tB} | r-1={r-1}  match={tB==([r-1] if r-1>1 else [])}")

def complete(n): return [[1 if i!=j else 0 for j in range(n)] for i in range(n)]
def prism(k):
    n=2*k; a=[[0]*n for _ in range(n)]
    for i in range(k):
        for (x,y) in [(i,(i+1)%k),(k+i,k+(i+1)%k),(i,k+i)]: a[x][y]=a[y][x]=1
    return a
def mobius(k):
    n=2*k; a=[[0]*n for _ in range(n)]
    for i in range(n): a[i][(i+1)%n]=a[(i+1)%n][i]=1
    for i in range(k): a[i][i+k]=a[i+k][i]=1
    return a
def octa():
    n=6; a=[[1]*n for _ in range(n)]
    for i in range(n): a[i][i]=0
    for i in range(3): a[i][i+3]=a[i+3][i]=0
    return a
def cyc(n,k):
    a=[[0]*n for _ in range(n)]
    for i in range(n):
        for d in range(1,k+1): a[i][(i+d)%n]=a[(i+d)%n][i]=1
    return a
def petersen():
    import itertools
    V=list(itertools.combinations(range(5),2)); n=len(V)
    a=[[0]*n for _ in range(n)]
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

rep("K_4 (q=2)",complete(4),2)
rep("K_5 (q=3)",complete(5),3)
rep("K_6 (q=4)",complete(6),4)
rep("prism C_4xK_2 (q=2)",prism(4),2)
rep("Wagner V_8 (q=2)",mobius(4),2)
rep("prism C_5xK_2 (q=2)",prism(5),2)
rep("Moebius V_10 (q=2)",mobius(5),2)
rep("Petersen (q=2)",petersen(),2)
rep("Heawood (q=2)",heawood(),2)
rep("octahedron (q=3)",octa(),3)
rep("C_7^2 (q=3)",cyc(7,2),3)
rep("C_8^2 (q=3)",cyc(8,2),3)
