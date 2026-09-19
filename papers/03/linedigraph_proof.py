# Verification companion to Paper III, Theorem 2.1 (line-digraph recursion):
# q=3 instance (complete digraph on 4 vertices) and the AB/BA spectral identity.
import sympy as sp
n=4; A=sp.ones(n,n)-sp.eye(n); q=3
kG=int((q*sp.eye(n)-A)[1:,1:].det())
arcs=[(u,v) for u in range(n) for v in range(n) if u!=v]; m=len(arcs)
AL=sp.zeros(m,m)
for i,(u,v) in enumerate(arcs):
    for j,(x,y) in enumerate(arcs):
        if x==v: AL[i,j]=1
kL=int((q*sp.eye(m)-AL)[1:,1:].det())
assert kL == q**(n*(q-1)-1)*kG == 34992
print("q=3 check passed: kappa(L) =", kL, "= 3^7 * 16")
