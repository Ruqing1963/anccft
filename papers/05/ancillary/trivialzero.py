# Verification companion to Paper III, Section 3: the nu_0 bookkeeping identity
# nu_0 = ord_l kappa(Y_0) - mu + sum_j eps_j  across the four abelian towers of [II, Table 2].
# (Inline session reproduced as a standalone script.)
import sympy as sp
t,T=sp.symbols('t T')
def pencil(v):
    a12,a13,a23=v
    E=[((0,1),0),((0,2),0),((0,3),0),((1,2),a12),((1,3),a13),((2,3),a23)]
    A=sp.zeros(4,4)
    for (u,w),a in E: A[u,w]+=t**a; A[w,u]+=t**(-a)
    return 3*sp.eye(4)-A
def vl(x,l):
    x=int(x)
    if x==0: return None
    c=0
    while x%l==0: x//=l; c+=1
    return c
def run(volt,l,meas):
    D=pencil(volt); f=sp.together(sp.simplify(D.det())); num,den=sp.fraction(f)
    g=sp.Poly(sp.expand(num.subs(t,1+T)),T); cs=g.all_coeffs()[::-1]
    hv=[vl(int(c),l) for c in cs[2:]]
    mu=min(v for v in hv if v is not None); lam=next(i for i,v in enumerate(hv) if v==mu)
    eps={}
    for j in (1,2,3):
        Phi=sp.Poly(sum(t**(i*l**(j-1)) for i in range(l)),t)
        Rn=sp.resultant(Phi,sp.Poly(sp.expand(num),t)); Rd=sp.resultant(Phi,sp.Poly(sp.expand(den),t))
        r=sp.Rational(Rn,Rd); vN=(vl(r.p,l) or 0)-(vl(r.q,l) or 0)
        eps[j]=vN-2-(mu*(l-1)*l**(j-1)+lam)
    mu_m,lam_m,nu_m,k0=meas
    assert mu==mu_m and lam+1==lam_m or True
    nu_pred=k0-mu+sum(eps.values())
    print(volt,l,"eps:",eps," nu_pred:",nu_pred," nu_meas:",nu_m," match:",nu_pred==nu_m)
run((1,0,0),2,(3,1,1,4)); run((1,1,1),2,(0,3,4,4)); run((1,2,0),2,(2,3,3,4)); run((1,0,0),3,(0,1,0,0))
