import numpy as np, sympy as sp
from sympy.matrices.normalforms import smith_normal_form
for k in (1,2):
    Do=np.load(f"Do_{k}.npy"); Di=np.load(f"Di_{k}.npy")
    dd=sp.Matrix((Do-Di).tolist())
    ker=dd.T.nullspace()
    Eker = len(ker)==1 and all(x==ker[0][0] for x in ker[0]) and ker[0][0]!=0
    S=smith_normal_form(dd)
    inv=[abs(S[i,i]) for i in range(min(S.rows,S.cols))]
    Ecok = inv.count(1)==dd.rows-1 and all(x in (0,1) for x in inv)
    print(f"k={k}: (E-ker) ker(tau^*-eta^*)=Z*1: {Eker}   (E-cok) SNF(tau_*-eta_*) units={inv.count(1)} of n_k-1={dd.rows-1}, torsion-free: {Ecok}")
