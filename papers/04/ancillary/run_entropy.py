import entropy as E, numpy as np, math
rows = []

# --- genuine LPS Ramanujan graphs X^{5,l}, degree 6 ---
for l in (13, 17):
    A = E.lps(5, l)
    rows.append(E.stats(A, 6, f"LPS X^(5,{l})"))
    print("done LPS", l, A.shape, flush=True)

# --- exactly-Ramanujan but NOT tree-like: incidence graph of PG(2,5), degree 6 ---
rows.append(E.stats(E.incidence_pg(5), 6, "PG(2,5) incidence"))

# --- random regular, degree 6 and 3 ---
for d in (6, 3):
    for n in (200, 500, 1000, 2000):
        if n*d % 2: continue
        A = E.random_regular(n, d, seed=n+d)
        if A is None: print("skip", n, d); continue
        rows.append(E.stats(A, d, f"random {d}-reg n={n}"))
        print("done random", d, n, flush=True)

hdr = f"{'family':22s} {'n':>5s} {'d':>2s} {'lam2':>7s} {'Ram':>4s} {'KS':>6s} {'(1/n)log|Jac|':>14s} {'h_d':>8s} {'bound':>8s} {'emp/h':>7s} {'emp/bd':>7s}"
print("\n"+hdr); print("-"*len(hdr))
for r in sorted(rows, key=lambda r: (r['d'], r['n'])):
    h, b = E.h_tree_closed(r['d']), E.bound(r['d'])
    print(f"{r['name']:22s} {r['n']:5d} {r['d']:2d} {r['lam2']:7.4f} "
          f"{'yes' if r['ram'] else 'NO':>4s} {r['ks']:6.4f} {r['ent']:14.6f} "
          f"{h:8.6f} {b:8.6f} {r['ent']/h:7.4f} {r['ent']/b:7.4f}")
