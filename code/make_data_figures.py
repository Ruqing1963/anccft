# -*- coding: utf-8 -*-
"""
make_data_figures.py -- write the repository's CSV data files and PDF figures from the
verified numbers of Volume I (every number below is reproduced by a paper script in this
repository, and is copied from that script's archived *_output.txt; this file only
tabulates and plots).

    python make_data_figures.py        # writes ../data/*.csv and ../figures/*.pdf
"""
import os, csv, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
FIGS = os.path.join(os.path.dirname(HERE), "figures")
os.makedirs(DATA, exist_ok=True)
os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({"font.size": 9, "figure.dpi": 150, "axes.spines.top": False, "axes.spines.right": False})


def write_csv(name, header, rows):
    with open(os.path.join(DATA, name), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print("  wrote", name, f"({len(rows)} rows)")


# ---------------------------------------------------------------- K4 Iwahori tower (Chapters 13, 7, 12 = Papers VIII, IX, XIII)
k4 = [
    # level, n_k (arcs), kappa(Y_k) (spanning trees), ord_2 kappa, Jac(Y_k), coker(Delta_k^0), det Delta_k^0
    (1, 12, 2**7 * 3, 7, "Z/2+Z/8+Z/24", "Z/8+(Z/24)^2", 4608),
    (2, 24, 2**18 * 3, 18, "(Z/2)^8+Z/4+Z/16+Z/48", "(Z/2)^9+Z/16+(Z/48)^2", 18874368),
    (3, 48, 2**41 * 3, 41, "(Z/2)^12+(Z/4)^8+Z/8+Z/32+Z/96", "(Z/2)^12+(Z/4)^9+Z/32+(Z/96)^2", 316659348799488),
]
write_csv("k4_iwahori_tower.csv",
          ["level", "n_k", "kappa", "ord2_kappa", "Jac", "coker_restricted_laplacian", "det_restricted_laplacian"], k4)

# abelian voltage towers over K4 (Chapter 7 = Paper IX): cleared pencil, L-invariant, Iwasawa (mu, lambda), defects
volt = [
    ("(1,0,0)", 2, "-8T^2", "-1/2", 3, 0, "0,0"),
    ("(1,1,1)", 2, "-T^4-16T^3-16T^2", "-1", 0, 2, "0,0"),
    ("(1,2,0)", 2, "-4T^4-24T^3-24T^2", "-3/2", 2, 2, "-1,2"),
    ("(1,0,0)", 3, "-8T^2", "-1/2", 0, 0, "0,0"),
]
write_csv("k4_voltage_towers.csv", ["voltage", "l", "cleared_pencil", "L_invariant", "mu", "lambda", "defects_eps1_eps2"], volt)

# ---------------------------------------------------------------- composite q (Chapter 14 = Paper XVII, composite_tail_norm.py)
comp = [
    # tower, q, level k+1, n_{k+1}, prime p | q, p-rank of Jac(Y_{k+1}), r_k = n_k(q-1)-1, min p-adic valuation of the
    # invariant factors divisible by p (= v_p(q), the divisibility statement)
    ("K4", 2, 2, 24, 2, 11, 11, 1), ("K4", 2, 3, 48, 2, 23, 23, 1),
    ("Petersen", 2, 2, 60, 2, 29, 29, 1), ("Petersen", 2, 3, 120, 2, 59, 59, 1),
    ("K5", 3, 2, 60, 3, 39, 39, 1), ("K5", 3, 3, 180, 3, 119, 119, 1),
    ("K4,4", 3, 2, 96, 3, 63, 63, 1),
    ("K6", 4, 2, 120, 2, 89, 89, 2), ("K6", 4, 3, 480, 2, 359, 359, 2),
    ("K5,5", 4, 2, 200, 2, 149, 149, 2),
    ("K7", 5, 2, 210, 5, 167, 167, 1),
    ("K8", 6, 2, 336, 2, 279, 279, 1), ("K8", 6, 2, 336, 3, 279, 279, 1),
]
write_csv("composite_q_towers.csv",
          ["tower", "q", "level", "n_level", "p", "p_rank_Jac", "r_k", "min_valuation_vp"], comp)

# ---------------------------------------------------------------- A~_2 complexes (Chapters 16-19 = Papers X-XII, XV)
a2 = [
    ("X_3", 3, "Z/2+Z/2+Z/14", "triangle presentation over PG(2,2), Singer cycle"),
    ("Y_21", 21, "(Z/2)^6", "abelian cover, Heawood links, K_{7,7,7}"),
    ("Y_24", 24, "(Z/2)^4+Z/14", "abelian cover, Heawood links"),
    ("Y_42", 42, "(Z/2)^3+(Z/4)^2", "abelian cover, Heawood links"),
    ("Y_84", 84, "", "abelian cover"),
    ("Y_168", 168, "", "abelian cover; edge flow carries all 498 nontrivial pencil roots, 678 remainder eigenvalues |.|=sqrt q"),
]
write_csv("a2_complexes.csv", ["complex", "vertices", "H1", "notes"], a2)

# the geodesic edge flow on the three thick complexes (Chapters 20, 21 = Papers XVIII, XX)
edge = [
    # complex, n, |E|, q, singular values equal to sqrt q (= |E|-n), rank Psi_2 (of 3n), corank, dim W_2 = |E|-3n+6
    ("Y_21", 21, 147, 2, 126, 57, 6, 90), ("Y_24", 24, 168, 2, 144, 66, 6, 102), ("Y_42", 42, 294, 2, 252, 120, 6, 174),
]
write_csv("a2_edge_flow.csv",
          ["complex", "n", "edges", "q", "singular_values_sqrt_q", "rank_Psi", "corank_Psi", "dim_W"], edge)

# PG(d,q): the point/hyperplane non-incidence identity (Chapter 20, check H) and the flag count (Chapter 21, check F)
pg = [(2, 2, 7, 4, 2), (2, 3, 13, 9, 6), (2, 5, 31, 25, 20), (2, 7, 57, 49, 42), (2, 11, 133, 121, 110),
      (3, 2, 15, 8, 4), (3, 3, 40, 27, 18), (3, 5, 156, 125, 100), (3, 7, 400, 343, 294),
      (4, 2, 31, 16, 8), (4, 3, 121, 81, 54), (4, 5, 781, 625, 500), (5, 2, 63, 32, 16), (5, 3, 364, 243, 162)]
write_csv("pg_unitarity.csv", ["d", "q", "points", "BBt_diagonal", "BBt_offdiagonal"], pg)
flags = [(2, 2, 0, 3, 0, 3), (2, 2, 1, 3, 1, 2), (2, 3, 0, 4, 0, 4), (2, 3, 1, 4, 1, 3),
         (3, 2, 0, 21, 0, 21), (3, 2, 1, 9, 3, 6), (3, 2, 2, 21, 9, 12),
         (3, 3, 0, 52, 0, 52), (3, 3, 1, 16, 4, 12), (3, 3, 2, 52, 16, 36),
         (4, 2, 0, 315, 0, 315), (4, 2, 1, 63, 21, 42), (4, 2, 2, 63, 27, 36), (4, 2, 3, 315, 147, 168)]
write_csv("pg_flag_counts.csv", ["d", "q", "i", "N_V_in_H", "N_V_not_in_H", "K_i"], flags)

# ---------------------------------------------------------------- the odd cyclic class (Chapter 9 = Paper XIX)
cyc = [("K4 (1,0,0)", 3, 16, "1/2", 8, 2), ("K4 (1,1,1)", 3, 16, "1", 16, 2), ("K4 (1,2,0)", 3, 16, "3/2", 24, 2),
       ("K5 single edge", 6, 125, "3/5", 75, 2), ("K5 triangle", 6, 125, "7/5", 175, 2)]
write_csv("cyclic_obstruction.csv",
          ["voltage_class", "b1", "kappa0", "h_norm_squared", "kappa0_h_norm_squared", "winding_number_of_det_D"], cyc)

# ---------------------------------------------------------------- Holo 1 / Net 1 small tables (Chapters 23, 24)
write_csv("holo_witnesses.csv", ["quantity", "witness_1", "witness_2"],
          [("K_0 of boundary crossed product", "Z^r+Z/(r-1)", "Z^r+Z/(r-1)"), ("type", "same", "same"),
           ("|Jac|", 256, 363), ("coprime", "yes", "")])
write_csv("net_lambda_signs.csv", ["tower", "lambda", "reading"],
          [("abelian voltage (1,1,1) l=2", 2, "normal, abelian"), ("dihedral", 3, "normal, non-abelian: sign does not detect commutativity"),
           ("generic non-normal", "<0", "certifies non-normal"), ("lambda = -1", -1, "telescopes from Knuth 1967")])

# ================================================================= figures
# 1. K4 Iwahori tower: 2-adic growth of kappa and of det(restricted Laplacian)
fig, ax = plt.subplots(figsize=(4.2, 2.8))
lv = [r[0] for r in k4]
ax.plot(lv, [r[3] for r in k4], "o-", color="#3182bd", label="$\\mathrm{ord}_2\\,\\kappa(Y_k)$")
ax.plot(lv, [r[1] - 3 for r in k4], "s--", color="#e6550d", label="$n_k-3=\\mathrm{ord}_2\\det\\Delta_k^0$")
ax.set_xlabel("level $k$")
ax.set_xticks(lv)
ax.set_title("$K_4$ Iwahori tower ($q=2$): pure arc growth")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "k4_tower_growth.pdf"))
print("  wrote k4_tower_growth.pdf")

# 2. composite q: the p-rank of Jac(Y_{k+1}) equals r_k = n_k(q-1)-1 for every prime dividing q
fig, ax = plt.subplots(figsize=(4.2, 2.8))
rk = [r[6] for r in comp]
pr = [r[5] for r in comp]
ax.loglog(rk, pr, "o", color="#3182bd")
ax.loglog([min(rk), max(rk)], [min(rk), max(rk)], "--", color="#e6550d", lw=0.8, label="$p$-rank $=r_k$")
for r in comp:
    if r[7] == 2:
        ax.annotate(f"{r[0]}, $q=4$", (r[6], r[5]), fontsize=6, xytext=(4, -8), textcoords="offset points")
ax.set_xlabel("$r_k=n_k(q-1)-1$")
ax.set_ylabel("$p$-rank of $\\mathrm{Jac}(Y_{k+1})$")
ax.set_title("Composite $q$: $\\ker\\tau_*=\\mathrm{Jac}(Y_{k+1})[q]$")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "composite_q_prank.pdf"))
print("  wrote composite_q_prank.pdf")

# 3. the edge flow of the thick A~_2 complexes: |E| = rank Psi (Hecke part) + dim W (unitary remainder)
fig, ax = plt.subplots(figsize=(4.2, 2.8))
names = [r[0] for r in edge]
hecke = [r[5] for r in edge]
unit = [r[7] for r in edge]
ax.bar(names, hecke, color="#3182bd", label="Hecke part $\\operatorname{rank}\\Psi_2$")
ax.bar(names, unit, bottom=hecke, color="#9ecae1", label="unitary remainder $\\dim W_2$")
for i, r in enumerate(edge):
    ax.text(i, r[2] + 4, f"$|E|={r[2]}$", ha="center", fontsize=7)
ax.set_ylabel("dimension")
ax.set_title("$\\widetilde A_2$ edge flow: Hecke part and critical-circle part")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "a2_edge_flow_split.pdf"))
print("  wrote a2_edge_flow_split.pdf")

# 4. Ctrl 1: the power-sum bound s >= ceil(max_k [-p_k]_+ / rho^k) on the two witnesses of the
#    note -- M_3(0.9) = {1} u 3{+-0.9i} (Thm 3.1) and the tutorial poles {1, +-0.9i} (Rem 4.3):
#    k = 1 gives nothing (Benvenuti's zeta = 0), k = 2 binds.
fig, ax = plt.subplots(figsize=(4.2, 2.8))
kmax = 8
for lam, lab, col, mk in ((([1.0] + [0.9j, -0.9j] * 3), "$\\mathcal{M}_3(0.9)$: bound $3.86$, so $s\\geq4$", "#3182bd", "o"),
                          (([1.0, 0.9j, -0.9j]), "poles $\\{1,\\pm0.9i\\}$: bound $0.62$, so $s\\geq1$", "#e6550d", "s")):
    rho = max(abs(l) for l in lam)
    bounds = [max(-sum(l ** k for l in lam).real, 0.0) / rho ** k for k in range(1, kmax + 1)]
    ax.plot(range(1, kmax + 1), bounds, mk + "-", color=col, label=lab)
ax.set_xlabel("$k$")
ax.set_ylabel("$[-p_k]_+/\\rho^k$")
ax.set_xticks(range(1, kmax + 1))
ax.set_title("Power-sum bound: $k=2$ binds, $k=1$ sees nothing")
ax.legend(frameon=False, fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "ctrl_powersum_bound.pdf"))
print("  wrote ctrl_powersum_bound.pdf")
print("done")
