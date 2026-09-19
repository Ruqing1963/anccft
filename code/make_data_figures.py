# -*- coding: utf-8 -*-
"""
make_data_figures.py -- write the repository's CSV data files and PDF figures from the
verified numbers of the ANCCFT series (every number below is reproduced by the paper
scripts in code/; this file only tabulates and plots them).

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


# ---------------------------------------------------------------- K4 Iwahori tower (Papers VIII, IX, XIII)
k4 = [
    # level, n_k (arcs), kappa(Y_k) (spanning trees), ord_2 kappa, Jac(Y_k), coker(Delta_k^0), det Delta_k^0
    (1, 12, 2**7 * 3, 7, "Z/2+Z/8+Z/24", "Z/8+(Z/24)^2", 4608),
    (2, 24, 2**18 * 3, 18, "(Z/2)^8+Z/4+Z/16+Z/48", "(Z/2)^9+Z/16+(Z/48)^2", 18874368),
    (3, 48, 2**41 * 3, 41, "(Z/2)^12+(Z/4)^8+Z/8+Z/32+Z/96", "(Z/2)^12+(Z/4)^9+Z/32+(Z/96)^2", 316659348799488),
]
write_csv("k4_iwahori_tower.csv",
          ["level", "n_k", "kappa", "ord2_kappa", "Jac", "coker_restricted_laplacian", "det_restricted_laplacian"], k4)

# abelian voltage towers over K4 (Paper IX): cleared pencil, L-invariant, Iwasawa (mu, lambda), defects
volt = [
    ("(1,0,0)", 2, "-8T^2", "-1/2", 3, 0, "0,0"),
    ("(1,1,1)", 2, "-T^4-16T^3-16T^2", "-1", 0, 2, "0,0"),
    ("(1,2,0)", 2, "-4T^4-24T^3-24T^2", "-3/2", 2, 2, "-1,2"),
    ("(1,0,0)", 3, "-8T^2", "-1/2", 0, 0, "0,0"),
]
write_csv("k4_voltage_towers.csv", ["voltage", "l", "cleared_pencil", "L_invariant", "mu", "lambda", "defects_eps1_eps2"], volt)

# ---------------------------------------------------------------- A~_2 complexes (Papers X-XII, XV)
a2 = [
    ("X_3", 3, "Z/2+Z/2+Z/14", "triangle presentation over PG(2,2), Singer cycle"),
    ("Y_21", 21, "(Z/2)^6", "abelian cover, Heawood links, K_{7,7,7}"),
    ("Y_24", 24, "(Z/2)^4+Z/14", "abelian cover, Heawood links"),
    ("Y_42", 42, "(Z/2)^3+(Z/4)^2", "abelian cover, Heawood links"),
    ("Y_84", 84, "", "abelian cover"),
    ("Y_168", 168, "", "abelian cover; edge flow carries all 498 nontrivial pencil roots, 678 remainder eigenvalues |.|=sqrt q"),
]
write_csv("a2_complexes.csv", ["complex", "vertices", "H1", "notes"], a2)

# ---------------------------------------------------------------- phiX174 sandpile spectrum (Bio 1-4)
phix = [
    # k, branch vertices of compacted graph, chains, skeleton vertices, K(skeleton), length-driven part, log10|K|
    (12, 2, 1, 1, "0", "Z/2", math.log10(2)),
    (11, 10, 2, 8, "Z/3+Z/11", "(Z/2)^2", math.log10(4 * 33)),
    (10, 47, 8, 37, "Z/4+Z/1076338579", "(Z/2)^10", math.log10(2**10 * 4 * 1076338579)),
    (9, 146, 28, 110, "(Z/2)^6+(Z/4)^4+Z/7+Z/523+Z/3779+Z/5297+Z/6561+Z/7669+Z/2707831823", "(Z/2)^36",
     math.log10(2**36 * 2**6 * 4**4 * 7 * 523 * 3779 * 5297 * 6561 * 7669 * 2707831823)),
]
write_csv("phix174_sandpile_spectrum.csv",
          ["k", "branch_vertices", "bundle_chains", "skeleton_vertices", "K_skeleton", "length_driven", "log10_order_K"],
          [(k, b, c, s, K, L, round(l, 4)) for k, b, c, s, K, L, l in phix])

# skeleton degree profile (Bio 4)
prof = [(12, 1, 0, "0:1"), (11, 8, 16, "2:8"), (10, 37, 74, "2:37"), (9, 110, 227, "2:104,3:5,4:1"),
        (8, 312, 679, "2:264,3:42,4:5,5:1"), (7, 826, 1945, "2:619,3:147,4:41,5:12,6:7")]
write_csv("phix174_skeleton_profile.csv", ["k", "skeleton_vertices", "arcs", "outdegree_multiset"], prof)

# phiX174 k=12: the two reconstructions (Bio 2)
write_csv("phix174_k12_reconstruction.csv",
          ["segment", "start", "end", "length_nt", "gene_context"],
          [("A", 457, 468, 12, "repeated 12-mer CTTCTGCCGTTT (also 3205..3216)"),
           ("B", 264, 275, 12, "repeated 12-mer CGTCAAGGACTG (also 2760..2771)"),
           ("X", 3217, 263, 2433, "wrapping"), ("U", 276, 456, 181, "spans C/D junction"),
           ("Z", 469, 2759, 2291, ""), ("Y", 2772, 3204, 433, "spans G/H spike junction")])

# ---------------------------------------------------------------- Bio 4 census
census = [
    ("(3,3)", 4, 3, 3, 3, "", 3, 0), ("(4,4)", 10, 4, 4, 4, "", 7, 0),
    ("(3,3,3)", 188, 12, 9, 9, "(Z/2)^2;(Z/3)^2", 15, 6),
    ("(3,3,2)", 70, 15, 7, 6, "(Z/2)^2", 10, 4), ("(2,2,3)", 30, 10, 5, 4, "(Z/2)^2", 6, 2),
    ("(4,4,4)", 2896, 26, 17, 16, "(Z/2)^2;(Z/2)^2+Z/3;(Z/3)^2;(Z/4)^2;Z/2+Z/4", 174, 86),
    ("(3,3,3,3)", 30804, 96, 28, 27, "(Z/2)^2;(Z/2)^2+Z/3;(Z/2)^3;(Z/3)^2;(Z/3)^3;(Z/4)^2;Z/2+(Z/3)^2;Z/2+Z/4", 180, 132),
]
write_csv("multicopy_census.csv",
          ["multiplicities", "circular_words", "transition_digraphs", "groups", "max_order", "noncyclic_groups",
           "pairwise_classes", "pairwise_classes_mixed"], census)

# chord-diagram invariants, r = 2 (Bio 3 check D)
write_csv("chord_invariants_r2.csv", ["m", "words", "interlace_graphs", "graph_classes_with_more_than_one_group"],
          [(3, 30, 4, 0), (4, 630, 11, 0), (5, 22680, 34, 0), (6, "MMN pair 123456132465 / 123654132564", "same", "Z/3+Z/6 vs Z/18")])

# ---------------------------------------------------------------- Holo 1 / Net 1 / Ctrl 1 small tables
write_csv("holo_witnesses.csv", ["quantity", "witness_1", "witness_2"],
          [("K_0 of boundary crossed product", "Z^r+Z/(r-1)", "Z^r+Z/(r-1)"), ("type", "same", "same"),
           ("|Jac|", 256, 363), ("coprime", "yes", "")])
write_csv("net_lambda_signs.csv", ["tower", "lambda", "reading"],
          [("abelian voltage (1,1,1) l=2", 2, "normal, abelian"), ("dihedral", 3, "normal, non-abelian: sign does not detect commutativity"),
           ("generic non-normal", "<0", "certifies non-normal"), ("lambda = -1", -1, "telescopes from Knuth 1967")])

# ================================================================= figures
# 1. phiX174 spectrum: log10 |K| split into length-driven and arrangement parts
ks = [k for k, *_ in phix]
ld = [math.log10({12: 2, 11: 4, 10: 2**10, 9: 2**36}[k]) for k in ks]
tot = [l for *_, l in phix]
arr = [t - d for t, d in zip(tot, ld)]
fig, ax = plt.subplots(figsize=(4.2, 2.8))
ax.bar(ks, ld, color="#9ecae1", label="length-driven  $\\bigoplus(\\mathbb{Z}/r)^{\\ell-1}$")
ax.bar(ks, arr, bottom=ld, color="#3182bd", label="arrangement  $\\mathcal{K}(\\mathrm{Sk})$")
ax.set_xlabel("$k$")
ax.set_ylabel("$\\log_{10}|\\mathcal{K}(G_k)|$")
ax.set_title("$\\varphi$X174: sandpile spectrum, split by repeats")
ax.invert_xaxis()
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "phix174_spectrum.pdf"))
print("  wrote phix174_spectrum.pdf")

# 2. phiX174 skeleton: size and multi-copy vertices vs k
fig, ax = plt.subplots(figsize=(4.2, 2.8))
kk = [p[0] for p in prof]
sz = [p[1] for p in prof]
multi = []
for p in prof:
    d = dict((int(a), int(b)) for a, b in (x.split(":") for x in p[3].split(",")))
    multi.append(sum(v for dd, v in d.items() if dd >= 3))
ax.semilogy(kk, sz, "o-", color="#3182bd", label="skeleton vertices")
ax.semilogy(kk, [max(m, 0.5) for m in multi], "s--", color="#e6550d", label="vertices of out-degree $\\geq 3$")
ax.set_xlabel("$k$")
ax.invert_xaxis()
ax.set_title("$\\varphi$X174: where multi-copy repeats appear")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "phix174_skeleton_profile.pdf"))
print("  wrote phix174_skeleton_profile.pdf")

# 3. K4 Iwahori tower: 2-adic growth of kappa and of det(restricted Laplacian)
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

# 4. multi-copy census: number of distinct groups vs number of words
fig, ax = plt.subplots(figsize=(4.2, 2.8))
labels = [c[0] for c in census]
ax.bar(range(len(census)), [c[3] for c in census], color="#3182bd")
ax.set_xticks(range(len(census)))
ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=7)
ax.set_ylabel("distinct arrangement groups")
ax.set_title("Multi-copy repeats: census of $\\mathcal{K}(D_w)$")
for i, c in enumerate(census):
    ax.text(i, c[3] + 0.4, f"max {c[4]}", ha="center", fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "multicopy_census.pdf"))
print("  wrote multicopy_census.pdf")

# 5. Ctrl 1: the power-sum bound s >= ceil(max_k [-p_k]_+ / rho^k) is illustrated on a generic spectrum
fig, ax = plt.subplots(figsize=(4.2, 2.8))
import cmath
lam = [1.0, 0.9 * cmath.exp(2j * math.pi / 3), 0.9 * cmath.exp(-2j * math.pi / 3)]
rho = 1.0
kmax = 12
bounds = []
for k in range(1, kmax + 1):
    pk = sum(l ** k for l in lam).real
    bounds.append(max(-pk, 0) / rho ** k)
ax.plot(range(1, kmax + 1), bounds, "o-", color="#3182bd")
ax.axhline(max(bounds), color="#e6550d", ls="--", lw=0.8)
ax.set_xlabel("$k$")
ax.set_ylabel("$[-p_k]_+/\\rho^k$")
ax.set_title("Hidden-eigenvalue bound: the $k$ that binds need not be $k=1$")
fig.tight_layout()
fig.savefig(os.path.join(FIGS, "ctrl_powersum_bound.pdf"))
print("  wrote ctrl_powersum_bound.pdf")
print("done")
