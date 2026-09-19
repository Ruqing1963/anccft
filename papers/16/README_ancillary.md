# Ancillary bundle for Paper XVI (`anccft16.tex`)

`head_correspondence.py` reproduces every computational claim of the paper on the Iwahori path tower over
K4 (q = 2), levels k = 1, 2, 3 (shadow graphs on 2 n_k = 24, 48, 96 vertices). It imports
`../Paper 10/pgl3_building.py` for the exact p-adic Smith form. Requires `numpy`, `sympy`. Runs in seconds.
The archived output is `head_correspondence_output.txt`.

Conventions of Paper VIII: (CK1) s_e* s_e = p_{r(e)}, (CK2) p_v = sum_{s(e)=v} s_e s_e*; K_0(C*(E)) = coker(1 - A_E^t),
K_1 = ker(1 - A_E^t). Shadow graph Y^# has vertices v, v'; the arcs of Y; one arc v -> v'; c = q-1 return arcs
v' -> v; two loops at v'. Its vertex matrix is B = [[A, I],[cI, 2I]]. tau_*, eta_* are the tail and head
pushforwards of arcs; tau^*, eta^* the transposes (pullbacks); H_k = diag(eta_*, eta_*), T^k = diag(tau^*, tau^*).

| key | content |
|-----|---------|
| (O) | eta^# is a graph morphism, bijective on OUT-arcs (out-covering), q-to-1 (not bijective) on in-arcs; tau^# is the mirror in-covering |
| (R) | path reversal rho conjugates A_k to A_k^t and satisfies rho_k tau_* = eta_* rho_{k+1}: the (eta, Delta^t) theory is the reversal of Paper VIII's (tau, Delta) theory |
| (K) | H_k descends on the K_0 side: H_k (I - B_{k+1}^t) = (I - B_k^t) H_k (whereas diag(tau_*,tau_*) does NOT); H_k 1 = q1 (so [X_eta] = q[1], no unital homomorphism); K_1-map = q |
| (H) | round-trip Kasparov products on K_0: H_k T^k = diag(A_k^t), T^k H_k = diag(A_{k+1}^t) (the two Hecke operators); explicit integer witness diag(A^t - q) = (I-B^t) V N that both act as q |
| (B) | opposite-algebra (Bowen-Franks) mirror: diag(tau_*,tau_*) and diag(eta^*,eta^*) descend on coker(I-B); composites give A_k, A_{k+1} |
| (N) | eta_* on Tor coker Delta^t (the in-module transitions): surjective, kernel = the 2-torsion of order 2^{n_k(q-1)-1}; group structures; rank_{F2} A_{k+1} = n_k |
| (X) | UCT/Ext bookkeeping: Ext(K_0(O_{k+1}), K_1(O_k)) = Jac(Y_{k+1}) (KK-lift ambiguity); Ext(K_0(O_k), K_1(O_k)) = Jac(Y_k) |

The K4 tower data matches Papers V, VII, VIII: Jac(Y_1) = Z/2+Z/8+Z/24, Jac(Y_2) = (Z/2)^8+Z/4+Z/16+Z/48,
Jac(Y_3) = (Z/2)^12+(Z/4)^8+Z/8+Z/32+Z/96.
