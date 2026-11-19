# Binding interaction graph, TACE-AS style (4 ligand x 6 receptor = 24 vertices)
import itertools, numpy as np
from sklearn.manifold import MDS
# ligand distances from Fig 1D (Angstrom). Assumption: Fig 2/4 ligand graph = D1, A1, H1, H2 subset
L = ["D1", "A1", "H1", "H2"]
dL = np.array([[0, 4.6, 9.1, 9.9], [4.6, 0, 8.1, 8.4], [9.1, 8.1, 0, 1.2], [9.9, 8.4, 1.2, 0]])
R = ["d1", "d2", "a1", "a2", "a3", "h1"]
# receptor distances are NOT in the main paper (supplement) -> made up: ligand 3D shape + noise + decoys
rng = np.random.default_rng(3)
pl = MDS(n_components=3, dissimilarity="precomputed", random_state=0).fit_transform(dL)
rec = np.array([pl[1] + rng.normal(0, .3, 3),   # d1 sits where A1 is  (donor facing acceptor)
                pl[1] + rng.normal(0, 1.5, 3),  # d2 decoy near A1
                pl[0] + rng.normal(0, .3, 3),   # a1 where D1 is
                pl[0] + rng.normal(0, 2.0, 3),  # a2 decoy
                pl[0] + rng.normal(0, 3.0, 3),  # a3 decoy
                (pl[2] + pl[3]) / 2 + rng.normal(0, .3, 3)])  # h1 between H1,H2
dR = np.linalg.norm(rec[:, None] - rec[None], axis=2)
# made-up contact potential: complementary labels strong, others weak
good = {("D1", "a"), ("A1", "d"), ("H1", "h"), ("H2", "h")}
def kappa(l, r): return 1.0 if (l, r[0]) in good else 0.2
V = [(i, j) for i in range(4) for j in range(6)]           # vertex = (ligand point, receptor point)
w = [kappa(L[i], R[j]) for i, j in V]
tol = 1.0   # stands in for tau + 2*epsilon
edges = set()
for a, b in itertools.combinations(range(len(V)), 2):
    (i, j), (k, l) = V[a], V[b]
    if i != k and j != l and abs(dL[i, k] - dR[j, l]) <= tol: edges.add((a, b))
if __name__ == "__main__":
    from clique_vqe import brute_force
    print("vertices", len(V), "edges", len(edges), "of", len(V)*(len(V)-1)//2, "possible")
    print("best clique (brute force):", brute_force(w, edges))
    best = brute_force(w, edges)[1]; print([(L[V[a][0]], R[V[a][1]]) for a in best])
