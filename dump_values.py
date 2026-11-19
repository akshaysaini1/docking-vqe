import warnings, json, itertools; warnings.filterwarnings("ignore")
import numpy as np
from docking_graph import *
from clique_vqe import brute_force
best = brute_force(w, edges)
# 8-vertex random toy graph (toy2.py)
rng = np.random.default_rng(0); n = 8
w8 = [round(float(x), 2) for x in rng.uniform(0.3, 1.5, n)]
e8 = set(); planted = (0, 2, 3, 5)
for a, b in itertools.combinations(planted, 2): e8.add((a, b))
for a, b in itertools.combinations(range(n), 2):
    if rng.random() < 0.35: e8.add((a, b))
out = dict(L=L, R=R, dL=dL.tolist(), dR=dR.round(2).tolist(), lig_xyz=pl.round(2).tolist(), rec_xyz=rec.round(2).tolist(),
           V=[(L[i], R[j]) for i, j in V], w=w, tol=tol, edges=sorted(edges), best=best,
           w8=w8, e8=sorted(e8), best8=brute_force(w8, e8))
json.dump(out, open("values.json", "w"))
print("edges", len(edges), "best", best, [out["V"][a] for a in best[1]])
print("rec_xyz", out["rec_xyz"]); print("dR", out["dR"])
print("w8", w8, "e8", sorted(e8), brute_force(w8, e8))
# distances of receptor points vs ligand points of the best clique, for the check
