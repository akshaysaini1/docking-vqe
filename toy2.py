from clique_vqe import *
rng = np.random.default_rng(0)
n = 8
w = [round(float(x), 2) for x in rng.uniform(0.3, 1.5, n)]
edges = set(); planted = (0, 2, 3, 5)
for a, b in itertools.combinations(planted, 2): edges.add((a, b))
for a, b in itertools.combinations(range(n), 2):
    if rng.random() < 0.35: edges.add((a, b))
opt = brute_force(w, edges); print("optimum", opt)
qp = build_qubo(w, edges, 2.0)
for reps in [2, 4]:
    hits = 0; res_list = []
    for seed in range(10):
        E, probs, dt, npar, nev = run_vqe(qp, reps=reps, maxiter=1000, seed=seed)
        b = max(probs, key=probs.get); S = decode(b)
        ok = all((a, c) in edges for a, c in itertools.combinations(S, 2))
        wt = sum(w[i] for i in S) if ok else 0
        hits += abs(wt - opt[0]) < 1e-9; res_list.append((S, round(wt, 2)))
    print(f"reps={reps}: found optimum in {hits}/10 starts; results {res_list}", flush=True)
