import warnings, time; warnings.filterwarnings("ignore")
from docking_graph import *
from clique_vqe import *
from qiskit import transpile
# (a) keep only complementary-label contacts (weight 1.0) -> 7 vertices
keep = [a for a in range(len(V)) if w[a] == 1.0]
idx = {v: i for i, v in enumerate(keep)}
w2 = [w[v] for v in keep]; e2 = {(idx[a], idx[b]) for a, b in edges if a in idx and b in idx}
opt = brute_force(w2, e2); print("7-vertex graph: edges", len(e2), "optimum", opt, flush=True)
qp = build_qubo(w2, e2, 2.0); hits = 0
for seed in range(10):
    E, probs, dt, npar, nev = run_vqe(qp, reps=2, maxiter=600, seed=seed)
    S = decode(max(probs, key=probs.get)); ok = all((a, c) in e2 for a, c in itertools.combinations(S, 2))
    wt = sum(w2[i] for i in S) if ok else 0; hits += abs(wt - opt[0]) < 1e-9
print(f"7 qubits: optimum found in {hits}/10 starts", flush=True)
# (b) cost of ONE energy evaluation vs number of qubits
est = AerEstimator(options={"backend_options": {"method": "statevector"}, "default_precision": 0.0})
for n in [10, 14, 18, 20, 22]:
    from qiskit.quantum_info import SparsePauliOp
    op = SparsePauliOp.from_sparse_list([("ZZ", [i, (i + 1) % n], 1.0) for i in range(n)] + [("Z", [i], 0.5) for i in range(n)], num_qubits=n)
    a = transpile(RealAmplitudes(n, reps=2), basis_gates=["u", "cx", "rz", "sx", "x", "ry"], optimization_level=1)
    p = np.random.rand(a.num_parameters); est.run([(a, op, p)]).result()
    t = time.time(); est.run([(a, op, p)]).result(); print(f"{n} qubits: {time.time()-t:.2f} s per energy evaluation", flush=True)
