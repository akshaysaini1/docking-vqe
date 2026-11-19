# Mapping: maximum weighted clique  ->  QUBO  ->  Ising Hamiltonian  ->  VQE
import itertools, time, numpy as np
from qiskit_optimization import QuadraticProgram
from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import Statevector
from qiskit_aer.primitives import EstimatorV2 as AerEstimator
from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA

def build_qubo(weights, edges, penalty):
    """x_i = 1 if vertex i is in the clique.
       minimise  -sum w_i x_i  +  P * sum_{i<j, (i,j) NOT an edge} x_i x_j
       (two vertices that are not connected must not both be chosen)"""
    n = len(weights)
    qp = QuadraticProgram()
    for i in range(n): qp.binary_var(f"x{i}")
    lin = {f"x{i}": -weights[i] for i in range(n)}
    quad = {(f"x{i}", f"x{j}"): penalty for i, j in itertools.combinations(range(n), 2) if (i, j) not in edges}
    qp.minimize(linear=lin, quadratic=quad)
    return qp

def brute_force(weights, edges):
    n = len(weights); best = (0, ())
    for r in range(1, n + 1):
        for S in itertools.combinations(range(n), r):
            if all((a, b) in edges for a, b in itertools.combinations(S, 2)):
                w = sum(weights[i] for i in S)
                if w > best[0]: best = (w, S)
    return best

def run_vqe(qp, reps=2, maxiter=400, seed=1):
    op, offset = qp.to_ising()
    ans = RealAmplitudes(op.num_qubits, reps=reps, entanglement="linear")
    est = AerEstimator(options={"backend_options": {"method": "statevector"}, "default_precision": 0.0})
    from qiskit import transpile
    ans_t = transpile(ans, basis_gates=["u", "cx", "rz", "sx", "x", "ry"], optimization_level=1)
    rng = np.random.default_rng(seed)
    x0 = rng.uniform(-np.pi, np.pi, ans.num_parameters)
    hist = []
    vqe = VQE(est, ans_t, COBYLA(maxiter=maxiter), initial_point=x0, callback=lambda n, p, v, m: hist.append(v + offset))
    t = time.time(); res = vqe.compute_minimum_eigenvalue(op); dt = time.time() - t
    # IMPORTANT: the answer is read from the final state's bitstrings, not from the energy number
    sv = Statevector(ans.assign_parameters(res.optimal_point))
    probs = sv.probabilities_dict()
    return res.eigenvalue.real + offset, probs, dt, ans.num_parameters, len(hist)

def decode(bitstr):  # qiskit bitstrings are little-endian
    return tuple(i for i, b in enumerate(reversed(bitstr)) if b == "1")
