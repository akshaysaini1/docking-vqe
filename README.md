# Molecular docking as a maximum weighted clique, mapped to VQE

Beginner-level attempt to map the docking problem of
Banchi et al., *Molecular docking with Gaussian Boson Sampling*, Science Advances 6, eaax1950 (2020)
to a VQE problem in Qiskit.

The code was written with help from Claude (AI assistant) and was run and checked by the author.

## Idea
1. Vertex = a possible contact (ligand point, receptor point). Ligand 4 points x receptor 6 points = 24 vertices.
2. Edge = two contacts that are geometrically compatible. Weight = how good the contact is.
3. Best docking pose = maximum weighted clique.
4. Mapping to VQE: one qubit per vertex, x_i = 1 if vertex i is chosen.
   Cost = - sum_i w_i x_i + P * sum_{(i,j) not an edge} x_i x_j, with P = 2.0.
   Substituting x_i = (1 - Z_i)/2 gives an Ising Hamiltonian. VQE (RealAmplitudes + COBYLA) finds its ground state.

## Files
| File | What it does |
|---|---|
| `clique_vqe.py` | The mapping (clique -> QUBO -> Ising) and the VQE run |
| `docking_graph.py` | Builds the 24-vertex toy docking graph |
| `toy2.py` | Random 8-vertex test graph, 10 random starts |
| `final_runs.py` | 7-vertex reduced graph (10 starts) and time of one energy evaluation vs number of qubits |
| `dump_values.py` | Saves all toy model values to `values.json` |
| `docs/Toy_Model_Values.pdf` | Every value of the toy model and its source |
| `docs/Docking_VQE_Findings.pptx` | Presentation of the findings |

## Run
```
pip install -r requirements.txt
python docking_graph.py     # builds the graph, prints the exact best clique
python toy2.py              # VQE on the 8-vertex test graph
python final_runs.py        # 7-vertex VQE + timing
```

## Which values come from the paper and which are made up
**From the paper:** 4 ligand points (D1, A1, H1, H2) and 6 receptor points (d1, d2, a1, a2, a3, h1) with their labels (p.5, Fig. 2);
24 vertices (p.2); ligand distances (Fig. 1D, p.2, assuming the Fig. 1 ligand is the 4-point ligand of Fig. 2/4); edge rule |d_ligand - d_receptor| <= tau + 2 epsilon (p.3);
allowed label pairs (p.5); goal = maximum weighted clique (p.3).

**MADE UP (not in the paper, real values are in the supplementary material):**
- receptor 3D positions (built from the ligand shape plus random noise, seed 3, decoy points added)
- contact weights (1.0 for complementary labels, 0.2 otherwise; real potential is table S1)
- tolerance tau + 2 epsilon = 1.0 Angstrom (real values in table S2)
- the 8-vertex random test graph (seed 0)

**My choices (the paper does not use VQE):** penalty P = 2.0, RealAmplitudes reps=2, COBYLA, random starting angles (seeds 0-9), Qiskit Aer exact statevector simulator.

Because of the made-up values this toy graph is **not** the paper's graph, and results must not be compared with the paper's numbers
(the paper's best clique has 8 vertices and weight about 3.99; the toy graph's best clique has 3 vertices and weight 3.0, and was planted on purpose).

## Results
| Test | Result |
|---|---|
| 8-vertex random graph, reps=2 | best clique found in 7/10 starts |
| 8-vertex random graph, reps=4 | best clique found in 5/10 starts |
| 7-vertex reduced docking graph | best clique found in 10/10 starts |
| 24-vertex docking graph | not solved with VQE: about 31 s per energy evaluation, a full run needs hours |

Time of one energy evaluation (simple test Hamiltonian for 10-22 qubits, docking Hamiltonian for 24 qubits): 10q 0.01 s, 14q 0.01 s, 18q 0.11 s, 20q 0.47 s, 22q 2.35 s, 24q 30.8 s.
A classical computer finds the exact answer for 24 vertices in seconds, so VQE gives no advantage at this size.

## Not done
- The paper's real graph (needs the supplementary material).
- A full 24-qubit VQE run.
- Comparison with QAOA or the paper's classical methods (greedy shrinking, local search).
