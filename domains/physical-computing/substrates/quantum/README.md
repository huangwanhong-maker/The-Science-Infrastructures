# quantum: Quantum Reservoirs

*Status: planned (ADR 0076).* Mode: reservoir. **Breadth:** the open dynamics of a few interacting spins, whose state space grows as $4^n$. Gate-based quantum processors are calx's backends and heuresis's samplers; this family holds quantum systems used as reservoirs.

## 1. The Substrates

| Substrate | What it shows | Materials and reuse |
|---|---|---|
| Disordered spin ensembles (Fujii and Nakajima 2017): nuclear spins in NMR, Rydberg arrays, superconducting qubits | a reservoir of $2^n$ amplitudes; memory from the injection of inputs | NMR at room temperature; superconducting qubits cryogenic; Rydberg arrays need lasers and vacuum |

## 2. The Simulators

1. **A transverse-field Ising reservoir**, of $n \le 8$ qubits. The Hamiltonian $H = \sum_{ij} J_{ij}X_iX_j + h\sum_i Z_i$ is real symmetric, and the state is a density matrix $\rho$. Each input is injected by replacing the first qubit's state with $\sqrt{1 - u}\,|0\rangle + \sqrt u\,|1\rangle$, and the state then evolves by $U = e^{-iH\tau}$, through virtual nodes at sub-steps. The readouts are $\langle Z_i\rangle$ and $\langle Z_iZ_j\rangle$, with the noise of finite shots; decoherence is by GKSL dephasing. Checks:
   - unitarity, and the trace;
   - a single spin's precession, in closed form;
   - the echo state property from the injection's contraction;
   - memory capacity as $n$ grows.

## 3. Tasks and their Lowering

**Filters:** the input is encoded in the injected state, and the schedule sets the virtual nodes.

## 4. What the Foundations Must Gain

physis `quantum`: registers of qubits, Hamiltonians of Pauli strings, density matrices with their unitary and GKSL evolution, partial traces. heuresis has asked for the same registers (the idea of 2026-11-05 on what the search case studies ask of the libraries).

## 5. References

- Fujii and Nakajima (2017): in the theory note.
- Mujal, P. et al. (2021), Opportunities in quantum reservoir computing and extreme learning machines, *Advanced Quantum Technologies* 4, 2100027.
- Nakajima, K., Fujii, K., Negoro, M., Mitarai, K. and Kitagawa, M. (2019), Boosting computational power through spatial multiplexing in quantum reservoir computing, *Physical Review Applied* 11, 034021.
