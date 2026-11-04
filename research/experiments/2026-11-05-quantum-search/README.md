# Quantum Search

*Status: draft (2026-11-05).* Research: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)) and its quantum module.

## 1. The Question

Quantum search is often described by its promise, a quadratic speed-up, and quantum-inspired algorithms by their vocabulary of qubits and rotation gates. heuresis states both more narrowly (decision 0003 §2.8; theory note §5 and §10). A device is a sampler of a law, and for amplitude amplification the law is known in closed form: after r Grover iterations with t of N points marked, a measurement returns a marked point with probability sin²((2r + 1)θ), sin²θ = t/N, uniformly within the marked and within the unmarked points (Boyer, Brassard, Høyer and Tapp 1998). A quantum-inspired evolutionary algorithm whose Q-bits are never entangled draws from a product Bernoulli law, and is a search distribution of that family (Proposition 5.3).

This study measures what those statements imply:
1. that heuresis's emulator draws from Grover's law;
2. how many oracle calls Grover adaptive search spends to find a minimum, against Dürr and Høyer's bound and against classical random search, and how its threshold descends on a QUBO;
3. that the quantum-inspired EA's observations are the product family's reference draws, and how its angles converge beside PBIL.

**What is and is not computed.** No quantum computer is simulated: there are no registers, gates or amplitudes. `quantum::EmulatedGrover` draws each measurement from the exact law above, and to know which points are marked it evaluates every point once, classically, a cost it reports apart. The speed-up is therefore measured in oracle calls (Grover iterations), the unit in which Dürr and Høyer state it, and not in seconds: on a classical machine the emulation costs N evaluations before the first measurement.

## 2. Method

**The law.** N = 4096 points with f(i) = i, the first t ∈ {1024, 256, 64, 16, 4, 1} marked. For every t and every even r from 0 to 60, 4000 measurements from the emulator, each from the stream keyed by (t, r, k); the frequency of a marked point against sin²((2r + 1)θ), with z = (frequency − p)/√(p(1 − p)/4000) and Σz² against χ² on the cells where 0 < p < 1.

**Grover adaptive search.** `quantum::GroverAdaptiveSearch` on the emulator (λ = 6/5, no cutoff), driven by `search::run` until the threshold holds the minimum (a target). Landscapes of N = 2⁶ … 2²⁰ points whose values are a random permutation of 0 … N − 1, so the minimum is unique. 200 runs at each N; classical random search (`search::RandomSearch`, uniform with replacement, expectation N evaluations) on the same landscapes, 100 runs at each N up to 2¹⁴. Compared with Dürr and Høyer's bound on the expected oracle calls, (45/4)√N + (7/10) lg² N (`quantum::dh_expected_bound`). On a QUBO of 20 bits (`benchmarks::Qubo::random(20, 7)`, Q_ij ~ N(0, 1) on and above the diagonal, its minimum −32.930581 found by enumerating its 2²⁰ points), the threshold's rank after every measurement in 12 runs, and the oracle calls to the minimum in 100.

**The quantum-inspired EA.** `quantum::QuantumInspired` (Han and Kim 2002: 10 individuals, the rotation table of their knapsack setting, global migration every 20 generations) on the same QUBO for 400 generations, 20 seeds, driven by hand. Before each generation its angles are copied, and each individual's observation is compared bit for bit with the product-Bernoulli family's reference draw (`distributions::ProductBernoulli::sample`) at pᵢ = sin²φᵢ from the same keyed stream. Beside it, PBIL (`distributions::DistributionSearch` with the `Pbil` rule: λ = 10, the best one learnt at rate 0.1, borders [1/n, 1 − 1/n]) with the same budget. The heat maps show, for seed 0, the probability of reading each bit as the optimum's.

## 3. Code Version

- The Science Infrastructures `4153cb9`; heuresis `6b3cc88`; calx `90aa710`, morphe `14b3fb8`, physis `e2b98d8`, tyche `b4034b8`.
- `quantum.cpp`, built by this directory's `CMakeLists.txt` (C++20, Release; GCC 11.5 in WSL) with no warnings.

```
cmake -S . -B ~/.cache/science-build/exp-quantum-search -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-quantum-search -j4
~/.cache/science-build/exp-quantum-search/quantum results 3     # in this directory; 13 s on three threads
.venv/Scripts/python.exe research/experiments/2026-11-05-quantum-search/figures.py
```

## 4. Parameters

As in §2. Seeds: the law's streams `Draws(2026)`; the landscapes `Draws(11)`; GAS runs `Draws(r)` for run r (the QUBO's `Draws(1000 + r)`); random search `Draws(r)`; QIEA and PBIL `Draws(seed)`, seeds 0–19.

## 5. Data

None: the landscapes are random permutations and heuresis's QUBO benchmark, generated from keyed streams.

## 6. Results

**The emulator keeps the law.** Over 176 cells of 4000 measurements where 0 < p < 1, Σz² = 157.9 on 176 degrees of freedom (p = 0.83), and the largest |z| is 3.30. Against the angle (2r + 1)θ every t/N falls on the one curve sin². The optimal numbers of iterations are r* = ⌊π/(4θ)⌋ = 6, 12, 25 and 50 for t/N = 1/64, 1/256, 1/1024 and 1/4096, where the probability of a marked point is within 0.4 % of 1.

**Grover adaptive search** (oracle calls until the threshold holds the minimum):

| N | 2⁶ | 2⁸ | 2¹⁰ | 2¹² | 2¹⁴ | 2¹⁶ | 2¹⁸ | 2²⁰ |
|---|---|---|---|---|---|---|---|---|
| mean oracle calls | 7.1 | 21.6 | 57.6 | 123.7 | 289.8 | 647.9 | 1349.1 | 2684.9 |
| Dürr and Høyer's bound | 115.2 | 224.8 | 430.0 | 820.8 | 1577.2 | 3059.2 | 5986.8 | 11 800.0 |
| mean measurements | 13.9 | 25.9 | 42.9 | 62.6 | 91.4 | 124.0 | 161.7 | 203.1 |
| random search, mean evaluations / N | 1.09 | 0.97 | 1.23 | 1.08 | 0.92 | | | |

- The mean stays below the bound at every size, at 6 % of it for N = 2⁶ and 23 % for N = 2²⁰. Its log–log slope is 0.60 over the whole range, 0.56 from 2¹² and 0.52 from 2¹⁶, approaching the ½ of the theory; random search's is 1, its mean within 25 % of N at every size measured (100 runs, whose standard error is N/10).
- The classical checks of the measured points (the measurements) grow far more slowly, from 13.9 to 203.1 over the range.
- **On the QUBO of 20 bits** the threshold's rank among the 2²⁰ points falls in steps (an improvement draws uniformly among the points below the threshold, so it halves the rank on average), and the minimum is reached in a median of 2758 oracle calls over 100 runs (mean 2793, range 5 to 7420), against the bound 11 800 and random search's expectation of 1 048 576 evaluations. The emulation itself evaluated all 1 048 576 points.

**The quantum-inspired EA is a product law.** Of its 1 600 000 observed bits (20 runs × 400 generations × 10 individuals × 20 bits), 0 differ from the product-Bernoulli family's reference draws at pᵢ = sin²φᵢ on the same streams. Both it and PBIL found the QUBO's minimum in all 20 runs, QIEA in a median of 950 evaluations (610 to 1940), PBIL in 395 (220 to 640). QIEA's angles keep turning after the minimum is found, to probabilities of 0.984–0.999 of reading the optimum's bits; PBIL's stop at its border 1 − 1/20 = 0.95. For seed 0, three of QIEA's bits (5, 6 and 7) lean away from the optimum's values at generation 50, and two of them (5 and 7) still do at generation 100, with probabilities 0.10 and 0.06 of reading the optimum's; by generation 150 all three have turned.

## 7. Figures

- `search-37-grover-law.jpg`: the law against r and against the angle, with the emulator's frequencies.
- `search-38-grover-adaptive-search.jpg`: oracle calls against N, beside the bound and random search; the threshold's descent on the QUBO.
- `search-39-quantum-inspired.jpg`: QIEA's and PBIL's probabilities bit by bit, their convergence, the identity.

The pictures and captions for the gallery are in [gallery.md](gallery.md).

## 8. Caveats

- **The device is emulated by its law, not simulated.** An emulator that knows the marked set cannot be faster than enumerating the space; the oracle calls count what a device keeping Grover's law would spend, not what this computation spent. A real device adds noise, decoherence and the cost of building the oracle as a circuit, none of which is modelled.
- **The speed-up is in oracle calls.** Classical random search is charged evaluations, GAS oracle calls; the comparison assumes one oracle call costs about one classical evaluation, which a quantum oracle for an arbitrary f does not guarantee.
- **The landscapes do not matter to GAS.** Grover adaptive search depends on the values only through their ranks, and the emulated law is uniform within the marked points, so its law is the same on any landscape with a unique minimum; the QUBO shows the same search on a structured problem, not a different one.
- **QIEA and PBIL are compared on one small problem** (2²⁰ points), which both solve in every run. The comparison shows what QIEA is, not which is better; the rotation table and PBIL's rate are their customary settings, not tuned.
