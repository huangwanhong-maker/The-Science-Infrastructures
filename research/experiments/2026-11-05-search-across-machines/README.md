# Search across Machines

*Status: draft (2026-11-05).* Research: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)) on distributed search.

## 1. The Question

A search spread over many machines faces three questions before any machine is bought. How fast does a good individual cross a structured population? What does it cost to keep a generational algorithm in step when the machines run at uneven speeds? And how should an archipelago of islands be joined, and how often should they exchange? heuresis answers each on its own core, without a cluster: cellular populations on graphs (`populations::CellularEA`), a discrete-event simulation of the master–worker model (`search::run_asynchronous`), and island models synchronous and asynchronous (`populations::run_islands`, `run_islands_asynchronous`).

Its protocol also makes a promise this study holds it to (decision 0003 §2.2 and §2.5): a generational algorithm's states do not depend on where or in what order its evaluations are made, because it updates once a generation is complete, reading it in ticket order, and every random draw is keyed by what it is for.

## 2. Method

### 2.1 Takeover

Selection alone on a graph: one cell holds the best value (0), the others the values 1 to N − 1 in a random order, and there is no variation. Each generation every cell draws a parent from its closed neighbourhood and keeps it when it is not worse than its occupant (synchronous update).
- **Best of the neighbourhood.** Cell i then holds after g generations the best of its ball of radius g, so the takeover time is the eccentricity of the starting cell, exactly (theory note §12; Rudolph 2000; Giacobini et al. 2005). Eccentricities are computed by breadth-first search (`morphe::combinatorics::distances`).
- **Binary tournaments** (two draws with replacement from the closed neighbourhood, the better kept). Two cases have exact expectations, computed in `figures.py`:
  - *Ring.* A cell beside a front takes the best unless both draws miss the best neighbour, with probability 1 − (2/3)² = 5/9; the last cell, between both fronts, with 1 − (1/3)² = 8/9. Over the gap g of cells not yet taken, E[g] = 1 + p²E[g − 2] + 2p(1 − p)E[g − 1] + (1 − p)²E[g] with p = 5/9, so a ring of k cells is taken in about (9/10)(k − 1) generations, 9/5 of its eccentricity.
  - *Complete graph.* With k copies of the best, every other cell takes it independently with probability 1 − (1 − k/N)², so the share without it follows q ← q³ in mean; the chain on k is solved exactly.

Graphs: four of 256 cells drawn as heat maps (a ring, a 16 × 16 torus with von Neumann's neighbourhood, a hypercube of dimension 8, the complete graph), and 25 in all (rings of 16 to 256, tori from 4 × 4 to 32 × 32 with von Neumann's and Moore's neighbourhoods, hypercubes of dimension 4 to 10, complete graphs of 16 to 1024); tournaments over 20 seeds.

### 2.2 Synchronisation

A real-coded GA (`genetic::GeneticAlgorithm` or `SteadyStateGA` with simulated binary crossover, η = 15, rate 0.9, and polynomial mutation, η = 20, rate 1/20; binary tournaments) on Rastrigin's function in 20 dimensions on [−5.12, 5.12]²⁰ (`benchmarks::rastrigin`, minimum 0), 32 000 evaluations, on 32 simulated workers. An evaluation lasts exp(σZ − σ²/2), Z standard normal from the stream keyed by its ticket: lognormal of mean 1 and log-spread σ ∈ {0, 0.25, …, 2}. In the straggler case every duration is 1 and worker 0 takes 10.

| Scheme | Driver | Population |
|---|---|---|
| Generational GA | `run_asynchronous`, 32 workers | μ = λ = 32, one elite: one generation per round of the workers |
| Steady-state GA | `run_asynchronous`, 32 workers | μ = 32, one offspring per ask (breadth 1), replacing the worst |
| Asynchronous islands | `run_islands_asynchronous`, one machine per island | 32 islands of a generational GA (μ = λ = 8) on a bidirectional ring, one emigrant every 5 generations, no latency; each island 1000 evaluations |

Measured over 10 seeds: the simulated wall time to spend the budget (the makespan), the best value at the budget, and the simulated time until the best is at most 10 and at most 1 (runs with the budget's target, which end there). For the generational GA the expected makespan is 1000 E[max of 32 durations] = 1000 ∫₀^∞ (1 − F(x)³²) dx, computed in `figures.py` by quadrature.

**The identity.** One generational GA (seed 0) is run on six placements: the serial driver, a pool of four threads (`search::Pool`), generations cut into parts of 7, and the workers at σ = 0.5, at σ = 2 and with the straggler. An adaptor that is itself an optimiser on the protocol records the best of the population after every generation; the six records are compared bit for bit. The steady-state GA is run on the same placements.

### 2.3 Islands on a Deceptive Problem

Eight concatenated traps of order 5 (`benchmarks::trap`, 40 bits; each block worth 5 when all ones, else 4 minus its ones, so that every block pulls towards all zeros; Deb and Goldberg 1993). 16 islands, each a generational GA of μ = λ = 20 (two-point crossover at rate 0.9, standard bit mutation at 1/40, binary tournaments, one elite), synchronous migration of the best along every edge replacing the receiver's worst (`run_islands`), 128 000 evaluations in all, stopping at the optimum. Topologies: isolated, ring, torus 4 × 4, hypercube of dimension 4, complete; intervals 1, 2, 5, 10, 25 and 50 generations; and one panmictic GA of μ = λ = 320 for reference. 40 seeds each; the success rate with Clopper and Pearson's 95 % interval and the expected running time (evaluations spent over successes) by `heuresis::analysis`.

## 3. Code Version

- The Science Infrastructures `4153cb9`; heuresis `6b3cc88`; calx `90aa710`, morphe `14b3fb8`, physis `e2b98d8`, tyche `b4034b8`.
- `machines.cpp`, built by this directory's `CMakeLists.txt` (C++20, Release; GCC 11.5 in WSL) with no warnings.

```
cmake -S . -B ~/.cache/science-build/exp-search-machines -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-search-machines -j4
~/.cache/science-build/exp-search-machines/machines results all 5     # in this directory
.venv/Scripts/python.exe research/experiments/2026-11-05-search-across-machines/figures.py
```

The three parts ran in 19 s, 41 s and 448 s on five threads of a machine shared with other studies.

## 4. Parameters

As in §2. Seeds: takeover heat maps seed 1, takeover times seeds 0–19; synchronisation seeds 0–9 (the duration streams keyed by the same seed); islands seeds 0–39, island i drawing from `Draws(seed, i)`. The problem sizes were fixed by pilots: Rastrigin in 10 dimensions was solved to 0.002 by every scheme and told nothing; traps of 60 bits defeated almost every archipelago (1 success in 8 runs at best).

## 5. Data

None: the problems are heuresis's benchmarks, whose optima are known.

## 6. Results

### 6.1 Takeover

- **Best of the neighbourhood: takeover equals the eccentricity on all 25 graphs**, from 1 generation (complete graphs) to 128 (the ring of 256). On the heat maps the wave of the best arrives at each cell exactly at its distance.
- **Binary tournaments are slower by a factor that depends on the graph:** about 1.8 on rings (the ring of 256 is taken in 230.7 ± 9.7 generations against the exact expectation 229.7), 1.35–1.6 on von Neumann tori, 2.0–2.2 on Moore tori, 1.3–1.45 on hypercubes. On complete graphs the eccentricity is 1 and says nothing: the take-over is the logistic q ← q³, measured 4.2, 5.9, 7.6 and 8.7 generations on 16, 64, 256 and 1024 cells against exact expectations of 4.22, 5.84, 7.36 and 8.81.

### 6.2 Synchronisation

| | σ = 0 | σ = 0.5 | σ = 1 | σ = 2 | straggler |
|---|---|---|---|---|---|
| **Makespan** (ideal 1000): generational | 1000 | 2561 | 5491 | 15 341 | 10 000 |
| — its expectation, 1000 E[max of 32] | 1000 | 2564 | 5490 | 15 451 | 10 000 |
| — steady-state | 1000 | 1001 | 1006 | 1163 | 1030 |
| — islands | 1000 | 1030 | 1081 | 1379 | 10 000 |
| **Best at the budget** (median): generational | 0.272 | 0.272 | 0.272 | 0.272 | 0.272 |
| — steady-state | 0.0099 | 0.0091 | 0.0059 | 0.0065 | 0.0072 |
| — islands | 5.86 | 4.54 | 5.44 | 5.92 | 5.04 |
| **Time until f ≤ 10** (median): generational | 262 | 662 | 1402 | 3868 | 2606 |
| — steady-state | 151 | 153 | 164 | 152 | 147 |
| — islands | 780 | 745 | 798 | 832 | 792 |

- **The generational GA pays the maximum.** Its makespan follows 1000 E[max of 32 durations] within 0.8 % at every σ (the median of 10 runs), 15 times the ideal at σ = 2, and the straggler makes it exactly 10 times the ideal, since every generation waits for it. Its quality at the budget does not change at all, because its runs are the same run (§6.3).
- **The steady-state GA keeps every worker busy:** its makespan stays within 17 % of the ideal even at σ = 2, and one slow worker of 32 costs it 3 %. It also reached f ≤ 1 in all 100 of its runs (every σ and the straggler), in medians of 329 to 386 time units, while the generational GA reached it in 6 of 10 runs, the same six seeds at every σ, in 852 time units at σ = 0 and 13 123 at σ = 2. On this problem its replacement of the worst is a stronger selection than the generational GA's, so part of that difference is the algorithm, not the waiting.
- **Asynchronous islands** never wait for each other, so their makespan is the slowest island's sum of 1000 durations (1379 at σ = 2). But with an equal share of the budget per island, the straggler's island must spend its share at a tenth of the speed, and the makespan is 10 000. With populations of 8 they stall near f = 5 and never reach f ≤ 1.

### 6.3 The Protocol's Identity

The generational GA's population best, generation by generation over all 1000 generations, is identical bit for bit on all six placements: serially, on four threads, in parts of 7, and on workers at σ = 0.5, at σ = 2 and with a straggler. The steady-state GA is identical on the first three (each of its candidates is told before the next is asked) and differs on the workers, where the order of arrival is the durations': its population best departs from the serial run's by up to 41.9, 83.7 and 44.6.

### 6.4 Islands on Deceptive Traps

Successes in 40 runs (the panmictic GA of 320: 10; isolated islands: 0):

| Topology | 1 | 2 | 5 | 10 | 25 | 50 |
|---|---|---|---|---|---|---|
| ring | 3 | 9 | 10 | 16 | 16 | 2 |
| torus 4 × 4 | 0 | 4 | 6 | 14 | 12 | 9 |
| hypercube | 0 | 3 | 9 | 17 | 15 | 16 |
| complete | 0 | 0 | 2 | 6 | 8 | 5 |

- **Isolated islands never succeed**: a GA of 20 settles each block into the deceptive all-zeros before it finds all-ones. **Frequent migration fails too**: every generation (interval 1) the islands become one population that converges as fast, and the complete graph with intervals 1 or 2 never succeeds.
- **Intermediate intervals are best**: at 10 generations the hypercube succeeds 17 times in 40 (Clopper–Pearson 95 % interval 0.27–0.59) and the ring 16 times, against 10 for the panmictic GA with the same evaluations; their expected running times are 198 400 and 239 000 evaluations against 477 088. Sparse topologies tolerate shorter intervals than the complete graph, as their takeover times in §6.1 predict; the ring fails again at 50 generations, when good blocks travel too slowly round it.

## 7. Figures

- `search-33-takeover-waves.jpg`: cell × generation on four graphs of 256 cells, under both rules.
- `search-34-takeover-times.jpg`: takeover against eccentricity and against size, with the exact expectations.
- `search-35-synchronisation.jpg`: makespan, quality and times to targets against σ and with a straggler; the identity.
- `search-36-island-topologies.jpg`: success rates by topology and interval.

The pictures and captions for the gallery are in [gallery.md](gallery.md).

## 8. Caveats

- **The workers are simulated.** Durations are independent of the genotype and of the machine's load, and there is no latency to the master; on real machines both matter. The asynchronous islands had no migration latency.
- **Structure and synchronisation are confounded in §6.2.** The islands' populations of 8 and the steady-state's replacement of the worst are different algorithms from the generational GA, so their quality differences are not only the cost of waiting. The makespans are comparable; the qualities are of three algorithms.
- **Equal shares of the budget per island** (the drivers' rule) make the slowest machine set the makespan. A wall-clock budget, which the asynchronous drivers do not offer, would let the other islands go on.
- **One problem each, fixed by pilots** (§4), and one setting of each algorithm. The intervals that suit the traps depend on the island's size, the problem's order and the budget.
- **Tournament takeover on tori and hypercubes** has no closed form here; only rings and complete graphs are compared with exact expectations.
