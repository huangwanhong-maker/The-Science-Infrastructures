# A Walk past Melbourne's Sensors

*Status: draft (2026-11-05).* Research: a case study of the search foundation heuresis ([ADR 0073](../../../development/decisions/0073-the-search-foundation-heuresis.md)) on a real travelling-salesman problem.

## 1. The Question

The City of Melbourne counts pedestrians at 134 places in and around its centre. A technician who visits every counter once and returns to the start walks a closed tour, and the shortest such tour is an instance of the symmetric travelling-salesman problem on 134 points. The question is how heuresis's searches on permutations compare on it: a genetic algorithm with three classical crossovers, the MAX–MIN ant system, iterated local search with 2-opt, simulated annealing, and a Plackett–Luce search distribution, against random tours, at equal budgets of evaluations.

It matters for two reasons. The problem is small enough that its best tour can be found with confidence and large enough (134!/268 ≈ 10²²⁶ tours) that no search can enumerate it, so the searches can be ranked by their excess over a known value. And each algorithm assumes a structure of the space (theory note §2 and §4: by the No Free Lunch theorems a search can only beat random search by exploiting one): the tour's value depends only on which sensors are adjacent, not on where in the sequence they stand, and the algorithms differ in whether their kernels respect that.

## 2. Method

**The problem.** The sensors' positions (latitude and longitude on WGS 84) are read from the data library's `sensors.csv`; nothing else in the file is used. The distance between two sensors is the length of the geodesic between them on the WGS 84 ellipsoid, by the geodesy domain's port of Karney's algorithm (`geodesy::geodesic::Geodesic::wgs84().inverse_degrees`). The objective is the closed tour's length in metres, minimised. The sensors are also placed in the local east–north–up frame of the data library's city model of Melbourne (origin 37.8085° S, 144.9525° E; `geodesy::figure::Place`), the frame in which the figures draw them.

*Remark.* The tour runs as the crow flies between counters. It is not a route along the streets, which are drawn beneath it only to place it in the city.

**The best known tour.** Iterated local search with white-box 2-opt (`heuresis::local::two_opt_descent`, which prices a move from four distances in O(1)) and the double-bridge kick (`heuresis::local::double_bridge`), accepting a kicked local optimum when it is not worse, from 32 random starts with 20 000 kicks each (Martin, Otto and Felten 1991; Lourenço, Martin and Stützle 2003). It is not part of the comparison: it reads the distances, which the compared algorithms do not.

**The comparison.** Every compared algorithm is a black box on the ask–tell protocol: it is told a tour's length and nothing else, except the ant colony, which by its definition also receives the heuristic η_ij = 1/d_ij. Each run is driven by `heuresis::search::run` for 10⁶ evaluations; its anytime curve (each improvement of the best) is read at 81 log-spaced budgets from 10² to 10⁶ by `heuresis::analysis::best_at`. Quality is the excess over the best known length. Runs are independent and run on a calx thread pool, one run per task.

| Algorithm | heuresis parts | Settings |
|---|---|---|
| Random tours | `search::RandomSearch` | uniform permutations, asked 1000 at a time |
| Plackett–Luce EDA | `distributions::DistributionSearch` with `PlackettLuce`, `PlackettLuceRule`, the Gumbel sampler through the sampler seam | λ = 200, the μ = 20 best fitted by 10 iterations of Hunter's MM, smoothing 0.3, floor 0.01 |
| GA (CX, PMX, OX) | `genetic::GeneticAlgorithm` with `variation::cycle_crossover`, `partially_mapped_crossover`, `order_crossover` | μ = λ = 100, tournaments of 3, crossover rate 0.9, 2 elites; each child inverted with probability 0.2, else an element moved with probability 0.2 (`inversion_mutation`, `insertion_mutation`) |
| MAX–MIN ant system | `swarms::AntColony` (`AcoVariant::max_min`) | 25 ants, α = 1, β = 3, ρ = 0.02, p_best = 0.05, the generation's best depositing |
| Simulated annealing | `annealing::SimulatedAnnealing` with `inversion_mutation` | geometric cooling from T₀ = 820.1 m (the mean \|Δ\| of 2000 random inversions of random tours) to T₀/10⁴ over the budget, 100 steps per level |
| Iterated local search | `local::IteratedLocalSearch` with `local::two_opt()` and `local::double_bridge` | first improvement, accepting kicked optima that are not worse; every neighbour is one evaluation |

## 3. Code Version

- The Science Infrastructures `4153cb9`; heuresis `6b3cc88`; calx `90aa710`, morphe `14b3fb8`, physis `e2b98d8`, tyche `b4034b8`; the geodesy domain as of the root's commit.
- `walk.cpp`, built by this directory's `CMakeLists.txt` (C++20, Release; GCC 11.5 in WSL) with no warnings.

```
cmake -S . -B ~/.cache/science-build/exp-melbourne-walk -DCMAKE_BUILD_TYPE=Release
cmake --build ~/.cache/science-build/exp-melbourne-walk -j4
# from the repository's root:
~/.cache/science-build/exp-melbourne-walk/walk data/cache/melbourne-pedestrian-counts/current/sensors.csv \
    research/experiments/2026-11-05-a-walk-past-melbournes-sensors/results 1000000 10 6 32 20000 all 5
.venv/Scripts/python.exe research/experiments/2026-11-05-a-walk-past-melbournes-sensors/figures.py
```

## 4. Parameters

- Budget: 10⁶ evaluations per run. Seeds 0–9 for every algorithm except the ant system, seeds 0–4 (its construction costs n² powers per ant, about 1500 s a run on the shared machine).
- Best known: 32 restarts (seeds 100 000–100 031) × 20 000 kicks.
- Threads: 6. The algorithms' own settings are in §2's table; none was tuned on this problem beyond a pilot at 10⁵ evaluations that fixed the budget.

## 5. Data

- **Sensor positions.** City of Melbourne, *Pedestrian Counting System – Sensor Locations*, Open Data Platform (`pedestrian-counting-system-sensor-locations`), CC BY 4.0, fetched by the data library (`data/registry/melbourne-pedestrian-counts.json`) on 2026-10-30 15:14 UTC, `sensors.csv` SHA-256 `b5d42cbe…f28edbf9`. 134 sensors, all active, none sharing a position. `results/sensors.csv` repeats their positions under the same licence, attributed to the City of Melbourne.
- **Streets (figures only).** The data library's city model of Melbourne (`data/cache/products/city-model-melbourne.json`, built 2026-10-31 from the `osm-city` areas, OpenStreetMap base of 2026-07-15), © OpenStreetMap contributors, ODbL 1.0. It is read by `figures.py` from the cache and is not copied into `results/`. Five sensors lie just outside the model's box (in Kensington, East Melbourne and the Fitzroy Gardens), where no streets are drawn.

## 6. Results

**The best known tour** is 26 306.31 m. 27 of the 32 restarts end on it, at kicks between 51 and 11 041 (median 824); the other five end on one tour of 26 376.08 m (+0.27 %). Measured by chords in the city model's local frame the same tour is 26 306.31 m too: over all 8911 pairs, chord and geodesic differ by at most 0.20 mm (9.0 × 10⁻⁸ of the distance), so the frame the figures use is exact for drawing.

**The comparison** (median excess over the best known length, in per cent; 10 runs, the ant system 5):

| Algorithm | 10³ | 10⁴ | 10⁵ | 10⁶ evaluations | best run at 10⁶ |
|---|---|---|---|---|---|
| Random tours | 548.5 | 530.8 | 514.1 | 497.7 | 490.1 |
| Plackett–Luce EDA | 543.4 | 513.3 | 392.3 | 388.6 | 363.2 |
| GA, cycle crossover | 508.4 | 309.0 | 88.9 | 5.1 | 1.9 |
| GA, partially mapped | 496.7 | 276.2 | 19.3 | 5.9 | 2.1 |
| GA, order crossover | 470.0 | 245.3 | 40.0 | 4.3 | 2.6 |
| MAX–MIN ant system | 44.7 | 5.9 | 1.7 | 1.7 | 1.1 |
| Simulated annealing | 386.6 | 308.5 | 151.8 | 2.3 | 2.2 |
| Iterated local search | 229.1 | 34.4 | 5.8 | 2.8 | 1.6 |

- **The ant system** is far ahead early, because its heuristic η = 1/d is the problem's own structure: at 10³ evaluations its tours are already within 45 % of the best. It then stagnates. Four of its five runs end on two tours of 26 759.6 m and 26 762.1 m (+1.7 %), unchanged after about 3 × 10⁴ evaluations, as the trails concentrate within their MAX–MIN bounds.
- **Simulated annealing** holds its excess above 100 % until about 3 × 10⁵ evaluations, while the temperature is high, and then falls to the best median of the black-box searches, 2.3 %, with the narrowest spread (2.2 % to 4.2 %).
- **Black-box iterated local search** improves steadily from the start but pays 8910 evaluations, one per 2-opt neighbour, to confirm each local optimum, which caps it at about a hundred kicks in the budget; it ends at 2.8 %. The white-box form of the same method, with the move priced from four distances, found the best known tour.
- **The GAs** end between 4.3 % and 5.9 %. Order crossover, which keeps relative order (an adjacency-like structure), leads at the end; partially mapped crossover leads at 10⁵; cycle crossover, which keeps absolute positions, is slowest until 10⁵.
- **The Plackett–Luce EDA** barely beats random tours (388.6 % against 497.7 %). Its law ranks items by position, while a tour's length is invariant under the 268 rotations and reflections of the sequence (the dihedral group D₁₃₄ acting on positions): the structure it assumes is not the problem's.

## 7. Figures

- `search-30-melbourne-best-tour.jpg`: the best known tour over the city's streets, with the 32 restarts' results.
- `search-31-melbourne-tours-by-algorithm.jpg`: each algorithm's best tour at 10⁶ evaluations, with its excess and its median's.
- `search-32-melbourne-convergence.jpg`: excess against evaluations, medians and quartiles.

The pictures and captions for the gallery are in [gallery.md](gallery.md).

## 8. Caveats

- **The best known tour is not proved optimal.** That 27 of 32 independent restarts reach the same length is evidence, not proof. Neither the foundations nor the domains hold a lower bound for the travelling-salesman problem (a Held–Karp 1-tree bound) or an exact solver, so the gap of the best known tour to the optimum is not bounded here.
- **One setting per algorithm.** The settings are standard ones, not tuned; a tuned GA, a larger ant colony or a different cooling schedule would change the ranking's details.
- **Evaluations are the unit, not seconds.** Every algorithm is charged one evaluation per tour, which favours white-box methods not at all and black-box local search least: its moves change two edges, yet each costs a full evaluation. The ant system receives the distances as its heuristic, which the others do not.
- **The ant system had five runs**, and its compute dominated: on a machine shared with other studies each of its runs took about 1500 s, four times the uncontended pace measured in a pilot (38 s per 10⁵ evaluations). The CPU times in `results/timing.json` are wall times of runs under that contention.
- **The walk is not a route along the streets**; walking distances on the street network would be longer and would change the tour.
