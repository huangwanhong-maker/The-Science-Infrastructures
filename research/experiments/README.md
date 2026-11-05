# experiments: One Record per Experiment

*Status: draft (2026-11-05).* Research.

## 1. Purpose

An experiment asks one question and keeps everything needed to answer it again: the code at a stated version, the parameters, the data and their provenance, the results, and what they do and do not show.

## 2. The Record

Each experiment is a directory `YYYY-MM-DD-slug/` holding:

- `README.md`, the record:
  1. the question, and why it matters;
  2. the method;
  3. the code version (the commits of the libraries used);
  4. the parameters;
  5. the data, with source, licence and retrieval date;
  6. the results, with their numbers;
  7. the figures;
  8. the caveats.
- `CMakeLists.txt` and the C++ of the study, built against the foundations and the domains as a project of its own, which writes its results to `results/`.
- `results/`: small tables (CSV or JSON) from which every figure and number of the record is drawn. They are committed; large or licensed inputs stay in the data library's cache.
- `figures.py`, which draws the figures from `results/` in the house style of [`style/search_style.py`](style/search_style.py) and saves them as `assets/gallery/search-NN-subject.jpg` (1600 × 1000).
- `gallery.md`: the record's pictures and captions for [docs/gallery.md](../../docs/gallery.md).

## 3. Index

| Date | Experiment | Question |
|---|---|---|
| 2026-11-05 | [Six Searches on One Landscape](2026-11-05-search-on-a-landscape/README.md) | How six algorithms of different kinds move over one two-dimensional landscape, at equal numbers of evaluations |
| 2026-11-05 | [Benchmarking Continuous Search](2026-11-05-benchmarking-continuous-search/README.md) | Which continuous searches solve which functions, in which dimensions, at what cost (a COCO-style comparison) |
| 2026-11-05 | [Kepler's Third Law Rediscovered](2026-11-05-keplers-third-law/README.md) | Whether genetic programming typed by physical dimensions finds T = 2π√(a³/GM) in JPL's DE440 where untyped search does not |
| 2026-11-05 | [A Body Learns to Move](2026-11-05-a-body-learns-to-move/README.md) | How far the robotics domain's rabbit on the Moon can be tuned, and the trade-off of speed against power |
| 2026-11-05 | [Fronts in Three Objectives](2026-11-05-fronts-in-three-objectives/README.md) | How five many-objective searches cover simplex, sphere and disconnected fronts |
| 2026-11-05 | [Repertoires of an Arm](2026-11-05-repertoires-of-an-arm/README.md) | How quality-diversity searches illuminate everything a twelve-jointed arm can reach |
| 2026-11-05 | [Rugged Landscapes](2026-11-05-rugged-landscapes/README.md) | How NK landscapes turn rugged with k, enumerated whole, and how searches fare as they do |
| 2026-11-05 | [A Walk past Melbourne's Sensors](2026-11-05-a-walk-past-melbournes-sensors/README.md) | The shortest closed walk past the city's 134 pedestrian counters, and how eight searches approach it |
| 2026-11-05 | [Search across Machines](2026-11-05-search-across-machines/README.md) | Takeover on topologies, the cost of synchronisation on heterogeneous workers, and islands on deceptive problems |
| 2026-11-05 | [Quantum Search](2026-11-05-quantum-search/README.md) | Grover adaptive search on the exact law of its device, against classical search, and the quantum-inspired EA as a product distribution |
| 2026-11-05 | [Fitting Neurons](2026-11-05-fitting-neurons/README.md) | Which of a spiking neuron's parameters its firing determines, and the repertoire of firing classes it can produce |
| 2026-11-05 | [The Genetic Algorithm as a Population](2026-11-05-the-genetic-algorithm-as-a-population/README.md) | A GA with fitness-proportional selection as a Wright–Fisher population: fixation, drift and balance |
| 2026-11-05 | [Tuning an Image Pipeline](2026-11-05-tuning-an-image-pipeline/README.md) | Tuning eikon's vessel reconstruction on synthetic phantoms, and how tuned chains overfit |
| 2026-11-05 | [Designing a District](2026-11-05-designing-a-district/README.md) | A district of Melbourne's Hoddle Grid redesigned under objectives from several domains |
| 2026-11-05 | [A District's Economy: the Problem](2026-11-05-a-district-economy-the-problem/README.md) | One planning problem for three models of people: the ground floors, their hours and the links of a district, valued by the spending they capture |
| 2026-11-05 | [A District's Economy without a Model](2026-11-05-a-district-economy-without-a-model/README.md) | What searches find when the people are interpolated from the counters and the evaluation is a black box, and where they are blind |
| 2026-11-05 | [A District's Economy as a Flow](2026-11-05-a-district-economy-as-a-flow/README.md) | People as conserved streams on the network, with their circulation and tides, and what a search that uses the model's linearity buys |
| 2026-11-05 | [A District's Economy as a Field Theory](2026-11-05-a-district-economy-as-a-field-theory/README.md) | People as a field with an action S on the district's space–time complex, solved in its modes, and what the action's locality buys a search |
| 2026-11-05 | [Walkers of the Hoddle Grid](2026-11-05-walkers-of-the-hoddle-grid/README.md) | Walkers as agents, calibrated to the counters: a week, a street closed, a festival's crowd, and the street's regime decided online |
| 2026-11-05 | [Households of the Centre](2026-11-05-households-of-the-centre/README.md) | Households choosing the centre's dwellings, calibrated to rents and the Census: supply shocks, Schelling's limit, and the planner's permits |
| 2026-11-05 | [Owners and Occupiers](2026-11-05-owners-and-occupiers/README.md) | The uses of the Hoddle Grid's buildings 2002–2024 as agents on a nested chain: why vacancy persists, and interventions under a budget |

The first experiments are case studies of the search foundation heuresis ([ADR 0073](../../development/decisions/0073-the-search-foundation-heuresis.md)): what its algorithms do, how they compare, and what they find on real problems of the domains. The next are the urban-science domain's agent-based models ([ADR 0075](../../development/decisions/0075-agents-in-the-city.md)).
