# Urban Science

*The street, the building and the city as objects of one science.* Status: **draft** (2026-10-30).

## 1. Purpose

This domain holds streets, their holdings and the flows through them as objects, with the models that describe them. It begins where the user's work begins, with **street oceanography**. W. Huang's *Street Oceanography: An Oceanographic Epistemology of Socio-Cultural Dynamical Systems* proposes a field organised around the street as an open system. Heterogeneous flows pass through a heterogeneous terrain, forced from outside on many timescales; each alters what passes through it and is altered in turn. The paper takes from oceanography a disciplinary form, not its equations. The domain keeps that standing ([theory/street-oceanography.md](theory/street-oceanography.md)):
- what the paper carries becomes a definition, computed from records;
- what holds of any street becomes a theorem;
- what proposes a mechanism becomes a representation, offered for exclusion and never for prediction.

The general methods it needs are in the foundations, and the domain writes none of them again:
- **Spectral analysis:** calx `signal` and tyche `time_series`.
- **Harmonic (tidal) analysis:** tyche.
- **Co-location and point patterns:** tyche `spatial`.
- **Kinematic waves and transport:** physis `continuum`.
- **Master equations, entropy production and mixing:** physis `statistical_mechanics`.
- **Reaction–diffusion and Turing's condition:** physis `chemistry`.
- **The Hodge decomposition:** physis `discretisation`.

The domain also holds the **landscape**: a district as a cell complex whose events rewrite it, with what it carries and how it is read ([theory/landscape-dynamics.md](theory/landscape-dynamics.md), [ADR 0062](../../development/decisions/0062-landscapes-as-cellular-configurations.md)). Structure, use and meaning are joined in one process: walkers and warmth flowing between events; stores opening and shutting; a rule coming into force, carried out as a lane that rewrites the cells, and a meaning adopted that changes what may happen next. Its general mathematics went to the foundations: regular complexes, their rewriting, traces and Möbius inversion (morphe decision 0008); composition, lumpability, jump and piecewise-deterministic processes, elimination, Doi–Peliti operators, Hodge Laplacians and plane faces (physis decision 0018). A building as a complex of rooms is the architecture domain's.

The cultural reading of streets and the analysis of human flow beyond counts are planned. They join as modules when their work begins. Buildings are a domain of their own, [architecture](../architecture/README.md), on which this domain will stand when a module needs them (ADR 0060 §2.8).

The Urban Building Workbench draws the domain ([applications/urban-building-workbench](../../applications/urban-building-workbench/README.md)), in its plug-ins *Street Oceanography* and *Landscape Dynamics*.

## 2. Modules

| Module | Content | Status |
|---|---|---|
| `street` | The street system (Definition 5.1): holdings with positions, uses and declared hours. Terrain as a field read through weightings, and its dominance order (Theorem 2.2: no weighting is invariant). The length of the terrain and the compression ratio. The departure from mixing, its partition across scales, and its erosion by exchange (Theorems 2.7–2.8). The effective terrain in time: the opening forcing's exact lines, the open fraction, ascent and descent. The visits of kinds of people (a representation) | draft |
| `rhythm` | Calendars as angles: the civil week, the synodic month, the tropical year. The lattice of their combinations, read in days and weeks, with the chance rate of an explanation of order d (Proposition 3.9). Records in UTC with civil offsets. Spectra against the integration null with their lines classified (forced, combination, undecided, own). Harmonic constants under the Rayleigh criterion. Selectivity. The civil share of the clock at changes of daylight saving. Bicoherence. The paper's Figure 1, with the case it misreads. The tide of a street: its prediction from the calendars in civil time, the surge, high and low water, the range and the form (Courtier's classes) | draft |
| `response` | Susceptibility to irregular and dated forcings (rain, temperature, holidays), with clustered errors. The comparison across a repeating alteration, by randomisation | draft |
| `affinity` | Co-location of uses against free and terrain-stratified nulls. Catalysis in space and time: the bridal store and the ring store, their coupling length, coupling window and creation without capture, and the plume of a visitor's state minute by minute. The spacing of clusters, predicted (Turing) and observed (Lomb–Scargle) | draft |
| `flow` | The crowd as a kinematic wave on a street of its own blocks and crossings, signalled or not, under a rule that opens the carriageway, with a window kept second by second (a representation). Currents past a counting line: transport, net and residual current, flood and ebb, the entropy production of the measured edges, and the sign test of handedness. Person-time as Little's identity. Movement among places: transit times (phase-type) and circulation, which only the non-exact part of the rates drives (Theorem 6.3) | draft |
| `classification` | The plane of compression and selectivity, and Clark and Evans's test of whether segments form groups (the field's first test, Claim 29.2) | draft |
| `landscape` | A district as a regular cell complex from its walking network (blocks with a dead end left as holes), its schema and the transfer of state by extensive and intensive fibres; the lane as three rewrites (two refinements, then a split), closed and planned again. The landscape running as a piecewise-deterministic process on the strata of its complex: walkers hopping each way in detailed balance, the blocks' warmth, stores, a rule in force, the rule carried out, a meaning adopted, with its record, causal order and registry of enabled events. The event skeleton composed from five parts: its law and irreversibility, what coarser descriptions forget (exactly, or within a bound and for how long), interventions' interactions, evidence. Resources carried and converted, and the lattice they conserve | draft |
| `culture` | The symbolic terrain: display, signage, description and its circulation (Street Oceanography §21) | planned |
| `movement` | Human flow beyond counts: trajectories, origin–destination, the instrument of §20 | planned |

## 3. Conventions

- **Positions.** Along a street: arc length s in metres from the axis's first point, with a signed offset (left positive). Around a sensor: x east and y north in metres, in a local tangent plane.
- **Time.**
  - Records are laid out in UTC, as Unix seconds, at equal steps.
  - Civil time is UTC plus the zone's offset at the instant, daylight saving included.
  - The week is the circle of 168 civil hours, with Monday 00:00 at 0.
- **Frequencies** are in cycles per hour internally, and cycles per day where shown.
- **Uses** are the data library's classes (`street-terrain-melbourne.json`), from ANZSIC codes and OSM tags. The bridal class is a heuristic on trading names.
- **Counts** are analysed as ln(1 + y) in spectra and harmonic constants; a tide is fitted to the counts themselves, so that its prediction and its surge are in people an hour.
- **Every scalar terrain carries its weighting.** No index ranks streets without naming its w (Theorem 2.2).

## 4. Dependencies

The foundations alone (`module.json`): calx, morphe, physis, tyche. The domain stands at height 0. Architecture (buildings), geography (places and maps) and economics (value) will be declared when a module needs them; with several of them beneath it, urban science becomes the bridge of the street and the city.

## 5. Checks

`node domains/urban-science/tests/run-tests.mjs` builds the C++ with CMake (in WSL on Windows) and runs one check program per module. Each module is held to:
- closed forms: the triangle correlation of periodic blocks; the lines of an opening interval; Little's identity; the exponential exit; the Rankine–Hugoniot spill;
- brute force: dominance and its transitive reduction; the closure of covers;
- its theorems: the chain rule and exchange for mixing; detailed balance, Schnakenberg's identity and the invariance of affinities under attraction for circulation; the edges' terms of the entropy production against physis's on a driven ring;
- Monte Carlo: the catalytic yield against simulated visitors;
- constructed records with known constants and verdicts: the paper's Figure 1, and a case it misreads; a street's tide predicted exactly, its high water, form and surge; a current meter's flood and ebb;
- the landscape (62 checks): regularity and homology of every configuration, extensive quantities kept by every rewrite, the two orders of composition agreeing, detailed balance exactly where Kolmogorov's criterion says, lumpability and its discrepancy, Möbius interactions against independent stores, walkers kept by the running process, each way's hops differing by the net flow.

## 6. Open Questions

- The instrument of carried properties at entry and departure (Street Oceanography §20) does not exist. Composition enters here only as a representation.
- The terrain's symbolic component is registered only through use classes.
- Whether the plane of compression and selectivity separates on a designed sample of segments, including the uninteresting one (Claim 26.1).
- When the domain should stand on geography and economics, and whether urban geography (now in geography's scope, README §4.5) moves here.

## 7. References

See [theory/street-oceanography.md](theory/street-oceanography.md), *References*.
