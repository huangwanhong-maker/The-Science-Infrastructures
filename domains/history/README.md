# history

*Status: draft. Field, height 0.* Histories as records of what sources say: the occasions at which beings meet, and the ties between them; their chronology across calendars; and the 2-complex a record makes, with its topology and its signed ties read as a gauge field, computed exactly on the foundations, and the history read as a field theory's Feynman graph and sample path ([ADR 0041](../../development/decisions/0041-the-history-domain.md), [ADR 0042](../../development/decisions/0042-the-history-as-a-field-theory.md)). Drawn by the [History Complex](../../applications/history-complex/). Serves research and learning.

## 1. Concepts

- **A record holds assertions, never facts.** Each assertion says what a source says, with a confidence, and with who drafted and who checked it; two sources that disagree are two assertions (the factoid model of prosopography). How a history is recorded is fixed by the [codebook](codebook.md), so that two coders reading the same sources write the same record.
- **Three kinds of thing are recorded:**
  - **occasions**, datable and attested, with the beings present in their roles: a letter, a lesson, a performance, a treaty, a battle, a death;
  - **beings**, which persist and attend occasions: persons, groups, institutions, states, works;
  - **ties**, dyadic, of a type and a colour (the discipline of their theme), with a sign where a source attests the relation's state at a meeting.
- **The complex is derived, not drawn.** From a record the construction makes, with no judgement, the CW 2-complex whose vertices are the occasions, whose edges are the beings between consecutive occasions they attend (the coupling of occasions is the being that attends both), and whose 2-cells are the lenses: each tie between each two consecutive meetings of its beings, bounded by their two worldlines. A 2-cell's cyclic order is time ([theory §2](theory/the-history-complex.md)).
- **Topology reads the record.**
  - The loops of H₁ are circles of encounters that no recorded tie fills: the record's silences.
  - The closed shells of H₂ are triangles of ties carried between two meetings.
  - A discipline's history is the sub-complex of its colour.
- **Signed ties are a ℤ₂ gauge field in temporal gauge** (theory §3):
  - a lens's holonomy is its tie's flip between meetings;
  - a triangle of ties at an occasion is frustrated when its product is −1 (structural balance);
  - on every shell the Bianchi identity binds the two;
  - switching a being, reversing all its ties, changes no reading;
  - the persistence of a colour's ties, tanh β, is estimated from its flips.
- **The signs live on the refined complex**, which opens each occasion into its presences, with a rung for each tie and a face for each triangle of ties there (theory §5): on the coarse complex a triangle that changes its frustration between two meetings would be curvature at a vertex.
- **The spin foam of a history** (theory §6) is the unified model of the idea note at the group the record measures, ℤ₂:
  - a Gibbs weight of persistence on the lenses and balance on the triangles: physis's ℤ₂ lattice gauge theory on the refined complex, its couplings the maximum a posteriori of the exact likelihood of the attested signs;
  - every mean, the posterior and the derivatives in the couplings computed exactly by physis's contraction along the order of time, whose states are the spin networks of the network of the moment (theory §6.6);
  - its character expansion is a sum over spin foams, the ℤ/2 2-cycles of the refined complex, in which every label of the social spin foam is derived: spins (whether a cell is in the foam), face amplitudes tanh β, the Gauss law at every edge as admissibility, trivial vertex amplitudes;
  - Wilson loops vanish round the record's silences and carry a tie's memory across its life;
  - the reconstruction conditions the weight on the attested signs and infers the others; power is the divergence of the reconstruction under an act;
  - SU(2) appears as kinematics: the channels by which a being holds its ties.
- **The history as a field theory** (theory §8): one Feynman graph of a group field theory over ℤ₂, whose quanta are the beings between occasions, each of rank the ties it holds, and whose vertices are the occasions, each a bubble; its amplitude is 2^(F + E) Z exactly (Proposition 8.1). In time it is one term of a Dyson series, whose rates of occasions, strands and binding are estimated and tested by their martingales.
- **Other complexes of a record** (theory §7): the history complex with its closed shells filled by 3-cells, the Dowker complex of who met whom with Atkin's Q-analysis, and the history filtered by time with its persistent homology.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `record` | The record read from JSON lines, its references resolved and checked against the codebook, every fault reported with where and why; the chronology: EDTF dates, the Gregorian and Julian calendars and Old Style, Julian day numbers ([`record/`](record/include/history/record/)) | draft |
| `construction` | The history complex of a record (morphe's worldline complex of the record's encounters), its worldlines, the lives and meetings of its ties, its lenses; the derivation of the bereaved; the loops of H₁ told by beings and occasions ([`construction/`](construction/include/history/construction/)) | draft |
| `readings` | Signed ties as a ℤ₂ field in temporal gauge: flips, frustrated triangles, closed shells and the Bianchi identity, persistence by colour, carried signs, the network of the moment; disciplines as colours, their carriers and meetings ([`readings/`](readings/include/history/readings/)) | draft |
| `foam` | The refined complex K̂ and the loops round a tie's life; the spin foam at ℤ₂ as physis's gauge theory on K̂: the theory and the order of time handed to physis's exact contraction, the couplings by maximum a posteriori with Laplace errors, physis's heat bath where the contraction is out of reach, the space of foams (morphe's cycles over ℤ/2), spins, Wilson loops, the divergence of marginals for power; the other complexes: K³, the Dowker complexes and Q-analysis (morphe's), persistence in time ([`foam/`](foam/include/history/foam/)) | draft |
| `field` | The history as a Feynman graph: its stranded graph, bubbles typed by morphe's canonical forms, the Feynman rules over ℤ₂ eliminated for each record; and as a sample path: rates of occasions, strands and binding, martingale checks by era ([`field/`](field/include/history/field/)) | draft |
| `canons` | The domain's records, six seed canons drafted by a language model and not yet checked: Athens 470–322 BC, Florence and Rome 1401–1564, the Scientific Revolution 1543–1727, Europe 1785–1828, the European alliances 1871–1914 and Vienna 1890–1938 ([`canons/`](canons/)) | seed |
| `gauge` | Other groups for other ties: ℤ₁₂ and D₁₂ for music, O(d) for concepts, U(1) for positions, with physis's gauge fields; SU(2) dynamics once a record measures it | planned |
| `generation` | The generating field theory of a record as a whole: an exponential family over complexes, with a hierarchical prior over its couplings (the idea note, §8.5); `field` holds its first part, the observed graph and the process | planned |

## 3. Conventions

- **The codebook** ([codebook.md](codebook.md), version 1) fixes the format, the vocabularies, the signs, the confidences and the life of ties; the checks hold the code to it.
- **Time.** Dates are EDTF as the source writes them, with their calendar; occasions are ordered by the Julian day numbers they may denote. Overlapping dates of one being's occasions are ordered by the file, and the construction says so.
- **Ties live from their first attestation to their last**, or to their recorded end; signs between attestations are carried and weighed by persistence, never treated as attested.

## 4. Dependencies

- **Domains:** none, so history has height 0. Its chronology stands alone. Dating by eclipses and other astronomical events would bring in astronomy and metrology, as the README's §4.4 foresees; the declaration will grow when that work begins.
- **Foundations:** morphe (CW and simplicial complexes, the worldline complex, homology with generators, π₁, persistence, cycles over ℤ/2, Q-analysis, canonical forms, exact integers), physis (the ℤ₂ lattice gauge theory: its exact contraction with derivatives, its heat bath, couplings by 2-cell and frozen links; the dimensions of SU(2) invariants; binned errors) and calx (random streams, jets, linear solves, and the report of the checks). The domain writes no general method of its own: what a foundation holds is used, and what it lacks is written there (ADR 0042, §2.5).
- **The data library** checks a canon's beings against Wikidata (`data/products/history_check.py`); the domain itself reads only the record whose path it is given.

## 5. Code and Checks

The C++ is header-only, C++20, built with CMake (`CMakeLists.txt`), each module a target `history::<module>`. The checks are run by

```
node domains/history/tests/run-tests.mjs
```

They hold the codebook's vocabularies to the C++, and every canon to the format. They then build and run the C++ checks:

| Check | What it verifies |
|---|---|
| chronology | Julian day numbers against Meeus; the reforms of 1582, 1700 and 1918; Old Style; the inverses over the whole range; the EDTF subset |
| record | A record read whole; each fault the codebook names, found where it is; the derivation of the bereaved |
| construction | A disc, a circle and its loop, a sphere; carried signs; every canon a valid CW complex with χ the alternating sum of its Betti numbers |
| readings | The Bianchi identity; frustration at an occasion and in the moment; a flip rate recovered within 3σ from ties drawn at a known rate; on every canon, the identity on every shell and invariance under switching each being |
| foam | The refined prism a sphere; physis's heat bath on K̂ against the exact contraction, prior and posterior; the model in temporal gauge as physis's theory on K̂ (log Z, and the flat connections at strong coupling); spins from the character expansion against foams enumerated one by one; the Gauss law; a Wilson loop round an unfilled circle vanishing; the closed forms tanh 2β and tanh β_Δ of the posterior; switching; K³, Dowker's theorem, Atkin's Q-vector; on every canon, attested signs kept, GKS, spins in [0, 1] and below t/(1 + t), the foams the 2-cycles of K̂, persistence ending at the Betti numbers |
| transfer | The contraction against brute force on a small history (log Z, every mean, two-point functions, the law at a cut), prior and posterior; the Wilson loops of ties' lives against the two-point functions; the jets' gradient and Hessian; tanh 2β and tᵏ; Proposition 6.4, the MAP of independent cells the KT estimate with its error; on every canon, Proposition 6.3, GKS and Corollary 6.2 exactly, the model against physis's theory on all of K̂, the heat bath in agreement, the MAP's gradient zero |
| field | The stranded graph's pairing complete and its strands the lenses; Proposition 8.1 (one free strand per lens, each kernel its lens's holonomy); bubble automorphisms; the rates, binding and martingale checks on a history known by hand; all of it on every canon |

## 6. Open Questions

1. **Pooling the couplings.** A colour with no fully attested cell stays at the prior's mode; a hierarchical prior across colours belongs to the inference foundation.
2. **The JSON reader** is a copy of the neuroscience domain's, now that a second library needs JSON; where its shared place should be is open (ADR 0041, §4).
3. **Two coders.** The canon's signs are one coder's; a second coder and Krippendorff's α come before any claim from them (codebook §9).
4. **Which group for which ties** beyond ℤ₂ (the idea note, §7), and how latent signs are sampled rather than carried.

## 7. References

- Bradley, J. and Short, H. (2005), Texts into databases: the evolving field of new-style prosopography, *Literary and Linguistic Computing* 20 suppl., 3–24.
- Breiger, R. L. (1974), The duality of persons and groups, *Social Forces* 53, 181–190.
- Cartwright, D. and Harary, F. (1956), Structural balance, *Psychological Review* 63, 277–293.
- Kogut, J. and Susskind, L. (1975), *Physical Review D* 11, 395–408.
- Meeus, J. (1998), *Astronomical Algorithms*, 2nd ed., ch. 7.
- *Explanatory Supplement to the Astronomical Almanac* (2013), 3rd ed., ch. 15.
- Zaslavsky, T. (1982), Signed graphs, *Discrete Applied Mathematics* 4, 47–74.
- Atkin, R. H. (1974), *Mathematical Structure in Human Affairs*; Dowker, C. H. (1952), Homology groups of relations, *Annals of Mathematics* 56, 84–95.
- Baez, J. C. (1998), Spin foam models, *Classical and Quantum Gravity* 15, 1827–1858; Reisenberger, M. and Rovelli, C. (1997), *Physical Review D* 56, 3490–3508.
- Drouffe, J.-M. and Zuber, J.-B. (1983), *Physics Reports* 102, 1–119; Gurau, R. (2016), *Random Tensors*, Oxford University Press; Andersen, P. K., Borgan, Ø., Gill, R. D. and Keiding, N. (1993), *Statistical Models Based on Counting Processes*, Springer.
