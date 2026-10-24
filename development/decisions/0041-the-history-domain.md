# 0041. The History Domain and the History Complex

- **Status:** Accepted; amended by [0042](0042-the-history-as-a-field-theory.md)
- **Date:** 2026-10-24
- **Decided by:** the user asked for a History Complex: histories of science, art, music, international relations and philosophy, and their integration, each represented as a 2-complex, to bridge social spacetime with theoretical physics and field models. They asked for the semantics to be designed first, and for Claude to form its own view rather than follow their discussion papers. They judged full historical data out of reach and proposed building from classical facts. Then they agreed: to build it properly, to make history a new domain, and to code the canon's signs by a codebook that Claude writes and the user reviews. While the work went on they asked for several things:

- histories of the several disciplines and integrated ones;
- a visualisation in three dimensions, and interactive design;
- that the meanings of spin, amplitude and the other labels of the spin foam be considered carefully, as the revised unified model of their papers the idea note had built;
- other complex models besides the CW complex, as a challenge.

Claude designed the domain, the construction, the spin foam and the application, and wrote this record.

## 1. Context

The user's papers place social spacetime in CW 2-complexes and spin foams: the event–coupling–primitive architecture; the social spin foam and its weights; the street spin foam's frame group; generative relational history. The idea note [The History Complex](../ideas/2026-10-23-history-complex.md) found one fault in them and proposed a carrier. The fault is that facts on faces make a hypergraph, not a CW complex: a 2-cell's cyclic order must mean something, and in a history it can only mean time. The carrier puts occasions on vertices and beings on edges, between consecutive occasions they attend, so that the coupling of two occasions is a being that attends both. The faces are ties between two beings across the interval between two of their meetings, and the themes (ideas, motifs, offices, disputes) are types of ties.

A test on classical facts (the idea note, §11.2) showed four things:

- the construction runs on a canon drafted by a language model;
- morphe's invariants read it meaningfully: loops where the record is silent, shells where triangles of ties persist;
- signed ties read as a ℤ₂ field give the frustration of the Tilsit system and of the crisis at the Congress of Vienna;
- three faults: carried signs go stale, ties never end, and the canon is lopsided.

## 2. Decision

1. **A domain `domains/history/`**, a field of height 0. It depends on no domain and uses morphe and calx. It is the first domain of the humanities, and its C++ is header-only, in three modules:
   - **`history::record`:** the record read from JSON lines and checked against the codebook, and its chronology (EDTF dates, the Gregorian and Julian calendars and Old Style, Julian day numbers);
   - **`history::construction`:** the history complex, a morphe `CWComplex`, with the loops of its H₁;
   - **`history::readings`:** the signed ties as a ℤ₂ field in temporal gauge (flips, frustration, closed shells and the Bianchi identity, persistence, the network of the moment), and the disciplines as colours.
2. **A codebook** ([`domains/history/codebook.md`](../../domains/history/codebook.md), version 1) fixes how a history is recorded:
   - the format;
   - the types of beings, occasions and ties;
   - six colours: science, art, music, letters, philosophy and international relations;
   - the sign of each type of tie;
   - confidences;
   - the life of a tie: from its first attestation to its last, or its recorded end;
   - provenance.

   The domain's checks hold the C++ to the codebook's vocabularies. A canon names the version it was coded by.
3. **Records hold assertions, never facts.** They follow the factoid model of prosopography: each assertion with its sources, a confidence, and who drafted and who checked it. Reading never throws on content: every fault becomes a problem, with where it is and why.
4. **Canons first.** The domain's first records are canons: the classical facts of a history, a coarse resolution that denser sources refine later. The first is Europe, 1785–1828: 33 beings, 92 occasions, 57 ties. It was drafted by a language model, says so, and is not yet checked against its sources.
5. **Checking identities is an edge.** `data/products/history_check.py` compares a canon's beings with Wikidata (CC0): identifiers, years of life and deaths. It reports agreement, disagreement and what it could not match; a person resolves disagreements in the record.
6. **The signs live on the refined complex K̂**, which opens each occasion into its presences, with a rung for each tie at each meeting and a face for each triangle of ties there. On K, a triangle whose frustration changes between two meetings would be curvature at a vertex, which no connection on K can carry (theory note, Proposition 5.1).
7. **The spin foam of a history is the unified model at the group the record measures, ℤ₂** (the module `history::foam`):
   - a Gibbs weight, persistence β_c on the lenses and balance β_Δ on the triangles, estimated from the record;
   - its character expansion is a sum over the ℤ/2 2-cycles of K̂, which derives every label of the social spin foam: spins (whether a cell is in the foam), face amplitudes tanh β, the Gauss law at every edge (admissibility, a theorem), and trivial vertex amplitudes;
   - spins are computed from the configuration ensemble by an exact identity, on the support found over GF(2);
   - Wilson loops vanish exactly round the record's silences;
   - the reconstruction ensemble conditions the weight on the attested signs, in place of carrying them forward;
   - power is the divergence of the reconstruction under an act;
   - SU(2) enters as kinematics only (dimensions of invariants), until a record measures it.
8. **Other complexes of a record:** K with its closed shells filled by 3-cells; the Dowker complex of the presence relation, with Atkin's Q-analysis; and K filtered by time, with persistent homology over ℤ/2 written in morphe (`topology/persistence.hpp`), being general mathematics.
9. **Records join.** `join_records` glues several canons along the beings, occasions and ties they share; the application offers all canons, joined, as one record. Canons that mean the same being use the same identifier for it (franz_joseph, wilhelm2, austria_hungary, germany, leibniz). An identifier identifies across canons only where they agree about what it names, by Wikidata item or else by name, as an authority file would; otherwise the later canon's is renamed `id@canon` (the theory note, Definition 1.3).
10. **The History Complex**, `applications/history-complex/`, on port 8870, in the pattern of the Complex Atlas:
   - an engine of its own on the domain and morphe, built by its recipe into `data/cache/engines/history/` and run by its server, one request per run;
   - a page with six views: the foam (three.js 0.170.0, the history in spacetime, closed or opened, with a cut at any date that can be played), the loom, the moment, the complex, the record and the readings.

   Its main colour is a faded mauve, provisionally, until the user chooses one.
11. **A reader's dates** (a cut, the end of a year, the axis of time) are read in the historians' convention: Julian before the Gregorian reform of 15 October 1582, Gregorian from then (Meeus 1998, ch. 7), so that the years of a canon dated in the Julian calendar, as those of antiquity and the Renaissance are, end where the reader's do (codebook §4).

## 3. Alternatives Considered

- **History inside the Complex Atlas.** That atlas is a work on complexes in general, and nothing may depend on an application. The record, its chronology and the construction name objects of the human world, so they belong to a domain (CLAUDE.md §3.7).
- **The facts-on-faces carrier of the social spin foam.** Rejected for the reason of §1: it gives the 2-cells no cyclic order.
- **Provisional JavaScript.** The construction feeds morphe's exact homology and π₁, which are C++. ADR 0003 asks for C++ wherever possible, and the user prefers a C++ core for the field theories this will reach.
- **Declaring metrology and astronomy now**, as the README's §4.4 foresaw. The chronology needs neither yet. A dependency is declared when it is used.
- **Wikidata as the record.** Its statements are mostly unsourced and its "influenced by" is not evidence. It serves to check identities, not to supply ties.

## 4. Consequences

- **The order of the domains:** history at height 0, beside metrology and neuroscience. The root `CMakeLists.txt` builds its C++ after neuroscience's. The domain uses morphe and calx (random streams, checks).
- **morphe** gains persistent homology over ℤ/2, with its own checks.
- **The applications are six.** CLAUDE.md §5, §7 and §9, the README and the domains' index name the new domain and application.
- **JSON in the C++.** The domain carries a copy of the neuroscience domain's JSON reader. That header says it moves to a shared place when a second library needs JSON, and this is the second. A domain may not use the API's JSON (an edge), and no foundation yet holds formats. Where JSON should live is open: a module of calx, a small formats library, or keeping the two copies identical. It is the user's decision.
- **What remains of the idea note:**
  - groups beyond ℤ₂ (D₁₂ for music, O(d) for concepts), and SU(2) dynamics;
  - the generating field theory fitted to a record;
  - two coders and their agreement;
  - denser records.

  The note keeps them, with this record as the place its first part went.
