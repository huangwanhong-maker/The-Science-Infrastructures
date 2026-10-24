---
status: seed
---

# Value in Space and Time: the Economics Domain and the Economy Atlas

*Status: seed, 2026-10-24.* The design of the infrastructure's economics. The user asked for it on 2026-10-24, after the History Complex: "the next domain: Economics, GIS and Spacetime (complex) … another single application; the former one mainly renders and investigates the value circulation; be creative; see my papers." The note is written before the code, as the contract the modules are built to, and is promoted into an ADR, the domain's README and its theory note as the work lands. It follows the user's rule for the repository: what the foundations hold is used, and what a domain needs that names no object in nature is written in a foundation (CLAUDE.md §2.1).

## 1. Principles

1. **Value is transported, not measured.** The user's paper on loop deficits (Huang 2026, *The Geometry of Value*) gives the one rigorous core to build on. A transfer network Γ carries, on each leg e = (u, v), a transport U_e > 0: the value registered at v per unit surrendered at u. Re-describing the value at each position by g_v > 0 changes U_e to g_v U_e g_u⁻¹. That is an abelian lattice gauge theory whose group is ℝ_{>0}, and the logarithm makes it (ℝ, +).
   - **The invariants** are the holonomies of closed walks, λ(γ) = Π U_e^{±1}, and the **loop deficit** is δ(γ) = −log λ(γ).
   - **Description classes** form Hom(H₁(Γ; ℤ), ℝ_{>0}) ≅ ℝ_{>0}^{b₁}.
   - **Every loop deficit vanishes** exactly when a potential φ exists with log U_e = φ(v) − φ(u), unique up to a constant.

   The mathematics is that of the gauge theory of currency exchange (Ilinski 2001; Young 1999), and the user's papers differ in what U means, not in the mathematics. *Remark.* On a finite graph the last statement is the fundamental theorem of asset pricing in its simplest form: there is no arbitrage round any loop exactly when a consistent system of prices exists.
2. **Flows and transports are different cochains.** A flow F (goods, money, people per unit time on each leg) is a 1-chain; its boundary is the net accumulation at each position, so stocks and flows obey ∂F = ds/dt + sources. A transport is a 1-cochain ω = log U. Each has a Hodge decomposition on a complex with 2-cells:
   - flow = gradient (driven by a potential) + curl (circulation round faces) + harmonic (circulation round the holes);
   - transport = exact (a price system) + co-exact (local deficits, a curvature dω on each face) + harmonic (deficits round unfilled loops), of dimension b₁.

   This is the precise form of the papers' "the gradient part carries no circulation". It has been used for rankings from pairwise comparisons (HodgeRank, Jiang, Lim, Yao and Ye 2011) and for games (Candogan, Menache, Ozdaglar and Parrilo 2011). The pairing ⟨F, ω⟩ is the value a flow loses or gains along its transports.
3. **Space gives the loops faces.** A transfer graph alone has holonomy but no curvature, as the papers concede. Geography supplies 2-cells:
   - regions bounded by their borders, or the faces of a triangulation of places on the sphere;
   - a deficit then has a place: by Stokes, the deficit round a loop is the sum of the curvatures of the faces it encloses;
   - the harmonic part lives on the holes the map leaves: seas, closed borders, missing links.

   This is what the GIS shows.
4. **Time is a direction of transport.** Holding value from t to t′ is a transport too: interest, inflation, depreciation, storage, spoilage. In the **spacetime complex** of a record of transactions, as in the History Complex:
   - agents have worldlines, and a transaction is an event on the worldlines of those taking part;
   - the 2-cells are the lenses between two agents' consecutive exchanges;
   - a loop of exchange closes through time, and its deficit is net of the time value of what was held;
   - curvature locates a deficit in an interval;
   - flatness is an intertemporal price system: state prices or discount factors, the no-arbitrage potential across agents and dates.

   The History Complex's construction is the same one with ℤ₂ in place of ℝ_{>0}. That construction names no object in nature, so it moves into morphe (§3.3), and history and economics both use it.
5. **Typed value, explicit conversion.** Value comes in kinds that need not share a unit: money, goods, labour time, care, attention, knowledge. Conversions between kinds are legs of their own, never silent. Stock–flow balances are kept for each kind (the user's typed lifecycle diagnostic, Huang 2026), with each conversion's rate a transport whose loops can be read.
6. **Checked against closed forms.** The checks use:
   - the papers' own worked numbers, with the triangle (1.2, 0.5, 1.0) as one of them;
   - the stock drift M = [[0.9, 0.6], [0.1, 0.4]] with π = (12/7, 2/7);
   - textbook input–output examples (Miller and Blair 2009);
   - the dimensions of the harmonic spaces against morphe's Betti numbers, and the orthogonality of the Hodge parts;
   - finite no-arbitrage examples, and brute force on small cases.

## 2. What the Papers Give, and What Is Not Yet Precise

Seven papers were read for this design (a summary is kept in the session's notes, and none of them is copied here).

- **Precise enough to implement and test:**
  - the ℝ_{>0} transport theory with its tree gauge, chord holonomies, cycle rank and potentials;
  - abelian transports in general, as Hom(H₁, G), with the pair (deficit, phase) for complex transports;
  - non-abelian matrix transports, read only through traces and spectra;
  - the retention-ratio identity Π r_i = λ under steady flow;
  - column-stochastic stock dynamics;
  - the typed temporal multigraph with its viability, reach, throughput and typed balances, once two choices are fixed: Kleene's strong three-valued logic for "?", and time-respecting paths;
  - partially ordered event complexes with provenance.
- **Not yet precise, and so not built as if it were:**
  - where a deficit "sits": any split of δ among positions whose product is λ can be realised, so the absorbing position is not an invariant;
  - a commutator of two holonomies taken at different basepoints, which is not gauge invariant (its spectrum is);
  - the "curvature" of preferences, and prices as Goldstone modes;
  - the functional of generative justice;
  - the closure of a conversion.

  The theory note will say so where it meets them, with the reason, as the History Complex's note did for the papers it drew on.

## 3. The Foundations

### 3.1 physis `discretisation`: Cochains and the Hodge Decomposition

Real cochains on a morphe `CWComplex`, written once for every science that needs discrete exterior calculus. Electromagnetism on complexes is planned to use them too (Desbrun, Hirani, Leok and Marsden 2005; Lim 2020).

| Header | Content |
|---|---|
| `cochains.hpp` | k-cochains; the coboundaries d₀, d₁ (sparse, from the attaching words); the pairing with chains and paths (holonomy, log form); closedness; integration along a spanning forest (the potential of a closed 1-cochain, and the tree gauge) |
| `hodge.hpp` | diagonal Hodge stars from weights; the Hodge Laplacians L₀, L₁; the decomposition of a 1-cochain into exact, co-exact and harmonic parts by least squares; harmonic bases, with dim = b₁ (morphe) |

The JavaScript stand-in for `discretisation` (physis ADR 0003) is untouched. These are the module's first C++.

### 3.2 calx `linalg`: Sparse Matrices and Conjugate Gradients

CSR matrices, matrix–vector products, the conjugate-gradient method for symmetric positive semidefinite systems (Hestenes and Stiefel 1952), including singular ones solved in the range (the Laplacians have kernels). The method is checked against the dense solve. calx is being extended by another session (eigenproblems, task graphs), so this is coordinated with it before it is written.

### 3.3 morphe `topology`: the Worldline Complex

The construction of the History Complex, made general: a sequence of events, the participants of each, and pairs of participants with the events at which they meet. From these:

- each participant's worldline runs through its events, giving the edges;
- each pair has a lens between consecutive meetings, giving the 2-cells;
- each event can be opened into its presences, with a rung for each pair meeting there and a face for each triangle of pairs: the refined complex.

The history domain's `construction` and `foam/refined.hpp` then use it, with the record's checks as they are, and so does the economics domain's spacetime complex.

## 4. The Domain

`domains/economics`, a field of height 0, using morphe, physis and calx. It is C++ from the start, header-only, each module a target `economics::<module>`:

| Module | Content |
|---|---|
| `transfer` | transfer networks: positions, legs with transports and flows, value kinds and conversions; loop deficits, the tree gauge and chord holonomies, potentials (no-arbitrage), retention ratios under steady flow; the papers' worked numbers as fixtures |
| `circulation` | the Hodge decomposition of flows and of transports on a network with faces (through physis); throughput, time-respecting reach with three-valued viability, typed stock–flow balances, column-stochastic stock dynamics |
| `production` | input–output tables: the Leontief inverse, output and value-added multipliers, the price dual p = pA + v, Perron–Frobenius of A (productivity, the maximal rate of profit), and the circulation of value round the cycles of the production network |
| `spacetime` | the spacetime complex of a record of transactions (morphe's worldline complex): holding transports on the worldlines, exchange transports on the rungs, the curvature of each lens, the no-arbitrage potential across agents and dates, Stokes's theorem checked on every loop |
| `places` | positions that are places: their coordinates and regions, the 2-cells geography supplies (borders, or a spherical triangulation of the places), and the geodesic arcs between them (from the geodesy domain's JavaScript in the page; the domain itself stays pure) |

Data are inputs: a record of positions, legs, flows and transports (JSON lines, with a codebook like history's), an input–output table, a record of transactions. The domain reads only what it is given.

## 5. The Economy Atlas

`applications/economy-atlas`, port 8910, in the pattern of the History Complex: an engine of its own on the domain, run by its server, and a page with three views.

1. **The map** (GIS), on the geography domain's projections (Equal Earth), with the data library's country geometry:
   - places as vertices, flows as arcs, the map's faces as 2-cells;
   - the potential of the flows as a field over the map (where value is drawn from and to);
   - circulation on the faces and round the holes, drawn as eddies;
   - the transports' deficits per face as curvature;
   - a year slider; the pairing of flow and transport per country.
2. **Spacetime**: the spacetime complex of a record of transactions in three dimensions, time rising as in the History Complex's foam (the drawing shared through `visualization/js`, not copied); the lenses coloured by curvature; the no-arbitrage potential as a height; a cut at any date.
3. **Production**: the input–output network of an economy, its multipliers and prices, and the circulation of value round its cycles.

Its colour is the user's to choose; provisionally a warm ochre within the family the user prefers.

## 6. Data

Candidates, each with its licence checked before it is fetched into `data/cache/` (never committed):

- bilateral trade (IMF Direction of Trade Statistics; CEPII BACI; UN Comtrade);
- input–output tables (OECD ICIO; WIOD);
- exchange rates and prices (ECB reference rates; World Bank WDI, CC BY 4.0);
- country geometry (Natural Earth, public domain, already used by the Earth Atlas).

Seed records of the kind the papers describe (a fee chain, a household's circuits, a platform's legs) are written by hand and kept small, as canons are, saying what they are.

## 7. The Order of the Work

1. physis `discretisation` (cochains, Hodge) and calx sparse CG, with their checks.
2. morphe's worldline complex; the history domain moved onto it, its checks unchanged.
3. The domain: `transfer`, `circulation`, `production`, `spacetime`, `places`, each with its checks.
4. The data products, and the application's engine, server and page.
5. The ADR, the theory note (with the papers' imprecise points stated), the gallery.

## 8. Open Questions

1. **The user's choice of name and colour** for the application.
2. **Which data first:** trade flows give the map its circulation, input–output tables the production view, and transactions the spacetime view. Licences decide part of it.
3. **The inference foundation.** Fitting transports and flows to noisy records needs it, as the history domain's couplings will.
4. **Non-abelian transports** (bundles of goods converted by matrices): read through traces and spectra (physis `symmetry`), when a record measures them.
