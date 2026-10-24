# 0044. Value in Space and Time: the Economics Domain and the Economy Atlas

- **Status:** Accepted, provisionally (the application's name and colour are the user's to choose; the inference foundation will take over the fitting of transports)
- **Date:** 2026-10-24
- **Decided by:** the user asked for the next domain after history: "Economics, GIS and Spacetime (complex) … another single application; the former one mainly renders and investigates the value circulation", to be designed with care and creativity, drawing on their papers if useful, and under their rule that a domain reuses the foundations. Claude read the papers, designed the domain, its foundations and the application, and wrote this record. The design was first written as the idea note [Value in space and time](../ideas/2026-10-24-value-in-space-and-time.md).

## 1. Context

The user's papers on the economics of generative relations offer one rigorous core: loop deficits (Huang 2026, *The Geometry of Value*). A transfer network carries on each leg a transport U > 0, and re-describing value at the positions changes it by a gauge transformation. The invariants are the holonomies of loops; they all vanish exactly when a price system exists. The papers also describe value circulating in space and in time, typed stock–flow accounts, and partially ordered event complexes. Some of their wider claims, which the theory note lists, are not yet definitions.

The infrastructure had no social science domain beside history. Its foundations had no discrete exterior calculus beyond a plan, no computational geometry, and no sparse solver.

## 2. Decision

### 2.1 The Domain

`domains/economics`, a field of height 0, C++ from the start, header-only:

| Module | Content |
|---|---|
| `transfer` | transfer networks: holonomies, loop deficits, re-description, circuits in the tree gauge, the price potential, retention |
| `circulation` | the Hodge decomposition of flows (driven, local, global), reciprocal exchange and net transfer, the pairing of flows and transports, balances by kind, stocks under a stochastic rule |
| `production` | input–output tables, the Leontief model and its dual, the circulation of deliveries |
| `spacetime` | the spacetime complex of a record of transactions: holding and exchange transports, curvature, three-valued reach |
| `places` | the faces a map gives a network of places |

Its theory note states each result with its proof:

- **Propositions 1.3–1.5:** the gauge invariants; prices iff flat, the fundamental theorem of asset pricing in its simplest form; retention.
- **Proposition 3.1:** only circulating flow meets a loop deficit.
- **Proposition 5.1 and Corollary 5.2:** a lens's curvature in closed form, and flatness as an intertemporal price system.
- **Proposition 5.3:** reach is monotone.

It also says where the papers are not yet precise: where a deficit sits, non-abelian commutators across basepoints, the curvature of preferences, generative justice.

### 2.2 The Foundations

Following the user's rule, the general methods were written in the foundations:

| Method | Where |
|---|---|
| real cochains, potentials, the Hodge decomposition (dense, and by conjugate gradients for large complexes) | physis `discretisation`, its first C++ ([decision 0008](https://github.com/huangwanhong-maker/physis/blob/main/development/decisions/0008-cochains-hodge-and-the-sphere.md)) |
| the spherical Delaunay triangulation, angles and areas on the sphere | physis `computational_geometry`, its first C++ (decision 0008) |
| conjugate gradients on compressed sparse rows | calx `linalg` ([decision 0007](https://github.com/huangwanhong-maker/calx/blob/main/development/decisions/0007-conjugate-gradients.md)) |
| the worldline complex of a record of encounters | morphe `topology` (`worldlines.hpp`) |

The history domain's construction now uses morphe's worldline complex, its answers unchanged to the fixtures' precision.

### 2.3 Space and Time

- **Space.** The faces of a network of places are the spherical Delaunay triangles of the places. Curvature locates deficits, and circulation splits into eddies among neighbours and circulation through distant links.
- **Flows have conductances.** A flow is a chain, and its decomposition needs a metric: each link conducts in proportion to its gross volume (Ohm's and Kirchhoff's laws, physis `flow_decomposition`). With unit weights, a tiny link between a small and a large economy carried billions of dollars of drive cancelled by as much circulation, and an economy's circulating share exceeded 1000. With conductances the attraction is a share of gross trade, without unit.
- **Time.** Holding is a transport. In the spacetime complex, a lens's curvature is the excess return round it; for currencies, the departure from uncovered interest parity.
- **Gauge-invariant measures.** The share of the transports in their exact part changes with the units in which currencies are quoted, a re-description of value, and is not reported. The curvature is invariant, and the atlas measures it against the exchange rates' alone.

### 2.4 The Data

Three new manifests and two products in the data library, each source's terms checked:

| Source | Manifest | Licence |
|---|---|---|
| the OECD's Balanced International Merchandise Trade Statistics, adjusted for re-exports, 2010–2024 | `oecd-bimts` | CC BY 4.0 |
| the ECB's monthly euro reference rates and its deposit facility rate | `ecb-exr` | reuse with the source acknowledged |
| the BIS's central bank policy rates | `bis-policy-rates` | unrestricted with the BIS cited |

The products:

- `trade-flows.json`: the flows by year, with Natural Earth's label points as the economies' places;
- `currency-spacetime.json`: currencies, months, rates and policy rates. The euro's is the ECB's deposit facility rate, not the BIS's euro-area series: that series follows the main refinancing rate until August 2024, 25 to 50 basis points above what euro money earned under the excess liquidity since 2015.

The IMF's trade statistics were set aside because its terms forbid bulk automated download without permission.

### 2.5 The Economy Atlas

`applications/economy-atlas`, port 8910, an engine of its own on the domain (with the API's JSON, so the domain needs no JSON reader) and four views:

- **the map:** the world's trade split into exchange, drive and circulation, on Equal Earth (the geography domain's projection);
- **spacetime:** the currencies' lenses and their cumulative curvature;
- **production:** an editable table and its Leontief model;
- **circuits:** an editable transfer network with the paper's examples.

Its colour is a warm ochre, provisionally.

## 3. Alternatives Considered

- **A domain of prices and markets first** (textbook microeconomics). The user asked for circulation in space and time; markets are planned as a module.
- **Real input–output tables now** (OECD ICIO). The files are large and not yet in the data library; the production view starts from tables given to it.
- **The IMF's trade statistics.** Excluded by their terms on automated download.
- **The domain's own least squares and triangulation.** Excluded by the user's rule.

## 4. Consequences

- **The order of the domains:** economics at height 0. The root `CMakeLists.txt` builds its C++ after history's.
- **The applications are eight:** the Genome Atlas (ADR 0043) and the Economy Atlas join them. CLAUDE.md §§1, 3, 5, 7 and 9, the README and the indexes name them.
- **Checks:** the domain's 68; physis's 107 for discretisation and 11 for the sphere; calx's 7 for conjugate gradients; morphe's 12 for worldlines; the data library's 14 for the products; the application's 24.
- **What the data say** is in the theory note, §6: the balances drive 36–44% of the net flows' energy in 2010–2024, and the rates of holding remove 2.4% of the exchange rates' curvature in 2015–2026.
- **Engines judged by their includes.** The Economy Atlas's builder judges its engine stale by the headers the engine includes, followed transitively, not by every header under its include directories; the History Complex's and the Complex Atlas's builders now do the same, so that another library's unrelated work no longer makes them stale. The builders are near copies of one another; a module they share is an open question, since it would touch every application.
- **What remains:** markets and distribution; real input–output tables; fitting transports to noisy records (the inference foundation); non-abelian transports.
