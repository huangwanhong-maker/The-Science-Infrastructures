# economics

*Status: draft. Field, height 0.* Value carried between positions and through time: its transports and the invariants of carrying it round (loop deficits, the existence of prices), the circulation of its flows on a network with faces, production as a circular flow, and the spacetime complex of a record of transactions ([ADR 0044](../../development/decisions/0044-value-in-space-and-time.md)). Drawn by the [Economy Atlas](../../applications/economy-atlas/). Serves research and learning.

## 1. Concepts

- **Value is transported, not measured.** A leg from one position to another carries a transport U > 0, the value registered at its head for each unit surrendered at its tail. Re-describing value at each position by a factor g changes U to g_head U / g_tail: an abelian lattice gauge theory with group ℝ_{>0}, which the logarithm makes (ℝ, +) and physis's discrete exterior calculus computes ([theory §2](theory/value-in-space-and-time.md)).
  - The invariants are the holonomies of loops, λ = Π U^{±1}, and their **loop deficits** δ = −ln λ (the user's paper on loop deficits, Huang 2026).
  - Every deficit vanishes exactly when a **price system** exists, a potential φ with ln U = φ(head) − φ(tail): on a finite network, the absence of arbitrage.
  - In the tree gauge, the chords of a spanning tree carry the holonomies of their cycles, complete invariants.
- **Flows circulate.** A flow, the units surrendered on each leg per unit of time, is a 1-chain. With a conductance on each leg (for trade, the link's gross volume), it has physis's decomposition by Ohm's and Kirchhoff's laws, on a network whose faces are given ([theory §3](theory/value-in-space-and-time.md)):
  - a driven part, up a single potential (the attraction), each leg carrying it in proportion to its conductance;
  - local circulation round the faces;
  - global circulation round the holes.

  What two positions send each other both ways is exchange, a round trip; what is left, one net flow for each pair, is driven or circulates. Only circulating flow meets a loop deficit: the pairing Σ F ln U splits into the pairings of like parts, whatever the conductances.
- **Geography gives loops faces.** The spherical Delaunay triangles of the places (physis computational geometry) are the faces of a network of places. Curvature then locates a deficit (Stokes), and circulation splits into eddies round neighbouring places and circulation through distant links.
- **Production is a circular flow** ([theory §4](theory/value-in-space-and-time.md)): an input–output table, the Leontief inverse, multipliers, the price dual, and the decomposition of the deliveries between sectors.
- **Time is a direction of transport** ([theory §5](theory/value-in-space-and-time.md)). In the spacetime complex of a record of transactions, morphe's worldline complex, agents hold value between events at their rates of holding and exchange it at events.
  - A lens's curvature is the excess return round it: for currencies, the departure from uncovered interest parity.
  - A triangle's curvature at an event is a triangular arbitrage.
  - Flatness is an intertemporal price system.
  - Reach is time-respecting, with three-valued viability.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `transfer` | Transfer networks: positions, legs with transports and flows, value kinds; holonomies, loop deficits, re-description, circuits in the tree gauge, the price potential, retention round directed cycles ([`transfer/`](transfer/include/economics/transfer/)) | draft |
| `circulation` | The decomposition of flows with conductances and the attraction; reciprocal exchange and net transfer, with each link's gross volume; the pairing of flows with transports; balances by kind; stocks under a column-stochastic rule ([`circulation/`](circulation/include/economics/circulation/)) | draft |
| `production` | Input–output tables, the Leontief model and its price dual, multipliers and linkages, the network of deliveries ([`production/`](production/include/economics/production/)) | draft |
| `spacetime` | The spacetime complex of a record of transactions: holding and exchange transports, curvature of lenses and triangles, reach with Kleene's three-valued viability ([`spacetime/`](spacetime/include/economics/spacetime/)) | draft |
| `places` | Places on the sphere, and the faces their spherical Delaunay triangulation gives a network ([`places/`](places/include/economics/places/)) | draft |
| `markets` | Prices and demand in the textbook sense: revealed preference (Afriat), Slutsky, general equilibrium | planned |
| `distribution` | Inequality and its measures (Lorenz, Gini), income and wealth | planned |

## 3. Conventions

- **Transports** are dimensionless ratios of value registered to value surrendered, and enter the computation as ln U. Flows carry their units (millions of US dollars for trade). Holding rates are per year, continuously compounded: ln(1 + r/100) for a rate of r per cent a year.
- **Orientation.** A leg runs from its tail, which surrenders, to its head, which registers. A face's word runs anticlockwise seen from outside the sphere. A lens runs a[p, q] · rung(q) · b[p, q]⁻¹ · rung(p)⁻¹ (morphe's worldline complex).
- **Data are inputs.** The domain takes C++ structures and reads no files: an application reads the data library's products and hands them over.

## 4. Dependencies

- **Domains:** none, so economics has height 0.
- **Foundations:**
  - physis: discrete exterior calculus (cochains, the Hodge decomposition, dense and by conjugate gradients) and the spherical Delaunay triangulation;
  - morphe: CW complexes, spanning forests, the worldline complex;
  - calx: dense linear algebra and conjugate gradients, random streams for its checks, the report of its checks.

  The domain writes no general method of its own (CLAUDE.md §2.1).

## 5. Code and Checks

The C++ is header-only, C++20, built with CMake (`CMakeLists.txt`), each module a target `economics::<module>`. The checks run by `node domains/economics/tests/run-tests.mjs`:

| Check | What it verifies |
|---|---|
| transfer | The paper's triangle: holonomy 0.60, deficit 0.511, unchanged under its re-description to legs 3.0, 0.08, 2.5. Its five positions: chord holonomies 0.975 and 1.0725. Deficits reverse with the loop and add under concatenation. Chord holonomies are invariant under random re-description while the log prices shift. A network built from prices has no deficit and gives its prices back. Retentions multiply to the holonomy |
| circulation | Flows round an unfilled triangle are wholly global, and wholly local once it is a face. A flow d₀α is driven, with attraction α. A round trip is global circulation. Shares add to 1. Balances by kind add up to the net inflow. Reciprocal plus net equals the whole flow. The pairing's terms add up, and a driven flow meets no deficit. With the gross volumes as conductances, the driven flow is gross × the difference in attraction and carries every net inflow; with any conductances a driven flow meets no deficit. The paper's stock drift follows π − (5/7)(0.3)ᵗ(1, −1) |
| production | Miller and Blair's example: A, L, L y = x, unit prices, multipliers, cost-push. An unproductive table is detected. The driven share of a cycle of deliveries, 250/5425, matches least squares by hand |
| spacetime | A market built from prices is flat. The lens closed form (r_a − r_b)Δt + ln U_q − ln U_p holds, with rates changing at events. A triangle has curvature ln(2·3/5). Stokes holds when one exchange changes. Reach is viable along a chain, unknown at an unknown link, lost at a non-viable one, never backwards in time, and monotone |
| places | On the octahedron's places: eight faces and no global circulation. Stokes holds round the equator. A leg between antipodes opens a hole. Neutral legs are added between neighbours. A region gives no faces across the empty side |

## 6. Open Questions

1. **Fitting transports and flows to noisy records** needs the inference foundation (planned).
2. **Non-abelian transports**: bundles of goods converted by matrices, read through traces and spectra, once a record measures them.
3. **What the user's papers leave imprecise**: where a deficit "sits", the curvature of preferences, the functional of generative justice. The theory note says why each is not computed as if it were.
4. **Real input–output tables** (OECD ICIO) for the production view: large files, not yet in the data library.

## 7. References

- Huang, W. (2026), *The Geometry of Value: Loop Deficits as Invariant Measures of Structural Exploitation*; *Toward a Pathology of Value Circulation: A Typed Lifecycle Diagnostic*.
- Ilinski, K. (2001), *Physics of Finance: Gauge Modelling in Non-Equilibrium Pricing*, Wiley; Young, K. (1999), Foreign exchange market as a lattice gauge theory, *American Journal of Physics* 67, 862–868.
- Jiang, X., Lim, L.-H., Yao, Y. and Ye, Y. (2011), Statistical ranking and combinatorial Hodge theory, *Mathematical Programming* 127, 203–244; Lim, L.-H. (2020), Hodge Laplacians on graphs, *SIAM Review* 62, 685–715.
- Leontief, W. (1986), *Input–Output Economics*, 2nd ed., Oxford University Press; Miller, R. E. and Blair, P. D. (2009), *Input–Output Analysis*, 2nd ed., Cambridge University Press.
- Fama, E. F. (1984), Forward and spot exchange rates, *Journal of Monetary Economics* 14, 319–338 (uncovered interest parity and its failure).
