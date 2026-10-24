# Value in Space and Time

*Status: draft, 2026-10-24.* The theory of the economics domain: value carried between positions, the invariants of carrying it round, the circulation of its flows, production as a circular flow, and the spacetime of transactions. It realises the precise core of the user's papers on the geometry of value (Huang 2026) with the foundations' mathematics, and says where the papers' wider claims are not yet precise. Serves research and learning.

## 1. Positions, Legs and Transports

**Definition 1.1.** A *transfer network* N = (V, E, U) is a finite directed multigraph with positions V and legs E, each leg e = (u, v) carrying a *transport* U_e > 0: the value registered at v for each unit surrendered at u. Its *log transport* is the 1-cochain ω = ln U on the network as a CW complex (physis `discretisation`).

A position may describe value in its own unit. Changing the unit at each position by g_v > 0 changes every transport to

U′_e = g_v U_e g_u⁻¹,   ω′ = ω + d₀(ln g),

a gauge transformation of a lattice gauge theory whose group is ℝ_{>0} ≅ (ℝ, +). Nothing in a transfer network that a position cannot see from its own books is changed by it.

**Definition 1.2.** The *holonomy* of a closed walk γ = e₁^{ε₁} ⋯ e_k^{ε_k} is λ(γ) = Π U_{e_i}^{ε_i} = exp ∫_γ ω, and its *loop deficit* is δ(γ) = −ln λ(γ) = −∫_γ ω.

**Proposition 1.3.** Holonomies are gauge invariant; δ(γ⁻¹) = −δ(γ); deficits add under concatenation at a shared position. The holonomies of the fundamental cycles of a spanning forest (one for each leg outside it, a chord) determine every holonomy. The description classes, transport systems up to gauge, form Hom(H₁(N; ℤ), ℝ_{>0}) ≅ ℝ_{>0}^{b₁}, b₁ = |E| − |V| + c.

*Proof.* ∫_γ d₀ f = 0 round any closed walk, which gives invariance. Reversal and concatenation are properties of the integral. In the gauge where the forest's legs carry 1 (integrate ω along the forest), each chord carries the integral round its fundamental cycle, and these cycles are a basis of H₁. Two systems with the same chord values differ by an exact cochain, and every assignment to the chords is realised. ∎

**Proposition 1.4 (prices).** Every loop deficit of N vanishes iff there is a potential φ: V → ℝ with ln U_e = φ(head) − φ(tail) for every leg; φ is unique up to a constant on each component.

*Proof.* If ω = d₀φ, its integral round every closed walk is 0. Conversely, integrate ω along a spanning forest from each root; a closed ω (all deficits zero) agrees with d₀ of that integral on the chords too. Uniqueness: d₀ψ = 0 iff ψ is locally constant. ∎

*Remark.* φ(v) is the logarithm of the units of value at v for a common unit of account: a price system, with transports the ratios of prices. On a finite network, Proposition 1.4 is the fundamental theorem of asset pricing in its simplest form: no arbitrage round any loop iff a consistent system of prices. The same mathematics describes currency exchange (Ilinski 2001; Young 1999); what the user's papers add is the reading of U as any transfer of value, wages, fees, gifts and conversions among them.

**Proposition 1.5 (retention).** Let γ be a directed cycle of legs e₁, …, e_k, with positive flows x_i surrendered on e_i at its tail. The position entered by e_{i−1} and left by e_i receives y_i = U_{e_{i−1}} x_{i−1} and passes on x_i. The retentions r_i = y_i / x_i satisfy Π r_i = λ(γ).

*Proof.* Π r_i = Π U_{e_{i−1}} · Π x_{i−1} / Π x_i = λ(γ). ∎

## 2. What the Papers Leave Imprecise

The domain computes what is precise and says so where a claim is not yet a definition.

- **Where a deficit sits.** By Proposition 1.5, any positive retentions whose product is λ are realised by some flows: the split of a deficit among positions is a property of the flows, not of the network. The papers' "position that absorbs the deficit" is therefore not an invariant, and is not computed as one.
- **Non-abelian transports.** For bundles of goods converted by matrices, only conjugation invariants (traces, spectra) of holonomies at one basepoint are invariant. A commutator of holonomies at two basepoints is not, and is not computed.
- **Preferences.** The "curvature of preferences" and "prices as Goldstone modes" of the field-theoretic paper are, respectively, a quantity that depends on the units of goods and the first-order condition of optimisation. Textbook microeconomics is planned (`markets`), without those readings.
- **Generative justice and closure.** The functional of generative justice and the closure of a conversion are not specified enough to compute; they remain the user's open programme.

## 3. Flows and Their Circulation

### 3.1 The Hodge Decomposition of Flows

A flow F, units surrendered per unit of time on each leg, is a 1-chain. Its boundary ∂F = d₀ᵀF is the net inflow at each position (inflow − outflow). A chain is not a cochain, and to decompose it one needs a metric. For flows the natural one is a *conductance* c_e > 0 on each leg, its capacity to carry flow (Ohm's law; Kirchhoff 1847). The *drop* E = F / c is a 1-cochain. It decomposes orthogonally in ⟨a, b⟩_c = Σ_e c_e a_e b_e on the network's complex (physis `flow_decomposition`; Lim 2020, Thm 5.3), and multiplying back by c splits the flow:

F = c d₀α + c δ₁γ + c h,

orthogonal in the energy Σ_e F_e G_e / c_e, whose shares the three parts take.

- **The driven part c d₀α.** Along each leg it is c_e (α(head) − α(tail)): value runs up the attraction α, each leg carrying it in proportion to its conductance. It carries all of the net inflow, since the other two parts are divergence-free (d₀ᵀ W₁ δ₁γ = d₀ᵀ W₁ h = 0, Kirchhoff's current law): ∂F = d₀ᵀ c d₀α. It is the flow the balances alone would drive through the conductances.
- **Local circulation c δ₁γ = d₁ᵀγ**, round the faces: the same amount on every leg of a face (face weights 1).
- **Global circulation c h**, round the holes no face fills (of dimension b₁ of the complex).

A flow is pure circulation iff every position is balanced. With every c = 1 this is the decomposition of HodgeRank (Jiang et al. 2011).

*Remark.* Without conductances, the same difference in α drives the same flow along every leg, whatever its size. A small link between a small and a large position then carries a large driven flow, cancelled by an equal and opposite circulation, and the two parts each exceed the flow they make up. With a link's gross volume as its conductance, α is a share of gross flow without unit: the driven net flow between two positions is their gross flow times their difference in attraction.

**Reciprocal exchange and net transfer.** Flows both ways between two positions are, up to the smaller of the two, a round trip, circulation round a 2-cycle. Writing F_uv for the total from u to v,

Σ (F_uv + F_vu) = 2 Σ min(F_uv, F_vu) + Σ |F_uv − F_vu|,

reciprocal exchange plus net transfer. The decomposition is then applied to the net network, one leg for each pair along its net flow, with the same net inflows, each leg conducting its pair's gross flow F_uv + F_vu.

**Proposition 3.1 (only circulation meets a deficit).** Let the drop of F and ω = ln U be decomposed in the same inner product ⟨·,·⟩_c. Then the pairing of the chain F with the cochain ω, which needs no metric, splits into the pairings of like parts:

Σ_e F_e ω_e = Σ_e (c d₀α_F)_e (d₀α_ω)_e + Σ_e (c δ₁γ_F)_e (δ₁γ_ω)_e + Σ_e (c h_F)_e (h_ω)_e.

In particular a driven flow (F = c d₀α_F) meets only the price system of ω, never its loop deficits, whatever the conductances.

*Proof.* Σ_e F_e ω_e = ⟨F/c, ω⟩_c, and the three subspaces are mutually orthogonal in ⟨·,·⟩_c, so the cross terms vanish. The deficits of ω are its integrals round cycles, which only its co-exact and harmonic parts carry (Proposition 1.4). ∎

So value is extracted round loops only by what circulates: a pattern of trade fully explained by balances, however unequal, moves value along a price system and realises no deficit.

### 3.2 Places and Faces

A transfer network among places has loops but no faces. Geography supplies them: the spherical Delaunay triangulation of the places (physis `computational_geometry`, the convex hull of their unit vectors), each triangle a 2-cell bounded by legs between its corners. A neutral leg (U = 1, no flow) is added where neighbours have none; for flows it needs a conductance, and the atlas gives it the least gross flow recorded, joining such neighbours at the resolution of the record, so that almost nothing circulates through it.

- **Curvature.** The curvature d₁ω on a triangle is the integral of ω round it. By Stokes, the deficit round any loop of neighbours is minus the sum of the curvatures it encloses.
- **Circulation** splits into eddies round neighbouring places and global circulation through the loops of distant links, which no triangle fills.
- **Regions.** When the places do not surround the sphere's centre, the hull triangles whose planes leave it outside cover the empty side, not the region, and are not faces.

### 3.3 Stocks

When each position passes shares M_ij of its stock in each period (columns summing to 1), stocks follow k(t + 1) = M k(t). The total is conserved, and for irreducible aperiodic M the stocks tend to the stationary π, M π = π. With M = [[0.9, 0.6], [0.1, 0.4]] and k(0) = (1, 1),

k(t) = π − (5/7)(0.3)ᵗ (1, −1),   π = (12/7, 2/7):

stocks drift with no loop and no deficit. This is the paper's own warning (its Counterexample 8.2) that accumulation is not exploitation.

## 4. Production

An input–output table records deliveries z_ij from sector i to sector j, final demand y and value added v. Gross output x satisfies x_i = Σ_j z_ij + y_i = Σ_j z_ji + v_i. Write a_ij = z_ij / x_j and L = (I − A)⁻¹ (Leontief 1936; Miller and Blair 2009, ch. 2).

- **Quantities.** x = L y. Output multipliers are the column sums of L.
- **Prices.** p = p A + c with c_j = v_j / x_j, so p = c L; every price is 1 in a table in money, and a change in value added propagates as cost-push.
- **Productivity.** For non-negative A, L ≥ 0 iff every leading principal minor of I − A is positive (Hawkins and Simon 1949).

The deliveries are a flow among sectors. Their split into reciprocal exchange, a driven part and circulation measures Quesnay's and Leontief's circular flow: how much of production runs round its own cycles.

## 5. Spacetime

### 5.1 The Spacetime Complex

A record of transactions has agents, events in time with the agents present, and exchanges between two agents at an event. Its encounters make morphe's worldline complex W and its refinement Ŵ:

- **vertices:** the presences (x, v);
- **edges:** each agent's holding between consecutive events, and each exchange (a rung, from the pair's first agent to its second);
- **2-cells:** the lenses (two agents between consecutive exchanges, a[p, q] · r(q) · b[p, q]⁻¹ · r(p)⁻¹) and the triangles of exchanges at an event.

Every edge carries a transport:

- a holding from p to q, exp ∫_p^q r_x(t) dt, where r_x is x's continuously compounded rate of holding (interest, inflation, depreciation);
- an exchange, its rate.

Ŵ is then a transfer network on the presences, with faces, and §§1–3 apply to it.

**Proposition 5.1 (the lens).** If a and b exchange a → b at events p and q with rates U_p and U_q, the curvature of their lens is κ = ∫_p^q (r_a − r_b) dt + ln U_q − ln U_p.

*Proof.* The integral of ω round a[p, q] · r(q) · b[p, q]⁻¹ · r(p)⁻¹. ∎

κ is the log excess return of holding with a and exchanging at q over exchanging at p and holding with b. A triangle of exchanges at an event has the curvature of its triangular arbitrage.

**Corollary 5.2.** Ŵ is flat iff there is an intertemporal price system: φ(x, v), the log price of x's value at v in a common numeraire, with every holding and every exchange its difference. On a finite record, flatness is the absence of arbitrage across agents and dates.

### 5.2 Currencies

Take the agents to be currencies, each held at its central bank's policy rate, and exchanged with the euro each month at the month's average rate S (units per euro). Proposition 5.1 gives the curvature of the lens of the euro and currency X over a month:

κ = (r_EUR − r_X) Δt + ln S_q − ln S_p.

This is the departure from uncovered interest parity, which predicts E[ln S_q − ln S_p] = (r_X − r_EUR) Δt. The curvature's mean is zero under parity. Summed over a currency's lenses it is the cumulative excess return of the euro over X, by Stokes the integral round the loop of all its lenses. Its persistent sign for high-yielding currencies is the forward premium puzzle (Fama 1984).

The exact part of the transports' Hodge decomposition is the price system closest to the record, and the norm of the rest, ‖ω − d₀α‖, its distance from the record.

*Caveat.* The share of ‖ω‖² in the exact part is no measure of how much a price system explains. A change of a currency's unit adds d₀β to ω (§1) and changes the exact part alone, so the share depends on the units in which the currencies are quoted; with rungs ln S of order 1 to 7 and holdings of order 10⁻³, it is near 1 for any record. Only gauge invariants measure the record: the curvatures, the distance (ω and the potential change together), and their ratios. The atlas measures the curvature against that of the exchange rates alone, κ₀ = ln S_q − ln S_p, the record's curvature were nothing earned by holding. The share of it that the rates of holding remove is 1 − Σκ² / Σκ₀²; parity with perfect foresight would make it 1.

### 5.3 Reach

Value surrendered at a presence travels forwards along its holder's worldline and across exchanges in their direction. With each exchange viable, unknown or not, a route's viability is the least of its steps and a presence's the greatest over its routes, in Kleene's strong three-valued logic (no < unknown < yes). This makes the user's regenerative reach precise, with its two open choices fixed.

**Proposition 5.3.** Lowering one exchange's viability never raises any presence's.

*Proof.* Every route's viability is a minimum of its steps, monotone in each; the reach is a maximum over routes, monotone too. ∎

## 6. What the Data Say

Computed by the Economy Atlas's engine on the data library's products, 2026-10-24.

### 6.1 The World's Trade

The OECD's balanced trade in goods, adjusted for re-exports, among 199 economies, each link conducting its gross trade. The economies' places give 394 spherical Delaunay triangles, which cover the sphere; 56 of their sides join neighbours with no recorded trade.

| Year | Exports | Reciprocal | Driven | Local | Global | α: USA | China | Germany |
|---|---|---|---|---|---|---|---|---|
| 2010 | 14.5 tn $ | 68.4% | 39.3% | 1.3% | 59.4% | −0.04 | −0.38 | −0.31 |
| 2014 | 18.0 tn $ | 67.8% | 35.7% | 1.2% | 63.1% | −0.02 | −0.35 | −0.34 |
| 2019 | 18.0 tn $ | 69.8% | 38.3% | 1.2% | 60.5% | +0.02 | −0.36 | −0.30 |
| 2020 | 16.7 tn $ | 69.9% | 39.8% | 1.3% | 58.9% | +0.05 | −0.37 | −0.28 |
| 2022 | 23.8 tn $ | 66.9% | 44.0% | 1.4% | 54.6% | +0.05 | −0.40 | −0.20 |
| 2024 | 22.8 tn $ | 67.7% | 42.3% | 1.3% | 56.4% | +0.04 | −0.43 | −0.26 |

Reciprocal is the share of exports exchanged both ways; driven, local and global are the shares of the net flows' energy.

- **Exchange.** About two thirds of the world's exports of goods are exchanged both ways between the same two economies.
- **Drive.** The balances drive between 36% and 44% of the net flows' energy. The rest circulates, almost all of it through distant links: the balances explain less than half of who sends value to whom.
- **Attraction.** Among the economies trading over 20 bn $, the exporters of oil, gas and ores have the lowest attraction every year (Qatar, Kuwait, Angola, Saudi Arabia, Guinea: −0.55 to −0.84). The highest belong to Panama and Liberia, whose registries record ships as merchandise imported by the flag state, and Lebanon.
- **The United States.** Its attraction is near the world's middle. Its large driven inflow comes from its partners' low attraction (China's, Germany's) times their large gross trade with it, not from a high attraction of its own.
- **Eddies.** The largest circulation round a triangle of neighbours is, in most years, Korea → Taiwan → China → Korea: 10 bn $ in 2010, 26 bn $ in 2024. Belgium, Germany and the Netherlands (12 to 15 bn $ since 2022) and Denmark, the United Kingdom and the Netherlands follow. Local eddies are only about 1.3% of the energy, since 394 triangles span few of the net network's cycles, but they are where neighbours' trade runs round.

### 6.2 Currencies

Twenty currencies, 2015-01 to 2026-06 (138 months, 137 monthly lenses; the yen is left out, its policy rate missing for twenty months). Each is held at its central bank's policy rate and the euro at the ECB's deposit facility rate.

- **Month by month.** The lenses curve by 2.15% a month (root mean square); the exchange rates alone would curve them by 2.17%. The rates of holding remove 2.4% of the exchange rates' curvature: month by month, interest differentials offset almost none of the movement of exchange rates. The lira is the exception (11%).
- **On average.** The euro's mean excess return a year is negative for the currencies of higher rates: the Mexican peso −4.8%, the real −3.2%, the koruna −3.0%, the yuan −2.5%, the rand −2.5%, the zloty −2.2%, the forint −2.2%, the dollar −1.5%. Holding them at their rates beat holding euros, after exchange. It is positive for the lira (+9.2%, its fall outrunning its rates), the won (+1.7%) and the krona (+1.2%). This is the forward premium puzzle (Fama 1984).
- **Significance.** Over 11.4 years the standard error of a mean is the volatility over √11.4 ≈ 3.4. Only the koruna (2.9 standard errors) and the lira (2.0) depart from zero by two or more. So does the krone (+0.31% a year, at a volatility of 0.2%): pegged to the euro, it bears the spread by which Danmarks Nationalbank kept its rate below the ECB's. The carry of the others is suggestive, not established, over this span.
- **Which rate.** The BIS's series for the euro area follows the main refinancing rate until August 2024, and the deposit facility rate from September 2024. In 116 of the 138 months it is 25 to 50 basis points above the deposit facility rate, at which euro money has been remunerated since the Eurosystem's excess liquidity of 2015. With it, every excess return of the euro rises by about 0.37% a year, and the krone's more than doubles. A rate of holding must be the rate holders earn.

## 7. Notation and Code

| Symbol | Meaning | Code |
|---|---|---|
| N = (V, E, U) | a transfer network | `economics::TransferNetwork`, `Leg` |
| ω = ln U | the log transport | `TransferNetwork::transport` |
| λ(γ), δ(γ) | holonomy, loop deficit | `holonomy`, `deficit` |
| chords, φ | the tree gauge, the price potential | `circuits`, `Circuits` |
| r_i | retention | `retention` |
| c, E = F / c | conductance, drop | `circulation`, `NetTransfer::gross` |
| F = c d₀α + c δ₁γ + c h | the circulation of a flow | `circulation`, `FlowCirculation` (physis `flow_decomposition`) |
| reciprocal, net | exchange and net transfer | `net_transfer` |
| Σ F ω | the pairing | `pairing` |
| M, π | stocks under a stochastic rule | `stock_path`, `stationary` |
| A, L, p | the Leontief model | `leontief`, `InputOutputTable`, `LeontiefModel` |
| Ŵ, κ | the spacetime complex, curvature | `spacetime`, `Spacetime`, `Transactions` |
| reach | three-valued reach | `reach`, `Viability` |
| faces of places | the map's triangles | `with_map_faces`, `Place` (physis `spherical_delaunay`) |

## 8. References

- Candogan, O., Menache, I., Ozdaglar, A. and Parrilo, P. A. (2011), Flows and decompositions of games: harmonic and potential games, *Mathematics of Operations Research* 36, 474–503.
- Fama, E. F. (1984), Forward and spot exchange rates, *Journal of Monetary Economics* 14, 319–338.
- Hawkins, D. and Simon, H. A. (1949), Some conditions of macroeconomic stability, *Econometrica* 17, 245–248.
- Huang, W. (2026), *The Geometry of Value: Loop Deficits as Invariant Measures of Structural Exploitation*; *Toward a Pathology of Value Circulation: A Typed Lifecycle Diagnostic*; *Symmetry, Dynamics, and Emergent Justice: A Field-Theoretic Framework for Economics*; *Generative Relational Economics*, ch. 15.
- Ilinski, K. (2001), *Physics of Finance: Gauge Modelling in Non-Equilibrium Pricing*, Wiley.
- Jiang, X., Lim, L.-H., Yao, Y. and Ye, Y. (2011), Statistical ranking and combinatorial Hodge theory, *Mathematical Programming* 127, 203–244.
- Kirchhoff, G. (1847), Ueber die Auflösung der Gleichungen, auf welche man bei der Untersuchung der linearen Vertheilung galvanischer Ströme geführt wird, *Annalen der Physik und Chemie* 72, 497–508.
- Kleene, S. C. (1952), *Introduction to Metamathematics*, North-Holland, §64.
- Leontief, W. (1936), Quantitative input and output relations in the economic system of the United States, *Review of Economics and Statistics* 18, 105–125.
- Lim, L.-H. (2020), Hodge Laplacians on graphs, *SIAM Review* 62, 685–715.
- Miller, R. E. and Blair, P. D. (2009), *Input–Output Analysis: Foundations and Extensions*, 2nd ed., Cambridge University Press.
- Young, K. (1999), Foreign exchange market as a lattice gauge theory, *American Journal of Physics* 67, 862–868.
