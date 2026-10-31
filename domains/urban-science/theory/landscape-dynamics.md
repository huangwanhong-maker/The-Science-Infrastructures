# Landscape Dynamics: Districts as Cellular Configurations, Their Events and Their Spacetime

- **Status:** draft
- **Activity:** research
- **Modules:** `landscape` (`landscape/include/urban/landscape/`), checked by `landscape/tests/check_landscape.cpp`; the building of §7 is the architecture domain's `space` and `environment` (`domains/architecture/`)
- **Decision:** [ADR 0062](../../../development/decisions/0062-landscapes-as-cellular-configurations.md)
- **Foundations:** morphe (regular complexes, their rewriting, traces, Möbius inversion; its decision 0008), physis (composition, lumping, jump and piecewise-deterministic processes, elimination, Doi–Peliti operators, Hodge Laplacians, planar faces; its decision 0018), calx (the θ-method, matrix exponentials, random streams)

This note states the formalism the module implements, with its definitions, its results and their proofs or sources, and says where a result is a theorem about any district and where it is a property of a model with declared rates. Notation follows the code: names in `code` are those of the headers.

## 1. The Landscape: Structure, Use and Meaning

A landscape, here, is a place described at three levels that change one another.

- **Structure** is what is built and how it is joined: junctions, streets and blocks; rooms, walls and floors. It is a cell complex, and it changes by *events* that rewrite cells.
- **Use** is what the structure carries and how it moves: walkers at junctions and along streets, the warmth of blocks and rooms, whether a store is open. It is a state on the cells, evolving continuously between events and jumping at them.
- **Meaning** is how a part of the place is read: a rule in force (a permit), the rule carried out (a lane built under it), a meaning adopted (the lane recognised as a passage). It is also state, but of a kind that changes *which events can happen*: a recognised lane can no longer be removed.

The formalism makes each level exact and joins them in one process. Three kinds of emergence can then be stated, not merely named.

1. **Dynamical emergence.** A pattern of use belongs to the whole and to no single cell: the day's tide of walkers across a district (§5, §10.3).
2. **Structural emergence.** The cells themselves have histories: a lane begins at an instant, the street it cuts ends there and two new segments begin (§3, §10.3).
3. **Interpretive emergence.** A reading becomes a state that changes the set of possible futures, not only their rates (§6.4, §10.4).

*Rationale.* Each level has a mature mathematics of its own: cell complexes and their rewriting for structure; Markov and hybrid processes for use; guarded transition systems for rules. The contribution is the joint, with the laws that hold across it: what a rewrite must carry and conserve (§2.3), what a coarser description can forget exactly (§6.3), what parts compose to (§4).

## 2. The District as a Typed Cellular Configuration

### 2.1 Regular Complexes

**Definition 2.1** (regular complex). A finite *regular cell complex* is presented by its face poset: cells $\sigma$ with dimensions $\dim\sigma$, each cell's facets (cells of one dimension less on its boundary) with incidence signs $[\sigma:\tau]\in\{\pm1\}$, such that the closure of every cell is homeomorphic to a closed ball whose boundary is the closure of its facets. Then the complex is determined up to homeomorphism by its face poset, which is a CW poset (Björner 1984, *Eur. J. Combin.* 5, 7).

`morphe::topology::admissible` decides regularity exactly up to dimension 3. A 0-cell has no facets. A 1-cell has exactly two distinct 0-facets. The facets of a 2-cell form a single cycle, in which every vertex lies on exactly two of its edges. The facets of a 3-cell form a closed connected surface (every edge on exactly two of its faces, every vertex link one cycle) of Euler characteristic 2, a sphere. In dimension 4 and above the test is the homology mod 2 of the boundary and is marked inexact. Orientations are propagated across the diamonds of the face poset (`orient`), and the cellular chain complex follows (`chain_complex`).

### 2.2 The District's Complex

A district is cut from the walking network of a city model: the ways where people walk, split where two ways share a point, in the city's exact local frame (`shared.js`, `networkOf`). The segments form a plane graph when no two cross without a shared point. Ways that cross at different levels (a ramp, a bridge) are found pairwise and the shorter of each crossing pair is dropped and named (`planarised`). The bounded faces of the plane graph are traced from the rotation system at each junction (`physis::computational_geometry::plane_faces`). Euler's relation over the face orbits of the rotation system, each component's outer face among them, $V - E + F = 2c$ for $c$ components, is checked, and a network that fails it (one not plane) is refused.

**Definition 2.2** (district). The *district* of a plane walking network is the 2-complex whose 0-cells are the junctions, whose 1-cells are the segments, and whose 2-cells are the bounded faces bounded by a single circuit, each oriented counter-clockwise. A bounded face whose boundary is not a circuit (a dead end or a loop of lanes inside it) is not a 2-cell; it is a *hole*.

**Proposition 2.3.** The district is a regular complex, and $b_1$ equals the number of holes; $b_0$ is the number of pieces of the network; $\chi = b_0 - b_1$.

*Proof.* Each 1-cell joins two distinct junctions (a closed way is split once more, `networkOf`), and each 2-cell's facets form a circuit: regular by Definition 2.1. A 2-cycle would be a combination of blocks with no boundary; but the outermost edge of its support borders one of its blocks only, so $b_2 = 0$. Euler's relation gives $V - E + F_\text{bounded} = c$, so $\chi = V - E + F_\text{cells} = c - (\text{holes})$, and $b_1 = b_0 - \chi = \text{holes}$. ∎

The box Elizabeth–Russell and Flinders–Lonsdale of Melbourne's Hoddle Grid gives 356 junctions, 399 segments and 37 blocks in 6 pieces with 12 holes. The pieces are ways the box cuts off, or ways joined to the streets only by the footways that are left out so that no street counts twice.

### 2.3 The Schema and the Transfer of State

**Definition 2.4** (schema). A *schema* assigns to each type of cell its fibres, each extensive or intensive. Junction: walkers (persons, extensive). Segment: flow (persons an hour, intensive along the segment). Block: temperature (°C, intensive) and heat capacity (J/K, extensive), whose product, heat, is extensive.

**Rule 2.5** (transfer). A rewrite that replaces cells carries each fibre by its kind. An extensive quantity is split in proportion to a measure (area, length) or summed. An intensive one is copied, or averaged with the weights that conserve its extensive partner. A new cell that nothing gives a quantity starts at zero.

**Proposition 2.6.** Under Rule 2.5 every extensive quantity is conserved by every rewrite of §3.4 and §7.3, and every intensive one is kept wherever its cell persists.

*Proof.* Splitting in proportion and summing preserve totals; averaging with weights $C_i$ the intensive $T_i$ gives $\bar T = \sum C_i T_i/\sum C_i$, so $\bar T\sum C_i = \sum C_i T_i$. ∎

The checks verify it on every rewrite: the blocks' areas and capacities summed across a split, temperatures copied, walkers kept (a new junction starts with none), heat kept by every merge and parting of rooms.

## 3. Change as Rewriting

### 3.1 Rewriting a Complex

A complex is rewritten as its *incidence graph*: a vertex for each cell, typed by dimension and label, an edge from each cell to each facet, labelled by the sign. A rule is a span $L \leftarrow K \rightarrow R$ of such graphs, applied at a match $m: L \to G$ by the double pushout (Ehrig, Pfender and Schneider 1973; Ehrig, Ehrig, Prange and Taentzer 2006, *Fundamentals of Algebraic Graph Transformation*, ch. 3). The pushout complement exists exactly when the *dangling condition* holds: no cell deleted is a facet of a cell kept. Graphs with injective matches form an adhesive category (Lack and Sobociński 2005, *RAIRO-ITA* 39, 511), so independent rewrites commute (the local Church–Rosser theorem). The result is read back as a complex and tested for regularity (`morphe::topology::rewrite`).

*Remark.* The automorphisms of $L$ act freely on the injective matches $L \to G$, so matches up to $\mathrm{Aut}(L)$ are orbits of a free action; no stabiliser enters.

The rules the module uses are built at a site rather than written by hand (`complex_rewriting.hpp`): `attach_cell`, `attach_edge`, `remove_cell`, `subdivide_edge`, `merge_at_vertex`, `partition_cell`, `split_face`, `merge_cells`, `relabel`. Cells keep persistent identifiers across rewrites, and each rewrite returns the cells it created and deleted and its *footprint*: what it read and what it wrote.

### 3.2 Refinement and Intervention

**Definition 3.1.** A rewrite is a *refinement* when it describes the same place more finely, and every quantity that held of the cells it replaces holds of their union: subdividing a segment at a point (its cofaces keep $e \mapsto e_1 + e_2$), merging two segments at a junction of degree two. It is an *intervention* when it changes the place: splitting a block by a new passage, merging rooms across a removed wall.

A refinement leaves the homeomorphism type and every cochain's integral unchanged; an intervention need not. The engine labels each step (`Step::refinement`), and the page shows which is which.

### 3.3 Traces and the Order of Events

Two events are *dependent* when one writes what the other reads or writes (Mazurkiewicz 1977; Diekert and Rozenberg 1995, *The Book of Traces*). A record of events $e_1, \dots, e_n$ with footprints $(R_i, W_i)$ has the causal order: $i \prec j$ when $i < j$ and the events are dependent, closed transitively (`morphe::formal::traces::causal_order`). The record is one linearisation of its trace; the trace is the history, and every linearisation of it is the same history. Its Foata normal form groups the events into steps of mutually independent events (Cartier and Foata 1969).

### 3.4 The Lane

A lane between two stores, through the block both face, is three rewrites: the street each store faces subdivided at the point facing it, and the block split between the two new junctions. The two subdivisions are independent; the split reads both. The trace has two linearisations, and the history is one. Closing the lane merges the two blocks across it and the streets at its ends, returning a complex isomorphic to the first. The closed streets get new identifiers, so a plan names the cells of the configuration it was made in; the process makes the plan again after every rewrite (`process.hpp`), and `build_lane` refuses a plan made for another configuration.

## 4. Composition of Open Systems

**Definition 4.1** (open system). An *open system* $S$ has a finite state space $X_S$, ports $p: X_S \to P$ (a port reads a function of the state), and moves: transitions $x \to y$ with rates, grouped into named mechanisms (`physis::statistical_mechanics::OpenSystem`).

**Definition 4.2** (composition). Given $S$, $T$, ports $p_S$, $p_T$ and a relation $\rho \subseteq P_S \times P_T$, the composite has the states $\{(x, y) : (p_S(x), p_T(y)) \in \rho\}$, the fibre product over $\rho$. Each component's moves are lifted where both ends stay in the fibre product; a lifted move that would leave it is *refused at the interface*. *Joint mechanisms*, declared on the components they read, are added at the first composition that holds all their components. With no ports the composite is the free product (juxtaposition).

**Proposition 4.3.** Composition is associative up to the canonical bijection of states, provided each joint mechanism is added once, at the first composition holding all its components.

*Proof.* The fibre product over relations is associative (it is the limit of the diagram of components and relations, which does not depend on the order of the pullbacks). Local moves are lifted identically in either order. A joint mechanism added at the first composition holding its components contributes the same rates whichever composition that is. ∎

The check composes the event skeleton (§6) in two orders and finds the same rate for every mechanism between every two joint states.

*Remark.* This is not the composition of open Markov processes by Baez, Fong and Pollard (2016, *J. Math. Phys.* 57, 033301), where processes are glued by identifying boundary states (cospans of finite sets). Here the parts keep their own state spaces, and the interfaces are relations between them: a synchronised product, closer to the interconnection of behaviours (Willems 2007, *IEEE Control Syst. Mag.* 27, 46).

## 5. The Landscape as a Hybrid Process

The landscape's state is $z = (K, x, \gamma, h)$: the complex $K$; the continuous and discrete state $x$ on its cells; the event kinds $\gamma(z)$ enabled at $z$; the record $h$. It is a piecewise-deterministic Markov process (Davis 1984, *J. R. Stat. Soc. B* 46, 353) on the strata of $K$: a flow between events, events at rates that read the state, and forced jumps where the flow meets the boundary of its domain (`physis::dynamics::PiecewiseDeterministic`).

### 5.1 Between Events

Walkers hop along the streets. A walker at junction $v$ crosses segment $e = vw$ at rate $(G_e/c)\, e^{(u_w - u_v)/2}$, where $G_e$ is the segment's conductance (the walking rate over its length), $c$ a junction's capacity, and $u$ the attraction at the ends of an open store's segment. So

$$\dot n_v = \sum_{e \ni v} (\text{hops in} - \text{hops out}) + e_v(t) - n_v/\tau,$$

with walkers entering at the district's gates at the hour's counted rate $e(t)$ and leaving at $1/\tau$.

**Proposition 5.1.** Without attraction the hops are $-\Delta_0 (n/c)$, the Hodge Laplacian with stars $c$ on the junctions and $G$ on the segments. With it they are a master equation in detailed balance with the law $n \propto c\, e^{u}$, and $n$ stays non-negative.

*Proof.* The hop rates are those of a reversible random walk with weights $c\,e^{u}$: $(G/c) e^{(u_w - u_v)/2} \cdot c\,e^{u_v} = G\, e^{(u_v+u_w)/2}$, symmetric in $v, w$. A linear master equation preserves the non-negative orthant. ∎

A potential-driven flux $G(n_v/c - u_v - n_w/c + u_w)$ would not preserve it: it can draw walkers out of an empty junction. Each segment's hops *each way* are recorded (`Frame::forth`, `Frame::back`); their difference is the net flow, their sum the traffic a counter on the segment would see. Each block warms toward the air and the Sun, $\dot T_f = (T_\text{air} + s\,\sin^+\!a(t) + h\,n_f - T_f)/\tau_f$.

### 5.2 Events

A store opens at $\alpha(1 + \eta[\text{the other open}] + \xi[\text{the lane there}])\,e^{-\kappa(T - T^*)^+}$ during its hours and shuts at $\beta$. The permit is granted at $\kappa_p$. Under it the lane is built at $\kappa_+$ (three rewrites, §3.4). While it is unrecognised it is removed at $\kappa_-$. Recognition comes at $\kappa_d$ times the lane's use, the walkers at its two ends against a reference. At the closing hour the stores close by a forced jump: the clock, part of the state, reaches the boundary of the stratum's domain, and the reset carries the next closing hour. The forced jump is a boundary condition of the generator, not a rate.

### 5.3 The Registry of Enabled Events

The set $\gamma(z)$ of event kinds that can happen is a function of the state: the lane can be built only under a permit and only when it does not stand; it can be removed or recognised only while it stands and is unrecognised. Each frame of a run carries $\gamma$ (`Frame::enabled`). The page shades, for each kind, the hours it was enabled behind the dots of the hours it happened. The check verifies that every frame's registry obeys these guards.

The continuous part keeps its walkers exactly: $\frac{d}{dt}\sum n = \sum e - \sum n/\tau$, and every rewrite carries them across (Proposition 2.6). The check finds the balance at midnight to rounding.

## 6. The Event Skeleton: What Follows Exactly

The discrete part of the landscape, with the walkers held at a reference, is a finite Markov process that can be solved without simulation (`skeleton.hpp`).

### 6.1 The Parts and the Composite

The parts are: $A$ and $B$, the stores (shut 0, open 1); $P$, a rule in force (the permit); $W$, the rule carried out (the lane); $R$, a meaning adopted (the lane recognised). The interfaces are the relations $w \le p$ and $r \le w$: the lane only under a permit, its recognition only while it stands. So a permit cannot be withdrawn under a standing lane, and a recognised lane cannot be removed. Those moves are refused at the interface (Definition 4.2), not forbidden by hand. The mechanisms are declared once on the parts they read:
- a store opens, on $\{A, B, W\}$, at $\alpha(1 + \eta b + \xi w)\pi$, with $\pi$ a promotion;
- a store shuts, on itself, at $\beta$;
- the permit is granted at $\kappa_g$ and withdrawn at $\kappa_r$;
- the lane is built, on $\{P, W\}$, at $\kappa_b$ under a permit, and removed at $\kappa_d$;
- the lane is recognised, on $\{A, B, W, R\}$, at $\kappa_u u(1 + \nu(a+b))$ while it stands, where $u$ is its use with both stores shut and $\nu$ what each open store adds to it; recognition lapses at $\kappa_l$.

There are $2\cdot2\cdot4 = 16$ joint states.

### 6.2 Law and Irreversibility

**Proposition 6.1.** The skeleton is in detailed balance if and only if $\xi_A = \xi_B = \nu = 0$ and $\eta_A = \eta_B$. Its stationary law is then $\pi(a,b,p,w,r) \propto K_A^a K_B^b (1+\eta)^{ab} \cdot \varpi(p,w,r)$, with $K = \alpha/\beta$ and $\varpi$ the institutions' path law $1 : \kappa_g/\kappa_r : \kappa_g\kappa_b/(\kappa_r\kappa_d) : \kappa_g\kappa_b\kappa_u u/(\kappa_r\kappa_d\kappa_l)$.

*Proof.* Kolmogorov's criterion: the rates are reversible pairwise, so detailed balance holds exactly when every cycle's forward and backward products agree (Kelly 1979, *Reversibility and Stochastic Networks*, §1.5). The institutions' states form a path, which has no cycle. Three squares remain. The square of a store's opening and the lane's building has ratio $1+\xi$. The square of a store's opening and the recognition has ratio $1+\nu$. The square of the two stores has ratio $(1+\eta_B)/(1+\eta_A)$. So the conditions are necessary. Under them the product law satisfies detailed balance on every move, by direct substitution. ∎

When the conditions fail, the process runs round its cycles. Its entropy production $\sigma = \tfrac12\sum_{ij}(\pi_i k_{ij} - \pi_j k_{ji})\ln\frac{\pi_i k_{ij}}{\pi_j k_{ji}} > 0$ is computed with the currents and affinities of a spanning tree's chords (Schnakenberg 1976, *Rev. Mod. Phys.* 48, 571). *A meaning adopted through use is a driven process*: with $\nu > 0$ no setting of the other rates restores balance.

### 6.3 What a Coarser Description Forgets

A *description* keeps some parts and forgets others, or keeps only how many stores are open; it partitions the joint states into blocks.

**Definition 6.2.** A partition is *strongly lumpable* when the rate $R(i, B) = \sum_{j \in B} k_{ij}$ into each other block is the same for all $i$ in a block (Kemeny and Snell 1960, *Finite Markov Chains*, §6.3); *eventwise lumpable* when this holds for every mechanism separately.

**Proposition 6.3.** The description by the stores alone is eventwise lumpable if and only if $\xi_A = \xi_B = 0$. The descriptions without the permit or without the recognition are never lumpable.

*Proof.* The institutions' moves stay inside the stores' blocks. A store's opening from $(a{=}0, b, w)$ goes to the block $(1, b)$ at a rate that is independent of $w$ exactly when $\xi = 0$. Building reads the permit, so forgetting the permit confuses states from which building has rate $\kappa_b$ and $0$. Removal is refused while recognised, so forgetting the recognition confuses states from which removal has rate $\kappa_d$ and $0$. ∎

The engine reports, for each description that fails, the mechanism and the two states it confuses (`ViewVerdict`); they differ exactly in what the description forgot. Counting the open stores is lumpable as a whole when the stores are alike, but never mechanism by mechanism, for the count cannot say which store opened.

**Proposition 6.4** (how close a coarse chain stays). Aggregate by the stationary law, $\bar k_{AB} = \sum_{i \in A}\pi_i R(i,B)/\sum_{i\in A}\pi_i$, and let $\delta = \max_i \sum_{B \ne C(i)} |R(i,B) - \bar k_{C(i)B}|$, the *outward discrepancy*. Then over $[0, T]$ the law of the projected fine path and that of the aggregated chain started alike differ in total variation by at most $1 - e^{-\delta T}$. The aggregated chain keeps the stationary law of the blocks.

*Proof.* Couple the two maximally (Lindvall 1992, *Lectures on the Coupling Method*, ch. 1). While both are in the same block, with the fine state at $i$, jump to each other block $B$ together at rate $\min(R(i,B), \bar k_{C(i)B})$, and alone at rate $|R(i,B) - \bar k_{C(i)B}|$. Moves within a block do not separate them. They separate at total rate at most $\delta$, so they are still together at $T$ with probability at least $e^{-\delta T}$. The stationary claim is the balance $\sum_A \bar\pi_A \bar k_{AB} = \sum_i \pi_i R(i,B)$. ∎

So a description stays within 10% for $T \le -\ln 0.9/\delta$, which the page states as a duration. What an inexact description remembers shows in its dwell times. The time a block lasts, entered with the stationary law on its states, is phase-type (Neuts 1981, *Matrix-Geometric Solutions*, ch. 2). It is exponential when the description is exact; otherwise its squared coefficient of variation differs from 1.

### 6.4 Interventions and Their Interactions

Three interventions are considered: promoting $A$, promoting $B$ (each opening multiplied by $1+g$), and granting permits faster. For an outcome $Y$ at stationarity, the *anchored interactions* are the Möbius inverse over the Boolean lattice of the subsets of interventions done, $g(S) = \sum_{T \subseteq S} (-1)^{|S \setminus T|} Y(T)$ (Rota 1964, *Z. Wahrsch.* 2, 340), computed by `morphe::combinatorics::anchored_interactions`.

**Proposition 6.5.** An interaction belongs to the outcome, not to the parts. For stores that answer nothing ($\eta = \xi = 0$), the number open has no interaction of two or three interventions. The chance that both are open has one between the two promotions.

*Proof.* With $\eta = \xi = 0$ the stores are independent and $E[a+b] = E[a] + E[b]$, where each term depends on one promotion only; Möbius inversion of an additive function vanishes on sets of two or more. $P(a=b=1) = P(a=1)P(b=1)$ is multiplicative, and its two-way interaction $(p_A' - p_A)(p_B' - p_B) \ne 0$. ∎

Where the stores answer the lane, permits granted faster interact with the promotions through the lane they answer. Where the lane's recognition is drawn by the open stores, the promotions reach the lane's meaning.

### 6.5 Evidence

A record in which every change is observed has an exact path density under the generator: the product of the rates of the changes taken and of the survival over each holding time (`coarse_log_density`). Over a grid of common couplings $\eta$ with a prior, it gives the posterior. Propositions are weighed through the worlds in which they hold, so $P(\eta \ge 1) \le P(\eta > 0)$ always.

## 7. A Building as a Complex of Rooms

The building belongs to the architecture domain, a domain of its own (ADR 0060 §2.8). Its modules `space` (`rooms.hpp`) and `environment` (`heat.hpp`) implement this section, and the urban-science domain does not hold them.

### 7.1 The Complex

A building is fitted with a box of $n_x \times n_y$ rooms on each floor, the cubical grid of `morphe::topology::cubical_grid`. Rooms are 3-cells. The 2-cells are typed by where they lie: ground, roof and slabs between floors (horizontal), facades outside and walls between rooms (vertical). It is a ball: $b = (1, 0, 0, 0)$.

### 7.2 Heat as a Hodge Laplacian

With $Q_c = C_c T_c$ the heat of room $c$, and $G_f = U_f A_f$ the conductance of face $f$, conduction is $\dot Q = -d\,G\,d^{\mathsf T} T$, with $d$ the coboundary from faces to rooms. This is exactly $-\Delta_3 Q$ for the stars $w_f = 1/G_f$ and $w_c = 1/C_c$: *the constitutive law is the Hodge star* (Tonti 1975; Bossavit 1998, *Computational Electromagnetism*, ch. 5). An outer face conducts to the outside. $\Delta_3$ is positive definite (no harmonic 3-cochain on a ball), so each mode decays at its own rate: the building's time constants. Each hour is propagated exactly with its forcing held, by the exponential of the augmented matrix (Van Loan 1978, *IEEE TAC* 23, 395).

The Sun's direct normal irradiance on a clear day is $I = 1353 \cdot 0.7^{m^{0.678}}$ W/m², $m = 1/\sin a$ (Meinel and Meinel 1976, *Applied Solar Energy*, §2). What reaches a room through a face of outward normal $n$ is the sol-air gain (ASHRAE, *Handbook — Fundamentals* 2017, ch. 18): $A I \cos\theta\,[\,g\,\mathrm{SHGC} + (1-g)\,U\alpha/h_o\,]$, with $g$ the glazed share (facades only), $\alpha$ the absorptance and $h_o \approx 25$ W/m²K the outer surface's coefficient. Most of what opaque fabric absorbs is lost to the air outside; only $U\alpha/h_o$ of it, a few per cent, enters. Diffuse light and shading by neighbours are left out.

### 7.3 Events in a Building

A door or a stair is a wall or a floor relabelled: the same cell carrying another conductance, and a way through for people. Merging two rooms across a wall is `merge_cells`: one 3-cell with the facets of both but the wall, its heat and capacity summed. Parting it again along the first room's facets copies the merged temperature to both, each with its share of capacity. By Proposition 2.6 every event keeps the building's heat, and every derived complex is regular. The rooms one can reach from the entrance are found by a search through doors and stairs.

## 8. Scales, Memory and Two Descriptions

### 8.1 The District's Frequencies and Its Flows Round the Holes

With the walking stars (weight 1 on a junction, $G_e$ on a segment), $-\Delta_0$ is the generator of one walker hopping along the streets. Its eigenvalues are the rates at which unevenness of crowding relaxes. Its slowest eigenvector on the largest piece of the network divides the district along its weakest joins, as the Fiedler vector of a weighted graph does.

**Proposition 8.1** (Hodge). In each degree $p$, $\dim\ker\Delta_p = b_p$ whatever the positive stars (Eckmann 1944; Lim 2020, *SIAM Rev.* 62, 685).

The harmonic 1-cochains $\varphi$ satisfy $d_1\varphi = 0$ and $d_0^{*}\varphi = 0$. Their fluxes $G\varphi$ are flows that neither gather nor lose walkers at any junction and circulate round no block. They circulate round the holes, one independent flow for each. No potential (an attraction, a gradient of crowding) can drive them, and nothing at a single junction can undo them.

### 8.2 What the Streets Remember of the Shops

Keep the walkers at the junctions, $x$, and eliminate those inside the shops, $y$:

$$\dot x = A x + B y,\quad \dot y = C x + D y \;\Longrightarrow\; \dot x(t) = A x(t) + B e^{Dt} y(0) + \int_0^t K(t-s)\,x(s)\,ds,\quad K(\tau) = B e^{D\tau} C,$$

exactly (Zwanzig 1961, *Phys. Rev.* 124, 983; Mori 1965, *Prog. Theor. Phys.* 33, 423; `physis::dynamics::LinearElimination`). Each shop opens on its nearest junction. A passer-by enters at a declared rate, stays a declared time by the shop's use, and comes out where they went in. Shops alike at one junction are one compartment, exactly: the sum of identical linear compartments obeys the same law. The *Markovian limit* $A - BD^{-1}C$ returns every visitor at once and keeps all walkers on the streets.

The district with its shops is a master equation in detailed balance, with weight $1$ at a junction and rate × stay in a compartment. In the variables $u = p/\pi$ it is $\Pi\dot u = -\mathcal A u$, with $\mathcal A = -L\Pi$ symmetric positive semidefinite. Crank–Nicolson (calx's θ-method, conjugate gradients on the sparse rows) propagates it and keeps $\sum\pi u$, the walkers, to the solver's tolerance.

*Result (a model).* In the district of §2.2, with its 151 holdings that passers-by walk into, a hundred walkers set down at the junction with most shops spread over the district within minutes. The shops take in up to about a sixth of them within ten minutes, and the streets go on changing for an hour or more as visitors go in and out: two time scales, where the forgetful description has one.

### 8.3 Walkers as Particles and as a Flow

Walkers hopping among junctions form a reaction network of first-order reactions $X_v \to X_w$. Its Doi–Peliti operator $\sum_r \kappa_r (a^\dagger_w - a^\dagger_v) a_v$ on the Fock space of counts (Doi 1976, *J. Phys. A* 9, 1465; Peliti 1985, *J. Phys.* 46, 1469) equals the generator of the master equation on the counts reachable, entry by entry (`physis::chemistry::network_operator`).

**Proposition 8.2.** For first-order reactions the mean counts obey the rate equation exactly.

*Proof.* The generator maps the linear functions of the counts to linear functions, so the first moments close. ∎

The flows drawn on the map are therefore the exact means of walkers counted one by one, and the checks find the operator equal to the master equation with deviation 0.

## 9. Resources Carried and Converted

The balance $\dot m = (B \otimes I_s)\,j + N r + u$ applies to $s$ typed stocks at $n$ places. Here $j$ is the flows along the links, $N$ the stoichiometric matrix of the conversion channels, and $u$ the inputs from outside. A weighting $w \in \mathbb Z^{ns}$ is conserved by every transport and conversion exactly when $w^{\mathsf T}[B \otimes I_s \mid N] = 0$; then $\frac{d}{dt}w^{\mathsf T}m = w^{\mathsf T}u$. The conserved weightings form a lattice, the left kernel over $\mathbb Z$, found from the Smith normal form (`value.hpp`, morphe).

The kinds of change differ in what they conserve. Moving a good conserves the good. A conversion conserves what its column is orthogonal to. Copying what is known (the source keeps it, the recipient gains it) conserves no positive count of holders: knowledge is not a stock. A sale entered in a ledger conserves claims less debts, not their number.

## 10. The Record and Its Spacetime

### 10.1 The Record

A run records the initial state, then each event's time, mechanism, match, footprint and, for a lane, its rewrites. With the flow between events, the record determines the path. It is anchored, so a path can be replayed from its record and checked.

### 10.2 Causal Order

The record's events are ordered by §3.3: a store's opening after the lane it answered to was built, the lane after the permit, the recognition after the lane. Events that touch nothing in common are unordered. The record is one of the linearisations of its trace, and the page counts them when they are few.

### 10.3 The Spacetime of a District

The page draws the day with time as height above the district, 10 m an hour from midnight at the ground (`layers/spacetime.js`).
- **The flow.** At every hour the street network is a thin slice, each street coloured by its traffic (both ways). Up the stack, the day's tide of walkers is a pattern of the whole.
- **The structure.** What has a history stands over the hours it exists. The two stores rise as worldlines, thick while open. The lane's sheet begins in mid-air at the instant it is built, with the worldlines of the junctions it made, and the blocks are outlined at that instant, the new ones in sand.
- **The meaning.** The lane's sheet deepens from the instant it is recognised.
- **The order.** Events are beads at their place and hour, joined by threads to the events that read what they wrote: the causal order drawn in spacetime. The past of an event is everything that can be reached down its threads.
- **The present.** A plane marks the clock's hour; beneath it the day as it went, above it, fainter, what the model goes on to do.

### 10.4 Meaning as State

A recognised lane cannot be removed: in the skeleton, the move is refused at the interface $r \le w$, and in the running process, removal is not enabled while it is recognised. An interpretation, once adopted, is a state that changes the set of possible futures. This is why no coarser description can forget it (Proposition 6.3), and why, when it is drawn by use, the whole is driven (Proposition 6.1).

## 11. What Is Proved and What Is Modelled

**Proved, for any district or building.** Proposition 2.3, Proposition 2.6, the commuting of independent rewrites (§3.1), Proposition 4.3, Proposition 5.1, Propositions 6.1–6.5 (for the skeleton as declared), Proposition 8.1, the exactness of the elimination and of the grouping (§8.2), Proposition 8.2, and the conserved lattice (§9).

**Modelled, with declared rates.** The stores' openings and the institutions' rates, the walkers' attraction and stay, the blocks' warming, the shops' stays and the rate of entering, the building's fabric. These explore how structure, use and meaning change one another. They do not predict a particular district.

**Measured.** The walking network (OpenStreetMap), the holdings and their uses (Street Oceanography's), the walkers entering each hour (half the City of Melbourne's counts at the sensors inside the district), and the Sun (the astronomy domain).

**Open questions.**
1. Learning the rates from records in which some changes are unseen (hidden-state inference on the skeleton).
2. The lane's recognition as a field over many lanes, with recognition spreading along the causal order.
3. Districts whose boundary moves: walkers leaving and entering through a box that is itself rewritten.
4. Composition of districts (Definition 4.2) at shared streets, and the lumpability of a city into its districts.
