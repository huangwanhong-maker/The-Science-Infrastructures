---
title: Event-Based 3D GIS on a 2-Complex: A Spin-Foam Approach. The formalism and its engineering
status: seed
date: 2026-11-04
related: [urban-science, architecture, the Urban Science Workbench (now the Urban Building Workbench), morphe topology (worldlines, regular complexes, rewriting), morphe formal (traces), morphe algebra (characters, group algebras), physis gauge (contraction), physis statistical_mechanics, physis quantum_gravity, the history domain, the economics domain, the data library]
---

# Event-Based 3D GIS on a 2-Complex: A Spin-Foam Approach

*This note designs two things before anything is built: the formalism of the planned paper, and what the Urban Science Workbench will implement of it.*

- **Part I** states the formalism in full, as the paper is to present it.
- **Part II** says what is implemented and with which real data. It also says where each piece sits in the infrastructure, what it costs, and how it is drawn.
- **Part III** outlines the paper.

The note starts from two sources:
- the user's discussion of the paper, kept [beside this note](2026-11-04-event-based-gis-on-a-2-complex-discussion.md) (here *the discussion*);
- the history domain ([theory note](../../domains/history/theory/the-history-complex.md)), whose carrier and spin foam this note generalises.

It takes its own position where it finds reason to; §0.2 lists what it changes in the discussion.

## 0. The Position in Brief

### 0.1 Eight Claims

1. **Three representations.** A city is held in three ways:
   - *geographically*, as its states $S_t$ in space with persistent identities: what a 3D GIS shows;
   - as its *event complex*, a 2-complex whose vertices are events, edges couplings and faces primitives: the user's spin-foam representation;
   - as its *foam*, the same complex labelled by modes and amplitudes.

   Exact maps join them: the construction $H$ from a record, the collapse of the refined complex, the cut at a date, the realization $\mathcal R$, and the transform to modes (§2, §4, §5).
2. **The user's ontology is the spin foam's own, read literally.** A face is the history of a link of the spin network of the moment, an edge the history of a node, a vertex an event. The result is a CW complex, and not a hypergraph, exactly when every primitive joins two couplings at each instant (§1.2). Three rules follow:
   - binary relations and interfaces are primitives;
   - a relation among three or more is held by a coupling of its own;
   - a property of one thing is a state of its coupling.
3. **For a city, the spin network of the moment is the dual of its decomposition of space.**
   - Volumes are the nodes and couplings: buildings, rooms, the spaces of streets and junctions.
   - Interfaces are the links and primitives: party walls, frontages, the cross-sections of streets, doors.

   This matches loop quantum gravity, where nodes are quanta of volume and links quanta of area. The discussion's construction A (the geometric dual of a cellular spacetime) and its construction B (the event-history complex) are then the same complex, up to refinement and with the plaquettes kept only where an event constrains them (Theorem 2.8).
4. **Spins are modes, not categories.**
   - An event's category is a rule: its footprint, its kernel, its rate.
   - A spin labels a mode of a primitive's state under the symmetry of its dynamics: an isotypic component of its state space or, with no symmetry, an eigenmode of its generator.

   The spin of a land-use class is the level of the classification at which the class is still remembered (Proposition 3.2).
5. **The amplitude is a likelihood, factorised over the complex.** Under the event-coupled Markov model, where cells interact only at events, the probability of a record is the contraction of a tensor network on the event complex (Theorem 4.2). Persistence sits on faces and edges, and event kernels sit at vertices. In the mode basis it is a sum over foams (Theorem 4.3), with:
   - face amplitude $e^{-\mu_k\tau_f}$, the survival of a mode over the face's duration;
   - vertex amplitude the matrix element of the event's kernel between modes, the discussion's $\langle\psi^+|T_\tau|\psi^-\rangle$ made exact;
   - selection rules wherever the dynamics is symmetric (Theorem 4.4).

   It does not depend on how time is sliced (Theorem 4.5). For a random walk on a group, the face weight is the heat kernel of two-dimensional Yang–Mills theory with duration in place of area.
6. **Foams carry evidence.** Without observations only the trivial foam survives and the total is 1. Every non-trivial foam reaches an observation, and its amplitude decays with its area, measured in time (Proposition 4.6). The amplitude of a single foam may be negative, or complex where a state is a phase (the hour of a street's daily crest); their sum is still a probability. This is interference without quantum mechanics (§4.8).
7. **Events are fundamental on both levels.**
   - In the record they are primary, and states are their realization.
   - In the theory they are the interaction vertices of a field theory whose Feynman graphs are event complexes: recorded vertices are fixed, unrecorded ones summed over.

   The engine implements the first and estimates the rates of the second. It never generates complexes to explain data, which removes the engineering difficulty the user noted (§1.6).
8. **What a 2-complex carries of three dimensions.** Volume lives in the couplings, not in 3-cells.
   - The areas and normals of interfaces form a classical spin network.
   - Closure, $\sum_f A_f n_f = 0$, holds for every closed volume.
   - By Minkowski's theorem, areas and normals determine each convex volume up to translation; a non-convex volume needs its face polygons (§5.4).

### 0.2 What This Note Changes in the Discussion

| The discussion | This note |
|---|---|
| $\mathcal K = (V, E, F, \partial_1, \partial_2, \prec, \lambda)$ with faces "relational primitives and higher-order compatibility structures" | A face is the life of a *binary* primitive between two events. Higher-order compatibility is held by a coupling (an edge) and its constraint, never by a face. Otherwise the cyclic order of a face's boundary would mean nothing (§1.2) |
| Weights $W(H\mid D) \propto \mathbf 1_{\text{valid}} \prod_f w_f \prod_e w_e \prod_v w_v$, stipulated | The factorisation is a theorem of a declared model (Theorem 4.2): faces carry persistence, edges the coupling's own persistence, vertices the events' kernels |
| $A_v = \langle\psi_v^+\mid T_{\tau_v}\mid\psi_v^-\rangle$ | Kept, exactly: the vertex amplitude in the mode basis, with $\psi^\pm$ the biorthogonal modes of the cells beginning and ending at $v$ (Theorem 4.3). These are not states of a Hilbert space |
| "Spin does not necessarily encode event category" | Agreed, and made constructive: categories are rules; spins are modes of states under the symmetry of their dynamics, found from the rates and not posited (§3) |
| Constructions A and B "are not automatically equivalent" | They are equivalent up to refinement, for spacetimes built of cylinders and rewrites (Theorem 2.8). The link of an event generalises the boundary spin network of a 4-simplex (Proposition 2.6) |
| A 2-complex may not carry enough to reconstruct 3D | It carries volume in its couplings. The class of configurations it determines is stated: with areas and normals, convex volumes up to position; with polygons, all (§5.4) |
| Fock space; $\hat L_{A\to B} = \lambda(\hat a^\dagger_B - \hat a^\dagger_A)\hat a_A$ | Two cases are kept apart. Identified objects (a building with its BIN) are hard-core, occupied 0 or 1. Anonymous populations (walkers) are Doi–Peliti. The event operator acts on the configuration and its complex, and keeps a demolished building's record (§4.10) |
| Noncommutative event algebra | Made precise: independent footprints commute. Events acting by group elements commute exactly when the elements do. Dependent stochastic events need not commute, and their commutator measures order-sensitivity (Proposition 3.5) |
| Complex amplitudes are to be avoided until a phenomenon needs them | Found: Fourier modes over cyclic states. The hour of a crest, the weekday, and a shift of civil time by daylight saving are vertices of definite phase. $Z$ stays a probability (§4.8) |
| "Should events be more fundamental?" (the user: an engineering difficulty) | Both, on two levels. The engine stays with the record, and the generative level is estimated, not simulated (§1.6, §6) |
| Title: "…Through Two-Complex Histories: A Spin-Foam-*Inspired* Formalism" | The user's title stands: *Event-Based 3D GIS on a 2-Complex: A Spin-Foam Approach*. The spin-foam form is derived (Theorem 4.3), not borrowed, so "approach" rather than "inspired" is earned |

# Part I. The Formalism

## 1. Ontology

### 1.1 Events, Couplings and Primitives

**Definition 1.1 (coupling).** A *coupling* $c$ is something that persists and holds primitives together. In a city these are:
- a building, a room, a plaza, the outside;
- the space of a street between two junctions;
- the space of a junction;
- in the network dictionary of §1.3, a junction.

A coupling has:
- a type;
- a life $[b_c, d_c]$, from the event of its birth to that of its death (either may lie outside the record);
- a geometry: its cell in space, and an anchor point;
- a state $x_c \in Q_c$: its unary properties, such as its use, its height class or its protection.

**Definition 1.2 (primitive).** A *primitive* $p$ is a relation or interface that persists between exactly two distinct couplings $a(p) \neq b(p)$, over a life $[b_p, d_p]$ within both of theirs. Examples:
- a party wall between two buildings;
- a frontage between a building and a street space;
- a cross-section between a street space and a junction space;
- a door between two rooms;
- in the network dictionary, a street segment between two junctions.

A primitive has:
- a type;
- a geometry: its cell in space, the polygon of a wall or the line of a segment;
- a state $x_p \in Q_p$: its permeability to walkers or vehicles, or whether a frontage is active.

**Definition 1.3 (event).** An *event* $v$ has:
- a time label $t_v$: an interval, with its time scale;
- a type $\tau_v$;
- a *match*: the couplings and primitives it is applied to;
- a *footprint* $(R_v, W_v)$: what it reads and what it writes;
- a *provenance*: recorded, inferred, modelled or synthetic.

Two derived relations:
- A coupling *attends* $v$ when it lies in the footprint, or when it is an end of a primitive in the footprint.
- A primitive is *met* at $v$ when it lies in the footprint.

*Rationale.* This is the architecture of the user's discussion papers (events, couplings, primitives), and the history domain's carrier with its terms generalised. Occasions are events, beings are couplings, and ties are primitives. A coupling *is* what holds its primitives together, as a being holds its ties, as a node of a spin network holds its links by an intertwiner.

### 1.2 Why Primitives Are Binary

**Proposition 1.4.** Suppose every primitive has two distinct couplings throughout its life, and both attend every event at which it is met. Let a face be the life of a primitive between two consecutive events at which it is met. Then:
- each face is attached along a closed path, up one coupling's worldline from the earlier event to the later and down the other's;
- the cyclic order of its boundary is time;
- the cells form a CW complex.

*Proof.* As for the history complex (its theory note, Proposition 2.2), whose construction is morphe's `worldline_complex`. Both couplings are present at both events, so both paths run between them, and the word $a[v,w]\cdot b[v,w]^{-1}$ returns to its base. ∎

The hypothesis cannot be weakened without losing the carrier.
- A primitive attached to one coupling would be attached along a path run forth and back, $w\,w^{-1}$. Its boundary is null in the chain complex, so each such face would be a 2-cycle on its own: a foam the data did not make.
- A primitive among three or more couplings at an instant does not sweep a surface whose time slices are arcs. No cyclic order of its couplings has a meaning.

This is the fault the history note found in the social spin foam's "facts on faces" (its idea note, §2). Here it fixes three rules:
- **A property of one thing is a state of its coupling.**
- **A relation among k ≥ 3 is held by a coupling of its own**, with a primitive to each party (Breiger 1974). For example, a rule in force over several streets is a coupling with its own life, enacted and repealed, joined to each street it governs by a primitive, its application there.
- **A loop is subdivided.** A street that leaves a junction and returns to it is split once more, as the landscape module already does (`networkOf`).

### 1.3 Two Dictionaries for a City

| Cell | Volumes and interfaces (the dual of space) | The walking network (the 1-skeleton of the district) |
|---|---|---|
| node of $\Gamma_t$, edge of $K$: **coupling** | a volume: a building, a room, a street space, a junction space, a plaza, the outside | a junction |
| link of $\Gamma_t$, face of $K$: **primitive** | an interface between two volumes: a party wall, a frontage, a cross-section of a street, a door, a stair | a street segment |
| plaquette of $\Gamma_t$, a 2-cell at an instant in $\hat K$ | the volumes round a hinge, for instance two buildings and the street at the corner where their party wall meets the facade line | a block |
| state of a node | use, height, protection | — |
| state of a link | permeability (to walkers, vehicles, heat, sight); the activity of a frontage | the regime of the segment |

The two dictionaries are related by Poincaré duality in space. In a decomposition of a $d$-dimensional region, the dual exchanges $d$-cells for nodes and $(d-1)$-cells for links.
- The walking network is the 1-skeleton of the landscape module's district (its §2.2). Its dual graph has the blocks as nodes.
- The volumes-and-interfaces dictionary is that of loop quantum gravity. There the spin network is dual to a cellular decomposition of space: its nodes are quanta of volume and its links quanta of area (Rovelli and Smolin 1995; Bianchi, Doná and Speziale 2011).

**This note adopts the second dictionary.** It holds buildings and streets in one network, and its natural state of a primitive, permeability, is what the discussion's pedestrian street changes. An Open Street is a street space whose cross-sections are closed to vehicles at its two ends, which is where the barriers stand. A door is a wall relabelled, as in the architecture domain (landscape note, §7.3).

### 1.4 Heterogeneous Events as Rules

Every category of event is a rule schema in the sense of the landscape module (its §3.1). Each is three things:
- a span $L_\tau \leftarrow K_\tau \to R_\tau$ on the typed incidence graph of the network;
- a kernel $T_\tau$ on the states of the cells it matches;
- a rate law $r_\tau$.

The categories:

| Category | Examples in the data | Rule on the network | Kernel on states | Footprint |
|---|---|---|---|---|
| **birth** | a new building (DOB's NB, a year of construction); a new business | attach a node with its links: new party walls and frontages | creation: the new cells' states drawn from a law | writes the new cells; reads the neighbours |
| **death** | a demolition (DOB's DM); a business closed | remove a node with its links | annihilation: the states summed out | writes the cells removed |
| **alteration** | a change of use or occupancy; an office building converted to housing | relabel a node | a map of classes, or a draw within a family | writes the node |
| **subdivision, merger** | a lot divided or merged; two rooms joined | partition or merge | copy $\Delta$; merge $\mu$ | writes the parts |
| **regime** | an Open Street on or off; a door locked; a café licensed on a frontage | relabel a link's state | a deterministic action $g\in G$ (a flip), or a kernel | writes the link |
| **institution** | a landmark designated; a rule in force; a permit granted | birth of a rule coupling with its links of application | creation | writes |
| **occurrence** | a collision; a complaint; an inspection | none | the identity | reads only |
| **observation** | a count at a sensor; a census; an annual snapshot | none | the identity, with a likelihood $\ell$ | reads only |
| **tick** | the clock | none | the identity | nothing |

The table rests on GIS precedent: the change types of identity-based change (Hornsby and Egenhofer 2000: create, destroy, continue, split, merge, and others) are rules here. What the rule adds to those types is its kernel on states and its rate.

### 1.5 Axioms of a Record

- **A1 (no change without an event).** Every change of a cell's existence, attachment or state happens at an event whose footprint writes it. A change the record does not hold is still an event, an *unrecorded* one (§4.1).
- **A2 (locality).** An event's footprint is the closure of its match: the couplings it touches with their primitives, and a primitive with its two couplings.
- **A3 (identity).** Couplings and primitives carry persistent identifiers. A rewrite says which cells persist (its $K$), which end and which begin.
  - Identity is a claim with a source: a BIN, a BBL, a property id, an OpenStreetMap way.
  - It is never a coincidence of geometry, which is the history domain's rule for beings across canons.
- **A4 (time as a label, order as structure).**
  - The causal order is the trace order of the footprints (landscape note, §3.3), refined by the time labels where these do not overlap.
  - Clock time is a label on the vertices, with its scale and uncertainty, not a coordinate of the complex.
- **A5 (interaction at events).** Between two events at which they are met, cells change independently of one another. Primitives interact only at events.

*Rationale.* A5 is the content of "event-based": events, not states, carry interaction. It makes the amplitude factorise exactly (Theorem 4.2).

A rate that reads several cells continuously, such as a store that opens more readily when its neighbour is open, is an interaction spread over time. Under A5 it is represented by unrecorded events joining the cells, which the theory sums (§6.3) and the engine does not.

### 1.6 Whether Events or Primitives Are Fundamental

The discussion asks two alternatives:
- should an event be generated whenever a primitive changes state?
- or should events be more fundamental, with the transformation of primitives emerging from their interaction?

The user notes that the second causes engineering difficulty. Both hold, on different levels, and keeping the levels apart removes the difficulty.

- **In the record, events are primary.** A state is the realization of the record up to a date (§5.1). Every change of state is an event (A1). An event need not change a state: occurrences and observations change none.
- **In the theory, events are interactions.** They are the interaction vertices of a field theory whose quanta are couplings carrying strands for their primitives (§6). A primitive changes only at a vertex, where it meets others (A5): its transformations emerge from interactions. Recorded events are fixed vertices; unrecorded ones are summed over.
- **The engine reads the record and never generates it.**
  - Unrecorded changes of a single cell are summed exactly, by the cell's persistence semigroup (§4.1).
  - Unrecorded interactions among several cells are the perturbative part of the theory. Their couplings are estimated from the rates of the record's event types (§6.2), not summed.
  - The engine never searches for complexes that would explain data.

  This is zeroth order in the latent interactions, and it is exact for what it computes.

## 2. The Complexes

### 2.1 The Record

**Definition 2.1 (record).** A *record* over a window $[T_0, T_1]$ consists of:
- a set $C$ of couplings and a set $P$ of primitives, with their types, their ends $a, b: P \to C$, their identifiers and their geometry;
- a set $E$ of events, with their time labels, types, matches, footprints and provenance;
- a total order $\le$ on $E$: by time label, then by the order of the source;
- the causal order $\prec$, the trace order of the footprints;
- two boundary events $\partial_0$ and $\partial_1$, at $T_0$ and $T_1$, attended by every cell alive then.

The lives of cells are read from the events that write their births and deaths. A cell born before $T_0$ is born at $\partial_0$. The record is that of the history domain (its §1) with typed beings and ties and with footprints, and that of the landscape module (its §10.1) on a network instead of a fixed district.

### 2.2 The Refined Complex and the Event Complex

**Definition 2.2 (refined event complex).** $\hat K$ is morphe's refined worldline complex (`refined_worldlines`) of the encounters in which:
- the participants are the couplings;
- the couplings present at an event are those that attend it;
- the pairs are the primitives;
- a primitive's meetings are the events at which it is met, its birth and death among them.

Its cells:
- **vertices:** the presences $(c, v)$, a coupling at an event: points of spacetime;
- **temporal edges:** the segments of a coupling between two consecutive events it attends;
- **rungs:** $(p, v)$, a primitive at an event at which it is met, from $(a(p), v)$ to $(b(p), v)$;
- **2-cells:**
  - the lens squares $a[v,w]\; r_{p,w}\; b[v,w]^{-1}\; r_{p,v}^{-1}$ for consecutive meetings $v, w$ of $p$;
  - the plaquettes at events: cycles of rungs declared as 2-cells, either morphe's triangles or the hinge cycles of §1.3.

**Definition 2.3 (event complex).** $K$ is morphe's worldline complex (`worldline_complex`), which is $\hat K$ with the rungs of each event collapsed. Its cells:
- **vertices:** the events, $\partial_0$ and $\partial_1$ among them;
- **edges:** the segments of couplings;
- **faces:** the lenses $a[v,w]\,b[v,w]^{-1}$.

The links of $\partial_0$ and $\partial_1$ are the spin networks $\Gamma_{T_0}$ and $\Gamma_{T_1}$. So $K$ is the foam with its two boundaries coned off, and $\hat K$ keeps them as they are.

**Proposition 2.4.** $\hat K$ and $K$ are CW complexes, and
$$\chi(K) = |E| + 2 - \sum_{c} n_c + \sum_{p} m_p,$$
where $n_c$ is the number of segments of $c$ and $m_p$ the number of lenses of $p$.

The collapse $\hat K \to K$ is a cellular quotient. It is a homotopy equivalence when, at each event, the graph of its rungs together with its plaquettes is contractible. Otherwise $H_1(K)$ loses the cycles of those graphs.

*Proof.* The first part is morphe's construction; the count is that of the cells. Collapsing a contractible subcomplex is a homotopy equivalence, and collapsing a graph kills its cycles. ∎

*Remark.* This is the third use of one construction. The history complex and the economics domain's spacetime of transactions are the other two, all on morphe's `worldlines.hpp`. The city needs no new complex.

### 2.3 The Cut and the Spin Network of the Moment

**Definition 2.5 (spin network of the moment).** For a time $t$ in the window that is not an event time, $\Gamma_t$ is the section of $\hat K$ at $t$:
- its nodes are the couplings alive at $t$, one for each temporal edge crossing $t$;
- its links are the primitives alive at $t$, one for each lens crossing $t$, joining their ends;
- its plaquettes are the declared 2-cells among them.

It is labelled by the states, or the modes (§4.6), of the cells it crosses. The geographic state at $t$ is its realization, $S_t = \mathcal R(\Gamma_t, \lambda)$ (§5.1), where $\lambda$ is the cells' geometry.

### 2.4 Time-Change Vertices and Refinement

The user's vertices "include time change". A tick at $t$ is attended by the couplings whose worldlines it cuts, and it meets no primitive. It subdivides their temporal edges, and with them the sides of the lenses that run along them: a refinement, a cellular subdivision that leaves the homotopy type unchanged.
- **Global ticks** (a clock) give the *fully sliced* complex: every worldline cut at every tick.
- **Observations** are local ticks: they cut only the cells they read.

**Why a tick meets no primitive.** A tick that met the primitives too would cut their lenses (still a refinement), but it would also bring the plaquettes among them at $t$. Those close a new shell at every hinge, which is not a refinement. Plaquettes record constraints at an instant (§4.7), and a clock must not invent them. The module's checks hold the first case, ticks that meet no primitive.

Theorem 4.5 shows that ticks change no amplitude. The slicing of time is therefore a choice of description, as a gauge is, and not a modelling choice.

### 2.5 The Link of an Event

**Proposition 2.6.** In $K$, the link of an interior event $v$ is the following graph:
- **vertices:** the half-edges at $v$, the incoming and outgoing segments of each attending coupling (only one of the two if the coupling is born or dies at $v$);
- **edges:** the corners of the lenses at $v$, of three kinds:
  - *closing* corners join two incoming half-edges, where a lens of a met primitive ends;
  - *opening* corners join two outgoing half-edges, where its next lens begins;
  - *passing* corners join the incoming and outgoing half-edges of one coupling, where a lens of an unmet primitive of an attending coupling runs past $v$.

Its incoming part (incoming half-edges and closing corners) is the network of the footprint before the event, $\Gamma_L$. Its outgoing part is the network after, $\Gamma_R$. The passing corners join them along what persists. This is the history domain's bubble (its §8.1), which its checks verify on every canon.

*Remark (Pachner moves).* Let a primitive change one of its couplings at an event: a *transfer* corner, from an incoming half-edge to an outgoing half-edge of a different coupling. Then the link of a 1–4 Pachner move of a 3-dimensional triangulation is the complete graph on five nodes. That is the boundary spin network of the 4-simplex dual to the vertex, $L \cup_\partial R$. So the link of an event generalises the boundary of a 4-simplex.

morphe's worldline complex keeps each pair's ends fixed, so it would need this generalisation. Without it, a transfer is an ending and a beginning joined by the event's kernel, which carries the same information (§3.4).

### 2.6 The Geometric Dual: Construction A Is Construction B

**Lemma 2.7.** Let $X$ be a regular cell decomposition of a $d$-manifold. For a region, take the dual relative to its boundary, with the outside as one more cell. Let $I = [T_0, T_1]$ be subdivided at instants. Then:
- the dual decomposition of $X \times I$ is $X^* \times I^*$, the product of the duals;
- $I^*$ is again a path, with a vertex for each period and an edge for each interior instant.

*Proof.* The dual cell of $\sigma \times \tau$ in a product is $\sigma^* \times \tau^*$ for the barycentric duals, and the dual of a subdivided interval is the interval subdivided at the midpoints of its periods (Munkres 1984, §64–65). ∎

**Theorem 2.8.** Suppose a city's space is decomposed by regular complexes $X_k$ on the periods $(t_k, t_{k+1})$ between events, and each event is a rewrite $X_{k-1} \Rightarrow X_k$ by a span $L \leftarrow K \to R$ whose $K$ persists. Let $M$ be the cellular spacetime whose cells are the products $\sigma\times[t_k, t_{k+1}]$ and $\sigma\times\{t_k\}$, glued along the persisting cells. Then:

1. **The spatial dual.** Dualise $M$ in space only: $d$-cells become nodes and $(d-1)$-cells links, with time left as it is. Its 2-skeleton is the fully sliced refined complex $\hat K_{\mathrm{glob}}$ of the record whose couplings are the $d$-cells and whose primitives are the $(d-1)$-cells:
   - every event attended by every coupling and meeting every primitive;
   - rungs at every instant, lens squares on every period, plaquettes at every instant.
2. **Local events, plaquettes aside.** Without its plaquettes, the refined complex $\hat K$ of the same record with local events is obtained from $\hat K_{\mathrm{glob}}$ without its plaquettes by undoing the subdivisions at events that do not touch a cell. These are refinements, so the two have the same homotopy type.

   The plaquettes are where the two differ:
   - the spatial dual has one at every hinge at every instant, closing a shell for each hinge and period;
   - $\hat K$ keeps one only at an event that meets all its interfaces, where it records a constraint at that instant (§4.7).
3. **The full dual.** The full dual of $M$ exchanges instants and periods along each worldline (Lemma 2.7). Up to subdivision it is the same decomposition of $|X^*| \times [T_0, T_1]$.
4. **Events as cells.** Taking each event's cobordism as a single $(d+1)$-cell of $M$ corresponds to the collapse $\hat K \to K$. The link of its dual vertex is the bubble of Proposition 2.6.

*Proof.*
- (1) By Lemma 2.7 on each period, with the time factor kept primal. Persisting cells glue identically across events.
- (2) Take an event that touches neither a coupling $c$ nor any primitive of $c$. The vertex $(c, t_k)$ has two temporal edges, and, once the plaquettes are set aside, its rungs at $t_k$ bound only lens squares of untouched primitives. Removing it, with those rungs, merges two temporal edges and the pairs of squares on either side: an elementary unsubdivision.
- (3) This is the instant–period exchange of Lemma 2.7.
- (4) This is Proposition 2.6. ∎

So the event-history complex *is* the dual of the city's cellular spacetime: the spatial dual, up to refinement, with the plaquettes kept where an event constrains them. The history domain's refined complex is an instance; there the space is the network of beings and ties at each occasion. A 3D GIS with versions and its event record are one structure, seen with time kept primal.

### 2.7 Readings of the Topology

The readings carry over from the history note (§4, §7):
- **$H_1(K)$:** circuits of couplings met pairwise at events, with no interface recorded across some pair. It marks a place to look for an unrecorded interface or event: a silence of the record.
- **$H_2(K)$:** closed shells of primitives, for example three couplings pairwise joined by primitives met at the same two events (two buildings and their street whose interfaces two events both changed: a corner rebuilt twice).
- **Persistence in time:** each circuit is born when an event closes it and dies when a later interface fills it (morphe's `persistent_homology`).
- **For a $\mathbb Z_2$ gauge model:** $H_2(\hat K;\mathbb Z/2)$ is the space of foams (§4.7).

## 3. Labels: Where Spins Come From

### 3.1 A Category Is a Rule, Not a Spin

An event's category determines three things (§1.4):
- which cells it can match;
- what it does to them (its rewrite and its kernel);
- at what rate it occurs.

Encoding categories as spins would give a name the meaning of a representation: a dimension $2j+1$ and a law of recoupling. Nothing in a construction or a closure supplies that meaning. The discussion is right to refuse it.

What representation theory does supply is the decomposition of *states* by the symmetry of their *changes*. That is where spins belong.

### 3.2 States, Symmetry and Modes

**Definition 3.1.** Each type of cell has:
- a finite state space $Q$;
- a generator $L$ of its unrecorded changes: a rate matrix with rows summing to zero;
- a symmetry group $G$: a group acting on $Q$ that commutes with $L$.

The largest such group is *found* from the rates (the automorphism group of the weighted graph of transitions); it is not posited.

The **modes** are:
- the isotypic components of $\mathbb R^Q$ (or $\mathbb C^Q$) under $G$, refined into eigenspaces of $L$;
- with no symmetry, the eigenspaces of $L$ alone.

Theorem 4.3 needs only that $L$ be diagonalisable. Reversibility makes the modes orthogonal in $L^2(\pi)$.

| State | $Q$ | Symmetry of its dynamics | Modes (spins) | What a spin says |
|---|---|---|---|---|
| existence of an identified cell | $\{0, 1\}$, hard-core | none | 0 and 1 | — |
| permeability, opening and closing alike | $\mathbb Z_2$ | $\mathbb Z_2$ | $j \in \{0, 1\}$ | whether the regime is remembered |
| a use class, with no structure | $q$ classes | $S_q$ | trivial; standard (dimension $q-1$) | whether the class is remembered |
| a use class in a classification of families | $n$ families of $m$ classes | $S_m \wr S_n$ | trivial; between families ($n-1$); within families ($n(m-1)$) | the level of the classification still remembered |
| a nucleotide (the genomics domain) | a torsor of the Klein group | $\mathbb Z_2^2$ | four characters | Kimura's three-parameter model, a precedent |
| the hour of a daily crest | $\mathbb Z_{24}$, or $U(1)$ | rotations | harmonics $n$ | which harmonic of the day; complex |
| the day of the week | $\mathbb Z_7$ | rotations, where the dynamics has them | seven characters | the weekly rhythm |
| a volume's shape (§5.4) | areas and normals | $SO(3)$ | spins on links, intertwiners on nodes | areas, and the polyhedron (classical) |

**Proposition 3.2 (the memory of a classification).** Let a cell's class lie among $n$ families of $m$ classes each. It jumps:
- within its family at total rate $\alpha$, to a uniformly chosen other class of that family;
- across families at total rate $\beta$, to a uniformly chosen class of another family.

Then the symmetry is $S_m \wr S_n$, $\mathbb R^Q$ splits into the trivial, between-family and within-family components, and their rates of decay are
$$\mu_{\text{between}} = \frac{n\beta}{n-1}, \qquad \mu_{\text{within}} = \frac{m\alpha}{m-1} + \beta.$$
When changes within a family are faster than changes across ($\alpha \gg \beta$), the family is remembered long after the class is forgotten. For $m = 1$ this is the Jukes–Cantor rate $q\beta/(q-1)$.

*Proof.*
- Let $f$ be constant on families and sum to zero over them. A jump within a family leaves $f$ unchanged. A jump across families lands in another family, chosen uniformly, whose value averages to $-f(x)/(n-1)$. So $Lf = \beta(-f/(n-1) - f) = -\frac{n\beta}{n-1}f$.
- Let $f$ sum to zero within each family. A jump within the family averages to $-f(x)/(m-1)$. A jump across families lands uniformly in a family where $f$ averages to zero. So $Lf = -\big(\frac{m\alpha}{m-1} + \beta\big)f$. ∎

So the spin of a land-use class is the level of the classification at which it is still carried: an answer of the classification's own kind. It is the genomics domain's reading of the bases (IUPAC's codes as cosets of the Klein group, the classifications as its characters) transposed to land use. Hierarchical classifications such as ANZSIC, with its divisions, groups and classes, or New York's land-use codes in their families, have it by construction.

### 3.3 Couplings' States and Intertwiners

A coupling's state is the label of its edge. Two situations arise.

- **Relational states (gauge models, §4.7).** Sometimes the states of a coupling's primitives are defined only relative to frames: a price ratio relative to a numéraire, a lag relative to a clock, an angle relative to a building's axes. Changing frames at a coupling is then a redundancy, the theory is a gauge theory, and the coupling's labels are intertwiners: invariant tensors of its primitives' representations. The Gauss law holds at edges.
- **Absolute states.** Where states are absolute, as an Open Street's closure is, there is no redundancy and no Gauss law at edges. The selection rules move to the vertices, where symmetric events conserve their modes' charges (Proposition 3.3).

In loop quantum gravity the two meet: the intertwiner is the state of a volume (its shape), and its gauge invariance is the closure of the polyhedron (§5.4).

### 3.4 What Events Do to Modes

| Kernel | Its action on modes |
|---|---|
| tick (identity) | each mode passes unchanged, $\delta_{kk'}$ |
| observation, a likelihood $\ell(x)$ | couples modes through $\ell$'s coefficients: an observation injects every mode it can distinguish |
| annihilation (a death, a cell forgotten): contraction with $\mathbf 1$ | only the trivial mode passes |
| creation (a birth) with initial law $p_0$ | the new cell begins in the modes of $p_0$; a stationary $p_0$ gives the trivial mode alone |
| copy $\Delta(x) = (x, x)$ (a split) | fusion. For abelian $G$, $j = j_1 + j_2$. For $S_q$, the standard representation lies in its own square, so the trivial mode coming in can give (standard, standard) going out: the two parts of a split lot are correlated even when their parent's use is unknown |
| merge $\mu$ | the dual of fusion |
| a deterministic action of $g \in G$ (a scheduled flip; a clock moved by $\Delta$ hours) | a phase $\chi_j(g)/d_j$: for a flip $(-1)^j$, for a shift $e^{2\pi i n\Delta/24}$ |
| a relabelling by a map $\varphi: Q\to Q$ (a reclassification) | the representation of $\varphi$: an element if $\varphi \in G$; otherwise it mixes the modes |
| transfer (a primitive re-attached, or carried to its successor) | the identity on its state: every mode passes to the successor |

**Proposition 3.3 (selection rule).** Suppose an event's kernel $T$ commutes with the action of $G$ on its incoming and outgoing cells. Then its matrix element between isotypic components of types $\rho_{\text{in}}$ and $\rho_{\text{out}}$ vanishes unless $\mathrm{Hom}_G(\rho_{\text{in}}, \rho_{\text{out}}) \neq 0$. For abelian $G$ this is the conservation of charge, $\sum j_{\text{in}} = \sum j_{\text{out}}$. The non-zero elements factor as an intertwiner times a reduced amplitude.

*Proof.* Restricted to isotypic components, $T$ is a $G$-map, so Schur's lemma applies. The factorisation is the Wigner–Eckart theorem (Serre 1977, §2). ∎

### 3.5 The Algebra of Events

Each event, with its match, acts as a stochastic operator on $\bigotimes_\sigma \mathbb R^{Q_\sigma}$: its kernel on the cells of its footprint, the identity elsewhere.

**Proposition 3.5.**
1. Events with independent footprints, in the sense of traces, commute.
2. Two deterministic actions by elements $g, h$ of a group acting faithfully commute exactly when $gh = hg$. So the noncommutativity of their events is that of the group: two reclassifications by permutations of $S_q$ need not commute.
3. Dependent stochastic events need not commute. The norm of their commutator, in total variation on the joint law, measures the sensitivity to order that the record's linearisations would show.

The converse of (1) fails: two flips of one street are dependent and commute.

*Proof.*
- (1) Operators on disjoint tensor factors commute.
- (2) A faithful action is a faithful representation.
- (3) A kernel that reads what the other writes is an example. ∎

*Remark.* The landscape module's anchored interactions (its §6.4) measure how interventions interact. The commutator measures how much order matters. They are different questions, and both are exact on a declared model.

## 4. Amplitudes

### 4.1 The Event-Coupled Markov Model

**Definition 4.1.** On the event complex $K$ of a record, give each cell $\sigma$ (a segment of a coupling, or a lens of a primitive) the state space $Q_\sigma$ of its type. Then:

1. **Persistence.** Along $\sigma$, of duration $\tau_\sigma$, the state follows a Markov chain with generator $L_\sigma - H_\sigma$, independently of every other cell (A5).
   - $L_\sigma$ generates the cell's unrecorded changes.
   - $H_\sigma$ is the diagonal of the hazards of the recorded event types whose rate reads $\sigma$ alone. Their survival, the chance that no such event happened along the cell, lies in the cell.
   - Write $P_\sigma(\tau) = e^{\tau(L_\sigma - H_\sigma)}$.
2. **Events.** At each event $v$, the states of the cells ending at $v$ are carried to the states of the cells beginning there by $V_v$:
   - the kernel $T_v$ of its type and match;
   - times its rate density, for types whose times are modelled;
   - times its likelihood $\ell_v$, at an observation.

   The times of other recorded events are conditioned on.
3. **Boundary.** The cells alive at $T_0$ start from an initial law $\pi_0$. At $T_1$ the states are summed.

**What a face's amplitude is.** A face carries the sum over the *unrecorded* histories of its primitive between two recorded events: changes that no source records, but which snapshots and observations betray. A fully recorded primitive has $L_\sigma = 0$, and its face carries only survival.

### 4.2 Factorisation

**Theorem 4.2.** Under Definition 4.1, consider the probability (density) of three things together:
- the record's events of modelled types;
- the kernels taken at its events;
- its observations,

conditional on the times of the other events. It equals the contraction $Z$ of the tensor network on $K$ with:
- the matrix $P_\sigma(\tau_\sigma)$ on each cell, from its start state to its end state;
- the tensor $V_v$ on each vertex;
- $\pi_0$ on the cells alive at $T_0$, and $\mathbf 1$ on those alive at $T_1$.

*Proof.*
- **Product.** By A5 and the Markov property, the joint law of every cell's states at its two ends, together with the event data, is the product of these factors in the order of time.
- **Marginalising.** Summing out the states contracts the indices. Each cell's start index meets the tensor of the vertex where it begins, and its end index that of the vertex where it ends.
- **Pattern.** So the contraction pattern is the incidence of $K$. ∎

*Remark.* This is a hidden Markov model whose hidden chain lives on the cells of a 2-complex rather than on a line: a factor graph (Kschischang, Frey and Loeliger 2001) whose factors are the cells and vertices of $K$. The complex adds no expressive power to factor graphs; what it adds is said in §7.

### 4.3 The Spin-Foam Form

**Theorem 4.3.** Suppose each $L_\sigma - H_\sigma$ is diagonalisable, with:
- eigenvalues $-\mu_{\sigma,k}$;
- right eigenvectors $r_{\sigma,k}$ and left eigenvectors $l_{\sigma,k}$, normalised so that $l_{\sigma,k}^{\mathsf T} r_{\sigma,m} = \delta_{km}$.

Then
$$Z = \sum_{k} \;\prod_{\sigma} e^{-\mu_{\sigma,k_\sigma}\tau_\sigma}\;\prod_{v} A_v(k),$$
where the sum runs over the labellings $k$ of the cells by modes, and the vertex amplitude is
$$A_v(k) = \sum_{y,\,x} \Big[\bigotimes_{\sigma\,\in\,\mathrm{in}(v)} l_{\sigma,k_\sigma}\Big](y)\; V_v(y, x)\; \Big[\bigotimes_{\sigma'\,\in\,\mathrm{out}(v)} r_{\sigma',k_{\sigma'}}\Big](x).$$
Here $y$ runs over the end states of the cells ending at $v$ and $x$ over the start states of those beginning there. At the boundary:
- $A_{\partial_0}(k) = \sum_x \pi_0(x) \prod r_{\sigma,k_\sigma}(x_\sigma)$;
- $A_{\partial_1}(k) = \prod_\sigma l_{\sigma,k_\sigma}^{\mathsf T}\mathbf 1$.

*Proof.* Insert $P_\sigma(\tau) = \sum_k e^{-\mu_{\sigma,k}\tau}\, r_{\sigma,k}\, l_{\sigma,k}^{\mathsf T}$ into the contraction of Theorem 4.2, and exchange the finite sums. ∎

The theorem is trivial to prove; its content is its reading.

| Cell | Label | Amplitude | Meaning |
|---|---|---|---|
| face (a primitive's lens) | a mode $k_f$ of its state | $e^{-\mu_{k}\tau_f}$ | how much of that mode of the primitive's state survives the interval between two events |
| edge (a coupling's segment) | a mode of the coupling's state | $e^{-\mu_{k}\tau_e}$ | the same, for the coupling (its use, its protection) |
| vertex (an event) | — | $A_v$, the matrix element of the event between incoming and outgoing modes | what the event does to what is remembered: the discussion's $\langle\psi^+\mid T_\tau\mid\psi^-\rangle$ |
| a labelling of every cell | a foam | the product | one channel by which the record's states inform one another |

**Theorem 4.4 (the symmetric case).** Suppose a group $G$ acts on the states of every cell, every $L_\sigma - H_\sigma$ commutes with it, and every $V_v$ commutes with the diagonal action on its cells. Then:
1. **Labels.** The modes can be taken in isotypic components, labelled by irreducible representations $j_\sigma$ with indices of multiplicity.
2. **Selection rules.** $A_v$ vanishes unless the incoming and outgoing representations admit a $G$-map (Proposition 3.3).
3. **The heat kernel.** For a random walk on $G$ (states a $G$-torsor; jumps by a conjugacy class $C$, closed under inverses, at rate $\lambda$),
$$\mu_j = \lambda\Big(1 - \frac{\chi_j(C)}{d_j}\Big), \qquad P_\sigma(\tau)(g, h) = K_\tau(g^{-1}h), \quad K_\tau(g) = \frac{1}{|G|}\sum_j d_j\, e^{-\mu_j\tau}\,\chi_j(g).$$
   $K_\tau$ is the heat kernel on $G$.

*Proof.*
- (1) Commuting operators share invariant subspaces.
- (2) This is Proposition 3.3.
- (3) The average of $\rho_j$ over a conjugacy class is $\chi_j(C)/d_j$ times the identity, by Schur's lemma (Diaconis 1988, ch. 3). At $\tau = 0$, $K_0(g) = \delta_{g,e}$ by the orthogonality of the columns of the character table. ∎

### 4.4 Refinement Invariance and Two-Dimensional Yang–Mills

**Theorem 4.5.** Subdividing a cell at an instant inside it, by a tick with the identity as its kernel, does not change $Z$. In modes:
- the tick's amplitude is $\delta_{kk'}$, so the mode is carried across;
- the face amplitudes multiply: $e^{-\mu\tau_1}e^{-\mu\tau_2} = e^{-\mu(\tau_1+\tau_2)}$.

So $Z$ is the same on $K$, on $\hat K$, and on every slicing of either by clocks or by observations that observe nothing.

For a random walk on $G$, this is the invariance under subdivision of two-dimensional Yang–Mills theory with the heat-kernel action (Migdal 1975; Witten 1991). The duration of a face plays the part of its area.

*Proof.* The semigroup property $P(\tau_1)P(\tau_2) = P(\tau_1 + \tau_2)$ (Chapman–Kolmogorov), and the biorthogonality of the modes. ∎

*Remark (a revision to propose to the history domain).* The history domain gives every lens of a colour one coupling $\beta$, whatever its duration. Theorem 4.5 says that, for a theory independent of how often a tie happens to be attested, a lens's coupling should depend on its duration: $\tanh\beta_f = e^{-2\lambda\tau_f}$ for flips at rate $\lambda$. The model with constant $\beta$ is invariant under refinement only for lenses of equal length. This note proposes the change; it does not make it.

### 4.5 Foams Carry Evidence

**Proposition 4.6.** Let every cell be conservative ($H_\sigma = 0$, with the event times conditioned on) and every kernel stochastic, so that $r_{\sigma,0} = \mathbf 1$. Then:

1. **No observations.** Without observations, $Z = 1$, and the trivial foam (every cell in its trivial mode) is the only labelling with a non-zero product.
2. **Where foams end.** With observations, a labelling has a non-zero product only if every cell in a non-trivial mode reaches, through cells in non-trivial modes, an observation in its future. Non-trivial foams are surfaces that begin where states are set (the initial law, births with a known state, events that set a state) and end on the observations that read them.
3. **The area law of evidence.** A foam's amplitude is the product of $e^{-\mu\tau}$ over its non-trivial cells. It decays exponentially with the area of its non-trivial part, area measured in duration. What an observation says about a primitive's state a time $\tau$ away decays as $e^{-\mu\tau}$ along every face it must cross.

*Proof.* Argue backward in time.
- At $\partial_1$, $l_{\sigma,k}^{\mathsf T}\mathbf 1 = l_{\sigma,k}^{\mathsf T} r_{\sigma,0} = \delta_{k0}$.
- Take a vertex without observation whose outgoing cells all carry the trivial mode. There $\sum_x V_v(y,x)\prod r_0(x) = 1$, because $V_v$ is stochastic. So $A_v = \prod_{\sigma\in\mathrm{in}(v)} l_{\sigma,k_\sigma}^{\mathsf T}\mathbf 1 = \prod \delta_{k_\sigma 0}$.
- Induction over the vertices in reverse order gives (1); the amplitude of the trivial foam is then $\sum\pi_0 = 1$.
- The same induction stops only at observations, which gives (2).
- (3) is Theorem 4.3. ∎

*Remark.* In the history domain's Gibbs model, foams are closed 2-cycles, and Wilson loops are what give them a boundary. Here time is directed and the model stochastic, so foams are surfaces with their boundary on the observations: relative cycles.

### 4.6 The Boundary State at a Cut

**Proposition 4.7.** At a cut $t$,
$$Z = \sum_k \hat F_t(k)\,\hat B_t(k), \qquad \hat F_t(k) = \langle F_t, r_k\rangle,\quad \hat B_t(k) = \langle l_k, B_t\rangle,$$
where:
- $k$ labels the modes of the cells crossing $t$;
- $F_t$ is the forward message: the joint law of the crossing cells' states given everything before $t$;
- $B_t$ is the backward message: the likelihood of everything after $t$ given those states.

The weights $w_k = \hat F_t(k)\hat B_t(k)/Z$ say how the foams cross the cut. They are the record's *boundary spin networks* at the date, signed in general. With no observation after $t$, $B_t = \mathbf 1$ and only the trivial network crosses.

*Proof.* $Z = \langle F_t, B_t\rangle$, by splitting the contraction at the cut, then the biorthogonal expansion. ∎

This is the history domain's Proposition 6.3 for a directed model. The spin network of the moment, labelled by its $w_k$, is the GIS's state at a date as one object together with its uncertainty.

### 4.7 Compatibility Weights and Relational States: Gauge Models

The history domain's model is a model of the other kind:
- its variables lie on the edges of $\hat K$: the relational values of ties at meetings, and the frames of beings;
- its weights are functions of holonomies on lenses and plaquettes;
- it is undirected, its $Z$ normalises, its foams are 2-cycles, and the Gauss law holds at edges;
- physis's contraction computes it exactly at $\mathbb Z_2$.

The economics domain's transports are another instance, with the group $\mathbb R_{>0}$.

A city has relational quantities of the same kind:
- the price per square metre of one lot against another, from recorded sales;
- the lag between the tides at two counters;
- the alignment of a facade with its street.

On these the city has gauge models on the same complex, read through the same character expansion (history note, §6.2). Where the two families meet, and where they part:
- **They agree on a chain of lenses.** On a chain with no plaquettes, the symmetric $\mathbb Z_2$ persistence of §4.1 is the Gibbs lens weight $e^{\beta s s'}/2\cosh\beta$ with $\tanh\beta = e^{-2\lambda\tau}$.
- **They part where plaquettes constrain an instant.** Frustration round a triangle is a constraint at one moment, and a directed model under A5 cannot impose it without an event.
- **The division of labour:**
  - Gibbs (compatibility) models are for constraints at an instant;
  - stochastic (transition) models are for change through time;
  - the city's first model is stochastic;
  - its gauge models of relational states come later (§8).

### 4.8 Signed and Complex Amplitudes

$Z \ge 0$ always: it is a probability. A foam's amplitude is a product of coefficients of modes.
- **Signed.** Amplitudes are signed when a kernel or a likelihood has negative coefficients: a flip's $(-1)^j$, or a likelihood that favours a state against the stationary law.
- **Complex.** Amplitudes are complex for the modes of $\mathbb Z_N$ ($N > 2$) and $U(1)$, such as the hour of a street's daily crest or the weekday.

Here is the phenomenon the discussion asked for, in which phases compose with an interpretable role:
- **The phase.** Daylight saving moves civil time by an hour. On every cell whose state is a civil hour, that is a vertex of phase $e^{2\pi i n/24}$; on cells kept in solar time it is no phase at all. Street Oceanography's question of whose clock a street keeps (the `rhythm` module) is the question of which cells carry the phase.
- **Composition.** Shifts of a cycle compose by adding their phases, and foams of opposite phase cancel in $Z$.
- **A topological observable.** The phase of the daily tide at the counters is a $U(1)$ field. Its winding round a block, summed from the wrapped differences along the streets, is an integer; where it is not zero, it is an *amphidrome* of the pedestrian tide, as the cotidal charts of the ocean have.

This is interference of Fourier coefficients, not of quantum amplitudes. A quantum reading would need boundary Hilbert spaces with an operational meaning and a unitary composition, and a city supplies neither. The note keeps such a reading out of scope, as the discussion advises.

### 4.9 The Meaning of an Amplitude

| The discussion's reading | Here | Section |
|---|---|---|
| compatibility (a non-negative weight) | the Gibbs weight of a gauge model; the support of $A_v$ as admissibility | §4.7, §3.3 |
| transition (a stochastic kernel) | the event-coupled Markov model: semigroups on faces and edges, kernels at vertices; $Z$ the probability of the record | §4.1–4.6 |
| a quantum-style history (complex) | the complex coefficients of modes over $\mathbb Z_N$ and $U(1)$, with $Z$ still a probability | §4.8 |

So the answer to "what does an amplitude mean" has two parts:
- **The whole.** $Z$, the amplitude of a whole history, is the probability of its record under a declared model.
- **One foam.** The amplitude of a single foam is a Fourier coefficient of that probability: one channel by which what was set informs what was seen, and it decays with the channel's area in time.

### 4.10 Generation of Events: the Fock Space and the Event Operator

**Two kinds of population.**
- **Anonymous populations.** For walkers at junctions, or the establishments of a class in a block, counts and their Doi–Peliti operators are exact (physis `chemistry/doi_peliti.hpp`; landscape note, §8.3).
- **Identified objects.** A building with its BIN is occupied once or not at all. So the Fock space of identified cells is $\bigotimes_\alpha \mathbb C^2$, the hard-core case the discussion anticipates. Creation $\hat a^\dagger_\alpha$ is the birth of $\alpha$; annihilation $\hat a_\alpha$ is its death.

**The process with its history.** The state of a city together with its history is a pair $|C, K\rangle$:
- $C$, the configuration: the cells alive, with their states;
- $K$, the complex so far.

**Definition 4.8.** The *event operator* of a type $\tau$ is
$$\hat{\mathcal E}_\tau = \sum_{m} r_\tau(m)\,\big(\hat T_{\tau,m} \otimes \mathrm{attach}_{\tau,m} - \mathbf 1\big),$$
summed over the matches $m$ of $L_\tau$ in $C$, where:
- $\hat T_{\tau,m}$ acts on the states of the matched cells;
- $\mathrm{attach}_{\tau,m}$ adds the event's vertex to $K$ with its link: the segments of the matched couplings ending and beginning there, and the lenses of its met primitives closing and opening.

The generator of the city with its history is $\hat L = \sum_\tau \hat{\mathcal E}_\tau$. It conserves probability, $\langle\mathbf 1|\hat L = 0$, because each term is a stochastic map less the identity.
- **Forgetting the history.** Summing $K$ away, $\hat L$ reduces to the generator of the configuration alone.
- **Keeping it.** With $K$ kept, the process records its own history. An annihilation removes a cell from $C$ but adds to $K$ the vertex of its death, so the record of a demolished building stays, as the discussion asks.

**Proposition 4.9.**
1. `attach` keeps $K$ a CW complex at every step.
2. The attachments of independent events commute.

*Proof.*
- (1) Every lens closed at the new vertex is closed along both its couplings' worldlines, because both attend the event (A2). This is Proposition 1.4, inductively.
- (2) Independent events touch disjoint cells, and their rewrites commute by the local Church–Rosser theorem for adhesive categories (landscape note, §3.1). ∎

This is the discussion's first proposition: a local event operator preserves the validity of the complex.

## 5. Realization and Reconstruction

### 5.1 The Realization Map

**Definition 5.1.** The *realization* $\mathcal R(t)$ gives the geographic state $S_t$ at a date. For each cell crossing $t$ it takes:
- its geometry, from its labels (a building's footprint, base and top; a street's line and width);
- its state.

Its scene is the 3D city at $t$.

### 5.2 Consistency

**Proposition 5.2.**
1. **Nothing lost at a date.** Under A1 and A3, $\mathcal R \circ \mathrm{cut}_t \circ H$ returns the states $S_t$ from which the record was built.
2. **Order does not matter.** $S_t$ does not depend on which linearisation of the trace is realised. Independent events commute (Proposition 3.5(1) and the local Church–Rosser theorem), so every linearisation realises the same state at every date after them.

*Proof.* By induction over the events, in any linearisation. ∎

### 5.3 What Slices Cannot Recover

Slices without identities do not determine the complex.
- **The lane.** A lane built and then removed returns a configuration isomorphic to the first (landscape note, §3.4). Its complex is not the trivial one.
- **Exchanged uses.** If only the multiset of uses is kept, two buildings that exchange their uses look the same as two that change nothing.

Identities (A3) are what CityGML's versioning keeps. They make the construction $H$ a function of the slices.

**Proposition 5.3.**
- **Slices at every event.** Given identities and the slices at every event time, $\hat K$ is determined.
- **Slices at coarser dates.** Given slices only at coarser dates (snapshots), $\hat K$ is determined up to the events between snapshots. These are unrecorded events, and §4's inference reconstructs them as the posterior of their faces.

*Proof.* Consecutive slices with identities determine which cells began, ended or changed between them, and hence the footprints. Between two snapshots, these are attributed to the period, not to an instant. ∎

### 5.4 What a 2-Complex Carries of Three Dimensions

The discussion's difficulty: a 2-complex has no 3-cells, and full volumetric topology needs 3-cells or more information. In the volumes-and-interfaces dictionary, a volume is a coupling. What a volume is, is carried by the labels of the coupling's edge and of its faces, as a node of a spin network carries a polyhedron.

Label each interface $f$ by its area $A_f$ and its unit normal $n_f$.

**Proposition 5.4 (closure).** For every bounded volume with a piecewise-planar boundary, $\sum_f A_f n_f = 0$, where the sum runs over its interfaces with outward normals.
- It is the classical Gauss constraint of $SU(2)$ intertwiners (Bianchi, Doná and Speziale 2011).
- On real data it checks the 3D model: an omitted wall or a misplaced roof shows as a defect $\lvert\sum_f A_f n_f\rvert$.

*Proof.* The divergence theorem applied to a constant vector field. ∎

**Theorem 5.5 (Minkowski 1897).** Let unit normals span $\mathbb R^3$, and let positive areas satisfy $\sum A n = 0$. Then there is a convex polyhedron with exactly those face normals and areas, and it is unique up to translation (Alexandrov 2005, ch. 7).

So the class of 3D configurations a labelled event complex determines depends on its labels:

| Labels | What is determined |
|---|---|
| none | the adjacency: which volumes touch, through which interfaces, and when |
| areas and normals | each convex volume up to translation. The interfaces' positions, one label more, then fix the neighbours' relative positions |
| face polygons | every volume, convex or not. These are the shapes the twisted geometries of Freidel and Speziale (2010) parametrise |

Most building footprints are not convex (an L-shape, a court), so the realization layer carries polygons. Whether a volume is solid or hollow is a state of its node, not a property of the network. The discussion's example of a ball and a sphere with one 2-skeleton is answered there.

### 5.5 CityGML 3.0

**The correspondence.**
- CityGML 3.0's *Versioning* module (Versions, VersionTransitions with their transactions of insert, delete and replace) is the record of rewrites; the construction $H$ reads it.
- Its *Dynamizer* (time series of attributes, from sensors) gives the observations at vertices (OGC 2021; Kutzner, Chaturvedi and Kolbe 2020).

**What the event complex adds.**
1. The causal order across features.
2. Faces: the persistence of a relation between two features, as one cell whose cyclic order is time.
3. Homology (silences, shells), and the cut as one object.
4. The factorised likelihood, with its expansion in modes, its selection rules and its invariance under refinement.
5. The conditions of reconstruction (§5.2–5.4).

**Export.** The realization at a date can be written as CityGML at LOD1 or LOD2, and the record as Versioning's transitions.

## 6. The City as a Field Theory

### 6.1 Couplings as Quanta, Events as Typed Interactions

As in the history domain (its §8.1), each segment of a coupling is a field quantum with one strand for each primitive it holds; its rank is its valence. Each event is an interaction vertex whose bubble is its link (Proposition 2.6).

The difference from a history is the types.
- **A history.** In a history nearly every occasion is a type of its own: 110 types among Vienna's 120 occasions.
- **A city.** A city's events are typed by their rule schemas. So the theory has one interaction for each event type, finitely many:
  $$S[\varphi] = S_0[\varphi] + \sum_\tau g_\tau\, V_\tau[\varphi],$$
  where $V_\tau$ is the bubble of the rule's two sides. This answers the history domain's difficulty that the couplings of single bubbles are not identifiable: here the types come from the rules, not from isomorphism classes of bubbles.

A record is one Feynman graph of this theory; in time, it is one term of its Dyson series.

### 6.2 Rates and Their Checks

**The rate of a type.** $g_\tau$ is its base rate, estimated as
$$\hat\lambda_\tau = \frac{n_\tau}{E_\tau}, \qquad E_\tau = \int (\text{matches of } L_\tau \text{ alive at } t)\, dt,$$
with Poisson errors.

**Its law.** A rate law $r_\tau(\text{state})$ is fitted by a log-linear model: tyche's generalised linear models, or nomos's Poisson regression with absorbed effects.

**The checks.** The model is checked by martingales by era, as in the history note (its §8.4): the count of events in an era, less its compensator.

These are the rates of the generative layer. They are estimated from the record and never used to generate complexes to fit it.

### 6.3 Latent Interactions and the Sum over Complexes

Unrecorded events among several cells, such as two stores closing together with no record of it, are latent vertices. The likelihood of a record is the sum over its completions by latent vertices:
$$Z = \sum_{K' \supseteq K}\;\prod_{v\,\text{latent}} g_{\tau(v)}\; W(K'),$$
the perturbative expansion in the couplings.
- **What the engine computes.** Order zero: no latent interactions. Unrecorded changes of single cells are still summed exactly by their semigroups.
- **What the theory leaves open.**
  - the higher orders;
  - whether a city's record has a regime in which the sum over complexes concentrates: the large-N, or condensed, regime of physis decision 0006.

## 7. What the 2-Complex Earns

The discussion asks that the 2-complex earn its place against simpler structures. It does so against each in a different way.

- **Against an event graph** (a directed acyclic graph of events and their dependencies).
  - The graph holds the causal order.
  - It has no cell for persistence, and no joint life of two features.
  - Its cycles are not homology, and it has no cyclic order to make a holonomy of.
- **Against a factor graph or tensor network.** The likelihood is one (Theorem 4.2), and the complex adds no expressive power. What it adds is structure and reading:
  1. **A topological bound on cost.** The contraction pattern is the incidence of a 2-complex, so the cost of an exact contraction along time is set by the network of the moment (Proposition 4.7; history note, Proposition 6.3).
  2. **Topological selection rules.** The expansion in modes is a sum over labellings of cells, constrained at edges and vertices, so selection rules are charges conserved along surfaces. A factor graph hides this.
  3. **Slicing as description.** Ticks are subdivisions, so slicing time is a choice of description, not a choice of model (Theorem 4.5).
  4. **A state at a date.** The boundary state at a cut is a labelled graph: the GIS's state at a date as one object with its uncertainty.
  5. **One structure.** By Theorem 2.8, the same complex is the dual of the city's cellular spacetime, so a 3D GIS and its event record are one structure.
- **Against CityGML** (§5.5).

**The test.** The paper must state results an event graph cannot. Five candidates, each computed by the engine:
1. the area law of evidence, with its measured rates (Proposition 4.6);
2. a selection rule found in a record, such as the parts of a split lot correlated in their (standard, standard) mode;
3. the posterior of unrecorded events, located on faces (§4.1);
4. the invariance under refinement, shown on one record sliced by hours and by years;
5. the closure defects of a 3D model (Proposition 5.4).

# Part II. The Engineering

## 8. What Is Implemented, and What Is Not

*As built on 2026-11-04:* E0 (the rename), E1 (Melbourne's record), E2 (the module `events`, 26 checks) and the drawing of E5 (the layer and the analysis of the record, without the amplitudes). E3 and E4 are next.

| Formalism | Implemented | Stage |
|---|---|---|
| Ontology, the type system, the schema of a record (§1) | built: volumes and interfaces for Melbourne's district | E1–E2 |
| $K$, $\hat K$, the cut, homology (§2) | built, on morphe; the link of an event to come | E2 |
| Theorem 2.8 (the dual) | checked on records with global ticks | E2 |
| Modes for finite symmetric dynamics: $\mathbb Z_2$, $S_q$, $S_m\wr S_n$ (§3) | first; $\mathbb Z_N$ next | E3 |
| The event-coupled Markov model: contraction, Theorems 4.2–4.5, Propositions 4.6–4.7 (§4.1–4.6) | first, exact where the support allows, truncated otherwise | E3 |
| Fitting the memories $\mu$ and the event rates | first | E4 |
| Gauge models of relational states: value from sales, the phase of the tide, amphidromes (§4.7–4.8) | later | E6 |
| Closure and Minkowski's reconstruction (§5.4) | later, in architecture | E6 |
| The event operator and Proposition 4.9 (§4.10) | the construction's checks only | E2 |
| Rates of event types and martingales (§6.2) | first | E4 |
| Latent interactions and the sum over complexes (§6.3) | theory only | — |
| Complex amplitudes (§4.8) | with $\mathbb Z_{24}$ hours, later | E6 |
| Re-attachment of primitives, transfer corners (§2.5) | if a case needs it | — |

## 9. Real Data

Every dataset below was queried live on 2026-11-04 for its keys, its span and its rows in the study areas. These are:
- Melbourne: the Hoddle Grid box, −37.8205, 144.9535 to −37.8075, 144.9735;
- Lower Manhattan: the area and Financial District boxes of `areas.json`.

The criterion of choice is A3: an event must attach to the persistent identifier of the cell it changes.

### 9.1 What the Data Library Already Holds

The workbench's recipe already fetches:
- the City of Melbourne's pedestrian counts (102 sensors, hourly, 2009 → now, with a gap from 2022-11 to 2024-10);
- the CLUE business establishments (2002–2024, 413,550 rows);
- OpenStreetMap's buildings, ways, spaces and premises;
- Valdrada's compiled city of New York (buildings, roads, trees);
- the Brooklyn Bridge promenade's counts (2017–2019);
- weather and holidays.

None of these is yet a record of events. Two of them need a change:
- **The CLUE establishments' manifest omits `base_property_id` from its selection.** It must be added, because it is the identifier that persists across censuses.
- **Valdrada's buildings carry no BIN or BBL.** New York's building footprints must be joined to them in space.

### 9.2 Melbourne: a Panel of Uses and the Works That Change Them

| Dataset (City of Melbourne unless stated) | Cells and events it gives | Identifier | Span; rows in the box |
|---|---|---|---|
| CLUE Building information (`buildings-with-name-age-size-accessibility-and-bicycle-facilities`) | buildings as couplings: their births (`construction_year`), refurbishments, floors and predominant use at every census | `base_property_id`, `property_id`, `block_id` | censuses 2002–2024; 35,138 rows over 1,524 base properties |
| CLUE establishments (held; re-fetched with `base_property_id`) | businesses as couplings, their occupancies of properties as primitives; births and deaths by census; ANZSIC4 classes | `property_id`, `base_property_id`; trading name within a property | 2002–2024; 229,042 rows in the CBD |
| CLUE café and restaurant seating (`cafes-and-restaurants-with-seating-capacity`) | outdoor seats, a state of a frontage | `base_property_id` | 2002–2024; 11,129 outdoor rows |
| Development Activity Monitor (`development-activity-monitor`) | developments: completions by year, with floors, dwellings and floor space by use | `property_id` (up to five), `clue_block` | 2002–2026; 320 in the box |
| Building Permits (`building-permits`) | works dated to the day: permits, occupancy permits, final inspections | an address only, geocoded to `property-boundaries` | 1983 → now, daily; 87,532 in Melbourne 3000 |
| Road segments with surface type (`road-segments-with-surface-type`) | works on streets, dated by `installationdate` | `assetid`, joined to `road-corridors.seg_id` by the "X between Y and Z" description | 1910–2025; 9,912 in the box |
| The Victorian Heritage Register (the Heritage Council's API, CC BY 4.0) | registrations: institution events on buildings | `vhr_number`, overlay numbers; a place at an address | `registered_on`; 160 places on the register, 1,876 with the overlays |
| Event permits 2014–2018 | street events, Bourke Street's among them | free-text places | 2,827 rows; 276 name Bourke Street |
| Victorian liquor licences by location (Liquor Control Victoria) | licences granted and cancelled, by differencing monthly snapshots | `Licence Number` | 2022-01 → 2026-09, monthly |
| Pedestrian counts (held) | observations: walkers each way, hourly | sensor | 2009 → now |

**What Melbourne offers.** It offers what no other source here does: a yearly panel, over 23 years, of every building's use and every business's class, under identifiers that persist.
- **Modes.** ANZSIC is a classification of four levels (division, subdivision, group, class). The memory of a premises' use at each level is Proposition 3.2 on real data.
- **Splits.** A `base_property_id` whose `property_id`s divide between censuses is a split. That makes the selection rule of §3.4 testable: are the parts' uses correlated in the (standard, standard) mode?

### 9.3 Lower Manhattan: Identified Buildings, Streets and Their Closures

| Dataset (NYC Open Data unless stated) | Cells and events it gives | Identifier | Span; rows in the area |
|---|---|---|---|
| BUILDING (`5zhs-2jue`) and BUILDING_HISTORIC (`ipkp-snf6`) | buildings as couplings, with their births (`construction_year`) and, for those gone, their deaths (`demolition_year`) | `bin`, `base_bbl` | 1719 → 2025; 2,238 standing, 238 demolished since 2000 |
| DOB NOW Job Application Filings (`w9ak-ipjd`) and the legacy filings (`ic3t-wcy2`) | new buildings, demolitions, major alterations (A1, ALT-CO), dated by filing, approval, permit and sign-off | `bin`, `bbl`, job number | 1990s → now; 732 (DOB NOW) and 6,234 (legacy, Community District 1) |
| LPC Designated and Calendared Buildings and Sites (`ncre-qhxs`) | designations: a rule as a coupling (a landmark, a historic district), with primitives of application to buildings | `bin_number`, `bbl`, `lp_number` | 1965 → now; 848 |
| NYC Permitted Event Information, historical (`bkfu-528j`) | street closures for events, a street's regime from vehicles to people and back, with start and end times: the discussion's pedestrian street, on real data | the event; places as "X between Y and Z", matched to LION by name | 2008 → 2026; 4,995 with a closure in the 1st Precinct |
| Open Streets Locations (`uiay-nctu`) | scheduled closures: Peck Slip; Water Street to Pearl Street, weekdays 07:30–16:00 | LION segment ids (version 23C) | the 2024 season only |
| Bike routes (`mzxg-pwib`) | lanes installed and retired, a state of a street space, with a chain of versions | LION `segmentid`, `bikeid`, `prevbikeid` | 1966 → 2025; 699 |
| Storefronts Reported Vacant or Not (`92iy-9c3n`, Local Law 157) | the state of each storefront each year: vacant or not, and its primary business activity. These are observations of frontages | `bbl`, `bin`, unit | 2019/20 → 2026; 15,128 (230–300 a year in the Financial District, vacant on 31 December rising from 25 to 72) |
| Motor vehicle collisions (`h9gi-nx95`), 311 requests (`erm2-nwe9`) | occurrences, read only: at street spaces and junctions, and at buildings | `collision_id`, `unique_key`, `bbl` | 2012 → now; 2,477 and 37,307 in the Financial District |
| LION (City Planning) | the street spaces and junction spaces, and the identifiers the closures, Open Streets and bike routes key on | segment and node ids | versions 25c–26b; its versions must agree with each source's |

**What New York offers.** Heterogeneous events, each attached to an identifier: births and deaths of buildings, works, designations, closures of streets, lanes, and the yearly state of every storefront.
- **The pedestrian street.** The closures for permitted events and the Open Streets of 2024 are the discussion's case: the regime of a street space, a $\mathbb Z_2$ state of its cross-sections, flipped by recorded events at recorded hours.
- **Storefronts.** The storefront reports are the observations a model of frontages needs: vacancy and activity, year by year, in the years after 2020.

### 9.4 From Datasets to a Record

A product of the data library, `city_events.py`, compiles each district's record. It writes:
- **couplings:** buildings (by BIN or `base_property_id`), street spaces (by LION segment, or Melbourne's `seg_id`), junction spaces, businesses, rules;
- **primitives:**
  - party walls and frontages, from the massing and the network (architecture `space/massing.hpp`; urban-science `landscape`);
  - cross-sections, from the network's incidences;
  - occupancies and applications of rules, from the datasets;
- **events:** with their type (§1.4), their time label (a day, a year, an interval), their match by identifier and their source;
- **observations:** the censuses, the storefront reports and the counts, as observation events with the cells they read.

Joins are by identifier where the sources share one. The exceptions:
- New York's buildings are joined to Valdrada's rings in space.
- Street events given as "X between Y and Z" are matched to segments by name, with every match that is not exact reported.
- Melbourne's permits are joined to properties by geocoded address.

Each join's failures are counted in the product's report. They are never silently dropped.

### 9.5 Provenance, Licences and Cautions

- **Provenance.**
  - Every event is *recorded* (from a source), *inferred* (a change between snapshots with no recorded event, §5.3), *modelled* or *synthetic*.
  - Time labels keep their resolution: a year of construction is an interval of a year, not a day.
- **Licences.**
  - The City of Melbourne's and Victoria's data are CC BY 4.0.
  - New York's datasets carry no licence field, so the NYC Open Data terms of use apply, as for the data already held.
  - Valdrada and OpenStreetMap are ODbL.
  - Nothing is committed; everything is fetched into `data/cache/`.
- **Personal names.** DOB filings, PLUTO, SLA licences and Victoria's liquor licences name owners and licensees. The queries select their columns so that no name is fetched.
- **Known limits.**
  - Past versions of PLUTO were not found at their old address (unverified).
  - Open Streets holds only 2024.
  - New York's closures for construction keep only a live window.
  - The SLA's issue dates were reset in 2023.
  - Some DOB dates are text, and some BBLs are malformed.
  - Melbourne's road closures and works are no longer published.

### 9.6 The First Case: Melbourne

The user chose Melbourne (2026-11-04). It is the first case for both the fit (stage E4) and the first foam drawn (stage E5).

**Melbourne's Hoddle Grid, 2002–2024: what a city remembers of its uses.**
- **The district.** The district of Landscape Dynamics, Elizabeth–Russell and Flinders–Lonsdale (landscape note, §2.2, §7.4):
  - 356 junction spaces and 399 street spaces;
  - 440 buildings sharing 565 walls.
- **Couplings, states and observations.** CLUE's buildings and establishments: the buildings by `base_property_id`, the businesses with their occupancies, their uses and classes at every census.
- **Events.**
  - the development monitor's completions;
  - permits and road works;
  - registrations on the heritage register;
  - for the day's scale, the street events of 2014–2018 from the event permits (Bourke Street's among them, matched to segments by name).
- **Observations.** The pedestrian counts, which give the day's scale its walkers.
- **What it yields.**
  - the memories $\mu$ at each level of ANZSIC (Proposition 3.2);
  - splits of properties and their selection rule (§3.4);
  - the posterior of unrecorded changes between censuses;
  - the foam in both scales, decades for the uses and a day for a street event.

**The Financial District follows,** with what Melbourne lacks: closures of streets for events (2008–2026) and the Open Streets of 2024, keyed by LION; buildings by BIN with their filings, designations and demolitions; the storefront reports.

## 10. Placement and Reuse

| Piece | Where | New or reused |
|---|---|---|
| $K$ and $\hat K$ from a record | morphe `topology/worldlines.hpp` (`Encounters`, `worldline_complex`, `refined_worldlines`) | reused unchanged |
| Transfer corners, a pair whose ends change | morphe `topology/worldlines.hpp` | later, only if a case needs it |
| Causal order, Foata steps, linearisations | morphe `formal/traces.hpp` | reused |
| Rewrites of the spatial network | morphe `topology/complex_rewriting.hpp`; urban-science `landscape/district.hpp` | reused |
| Characters, central idempotents, actions, Fourier on finite groups | morphe `algebra` (`characters.hpp`, `group_algebra.hpp`, `action.hpp`, `abelian_fourier.hpp`) | reused; the isotypic projectors of a permutation module are added to morphe if absent |
| The order of a worldline complex's cells in time, for contraction | now the history domain's `time_order` (`foam/transfer.hpp`) | moved down to physis, beside the contraction: its second user (CLAUDE.md §2.2, an abstraction once its pattern has appeared twice) |
| Contraction of an event-coupled Markov model on a worldline complex: forward and backward messages along time, biorthogonal modes, truncation, the laws at cuts, gradients | physis `statistical_mechanics` (a new header, physis decision 0020) | new and general: it names no object in nature |
| The exact $\mathbb Z_2$ contraction | physis `gauge/contraction.hpp` | reused for gauge models |
| Optimisation, jets | calx | reused, as the history domain's fit is |
| A city's record, its schema and event types, the construction, the check of A1, the model, the realization | urban-science, a new module `events` | new |
| Buildings as a 3-complex with shared walls | architecture `space/massing.hpp` | reused; the closure check later |
| The walking network and its district | urban-science `landscape/district.hpp` | reused |
| The datasets and the record of a district | data library: registry entries, a product `city_events.py` | new |
| The plug-in | the workbench: `plugins/event-spacetime/` | new |
| The opacity of the city | the workbench's shell, `src/app.js` | a new control |

*Remark (the order of the domains).* The module `events` works on a typed record and names no building. The network of volumes and interfaces is assembled in the plug-in's engine, which may compose both domains because it is a work. If the user lets urban science declare architecture, urban science becomes their bridge, as CLAUDE.md foresees, and the assembly can move into the domain (§14).

## 11. Algorithms and Their Cost

- **Construction.** Linear in the cells, plus $O(\sum_v m_v^3)$ for the triangles.

  For Melbourne's district (§9.6) the network of volumes and interfaces has:
  - about 1,200 couplings of space: 440 buildings, 399 street spaces and 356 junction spaces;
  - about 2,000 primitives of space: 565 party walls, the frontages, and two cross-sections for each street space;
  - the businesses and their occupancies, a few thousand at each census.

  Over 2002–2024, with its events and yearly observations, $K$ has some $10^5$ cells.

  Most businesses touch no other business at any event, so under A5 their cells form small clusters.
- **Contraction.** A forward pass in the order of time, then a backward pass.
  - **The message.** The forward message is a joint law of the live cells' modes, stored sparsely by its support.
  - **Why the support stays small.** Under A5, cells become correlated only at events, and a correlation between two cells decays as $e^{-(\mu_a + \mu_b)\tau}$. So the support is a union of clusters of cells joined by recent events.
  - **Exact or truncated.** The contraction is exact where the support stays below a limit (as physis's `ContractionLimits`). Otherwise components below a relative weight $\varepsilon$ are dropped, and the bound on the mass dropped is reported with every result.
  - **Later.** Forward filtering with backward sampling, where truncation is not enough.
- **Fitting.** $\log Z$ and its gradient, by the forward–backward pass on calx's jets; Newton's method from moment estimates; Laplace's errors (as the history domain's `fit_couplings`).
- **Checks.**
  - Brute-force enumeration on small complexes, for Theorems 4.2–4.5 and Propositions 4.6–4.7.
  - The model simulated on the real complex, then fitted, its rates recovered.
  - Proposition 2.4's Euler characteristic and $\partial\partial = 0$ on every record built.

## 12. The Plug-in and Its Visual Expression

### 12.1 The City Made Transparent

The shell gains a control for every plug-in: the city's opacity, continuous from solid to a ghost, over buildings, spaces, streets and trees. It sits beside the existing forms (as built, lowered, of glass), which are whole-layer settings today (`setBuildingForm`, app.js).

The plug-in also offers the city *as realised*. It draws $\mathcal R(t)$ at the clock's date: buildings not yet built are absent, and demolished ones are gone. This makes the realization map visible.

### 12.2 The Foam in Spacetime

**Time.** Time is height, $z = s\,(t - t_0)$, in one of two scales:
- **decades:** about 10 m a year, for structure and institutions;
- **a day:** 10 m an hour, for regimes and counts, as in Landscape Dynamics.

The origin is the ground, or the skyline (Lower Manhattan's towers reach 541 m, through any foam that starts at the ground).

**Edges (couplings).** Worldlines rise from the volume's anchor: a building's centroid, a street space's middle, a junction. They are thick while the coupling lives, coloured by type and, along their length, by state.

**Faces (primitives).** Translucent ribbons span the two worldlines they join, cut into lenses at the events that meet them. They are coloured by type, by state, by the probability of a non-trivial mode, or by the memory $e^{-\mu\tau}$. Each primitive's own geometry (a wall's line, a cross-section) can be drawn as a faint sheet extruded in time.

**Vertices (events).** Beads at their place and time.
- They are coloured by category, with a glyph by provenance (recorded, inferred, modelled, synthetic), as Cosmopoiesis's provenance hues are.
- In $\hat K$, there is a bead at each attending coupling, the beads joined by rungs.
- In $K$, there is one bead, with the worldlines bending into it: the event as a vertex of a Feynman graph in spacetime.

**Also drawn.**
- plaquettes, as faint horizontal polygons at their events;
- the causal threads, as in Landscape Dynamics.

**The cut.** A plane at the clock's date carries the spin network of the moment: nodes at the live couplings, links for the live primitives, coloured by state or by the weights $w_k$. Below it is the past; above it, faded, the future.

### 12.3 Views of One Complex

- **Spacetime points or events.** A switch between the refined complex $\hat K$ and the event complex $K$: the user's dual representation, made one gesture.
- **The link of an event.** Its bubble (before, after, and what passes), drawn small beside the bead, with its rule, its kernel and its amplitude.
- **The foam.** A chosen mode's surfaces: the cells with $P(k_\sigma \ne 0)$, spanning from where states are set to the observations that read them.
- **The unrecorded.** The posterior rate of unrecorded change along each cell: where the record is silent and the data say something happened.

### 12.4 Tools

- **The record.** Its sources, its counts by category and by provenance, and its silent changes.
- **An event.** Its rule, its footprint, its link, its kernel, its amplitude, and its causal past and future.
- **The foam.** Modes and memories, foams across a cut, and the area law.
- **The fit.** The rates by event type, with their errors; the memories $\mu$ by type of cell; the invariance under refinement shown.
- **Order.** Commutators of chosen events, and linearisations.
- **Reconstruction** (later). Realise a date, and show the closure defects.

**Palette.** The workbench's sea blue and sand; categories in the visualization layer's validated colours; provenance in hues validated for colour vision deficiencies.

## 13. Stages and Their Checks

| Stage | Content | Checks |
|---|---|---|
| E0 | Rename the Urban Building Workbench the Urban Science Workbench: its directory, page title, READMEs, CLAUDE.md, the README's §6, the parent tables and the test runner. Old decision records keep their paths, as CLAUDE.md §6.3 says. A new decision record (ADR 0067) records the rename with this design | the workbench's suite under its new path |
| E1 | Data: the event datasets of §9 fetched into the cache; a product compiling each district's record (couplings, primitives, events, snapshots) with provenance; a small record as a fixture | the product's checks; fixtures with their sources and dates of retrieval |
| E2 | Domain: the `events` module (record, schema, types, construction on morphe, cut, link, homology, the check of A1) | ∂∂ = 0; Proposition 2.4's χ; the collapse; the cut equals the network of the moment; Theorem 2.8 on synthetic records; Proposition 4.9's attachments |
| E3 | physis: the contraction of event-coupled Markov models and the time order moved down; morphe's isotypic projectors if absent | Theorems 4.2–4.5 and Propositions 4.6–4.7 against brute force; Proposition 3.2's rates; Proposition 3.3's zeros |
| E4 | Domain: the first model on real data, the fit, the posteriors of unrecorded changes, the rates of types with their martingales | recovery of the rates from records simulated on the real complex |
| E5 | The plug-in: layers (the foam, the cut, the events), tools, the city's opacity | headless Chrome, on DevTools ports agreed with the other sessions |
| E6 | Later: gauge models of relational states, closure and reconstruction, $\mathbb Z_{24}$ phases, transfer corners | each against its closed form |

## 14. Choices for the User

Decided by the user on 2026-11-04:
1. **The plug-in's name:** *Event Spacetime*.
2. **The dictionary:** volumes and interfaces (§1.3).
3. **The first case:** Melbourne (§9.6).
4. **The discussion's home:** moved beside this note and named by the convention, [2026-11-04-event-based-gis-on-a-2-complex-discussion.md](2026-11-04-event-based-gis-on-a-2-complex-discussion.md).

Open:

5. **The order of the domains.** The network of volumes and interfaces joins architecture's buildings to urban science's streets. Two ways to do it:
   - **Keep the two domains independent.** The join stays in the plug-in's engine, which as a work may compose both.
   - **Let urban science declare architecture.** It becomes their bridge, at height 1, as CLAUDE.md foresees. The join then moves into its `events` module, where the domain's own checks hold it, and architecture still depends on nothing but the foundations.

# Part III. The Paper

## 15. Outline

*Event-Based 3D GIS on a 2-Complex: A Spin-Foam Approach.*

1. **Introduction.** A geographic state as the realization of a constrained event history; the three representations.
2. **Related work.**
   - Event-based and spatiotemporal GIS: Langran 1992; Peuquet 1994; Peuquet and Duan 1995; Hornsby and Egenhofer 2000; Worboys 2005; Galton and Worboys 2005; the event-driven updating of 3D geodatabases that the discussion cites (Lv et al. 2016, to be verified).
   - Four-dimensional topology: Arroyo Ohori, Ledoux and Stoter 2017.
   - CityGML 3.0: OGC 2021; Kutzner, Chaturvedi and Kolbe 2020.
   - Spin foams and group field theory: Reisenberger and Rovelli 1997; Baez 1998; Oriti 2009; Perez 2013.
   - Stochastic graph rewriting: Danos et al. 2007; Behr, Danos and Garnier 2016.
   - Group-based phylogenetics, the same algebra on trees: Evans and Speed 1993; Székely, Steel and Erdős 1993.
3. **Ontology** (§1).
4. **The complexes and their duality** (§2): Theorem 2.8 and Proposition 2.6.
5. **Labels** (§3): Propositions 3.2, 3.3 and 3.5.
6. **Amplitudes** (§4): Theorems 4.2–4.5 and Propositions 4.6–4.7; the gauge models; signs and phases; the meaning of an amplitude.
7. **Realization and reconstruction** (§5).
8. **Generation** (§6).
9. **Implementation**: Part II as built, with §8's table stating what is and is not implemented.
10. **Case studies**: §9's districts, results only from what is built.
11. **Discussion**: §7, what the complex earns, and the limits (A5, order zero in latent interactions, the classical limit only).

Appendices hold the proofs. The paper presents §1–§7 in full, marks what the engine implements, and reports numbers only from the engine.

## 16. References

- Alexandrov, A. D. (2005), *Convex Polyhedra*, Springer (from the Russian of 1950).
- Arroyo Ohori, K., Ledoux, H. and Stoter, J. (2017), Modelling and manipulating spacetime objects in a true 4D model, *Journal of Spatial Information Science* 14, 61–93.
- Baez, J. C. (1998), Spin foam models, *Classical and Quantum Gravity* 15, 1827–1858.
- Behr, N., Danos, V. and Garnier, I. (2016), Stochastic mechanics of graph rewriting, *Proceedings of LICS 2016*, 46–55.
- Bianchi, E., Doná, P. and Speziale, S. (2011), Polyhedra in loop quantum gravity, *Physical Review D* 83, 044035.
- Breiger, R. L. (1974), The duality of persons and groups, *Social Forces* 53, 181–190.
- Danos, V., Feret, J., Fontana, W., Harmer, R. and Krivine, J. (2007), Rule-based modelling of cellular signalling, *CONCUR 2007*, LNCS 4703, 17–41.
- Diaconis, P. (1988), *Group Representations in Probability and Statistics*, IMS Lecture Notes 11.
- Evans, S. N. and Speed, T. P. (1993), Invariants of some probability models used in phylogenetic inference, *Annals of Statistics* 21, 355–377.
- Freidel, L. and Speziale, S. (2010), Twisted geometries: a geometric parametrization of SU(2) phase space, *Physical Review D* 82, 084040.
- Galton, A. and Worboys, M. (2005), Processes and events in dynamic geo-networks, *GeoS 2005*, LNCS 3799, 45–59.
- Hornsby, K. and Egenhofer, M. J. (2000), Identity-based change: a foundation for spatio-temporal knowledge representation, *International Journal of Geographical Information Science* 14, 207–224.
- Kschischang, F. R., Frey, B. J. and Loeliger, H.-A. (2001), Factor graphs and the sum-product algorithm, *IEEE Transactions on Information Theory* 47, 498–519.
- Kutzner, T., Chaturvedi, K. and Kolbe, T. H. (2020), CityGML 3.0: new functions open up new applications, *PFG* 88, 43–61.
- Langran, G. (1992), *Time in Geographic Information Systems*, Taylor & Francis.
- Migdal, A. A. (1975), Recursion equations in gauge field theories, *Soviet Physics JETP* 42, 413–418.
- Minkowski, H. (1897), Allgemeine Lehrsätze über die convexen Polyeder, *Nachrichten der Gesellschaft der Wissenschaften zu Göttingen*, 198–219.
- Munkres, J. R. (1984), *Elements of Algebraic Topology*, Addison-Wesley.
- OGC (2021), *OGC City Geography Markup Language (CityGML) Part 1: Conceptual Model Standard*, version 3.0, OGC 20-010.
- Oriti, D. (2009), The group field theory approach to quantum gravity, in *Approaches to Quantum Gravity*, Cambridge University Press, 310–331.
- Perez, A. (2013), The spin-foam approach to quantum gravity, *Living Reviews in Relativity* 16, 3.
- Peuquet, D. J. (1994), It's about time: a conceptual framework for the representation of temporal dynamics in geographic information systems, *Annals of the Association of American Geographers* 84, 441–461.
- Peuquet, D. J. and Duan, N. (1995), An event-based spatiotemporal data model (ESTDM) for temporal analysis of geographical data, *International Journal of Geographical Information Systems* 9, 7–24.
- Reisenberger, M. P. and Rovelli, C. (1997), "Sum over surfaces" form of loop quantum gravity, *Physical Review D* 56, 3490–3508.
- Rovelli, C. and Smolin, L. (1995), Discreteness of area and volume in quantum gravity, *Nuclear Physics B* 442, 593–619.
- Serre, J.-P. (1977), *Linear Representations of Finite Groups*, Springer.
- Székely, L. A., Steel, M. A. and Erdős, P. L. (1993), Fourier calculus on evolutionary trees, *Advances in Applied Mathematics* 14, 200–216.
- Witten, E. (1991), On quantum gauge theories in two dimensions, *Communications in Mathematical Physics* 141, 153–209.
- Worboys, M. (2005), Event-oriented approaches to geographic phenomena, *International Journal of Geographical Information Science* 19, 1–28.
