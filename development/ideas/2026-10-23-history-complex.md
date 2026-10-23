---
title: The History Complex: histories as social spin foams, from the record to the formalism
status: seed
date: 2026-10-23
related: [morphe topology, physis gauge, physis quantum_gravity, physis symmetry, physis field_theory, the inference foundation, the Complex Atlas, categories]
---

# The History Complex: Histories as Social Spin Foams, from the Record to the Formalism

*Promoted in part, 24 October 2026:* the carrier, the record and its construction, and the spin foam at ℤ₂ are now the history domain ([codebook](../../domains/history/codebook.md), [theory note](../../domains/history/theory/the-history-complex.md)) and the History Complex ([ADR 0041](../decisions/0041-the-history-domain.md)). The other groups (D₁₂, O(d), SU(2) dynamics), the generating field theory and its condensate, and the tests of §8 against baselines remain ideas here.

*The user's papers are discussion papers, and this note takes its own position where it finds reason to. It names what it takes from them (§1), sets out the carrier and its semantics (§2–4), and then the empirical side, which the user asked to be thought through before anything is built: what is recorded (§5), how the complex is built from it (§6), how values become group elements (§7), what the formalism computes and tests on a record (§8), where records come from (§9), what is hard (§10), and a first history (§11). §12 lists what the note changes in the papers.*

## 1. The Point of Departure

- **The event–coupling–primitive architecture** (*Social Spacetime Topology and Generative Grammars in CW Complexes*, 2026): a 2-complex K_{Q,ρ}, relative to a research question and a resolution, with vertices as events, edges as typed couplings and faces as primitives; persistent beings as distributed patterns rather than cells (§2.5); higher cells for composition; a cellular generative grammar Γ = (T, R, A) of interface-preserving rewrite rules, and a generative history as a realized derivation, chronology being the derivation and not a coordinate on the carrier (§6–7, Fig. 1).
- **The social spin foam** (*A Social Spin Foam Model for Power Analysis and Governance Decision-Making*, 2026): facts of any type on faces, couplings on edges, the joint instantiation of couplings at vertices (Def. 6.1, Ass. 5.3); admissibility as the support of a non-negative coupling weight (Prop. 6.4); configuration weights as products of local weights (Def. 6.8); histories with boundary and transition weights as sums over histories (Def. 7.1–7.2); non-negative weights unless the record criterion is met (Cor. 7.4, Rem. 7.5); weights estimated or stipulated, no group and no action (§7.2, Table 1); power as a divergence of normalised weights under an act (Def. 10.3); reconstruction as an ensemble (Def. 8.5).
- **The street spin foam** (*Street Spin Foam*, 2026): the same carrier with a group: valuation frames with no view from nowhere, a connection on 1-cells, holonomy on 2-cells, an alignment weight that is a class function, and the character expansion deriving spins, intertwiners and vertex amplitudes (§5–6). Two orientations of the cells and a counting argument between them (§4.4, §8.6); four open questions (§10.1): which group, whether to sum over complexes, what boundary states are, whether the vertex amplitude has content.
- **Generative Relational History** (*A Theory of Historical Generativity under Constraints*, 2026): constraint as a macroscopic phase of the relational field (Claim 1), with GFT condensates among its sources; the limits: no independently motivated symmetry, its place to be taken by invariants preserved under rewriting of the grammar (§3.5.1); no action principle, since an extremised historical magnitude is teleology (§3.5.4); no order parameter exhibited (§8.1.3).
- **The GFT music spin foam** (2026): pitch class ↔ spin j of SU(2), octave ↔ a U(1) charge, harmonic closure ↔ coupling to the singlet, and a positive vertex weight, the number of recoupling channels (Eq. 16).

## 2. The Proposal: One Carrier for the Three Readings

**The fault.** The social spin foam's "two-complex" is a set of faces, edges and vertices with incidence relations (SSF §4.1): a hypergraph with a third level, not a CW complex. A 2-cell of a CW complex is a disc, and its boundary is a *cyclic sequence* of edges and vertices; that cyclic order is the whole of what a 2-complex carries beyond a hypergraph, and every topological quantity the papers borrow (holonomy, ker ∂₂, homology, the spin foam sum) is computed from it. In physics the order means something: the adjacency of the simplices dual to the face. In the facts-on-faces reading it means nothing, and the street paper's orientation II shows the symptom: to close the boundary of an object-face "an object of odd degree has one membership dropped" (§8.6). A representation in which data are discarded to satisfy the carrier has the carrier the wrong way round.

**The fix.** Give the cyclic order a meaning, and the only meaning available in a history is time. A 2-cell must be the lifetime of something that runs between two things that persist. The reading that does this is the one of spin foams themselves, the worldsheet of a spin network, applied to the right network:

| Cell | Is | Carries | The user's term |
|---|---|---|---|
| **vertex** | an **occasion**: an interaction at which beings meet and the ties between them are recoupled: a letter, a congress, a publication, a trial, a treaty, a performance, a death; dated | the vertex amplitude A_v | event |
| **edge** | a **being between two occasions**: a segment of the worldline of a person, a school, an institution, a work, a journal, a state, oriented by time; a being *is* the coupling of the occasions it attends and the intertwiner of the ties it holds | the intertwiner ι_e | coupling |
| **face** | a **tie** between two beings across the interval between two of their meetings: a lens bounded by the two worldlines and pinched at the meetings; a relation's life is a chain of lenses, open at the record's boundary | the spin j_f | primitive |
| **a theme** | an idea, a motif, a doctrine, a norm, an office, a dispute: a **type** of ties; its history is the sub-complex of its type, its chain of custody the connectivity of that sub-complex | the group G_τ and stiffness β_τ of its type | type τ |
| **the boundary at a date** | the **spin network of the moment**: beings as nodes, the live ties between them as links | a state | configuration S |
| **the whole 2-complex** | the **history** | Z, the weight | history h |

Its attaching words are meaningful. A lens of the tie between A and B, between their meetings o and o′, has the word A[o, o′] · B[o, o′]⁻¹, where A[o, o′] is the path of A's segments from o to o′ (A may attend other occasions between). A dyadic relation born at v₁ and ending at v₂ is the street paper's bigon (Eq. 3), subdivided at the meetings. An office passed from holder to holder is a chain of lenses all bounded on one side by the institution's worldline. An idea passed from teacher to student is a tie between them of the idea's type. A relation among k beings is held through a being of its own (a school, a work, a treaty) with a tie to each: Breiger's duality of persons and groups (1974). Nothing is dropped and nothing is spurious.

*Remark (a revision).* The first draft made a theme one face and its chain of custody the word. A chain of custody is a path; it closes into a disc only when the theme is held between two carriers. Working out the record (§6) showed that every 2-cell is a tie between two beings, and that the theme is the type of ties. The carrier is unchanged; its faces are finer, and they are where attested values sit (§6.1).

**Why this orientation and not the other two.**

1. It *is* the user's event–coupling–primitive architecture, read literally: vertices are events, edges couplings and faces primitives, with "coupling" taken as what an entity does (it couples the occasions it attends) rather than as a pairwise link. It is also the Feynman-graph reading of a group field theory: an interaction vertex is an occasion, a propagator a being carried between occasions, a strand a tie. So the generation of histories by a field theory (physis decision 0006) is native to it; the facts-on-faces reading cannot offer that.
2. Time is native: edges are oriented by date, the grammar's derivation orders the occasions, and the boundary state at a date is a spin network of beings and ties. This answers the street paper's third question: a boundary state is the relational state of the world at an instant, and the transition weight between two such states is the weight of the interval (the amplitude of a foam with boundary, SSF Def. 7.2).
3. The counting argument (dim ker ∂₂ ≥ n₂ − n₁) asks that lenses outnumber being-segments, and the face-degree condition asks that every being hold at least three live ties between any two of its occasions. A richly recorded history does both. A being holding one tie alone is forced to charge zero by the Gauss law, which is right: a charge needs something to be measured against.
4. Facts of any type (SSF Ass. 5.3) are kept: a relation, an idea, a norm, an attribute or an event-as-fact is the type of the ties that carry it. A relation between relations (SSF Eq. 6) is held by a being that holds the two ties it bears on, which is where such relations are held in practice. Entities are persistent patterns (SSF Def. 5.1): a being is the chain of its segments, identified by the continuity of its intertwiners.
5. The street paper's orientation I (entities on edges, relations on faces) is this reading with occasions as vertices. Its orientation II (objects on faces) is what the counting argument reached for when occasions were edges and pairwise relations the faces. With occasions as vertices both difficulties go.

**Two structural consequences, before any weight.** With no 3-cells, ker ∂₂ = H₂(C; ℤ): the configuration space of an abelian spin sector (Street §4.2, §6.1) is the second homology of the history, its closed shells of ties, which morphe computes with generators. H₁ and π₁ are the loops of beings between occasions that no tie fills, on which flat connections wind (Hom(π₁, G), physis `bf.hpp`). Torsion is a loop that closes only after several traversals.

## 3. Where the Group Comes From, and What Is Derived

The social spin foam assumes no group and stipulates its weights; the street spin foam assumes one and derives them; GRH holds the social case has no independently motivated symmetry. The group is in the user's own framework already, in the redundancy of description:

**The frame group G of a type of tie is the group of re-descriptions that the inquiry's observables and admissibility predicates cannot distinguish.** Valuation, attribution, priority, canon and the reading of a norm are described in frames chosen at each occasion. A statement about a history is meaningful when it survives an independent change of frame at every occasion. That is gauge invariance, and G is the gauge group, motivated by the record and not by physics (§7 shows it is already in the literature for three kinds of tie). The group is to be *found*, by the street paper's own test (§10.1): under an abelian G the order of a being's occasions leaves the invariants unchanged, and under a non-abelian one it does not, so ordered records decide. The invariants GRH §3.5.1 asks for are the gauge-invariant quantities (holonomies, class functions), together with the carrier's topological invariants under the grammar's moves.

Given G and, for each type, an alignment stiffness β_τ (how strongly a tie holds its two ends to one relative frame), the street paper's lattice gauge theory transposes cell for cell. The parts:

- a connection g_e on beings (the transport of a being's frame between two occasions);
- the holonomy h_f around a lens (the mismatch of the two beings' frames accumulated between their meetings: the tie's tension);
- the weight e^{−β_τ(1 − Re χ(h_f))}, a Gibbs measure and not a principle.

The character expansion (Street §6.1; physis `gauge/exact.hpp`, `symmetry/weyl.hpp`) then derives what SSF stipulates:

| SSF object (stipulated) | derived from (G, β) | meaning |
|---|---|---|
| face label j(f), face weight A_f | the irrep j_f of G, with weight d_j c_j(β_τ) | the lens's **mode**: the Fourier dual of the tie's change across it, its resolution (§7.3). Spin 0 alone: the tie's value at one meeting says nothing about the next. High spins: the tie is held sharply. U(1): a winding; ℤ₂: whether the lens lies on a closed surface of flips; SU(2): its variety 2j + 1 |
| coupling weight C_e, admissibility | dim Inv(⊗_{f ∋ e} V_{j_f}), the intertwiner space of the being-segment | C_e > 0 exactly when the charges of the ties a being holds can balance: the Gauss law, so SSF Prop. 6.4 becomes a theorem. ι_e is the being's **mode of holding its ties together**, and dim Inv its internal freedom. At large spins the being is a quantum polyhedron whose faces are its ties (its shape) |
| vertex weight A_v | the contraction of the intertwiners of the beings meeting at the occasion: a {3nj}; the 6j for an occasion of four beings. On the refined complex, the partition function of the occasion's ties as they stand there (§7.3) | the **agreement of the accounts** of the occasion from the sides of its participants. The 6j is the change of basis between two orders of combination. It is zero when the triangle fails, and then the occasion cannot be jointly instantiated (SSF Def. 6.7) |
| W(S), A(S_in → S_out) | the spin foam amplitude with boundary | the weight of a configuration and of a transition between two states (SSF Def. 6.8, 7.2) |
| the sign | positive for ℤ₂, U(1) and for the channel count of the music paper; alternating for SU(2)'s 6j | SSF's convention holds automatically in the first three; the signed 6j only within a fibre of the record map (Rem. 7.5); physis's sampler reports ⟨sign A⟩ |
| power P_x = D(Â_x ‖ Â) (Def. 10.3) | Z under an act: a β changed, a tie removed, a being's channel forbidden, a type added | computable: exactly on small complexes by enumeration, by Metropolis on larger ones |

**The classical limit, without an action principle.** At large spins the 6j symbol is cos(S_Regge + π/4)/√(12πV), so the *weight* concentrates on labellings whose Regge geometry is stationary. The deficit angle around a lens is then its accumulated tension, and the dihedral angles are those of the quantum polyhedra of its carriers. Well-attested, much-structured relations behave geometrically; sparsely attested ones are quantum, superposed. This is where a measure concentrates, in the sense of GRH's own mean-field argument (§4.2). No magnitude that history extremises is asserted, so GRH §3.5.4 is respected. The distinction the framework needs is between a functional that trajectories *extremise* (excluded) and a weight e^{−S} that configurations *carry* (admitted; SSF Def. 6.8 is one).

## 4. The Grammar as a Field Theory, and the Phase as a Condensate

- **The grammar is a set of Feynman rules.**
  - A rule of Γ attaches cells along a preserved interface; a group field theory's interaction attaches a bubble along propagators.
  - A quantum of the field is a being with its d ties (a d-valent node of the spin network). An interaction is an occasion that gathers several beings and pairs their ties. A realized derivation H is a Feynman graph with its order of generation.
  - The sum over histories with a fixed boundary (SSF Def. 7.2) is an n-point function. The sum over all complexes the grammar generates, weighted by the rules' couplings and the spin-foam amplitude, is the perturbative expansion of physis decision 0006.
  - The admissibility predicates A_r are the supports of the rule weights: SSF Prop. 6.4 at the level of generation.
  - This answers the street paper's second question inside the user's framework: the complex is summed over by the grammar, with controlled weights. §8.5 states how such a theory is fitted to a record.
- **The constraint phase is the condensate.** GRH's Claim 1 names a macroscopic phase. The condensate of the generating theory (physis `condensate.hpp`) is its mean field σ: the relational field coarse-grained, the stage. Three candidates follow:
  - its modulus by kind, |σ_j|², for the order parameter that GRH §8.1.3 leaves unexhibited;
  - the critical coupling of the expansion (λ_c, where generated complexes grow without bound) for the boundary of stability of §5.4;
  - the bounce, a least volume below which the stage does not contract, for the departure from a stable stage (§4.5): a stage reconstitutes, it does not begin from nothing.

  These are candidates that the record must test.
- **Invariants preserved under rewriting** (GRH §3.5.1; Topology Def. 8.1–8.2): the homology and π₁ of the carrier under the grammar's moves, and the transition weight under refinement of the record. Under Pachner moves Z is Turaev–Viro's invariant; under the forgetting of meetings, §8.4 shows it is exact in two dimensions. morphe computes the first; physis the second where it holds.

## 5. The Record: What Is Collected

### 5.1 Three Layers, Kept Apart

- **The record:** assertions, each from a source, with its locator, the person or program that read it, and a confidence. Assertions are never merged into facts; two sources that disagree are two assertions. This is the factoid model of prosopography (Bradley and Short 2005), used by the Prosopography of the Byzantine World and of Anglo-Saxon England.
- **The construction:** a deterministic map from a consistent choice of assertions, a **reading**, to a typed 2-complex (§6). It makes no judgement. Conflicting or uncertain assertions give several readings, and so an ensemble of complexes: the reconstruction ensemble of SSF Def. 8.5 is the image of the readings.
- **The model:** labels and weights on the complex. Spins, intertwiners and unattested values are latent, couplings are fitted, and observables are predicted and compared with what was held out.

The empirical content is chiefly the complex itself and the values attested on it. Spins are never coded by hand: a coder who writes "spin 3/2" on a friendship has made the model into the record.

### 5.2 Three Kinds of Assertion

1. **Occasion:** something datable at which beings meet, attested by a source: a letter, a meeting, a lesson, a publication, a performance, an exhibition, a treaty, a session, a battle, a trial, a death.
   - The test is a date and a list of participants.
   - A period, a movement or a war of years is not an occasion. It is a type of ties, a sub-complex, or an occasion at a coarser resolution (§6, item 4).
   - An occasion's participants include every being whose ties it changes, each in its role: the author and the recipient of a letter, the reviewer and the reviewed, the signatories and the parties bound, the dead and the bereaved.
2. **Presence:** (being, occasion, role). These are two-mode data: beings and the occasions they attend. They are the most abundant and most reliable kind of historical datum: membership rolls, signatories, concert programmes, exhibition catalogues, lists of delegates, correspondents, co-authors.
3. **Tie:** a dyadic relation between two beings, of a type, born and ended at occasions (or open), and attested at joint occasions. Where its type has one, the tie carries a **value** at each attestation (§7): friend or enemy, ally or opponent; the transformation between two statements of a motif; the alignment between two usages of a concept. A relation among more than two is held through a being of its own (Breiger 1974).

Beings have kinds: person, group, institution, state, work, journal. Themes are types of tie, with a discipline among their labels.

### 5.3 Conventions of the Record

- **Identity** comes from authority files where they exist (Wikidata, GND, VIAF, ULAN, MusicBrainz, ISNI), and from our own identifier otherwise.
- **Dates** are written in EDTF (ISO 8601-2:2019), which states uncertainty and approximation ("1676-06?", "1712~", "1675/1676"), together with the **calendar** of the source:
  - England dated Old Style until 1752: Julian, with the legal year beginning on 25 March.
  - The Protestant German states kept the Julian calendar until 1700.
  - Russia kept it until 1918, which is why the October Revolution fell in November.

  Conversion goes through the Julian day. Metrology holds the Julian day for time scales, but not yet the civil calendars.
- **Places** come from a gazetteer (GeoNames, World Historical Gazetteer, Pleiades), and **periods** from PeriodO.
- **CIDOC CRM** (ISO 21127): the three kinds map onto it.
  - An occasion is an E5 Event or E7 Activity.
  - A presence is its P11 *had participant*, with a role.
  - A tie is an attribute assignment (E13) relating two actors.
  - The CRM even has E10 *Transfer of Custody*.

  Our own record is a plain file of lines (JSON Lines), one assertion each. It is exported to CIDOC CRM or Linked Art when another project needs it.
- **Codebooks:** for each discipline, the kinds of occasion, the roles, the types of tie and their values, each with a definition and examples. Two readers code a sample, and their agreement is reported (Krippendorff's α). The codebook is versioned: a change of codebook is a change of the record.

## 6. The Construction: Coupling Is Derived

The coupling of occasions is not chosen by the coder. Given a reading, the complex follows:

1. **Worldlines.** The presences of a being, ordered by date, are its points, and consecutive presences are its segments: the edges. That *is* the coupling: two occasions are coupled by a being exactly when it attends both with none between. A being's presences must be totally ordered, since no being is at two occasions at once. Where the dates leave them partially ordered (a reply after its letter, both undated), the reading chooses a linear extension, and the ensemble samples them.
2. **Ties to lenses.** A tie between A and B that lives across their joint occasions o₀ < o₁ < … < o_m contributes the lenses A[o_i, o_{i+1}] · B[o_i, o_{i+1}]⁻¹: one 2-cell each, pinched at every meeting. A tie born or ended where the record does not reach has an open lens at the record's boundary.
3. **Themes to sub-complexes.** The history of an idea, a motif, an office or a dispute is the sub-complex of the lenses of its type, and its custody is that sub-complex's connectivity. The history of a discipline is the sub-complex of its colour, and the integrated history is the whole. An occasion where colours meet is interdisciplinary by construction, not by a coder's label.
4. **Resolution.** An occasion at one resolution is a sub-complex at a finer one: the Congress of Vienna is a vertex, or else its sessions, dinners and letters.
   - Coarse views are computed from fine ones by two operations: collapsing sub-complexes (a cellular quotient, and a homotopy equivalence when the collapsed part is contractible), and forgetting meetings, which merges consecutive lenses.
   - The resolution is stated, and the model must relate the views (§8.4).
   - An event in Sewell's sense (2005), a sequence of occurrences that durably transforms structure, then becomes measurable: an occasion, or a collapsed sub-complex, across which the spin network of the moment changes in many ties.

### 6.1 The Refined Complex, Where Values Are Attested

Values are attested at meetings, and a vertex has no room for them. The refined complex K̂ has the following cells:

- **vertices:** one for each presence (b, o);
- **edges:** the worldline segments, and at each occasion a **rung** for each live tie between two beings present, from (a, o) to (b, o), carrying the tie's value there (attested or latent);
- **2-cells:** the lenses, now quadrilaterals a[o, o′] · rung(o′) · b[o, o′]⁻¹ · rung(o)⁻¹, and the triangles of rungs at an occasion.

Collapsing the rungs of each occasion gives back the carrier K of §2.

K̂ is the lattice of Hamiltonian gauge theory (Kogut and Susskind 1975) laid along a history:

| Hamiltonian lattice | K̂ |
|---|---|
| temporal links | worldline segments |
| spatial links | ties at an occasion |
| temporal plaquettes | lenses |
| spatial plaquettes | triangles of ties at an occasion |

## 7. Measurement: Values in Groups

### 7.1 What the Gauge Symmetry Claims

Temporal links form a forest, since each worldline is a path. So every configuration is gauge-equivalent to one in which all temporal links are the identity: the temporal gauge, which leaves one frame per being free. In this gauge:

- a lens's holonomy is the tie's change between the two meetings, u(o′) u(o)⁻¹;
- a triangle at an occasion has the product of its ties around it;
- both are defined up to conjugation by a being's frame.

The record gives values in the frames of its sources: friend or enemy, as the letters say. The frame group is therefore a symmetry of the *weights*, not of the record. The claim that G is the symmetry of a type of tie is the claim that the weights depend on its values only through holonomies (changes and cycle products). A term that prefers one value over another, such as a preference for friendship over enmity, breaks it, and its size is measured. A being whose ties all flip at once has changed its own frame: that is a gauge transformation from that occasion on, and it creates no curvature. Ties that flip singly do create curvature.

### 7.2 Groups Already in the Literature

- **ℤ₂: signed ties** (friendship and enmity, alliance and hostility, agreement and dispute). Structural balance (Heider 1946; Cartwright and Harary 1956) is ℤ₂ gauge theory on the graph of the moment, term for term:

  | Balance theory | ℤ₂ gauge theory |
  |---|---|
  | a balanced network | a flat connection |
  | Harary's two camps, s_ab = σ_a σ_b | its pure-gauge form |
  | switching in signed graph theory (Zaslavsky 1982) | a gauge transformation |
  | frustration (Toulouse 1977) | curvature |
  | the energy landscape of Marvel, Strogatz and Kleinberg (2009), the number of frustrated triads | the Wilson action of ℤ₂ on spatial plaquettes |

  The dynamics of Antal, Krapivsky and Redner (2005) takes the alliances of Europe before the First World War as its example. Italy's change of sides in 1915 is, in this reading, a change of Italy's frame where its ties flip together, and curvature where they do not.
- **D₁₂: music.**
  - Transposition and inversion act simply transitively on the 24 major and minor triads. The neo-Riemannian group generated by P, L and R is the commuting dual D₁₂ (Lewin 1987; Crans, Fiore and Satyendra 2009).
  - The keys and the forms of a motif are torsors, so absolute pitch is a gauge choice. A quotation, a modulation or a row form relates two statements by a group element.
  - Lewin's transformation networks require path consistency, which is a flat connection. Curvature measures where chains of quotation or modulation fail to compose.
  - ℤ₁₂ (the transpositions) is the abelian subgroup. With inversion the group is non-abelian, and the order test becomes musical. ℤ, the octave, is added where register is recorded.
  - This refines the ℤ₁₂ × ℤ of §12 wherever inversion is measured.
- **O(d): concepts.**
  - A corpus's word embedding is defined only up to rotation: the frame redundancy in its purest form.
  - The alignment of two corpora (two authors, two decades, two languages) is an orthogonal Procrustes problem (Schönemann 1966; across decades, Hamilton, Leskovec and Jurafsky 2016; across languages, Conneau et al. 2018), and it is a connection.
  - Around a loop of translations, the composed alignment differs from the identity: the curvature of meaning, measured. The loop is Aristotle from Greek through Arabic into Latin (Michael Scot), and from Greek into Latin directly (James of Venice, William of Moerbeke), the two read side by side in the thirteenth century. Its eigen-angles on a concept's vocabulary are the concept's drift.
  - The group is non-abelian and continuous. The data are noisy, so every holonomy is compared with that of resampled corpora.
- **ℝ or U(1): positions.** Positions on a latent scale are identified only up to translation, a redundancy of ℝ, or of U(1) if the scale is periodic. Examples are ideal points, the text scaling of Wordfish (Slapin and Proksch 2008), and the Goldstein scale of cooperation and conflict.
- **No group.** For types without values, the weights are estimated as SSF §7.2 has them.

### 7.3 What the Labels Then Mean

- **Spin.** The character expansion of a lens's weight is in the irreps of G, and the spin of a lens is the Fourier dual of the tie's change across it: its **resolution**.
  - Spin 0 alone means the tie's value at one meeting says nothing about the next. High spins mean the tie is held sharply.
  - For ℤ₂, the excited lenses form closed surfaces (mod 2), weighted by tanh β to their area.
  - For U(1), the mode of a lens is a winding n, weighted by I_n(β)/I₀(β).
- **Intertwiner.** The intertwiner on a being-segment is the being's way of holding its live ties together: the Gauss law on the charges of its lenses.
- **Vertex amplitude.** The vertex amplitude of an occasion is the partition function of its rungs: the evaluation of the occasion's spin network, with its triangle weights. It measures how far the ties, as they stand at the occasion, agree.
- **Boundary state.** A boundary state at a date is the spin network of the moment. It can be written in the holonomy basis (the values as last attested) or in the spin basis (the modes); the two are related by Peter–Weyl.

## 8. From the Record to the Formalism: What Is Computed and Tested

Each test is run against a baseline, and the formalism earns its place only where it does better.

1. **Persistence.** In two dimensions, the correlation of a tie's value at two meetings decays as the product, over the lenses between them, of the normalised character coefficients c_j(β)/(d_j c_0(β)) (Migdal 1975; Witten 1991; physis `gauge/exact.hpp`). β per type is fitted from the record: how long a friendship, an alliance or a reading of a motif holds its value. Baseline: a Markov chain per tie.
2. **Frustration and locality.** The frustration of the moment is compared with its distribution under resampled signs; we measure how it falls after flips (the balance dynamics) and the weight of the breaking term of §7.1. Locality is also testable: triangles of rungs exist only where three beings meet. Whether frustration among three beings resolves at their meetings, or without them, decides whether the model may stay local.
3. **The order test** (Street §10.1): whether the invariants depend on the order of a being's occasions. It decides abelian against non-abelian, and it is a formal test of path dependence in the sense of historical sociology (Abbott 2001; Mahoney 2000; Sewell 2005).
4. **Resolution.** Forgetting a meeting merges two lenses, and in two dimensions the merged weight is exactly the convolution: the normalised coefficients multiply. So β fitted on a fine record predicts β on a coarse one with no free parameter. A failure means the history has a memory that the model lacks.
5. **The generating theory.**
   - A group field theory's expansion weights a complex by ∏_b λ_b^{n_b} · N^F · A / sym. That is an exponential family over complexes: an exponential random graph model (Frank and Strauss 1986; Robins et al. 2007) on 2-complexes. Its sufficient statistics are the counts n_b of occasion kinds (the isomorphism classes of the links of the vertices: the bubbles) and the count F of faces.
   - Fitting it places a record relative to the critical coupling.
   - With one colour per discipline it is a coloured theory. Its mixed vertices measure the coupling between disciplines, and the disciplinary histories are its single-colour sectors.
   - This is a hypothesis, to be compared with ERGMs and stochastic block models on the same record.
6. **Topology and dimension:** χ, H₁, H₂ (closed shells of ties) and π₁; the spectral and Hausdorff dimensions of the dual graph, compared with branched polymers (d_H = 2, d_s = 4/3) and the Brownian map (4, 2). These are exploratory, since a record is not homogeneous.
7. **Prediction of what was held out:** presences, ties and values removed from the record and then predicted. Baselines: temporal-network link prediction, Hawkes processes and stochastic block models.
8. **Power** (SSF Def. 10.3) and acts: these are properties of the fitted model, and are stated as such, never as claims of what history would have done.

## 9. Where the Records Come From

| Field | Occasions with participants | Ties and values | Licence |
|---|---|---|---|
| all | Wikidata: a scaffold of identities, dates, memberships, education, "student of". Its "influenced by" is unsourced, and is not a record | — | CC0 |
| correspondence (every field) | correspSearch (CMIF, from scholarly editions) and Early Modern Letters Online: a letter is an occasion of two, dated and edited | the stance of the letter, coded | per edition |
| science | OpenAlex: publications as occasions of their authors | citations as ties between works | CC0 |
| music | the performance histories of Carnegie Hall and the New York Philharmonic, and MusicBrainz: programmes as occasions | quotations, transpositions and row forms from public-domain scores (D₁₂) | CC0 |
| art | MoMA's exhibition history; the Database of Modern Exhibitions, 1905–1915 (University of Vienna): catalogues as occasions | teacher of, member of, exhibited with | CC0; to check |
| philosophy | edited correspondence and records of meetings (Waismann's of Wittgenstein and the Vienna Circle) | disputes (ℤ₂), concepts (O(d)) | per edition |
| international relations | Correlates of War and ATOP (alliances, disputes, diplomatic exchange); GDELT and ICEWS for recent decades | signed (ℤ₂) and valued (Goldstein) | terms to check |

Every licence is checked before a source is adopted (CLAUDE.md §7). The pipeline has three routes:

- A small history is coded by hand from a scholarly edition.
- Adapters import the structured sources.
- Assertions that a language model extracts from digitised text are *proposals*, each with its quoted passage, and a person accepts or rejects each one. The record keeps who proposed and who accepted.

## 10. The Difficulties

1. **Selection.** Records over-represent the literate, the powerful, the West and the winners. The complex is a model of the record, not of the past. So each history states its sampling frame. It prefers census-like sources (rolls, programmes, catalogues, complete editions of correspondence) to encyclopaedic ones. It treats a missing presence as unobserved, not as absent, and it compares only within a frame.
2. **Resolution.** There is no natural size of an event; the resolution is a parameter, and the views are related (§6, §8.4).
3. **The identity of themes.** The coder does not decide whether the ether of 1880 is the ether of 1905. The re-description is recorded as a transport, and the curvature measures where identity fails to compose: the ship of Theseus becomes a holonomy.
4. **The coder's circularity.** The codebook shapes the complex. Two coders code a sample, their agreement is reported, and the ensemble is taken over the codebooks where they disagree.
5. **Underdetermination.** Many models fit one record. The tests of §8 are stated before they are run, and each has a baseline.
6. **Reading.** Amplitudes are the weights of a statistical model, not forces of history, and nothing in the model is extremised (§3).

## 11. A First History

The first history should test the method, not the ambition. It needs a bounded population with census-like sources, dense ties with values, a size small enough to code by hand, and enough study by historians that the construction can be checked against them.

1. **The calculus priority dispute, 1669–1726.**
   - About twenty beings: Newton, Leibniz, Oldenburg, Collins, Wallis, Fatio, Keill, Johann Bernoulli, Conti, Caroline of Wales, the Royal Society, the *Acta Eruditorum*, the *Commercium Epistolicum*.
   - Some hundred and fifty occasions, chiefly letters, in the editions of Newton's and Leibniz's correspondence and in EMLO.
   - Signed ties (ℤ₂) whose signs change, and triangles that frustrate and resolve (Leibniz, the Royal Society, and its president, Newton).
   - Two calendars.
   - It is studied by Hall (1980) and Bertoloni Meli (1993), and every part of the pipeline would be exercised in about a week of coding.
2. **The Vienna Circle and its neighbours, 1924–1938.** Several disciplines in one record:
   - the protocol-sentence dispute (Neurath, Carnap, Schlick, Popper);
   - the conversations with Wittgenstein, recorded by Waismann;
   - Carnap's lectures at the Bauhaus in 1929, and Neurath's Isotype;
   - the political end in 1934–38, and the emigration, where the translation of the Circle's concepts into English is a loop for O(d) (§7.2).
3. **The Second Viennese School and Schoenberg's Society for Private Musical Performances, 1904–1921:** lessons, premieres and the Society's concerts, whose programmes survive, with D₁₂ taken from the scores.

The proposal is the first as the pilot of the method, and the second as the application's first history.

### 11.1 The Canon First

The whole record of a history is beyond any project; a nation's libraries hold it. The formalism does not need it, as statistical mechanics does not need every molecule. The classical facts of a history, its **canon**, are a record at a coarse resolution (§6, item 4). On it the construction runs, and denser sources later refine a region of it without rebuilding the rest. What is built from a canon is the spacetime of *remembered* history, and the canon's own selection (the famous, the diplomatic, the written) is visible in its structure and has to be stated. Truth is not decided: classical facts are the layer of assertions on which the standard references agree, and contested ones become readings (§5.1).

### 11.2 A Test on Classical Facts (24 October 2026)

A canon of Europe 1785–1828 was drafted from general knowledge by a language model. It has 28 beings, 78 occasions, 208 presences and 45 ties across music, letters and philosophy, science, art and international relations, with the friend/enemy ties signed. Every entry names its standard reference and has a confidence. The construction of §6 was run on it, and the complex was given to the Complex Atlas's engine unchanged. The prototype is in the session's scratchpad, not in the repository.

- **The complex.** V = 78, E = 180, F = 129, χ = 27; it is connected and does not collapse. H₁ = ℤ⁶, π₁ is free of rank 6, and H₂ = ℤ³².
- **H₂** is generated by closed shells: the ribbons of three (sometimes four) mutually tied beings, closed off at two occasions that all of them attend. 14 of the 32 basis shells are the minimal ones of three lenses. On each shell the Bianchi identity binds the change of the triangle's frustration between the two occasions to the flips of the shell's lenses.
- **H₁** is generated by circles of encounters that no recorded tie fills. Examples:
  - Haydn, Beethoven and Lichnowsky around Mozart's death, with no tie between Haydn and Lichnowsky recorded;
  - Goethe, Humboldt, Volta and Napoleon, the circle of galvanism;
  - Napoleon, Goethe and Hegel around Jena.

  Each is a silence of the canon, so the topology gives a list of what to look for next.
- **ℤ₂ in temporal gauge.** 26 of the 129 lenses carry a flip. Triangles are frustrated at Tilsit (1807: Prussia left by Russia), at Erfurt (1808: Talleyrand's double game) and at the opening of the Congress (1 November 1814: the Saxon–Polish crisis). The secret treaty of 3 January 1815 resolves the Congress triangles.

  The frustration of the network of the moment by year:

  | Years | Frustrated triangles | Period |
  |---|---|---|
  | to 1806 | 0 | |
  | 1807–1811 | 0.50–0.83 | the Tilsit system |
  | 1812 | 0.29 | |
  | 1813 | 0.09 | the coalition |
  | 1814 | 0.35 | the crisis at the Congress |
  | from 1815 | 0.06 | |

  The signs are the drafter's coding, so this shows that the formalism computes and reads, not that balance holds. That needs independent coders and resampled baselines (§8).
- **Checking against structured sources.** Wikidata resolved 25 of the 28 beings automatically. The three it missed (Archduke Rudolph, Lobkowitz, Kinsky) are the least famous, which is where a model's knowledge also thins. All seven deaths in the canon agree with Wikidata to the day.
- **What the test showed is missing:**
  1. Values carried forward go stale. The frustrations of 13 March 1815 come from signs carried since January, though the Saxon question was settled in February. Unattested values must be inferred through the persistence of §8.1, not carried.
  2. Open ties live for ever. The one frustrated triangle left after 1815 is David–Charles IV, a commission of 1800 never attested again. Ties need ends or decay.
  3. The canon is thin and lopsided. 89 of the 129 lenses are diplomacy, 4 science and 3 art. Five ties have a single meeting and make no lens. Napoleon lies on four of the six loops: the great-man structure of a canon is visible in its topology.
  4. Resolution appears at once. The congress of Erfurt lasts weeks and contains Goethe's audience.

## 12. What This Note Changes in the Papers

1. **The carrier** (§2): facts on faces with couplings as edges is a hypergraph wearing a 2-complex's clothes. The cyclic order must be time, so occasions are vertices, beings edges, ties faces, and themes the types of ties. This changes the social spin foam's Def. 6.1 and the street paper's orientation II, and reads the topology paper's architecture rather than changing it.
2. **The symmetry** (§3, §7): GRH's constitutional limit is weaker than stated. The frame group is independently motivated by the redundancy of description that the street paper itself names (no view from nowhere), and for signed ties it has been in the social sciences since Cartwright and Harary as switching and frustration. The invariants GRH asks for are the gauge-invariant and topological ones. What remains open is which group each type of tie selects, an empirical question with a stated test.
3. **The weights** (§3): SSF's stipulated weights and Street's derived ones are the same object. The derivation applies wherever a group is motivated, and the estimated or stipulated weights of SSF §7.2 are the general case where it is not.
4. **The action** (§3): GRH's exclusion of an action principle is right about teleology, and it does not exclude Gibbs weights, which the framework already uses; the Regge phase of the large-spin limit is such a weight.
5. **The music mapping:** pitch class as a spin of SU(2) (the twelve lowest irreps) is a labelling, not a symmetry. The triangle rule then imposes constraints on chords that no musical symmetry motivates. The honest group of pitch class is ℤ₁₂ (with ℤ for the octave), abelian, under which harmonic closure is a Gauss law; or D₁₂ when inversion is measured (§7.2). Which one the music history carries is a case of the order test.

## 13. What the Semantics Asks of the Foundations and the Application

- **physis `gauge`:** lattice gauge theory on K̂, with the attested links fixed and the latent ones sampled; the temporal gauge; terms that break the symmetry; the exact persistence and convolution of §8.1 and §8.4; finite groups from morphe tables (ℤ₂, D₁₂; `FiniteGroupAdapter` exists), and O(d).
- **physis `symmetry`, `quantum_gravity`:** evaluation of closed spin networks of any valence (the {3nj} of an occasion with any number of beings), the channel count as a positive vertex weight, spin foams on arbitrary 2-complexes with open faces, and the order test.
- **morphe `topology`:** complexes with types, kinds, disciplines and confidences; the quotient K̂ → K; derivations as filtrations, with their invariants stage by stage; `categories` begun for the composition of histories (states as objects, histories as morphisms, periodisation as composition).
- **The inference foundation (planned):** readings and ensembles; exponential families over complexes; held-out prediction and its baselines.
- **A domain for history** (a proposal, to be agreed under CLAUDE.md §6.5): the record format, the codebooks, the civil calendars (Julian, Gregorian, Old Style), the construction and the readings. It would depend on metrology (time) and geography (places).
- **The application, the History Complex:**
  - the record, with every assertion beside its source;
  - the loom: worldlines running upward in time, ties as ribbons of lenses between them, occasions as knots;
  - the complex in space;
  - the spin network of the moment;
  - the tests of §8, with their baselines;
  - disciplines as colours.

## 14. Open Choices

1. The pilot (§11).
2. Who codes: by hand (the user, collaborators, Claude); by a model's proposals accepted by a person; or both.
3. The placement of the record and the construction: a domain for history (§13).
4. The first measured group: ℤ₂ for the pilot, then D₁₂ or O(d).

## 15. References

- **The user's papers:** Social Spacetime Topology and Generative Grammars in CW Complexes; A Social Spin Foam Model for Power Analysis and Governance Decision-Making; Street Spin Foam: A Preliminary Formulation; A Theory of Historical Generativity under Constraints; the GFT music spin foam (2026, Serendip).
- **Spin foams and group field theory:** Reisenberger and Rovelli (1997); Baez (1998), Spin foam models; Ponzano and Regge (1968); Turaev and Viro (1992); Bianchi, Doná and Speziale (2011), Polyhedra in loop quantum gravity; Oriti (2009, 2017) on group field theory and its condensates; Etingof, Gelaki, Nikshych and Ostrik (2015), *Tensor Categories*.
- **Gauge theory:** Kogut and Susskind (1975), Phys. Rev. D 11, 395; Migdal (1975); Witten (1991), Commun. Math. Phys. 141, 153; Toulouse (1977), Commun. Phys. 2, 115.
- **Balance and signed graphs:** Heider (1946), J. Psychology 21, 107; Cartwright and Harary (1956), Psychol. Rev. 63, 277; Zaslavsky (1982), Discrete Appl. Math. 4, 47; Antal, Krapivsky and Redner (2005), Phys. Rev. E 72, 036121; Marvel, Strogatz and Kleinberg (2009), Phys. Rev. Lett. 103, 198701.
- **Social networks:** Breiger (1974), Social Forces 53, 181; Frank and Strauss (1986), JASA 81, 832; Robins, Pattison, Kalish and Lusher (2007), Social Networks 29, 173.
- **Historical sociology:** Abbott (2001), *Time Matters*; Mahoney (2000), Theory and Society 29, 507; Sewell (2005), *Logics of History*.
- **Music theory:** Lewin (1987), *Generalized Musical Intervals and Transformations*; Crans, Fiore and Satyendra (2009), Amer. Math. Monthly 116, 479.
- **Concepts, texts and scales:** Schönemann (1966), Psychometrika 31, 1; Hamilton, Leskovec and Jurafsky (2016), ACL; Conneau et al. (2018), ICLR; Slapin and Proksch (2008), AJPS 52, 705.
- **Records and standards:** Bradley and Short (2005), Literary and Linguistic Computing 20 suppl., 3; CIDOC CRM (ISO 21127); EDTF (ISO 8601-2:2019); PeriodO; Krippendorff (2018), *Content Analysis*.
- **The pilot:** Hall (1980), *Philosophers at War*; Bertoloni Meli (1993), *Equivalence and Priority*.
