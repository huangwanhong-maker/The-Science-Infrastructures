# The History Complex

*Status: draft. Serves learning and research.* The construction of a 2-complex from the record of a history, and the reading of its signed ties as a ℤ₂ gauge field. The code is `construction/` and `readings/`; the conventions of the record are the [codebook](../codebook.md); the decision is [ADR 0041](../../../development/decisions/0041-the-history-domain.md); the semantics it begins to realise is the idea note [The History Complex](../../../development/ideas/2026-10-23-history-complex.md).

## 1. Records

**Definition 1.1 (record).** A record is a tuple R = (B, O, ≤, P, T) of:

- a finite set B of beings and a finite set O of occasions;
- a total order ≤ on O, the order of time: by the first day an occasion's date may denote, then its last, then the file (codebook §4);
- a presence relation P ⊆ B × O, with O_b = {o : (b, o) ∈ P} the occasions of b;
- a finite set T of ties.

A tie t ∈ T has two distinct beings a(t), b(t), a colour c(t), and attestations: a non-empty set A_t ⊆ O_{a(t)} ∩ O_{b(t)} of joint occasions, with a sign σ_t : A_t → {±1}. It may also have an end e(t) ≥ max A_t.

**Definition 1.2 (life, meetings).** The *life* of t is the interval [min A_t, ℓ_t] of O, where ℓ_t = max(max A_t, e(t)). Its *meetings* are M_t = O_{a(t)} ∩ O_{b(t)} ∩ [min A_t, ℓ_t], written in order as o₀ < o₁ < … < o_{m_t}.

*Remark.* A tie lives from its first attestation to its last, or to its end; it is never extended beyond what the record says (codebook §8). When t ends at the death of one of its beings, the codebook's derivation makes the other present there (role `bereaved`), so that e(t) ∈ M_t.

**Definition 1.3 (join).** Let R₁, …, R_n be records. An identifier is a name within one record; two records *share* what one identifier names when they agree about it:

- a being, when both give the same Wikidata item, or, when either gives none, the same name;
- an occasion, when both give it the same date;
- a tie, when it joins the same two shared beings.

The *join* R₁ ∨ … ∨ R_n is the union of the records with what they share identified: the union of the beings, of the occasions and of the ties; the order of time on the union; the union of the presences; and, for a shared tie, the union of its attestations. An identifier that two records use for different things names two things in the join: the later record's is renamed `id@canon` throughout that record, and the join reports it (`history::join_records`). Alexander the Great (Q8409) and Alexander I of Russia (Q15193), both `alexander` in their canons, stay two beings.

*Rationale.* This is the authority file of prosopography: identity across sources is a claim, settled by an authority (here Wikidata) or by agreement, and never by a coincidence of labels.

*Remark.* The join is the colimit of the records over the parts they share: each record maps into it, and it is glued only where they meet. Its complex is *not* the union of their complexes. A shared being's life runs through the occasions of every record that names it, and is cut by all of them, so its segments are re-cut; a tie attested in two records lives from the first attestation in either to the last. The join of canons of different periods is disconnected wherever no being crosses between them, and H₀ counts the separate histories (§2). The *integrated history* of the infrastructure is the join of all its canons; it grows as they come to share beings.

## 2. The Construction

**Definition 2.1 (the history complex).** The history complex K(R) is the CW 2-complex with these cells:

- **vertices:** the occasions O;
- **edges (the worldlines):** for each being b with O_b = {p₀ < p₁ < … < p_{n_b}}, an edge e_{b,i} from p_{i−1} to p_i, for i = 1, …, n_b; write b[p, q] for the path of b's edges from p to q;
- **2-cells (the lenses):** for each tie t and each i = 1, …, m_t, a disc f_{t,i} attached along the word

  ∂f_{t,i} = a(t)[o_{i−1}, o_i] · b(t)[o_{i−1}, o_i]⁻¹.

**Proposition 2.2.** K(R) is a CW complex, and

χ(K) = |O| − Σ_b (n_b) + Σ_t (m_t),

where n_b = |O_b| − 1 for a being present somewhere, and m_t = |M_t| − 1.

*Proof.* The word of f_{t,i} is a closed path: a(t) and b(t) are both present at o_{i−1} and at o_i, so both paths run from o_{i−1} to o_i, and the word returns to its base. The counts are those of the cells. ∎

*Remark (what the cyclic order means).* A 2-cell's attaching word is a cyclic sequence, and here its order is time: forward along one life between two meetings, back along the other. This is why the record's faces are ties between two persisting beings, and not facts on their own. A chain of custody is a path, and closes into a disc only when it is held between two carriers.

**Definition 2.3 (colours).** The sub-complex K_c of a colour c is the union of the lenses of the ties of colour c, with their boundaries. The history of a discipline is K_c, and the integrated history is K. An occasion where lenses of two colours meet is interdisciplinary by construction.

**Definition 2.4 (the refined complex).** K̂(R) has these cells:

- **vertices:** the presences (b, o) ∈ P;
- **temporal edges:** the worldline segments, from (b, p_{i−1}) to (b, p_i);
- **spatial edges (the rungs):** for each tie t and each meeting o ∈ M_t, a rung r_{t,o} from (a(t), o) to (b(t), o);
- **2-cells:** the lenses, now the squares

  a(t)[o_{i−1}, o_i] · r_{t,o_i} · b(t)[o_{i−1}, o_i]⁻¹ · r_{t,o_{i−1}}⁻¹,

  and, at each occasion, the triangles of rungs among three beings pairwise tied there.

Collapsing the rungs of each occasion to a point gives K back (with the triangles collapsed). K̂ is the lattice of Hamiltonian gauge theory (Kogut and Susskind 1975) laid along a history: worldlines are its temporal links, ties at an occasion its spatial links, lenses its temporal plaquettes and triangles its spatial plaquettes.

## 3. The Signed Ties as a ℤ₂ Field

Let a configuration assign g ∈ ℤ₂ = {±1} to every edge of K̂: to a rung, the tie's sign at that meeting; to a temporal edge, the change of its being's own frame between two occasions. A gauge transformation is a function φ: P → ℤ₂, acting by g_{(x→y)} ↦ φ(x) g φ(y).

**Proposition 3.1 (temporal gauge).** Every configuration is gauge equivalent to one whose temporal edges are all +1. The transformations that preserve this are those constant along each worldline: a sign φ_b for each being.

*Proof.* The temporal edges form a forest, one path for each being, so φ can be fixed along each path from its first presence. The residual freedom is φ's value at each path's start. ∎

The residual transformation φ_b = −1 reverses every tie of b: this is the *switching* of signed graph theory (Zaslavsky 1982). In temporal gauge:

- the holonomy of a lens is h(f_{t,i}) = σ_t(o_i) · σ_t(o_{i−1}), the tie's change between two meetings; the lens carries a **flip** when h = −1;
- a triangle (x, y, z) of ties at an occasion o has holonomy σ_{xy}(o) σ_{yz}(o) σ_{xz}(o); it is **frustrated** when this is −1 (Toulouse 1977; Cartwright and Harary 1956).

Signs at meetings where a tie is not attested are carried from its last attestation (Definition 3.5 weighs them).

**Proposition 3.2 (invariance).** Flips, frustrations and the products on shells (Proposition 3.3) are unchanged by switching any being.

*Proof.* Switching b reverses σ_t for every tie t at b. A lens of such a tie has both its signs reversed, so its holonomy is unchanged. A triangle containing b has exactly two of its ties at b. ∎

**Proposition 3.3 (closed shells and the Bianchi identity).** Let x, y, z be pairwise tied by ties whose consecutive meetings include the same pair p < q. Then the three lenses between p and q are a 2-cycle of K:

∂(f_{xy} + f_{yz} − f_{xz}) = 0.

Their holonomies satisfy

P(q) = P(p) · h(f_{xy}) h(f_{yz}) h(f_{xz}),

where P(o) is the triangle's product at o; so the triangle's frustration changes exactly when an odd number of its ties flip.

*Proof.* With the orientations of Definition 2.1, ∂f_{xy} = x[p, q] − y[p, q], and likewise for the others; the alternating sum vanishes. The identity is the product of the three holonomies σ(q)σ(p), each sign at p and at q appearing once. ∎

*Remark.* Every minimal shell is a class of H₂(K; ℤ); the classes of H₂ are combinations of such shells along chains of meetings. The configuration space of an abelian spin sector on K is H₂ (the idea note, §2), and this is where it lives in a record: triangles of ties carried through time.

**Definition 3.4 (persistence).** Let a lens carry the weight e^{β σ(o_{i−1}) σ(o_i)}. The character expansion

e^{βs} = cosh β (1 + s tanh β)

shows that along one chain of lenses the correlation of a tie's signs across a lens is tanh β. Let p̂ be the fraction of lenses whose two meetings are both attested that flip. Then the **memory** m̂ = 1 − 2p̂ estimates tanh β, and β̂ = artanh m̂ when 0 < m̂ < 1. It is estimated for each colour, and for all colours together.

**Definition 3.5 (carried signs).** A sign carried k lenses from its last attestation σ is expected to be σ m̂^k, the memory of the tie's colour clamped to [0, 1] (a Markov chain of flips).

**Definition 3.6 (the moment).** The network of the moment at a day d is made of:

- the beings and the ties alive at d, each at its last sign σ and its expected sign σ m̂^k;
- its triangles, with their products;
- their expected frustration Σ (1 − E[product])/2, the product's expectation taken as the product of the three expectations.

*Caveat.* The expected frustration treats the three ties of a triangle as independent, and the carried signs as a Markov chain. Both are the simplest assumptions, to be tested against the record (the idea note, §8). The estimates are only as good as the coding of the signs: a second coder and Krippendorff's α come before any claim (codebook §9).

## 4. Loops and Silence

The generators of H₁(K; ℤ) are circles of encounters that no recorded tie fills. Their beings meet in a cycle, and the record holds no relation spanning it. The engine tells each generator as closed walks along the worldlines: the beings it follows, and the occasions where it passes from one to another. Two readings follow:

- **for the record,** a loop marks a silence, and so a place to look for a tie;
- **for the field,** a loop is where a flat connection may wind: Hom(π₁(K), ℤ₂) = H¹(K; ℤ₂) counts the ways frames can twist round the history's loops with no tie to measure it.

*Caveat.* A canon of classical facts is thin and centred on its famous beings, and its loops gather round them. A loop is evidence about the record before it is evidence about the past.

## 5. Why Signs Live on the Refined Complex

The history complex K collapses each occasion to a point. A connection on K puts a group element on each being-segment, and the holonomy of a lens is the mismatch of its two beings' transports between their meetings. Read the signed record on K, and the holonomies of the three lenses of a closed shell multiply to the identity for every connection, because the shell is closed and ℤ₂ is abelian.

The record does not obey this. At the opening of the Congress of Vienna the triangle Alexander–Frederick William–Metternich is frustrated; at the secret treaty of 3 January 1815 it is not. The frustration changed while one tie of the three flipped: an odd number of flips on a closed shell.

**Proposition 5.1.** On K, a frustrated triangle at an occasion o is curvature concentrated at the vertex o: no connection on K has lens holonomies equal to the record's flips when a triangle changes frustration between two meetings. On K̂, where the triangle is a 2-cell, the record's signs are exactly a configuration, and the Bianchi identity of Proposition 3.3 holds with the triangle's own holonomy.

*Proof.* On K the three lens holonomies of a minimal shell multiply to 1 for any connection (the shell is a 2-cycle). The record's flips multiply to P(p)P(q) by Proposition 3.3, which is −1 when the frustration changes. On K̂ the rungs are edges and the triangles faces, and in temporal gauge the record's signs are the rung variables. ∎

So K is the carrier of the topology of encounters and of the themes without signs, and **K̂ is the carrier of the signed record and of its spin foam**. This is the precise form of the idea note's choice of occasions as vertices: an occasion is a vertex at the coarse resolution, and its link, the triangles of ties among those present, is where its internal agreement or disagreement lives.

## 6. The Spin Foam of a History

This section realises, for the group the record measures, the unified model of the idea note (§3). That model is the social spin foam's carrier and weights with the street spin foam's frame group, its weights derived by the character expansion rather than stipulated. With signed data the measured group is ℤ₂, and everything below is exact for it. For other groups the same derivation holds with other characters (§6.5).

### 6.1 The Weight

The variables are the rungs of K̂ in temporal gauge, u_r = ±1. The **Gibbs weight** of a configuration is a product over the 2-cells of K̂:

W(u) = Π_lenses e^{β_c(t) u_{t,o} u_{t,o′}} · Π_triangles e^{β_Δ u u′ u″}.

- β_c is the **persistence** of the ties of colour c: how strongly a tie holds its sign from one meeting to the next.
- β_Δ is **balance**: how strongly three ties at an occasion hold their product to +1.

W depends on u only through holonomies, so it is a class function, unchanged by switching (Proposition 3.2): the frame group is a symmetry of the weights and not of the record, as the idea note asks (§7.1). Nothing is extremised: W is a measure on configurations, not a functional on trajectories, so the exclusion of action principles by generative relational history is respected (the idea note, §3).

**Estimation.** Each coupling is estimated from the cells whose signs are all attested, with the Krichevsky–Trofimov estimate p = (k + ½)/(n + 1) of the flip or frustration rate, and tanh β = 1 − 2p. These are moment estimates, each made as if its cells were independent. They are finite where nothing has flipped.

### 6.2 The Character Expansion and the Meaning of Its Labels

Write t_f = tanh β_f and expand each factor as e^{β h} = cosh β (1 + h t). Summing over configurations, a product of holonomies survives only if every rung occurs in it an even number of times; the temporal edges being a forest, this holds exactly for the 2-cycles of K̂ over ℤ/2. Hence

Z = Π_f cosh β_f · 2^{#rungs} · Σ_{S ∈ Z₂(K̂; ℤ/2)} Π_{f ∈ S} t_f.

The sum is over **spin foams**: assignments of an irrep j_f ∈ {0, 1} of ℤ₂ to each 2-cell, admissible at every edge. Cell by cell, the labels mean the following.

| Label | Value at G = ℤ₂ | Meaning in a history |
|---|---|---|
| spin j_f | 0 or 1: the cell is in the foam or not | whether a tie's change across a lens (or a triangle's balance) takes part in a closed configuration of changes, one constrained by the others around it |
| face amplitude c₁(β)/c₀(β) | tanh β_f | how strongly a tie carries its sign across one interval between meetings (persistence), or a triangle its balance at one occasion |
| intertwiner at an edge | dim Inv(⊗ V_{j_f}) = 1 if Σ j_f is even, else 0 | **the Gauss law**: a rung or a life-segment lies on an even number of the foam's cells. A being cannot end a configuration of changes: changes close up through the beings that carry them. This is the social spin foam's admissibility (its Prop. 6.4), here a theorem |
| vertex amplitude | 1 on admissible configurations | ℤ₂ is abelian: the agreement of the accounts at an occasion is already in its triangles' amplitudes, since the triangles are the faces of the occasion's link |
| amplitude of a foam | Π_{f ∈ S} tanh β_f | how strongly the history carries that closed configuration of changes |
| the space of foams | Z₂(K̂; ℤ/2) = H₂(K̂; ℤ/2), of dimension b₂ | the independent ways in which changes of ties can close up: in the Europe canon, 66 |

**Proposition 6.1 (spins from configurations).** In the prior,

P(f ∈ S) = t_f (⟨h_f⟩ − t_f) / (1 − t_f²),

where ⟨h_f⟩ is the mean holonomy of f. A cell outside the support of every 2-cycle has ⟨h_f⟩ = t_f and never belongs to a foam.

*Proof.* Differentiate log Z with respect to β_f: ⟨h_f⟩ = t_f + ((1 − t_f²)/t_f) · P(f ∈ S). If no cycle contains f, the sum over foams does not involve t_f. ∎

**Corollary 6.2.** If every coupling is positive, then 0 ≤ P(f ∈ S) ≤ t_f/(1 + t_f) for every cell f.

*Proof.* The GKS inequality gives ⟨h_f⟩ ≥ t_f, and |h_f| = 1 gives ⟨h_f⟩ ≤ 1; P(f ∈ S) increases with ⟨h_f⟩, and equals t_f/(1 + t_f) at ⟨h_f⟩ = 1. ∎

So no cell lies in more than half of the foams' weight. The bound is reached by a cell whose holonomy is all but certain, such as a triangle among many beings at one well-attested occasion, which every configuration of the theory keeps balanced. In the Europe canon, with tanh β_Δ = 0.83, the triangles of the great congresses sit at the bound, 0.45.

The engine finds the support of the foams by Gaussian elimination of ∂₂ of K̂ over GF(2). It estimates ⟨h_f⟩ for the cells inside it by heat bath, with Rao–Blackwell estimators, and checks the identity against foams enumerated one by one on small records.

*Remark (mixing).* Single-site moves do not mix where balance is strong: on the alliances of 1871–1914, with tanh β_Δ = 0.95, the triangles of the great congresses freeze. With every coupling positive, the GKS inequality requires ⟨h_f⟩ ≥ tanh β_f, so a mean below it betrays the sampler. The engine therefore adds a cluster move at each presence (b, o), which flips together all the rungs of b's ties at o. Every triangle through the presence holds two of these rungs, so no triangle changes, only the lenses that begin or end there. The move is accepted by Metropolis, and the means then satisfy the inequality on every canon.

*Remark (signed foams).* The configuration side is a probability measure for every real β. The foam side is signed when some β < 0: a colour whose ties flip more often than not, or triangles frustrated more often than not, make tanh β negative and some foams' amplitudes negative. This is the social spin foam's signed weight (its Rem. 7.5), and the record decides when it occurs.

### 6.3 Wilson Loops: Coherence Round a Circle of Encounters

For a closed path C in K̂, the Wilson loop W(C) = Π_{r ∈ C} u_r is the product of the signs met round it, and

⟨W(C)⟩ = Σ_{S : ∂S = C} Π_{f ∈ S} t_f / Σ_{S : ∂S = 0} Π_{f ∈ S} t_f.

- **On the record's silences:** if C is not a boundary over ℤ/2 (a class of H₁(K̂; ℤ/2)), no S has ∂S = C and **⟨W(C)⟩ = 0 exactly**. Frames twist freely round a circle of encounters that no tie measures.
- **Area law:** if C bounds, the leading terms are the spanning foams of least area, ⟨W(C)⟩ ≈ N_min Π t_f over a minimal spanning surface. Coherence round a circle decays with the number of ties' intervals needed to close it.
- **A tie's life:** the boundary of the ribbon of a tie from its first meeting to its last gives ⟨u_first u_last⟩, the memory of a relation across its life. For an isolated tie it is t_c^k over k lenses. Every closed shell the ribbon belongs to adds spanning foams, so in a dense history **balance carries memory**: in the Europe canon, Alexander I and Frederick William III keep their sign across ten lenses with ⟨W⟩ ≈ 0.26, where an isolated chain would keep 0.003, and Metternich and Alexander keep 0.84 across six lenses, against 0.03.

### 6.4 The Posterior: Reconstruction, and Power

**The reconstruction ensemble** (the social spin foam's Def. 8.5) is W conditioned on the attested signs. An unattested sign is inferred from all the others at once: its tie's persistence before and after, and the balance of every triangle it closes. This replaces carrying the last sign forward (§3.6), whose fault the feasibility test found: a sign carried through a settlement it did not see. The engine samples the posterior by heat bath and reports:

- each rung's mean, which is the inferred sign and its certainty;
- the frustration of the network of the moment at each year's end, with its mean and its 5% and 95% points.

**Power** (the social spin foam's Def. 10.3) is the divergence D(P_x ‖ P) of the reconstruction under an act x from the reconstruction without it.

- **The power of a being:** the act removes the being's ties, with the couplings kept, and the divergence is taken over the posterior marginals of every other unattested sign. It is approximate where the signs are dependent.
- **The power of balance itself:** the act sets β_Δ = 0.

In the Europe canon, Metternich holds 0.43 nats over the reconstruction, Alexander I 0.40 and Talleyrand 0.02, and balance holds 2.3 nats; in Vienna, the Circle holds 2.1 and Wittgenstein 1.2, and balance 7.1. These are properties of the fitted model, stated as such, never claims about what history would have done.

### 6.5 Beyond ℤ₂: SU(2) Kinematics

For a non-abelian group, the same expansion puts:

- an irrep j_f on each cell, with face amplitude d_j c_j(β);
- an intertwiner in Inv(⊗ V_{j_f}) on each edge;
- at each vertex the evaluation of the spin network of its link.

This is Reisenberger and Rovelli's spin foam on a 2-complex (Baez 1998, §3). The record does not yet measure SU(2) elements, so the dynamics is not fitted. Its **kinematics** is computed: with spin ½ on every live tie, a being holding n ties across a segment of its life has dim Inv(V_½^{⊗n}) channels, the Catalan number C_{n/2} for n even and none for n odd. This is the being's internal freedom, the social spin foam's coupling weight C_e derived (the idea note, §3). A being holding an odd number of ties cannot hold them all at spin ½; it needs a tie of integer spin, a relation of another kind. This is a selection rule a record with SU(2) data could test.

## 7. Other Complexes of a Record

The record makes more than one complex. Each answers its own question, and the engine computes them all.

| Complex | Cells | What its homology says | In the Europe canon |
|---|---|---|---|
| K, the history complex (§2) | occasions; lives between them; lenses of ties | H₁: circles of encounters no tie fills (silences); H₂: closed shells of ties | ℤ, ℤ⁹, ℤ⁴⁶ |
| K̂, the refined complex (§2.4, §5) | presences; lives and rungs; lenses and triangles | H₂ over ℤ/2: the space of spin foams | ℤ, ℤ⁸, ℤ⁶⁶ |
| K³, K with its shells filled | K and a 3-cell on each minimal closed shell | H₂: shells that are not filled episodes; H₃: quartets of mutually tied beings whose four shells close into a 3-sphere | ℤ, ℤ⁹, ℤ¹⁵, ℤ⁶ |
| the Dowker complex on the beings (Dowker 1952) | a simplex for the beings present together at each occasion | Atkin's Q-analysis (1974): the number of q-connected components for each q; by Dowker's theorem, homotopy equivalent to the complex on the occasions with a simplex for each being's occasions (Breiger's duality made topological) | ℤ, ℤ⁹; Q = (1, 18, 14, 5, 1, 3, 3) |
| K filtered by time | each cell at the rank of its last occasion | persistent homology over ℤ/2 (morphe `persistence.hpp`): each circle of encounters is born when a being's arrival closes it, and dies when a later meeting's lens fills it, or lives on | 10 circles closed, 9 never filled |

*Remark.* K and the Dowker complex answer different questions. A circle of encounters is filled in K only by a tie's lens, and in the Dowker complex by any co-presence. Their first Betti numbers need not agree; in this canon both are 9.

## 8. Notation and Code

| Symbol | Meaning | Code |
|---|---|---|
| R = (B, O, ≤, P, T) | a record | `history::Record` (`record.hpp`) |
| ≤ | the order of time | `Record::order` |
| A_t, σ_t, e(t) | a tie's attestations, signs and end | `Tie::attested`, `Tie::ends` |
| M_t | a tie's meetings within its life | `HistoryComplex::rungs[t]` |
| K(R) | the history complex | `HistoryComplex::cw` (morphe's `CWComplex`) |
| e_{b,i}, b[p, q] | a worldline's edge, its path | `Segment`, `detail::path` |
| f_{t,i} | a lens | `Lens` |
| K_c | the sub-complex of a colour | `colours_of` (`colours.hpp`) |
| h(f), P(o) | a lens's holonomy, a triangle's product | `HistoryComplex::flip`, `Triangle::product` |
| p̂, m̂, β̂ | flip rate, memory, coupling | `Persistence::p`, `memory`, `beta` |
| σ m̂^k | a carried sign's expectation | `expected(Rung, m)` |
| K̂ | the refined complex | `RefinedComplex` (`foam/refined.hpp`) |
| β_c, β_Δ | persistence of a colour, balance | `Couplings::lens`, `Couplings::balance` (`foam/model.hpp`) |
| W(u), u_r | the Gibbs weight, a rung's sign | `HistoryModel`, `sample` |
| Z₂(K̂; ℤ/2), j_f | the space of foams, a cell's spin | `FoamSpace`, `spins` |
| ⟨W(C)⟩ | a Wilson loop | `Ensemble::loop` |
| D(P_x ‖ P) | power | `marginal_divergence`; the engine's `power` |
| dim Inv(⊗ V_j) | SU(2) channels | `invariant_dimension` |
| K³, Dowker, Q | the other complexes | `shells_filled`, `dowker_beings`, `q_analysis` (`foam/models.hpp`) |
| [b, d) | an interval of persistence | `persistence_in_time`; morphe `persistent_homology` |

## 9. References

- Atkin, R. H. (1974), *Mathematical Structure in Human Affairs*, London: Heinemann.
- Baez, J. C. (1998), Spin foam models, *Classical and Quantum Gravity* 15, 1827–1858.
- Dowker, C. H. (1952), Homology groups of relations, *Annals of Mathematics* 56, 84–95.
- Edelsbrunner, H. and Harer, J. (2010), *Computational Topology: An Introduction*, AMS, §VII.1.
- Oeckl, R. (2005), *Discrete Gauge Theory: From Lattices to TQFT*, Imperial College Press.
- Reisenberger, M. and Rovelli, C. (1997), "Sum over surfaces" form of loop quantum gravity, *Physical Review D* 56, 3490–3508.

- Cartwright, D. and Harary, F. (1956), Structural balance: a generalization of Heider's theory, *Psychological Review* 63, 277–293.
- Kogut, J. and Susskind, L. (1975), Hamiltonian formulation of Wilson's lattice gauge theories, *Physical Review D* 11, 395–408.
- Toulouse, G. (1977), Theory of the frustration effect in spin glasses: I, *Communications on Physics* 2, 115–119.
- Zaslavsky, T. (1982), Signed graphs, *Discrete Applied Mathematics* 4, 47–74.
- Hatcher, A. (2002), *Algebraic Topology*, ch. 2: cellular homology.
