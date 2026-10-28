# The Formalisms of *Chemical Structure Generation through Symmetry Breaking*, Made Exact

*Status: draft (2026-10-29). Activity: research.* A verification of the manuscript *An Inquiry into Chemical Structure Generation through Symmetry Breaking: Two Formalisms Drawing on Group-Equivariant Grammar and Group Field Theory* (W. Huang, v1.4, 2026-10-28), claim by claim, with the computations that the plug-in **Symmetry and generation** of the Chemistry Workbench runs (its engine's `checks`, which the workbench's tests require to pass). Section, equation and statement numbers refer to the manuscript as compiled: sections 1–12, appendices A–L, statements numbered within their section (Definition 4.2, Proposition C.1). The statements this note proves are numbered Result 1 to Result 11, apart from the manuscript's.

## 1. Purpose and Scope

The manuscript separates four levels of claim: representational admissibility, mathematical validity, physical adequacy and ontological attribution (§3). This note examines the second level only: whether each definition is well-formed, each proposition true under its hypotheses, and each computed value correct. It does not judge the conjectures' physical truth, though §6 notes where known experiments bear on them.

The verification has two halves. Every proposition was read with its proof. Every claim that admits computation was then executed on the foundations of The Science Infrastructures:

- **Graph rewriting.** morphe's double-pushout graph rewriting with relabelling (`formal/graph_rewriting.hpp`, new).
- **Symmetry.** physis's point groups (`symmetry/point_group.hpp`, new).
- **The field.** physis's stranded field theories on finite groups (`quantum_gravity/stranded.hpp`, new).
- **Kinetics.** physis's reaction kinetics (`chemistry/kinetics.hpp`, the first content of physis's planned `chemistry` module).
- **Existing methods.** morphe's homology, characters and canonical forms; calx's integrators and random streams.

The general methods were written in the foundations, each with checks against independent results. The plug-in holds only what names molecules.

Of the 34 claims examined:

| Status | Count | What it means |
|---|---|---|
| verified or proved | 20 | correct as stated |
| stronger than stated | 3 | correct; a stronger statement holds |
| to amend | 7 | correct only with a change of definition, hypothesis or wording |
| errors | 4 | wrong as written |

All 55 computational checks behind these claims pass. The errors and amendments are listed for the revision in §7.

## 2. Formalism I: Graph Grammars

### 2.1 Rules Need Relabelling

**Definition 4.2** takes a rule to be a span L ← K → R of monomorphisms in a category **ChemG** whose morphisms preserve atomic identities and the selected attributes.

**Result 1.** In a category of labelled graphs whose morphisms preserve labels, no span L ← K → R of monomorphisms changes the label of a vertex it keeps.

*Proof.* A kept vertex is v_L = l(k) in L and v_R = r(k) in R for one k ∈ K. Since l and r preserve labels, label(v_L) = label(k) = label(v_R). ∎

The acid–base edit of Fig. proton-graph-edit keeps every atom: the transferred H* is retained, and the acid and water oxygens change charge, O → O⁻ and O → O⁺. Charge is an attribute declared relevant, since §4.2 requires charge accounting. So the edit is not a rule of Definition 4.2. The text recognizes the difficulty for bond orders ("the choice of attributed-graph category must make the proposed relabeling operation valid") but not for atoms.

**Amendment.** Adopt double pushouts with relabelling (Habel and Plump 2002, *Relabelling in graph transformation*, ICGT, LNCS 2505, 135–147). K is partially labelled, and a kept item whose label differs between L and R is relabelled by the double pushout. This is the semantics of morphe's `graph_rewriting.hpp`.

In it, the edit is one rule (the Grammar view, "Acid–base edit"):
- it keeps all 7 matched atoms;
- it deletes O_acid–H* and creates O_water–H*;
- it relabels the two oxygens.

Its derivation conserves every element and the total charge, and the species-level balance AN = 0, zᵀN = 0 holds for CH₃COOH + H₂O → CH₃COO⁻ + H₃O⁺.

Hyperedge replacement (Definition 4.1, Eq. (production-schema)) is the special case in which L is a nonterminal vertex with its port edges, K its attachment vertices, and R the replacement graph. The ethanol derivation of Fig. grammar-syntax is executed this way: one terminal graph, by one derivation. As a check of the engine itself, one step of string rewriting (morphe's `grammar.hpp`) and one step of graph rewriting on the words drawn as paths give the same forms on 300 random forms of a context-sensitive grammar for aⁿbⁿcⁿ.

A precision on Eq. (production-schema): the condition χ_p(x, C, m) of a *syntactic* production reads the physical state x, which the sentential graph does not hold. In the syntactic language L_term(Γ; C) the condition should read the sentential graph and the context; conditions on x belong to Definitions 4.2–4.3.

### 2.2 The Stereocentre: τ, Its CIP Descriptor, and the Terminal Language

Case Study I's grammar has:
- a start graph S, the nonterminal A_tet with ports 1–4 to H, F, Cl, Br;
- two productions p± that relabel A_tet to a carbon carrying τ = ±1.

Executed, it has two terminal classes forming one orbit of the reflection p, and one constitutional class once τ is forgotten (**Proposition 9.1**, verified). Embedded on a regular tetrahedron, the two configurations are congruent by exactly one isometry, which is improper (det = −1). Their common point group is C1, so each is chiral.

**Type of L_term.** Eq. (chemical-terminal-language) defines L_term as a set of typed-graph isomorphism classes. Eq. (case-chiral-language) writes L_term(Γ) = {[X₊]_{SE(3)}, [X₋]_{SE(3)}}, classes of embeddings. The grammar produces τ-labelled graphs; the embeddings are their realization.

**Amendment.** Write L_term = {[X₊]_≅, [X₋]_≅}, with ≅ preserving τ, and give the embedding as the realization map.

**τ is relative to an order.** τ(X) = sgn det[r_F − r_H, r_Cl − r_H, r_Br − r_H] depends on the order of the ligands. Exchanging two ligand labels exchanges two columns (or moves the base point), so an odd permutation σ of the ports gives

  τ(σ, Q·X) = det Q · sgn σ · τ(X).

A relabelling therefore flips τ exactly as a reflection does; the engine checks this for F ↔ Cl. If ligand permutations are among the symmetries considered, G is Z₂ × S₄ acting on τ by the character (Q, σ) ↦ det Q · sgn σ, not Z₂ alone.

**Result 2 (τ and the CIP descriptor).** Let the ligands of a tetrahedral centre be single atoms of distinct atomic numbers, with priorities a > b > c > d (IUPAC 2013 Recommendations P-92.1.4, Rule 1a). The centre is R exactly when det[a − d, b − d, c − d] < 0.

*Proof.* Seen with d pointing away, the viewer is on the side of the plane through a, b, c away from d, so a − d points toward the viewer. The cycle a → b → c runs counterclockwise to the viewer exactly when ((b − a) × (c − a)) · (a − d) > 0, that is det[b − a, c − a, a − d] > 0. Subtracting the first column from the others and permuting cyclically,

  det[a − d, b − d, c − d] = det[a − d, b − a, c − a] = det[b − a, c − a, a − d].

So counterclockwise (S) is det > 0, and clockwise (R) is det < 0. ∎

The manuscript's order (H, F, Cl, Br) is increasing priority, d < c < b < a. Hence τ = sgn det[c − d, b − d, a − d] = −sgn det[a − d, b − d, c − d], and **τ = +1 is R**: X₊ is (R)-bromochlorofluoromethane and X₋ is (S). The engine confirms this on its embeddings (the Stereocentre view). Independently, RDKit's CIP labeller (Hanson et al. 2018, *J. Chem. Inf. Model.* 58, 1755), given the engine's coordinates, assigns R to X₊ and S to X₋ (`tests/fixtures/cip.json`, made by `tests/make_fixtures.py`).

### 2.3 Equivariance, Closure and Counting

**Proposition 4.1** (equivariance of derived rewritings) is correct, and close to tautological under its hypotheses. Its statement does not name the hypothesis that the rule family is closed under G, which sits in Definition 4.3. The hypothesis is needed. The engine checks equivariance at every derivable graph:
- the stereocentre's rules {p₊, p₋} are closed under p, and their derivations commute with it;
- the glycine coupling rule commutes with S₃ relabelling;
- {p₊} alone is not closed, and equivariance fails at the start graph.

**Result 3 (matches and occurrences).** For injective matches,

  #matches(L, H) = |Aut L| · #occurrences(L, H).

*Proof.* Aut L acts on the injective morphisms m: L → H by m ↦ m ∘ α. The action is free, because m is a monomorphism, so m ∘ α = m forces α = id. Two matches with the same image differ by the automorphism m⁻¹ ∘ m′ of L. So the orbits are the images. ∎

This makes exact the Remark of App. B.2 on state counting. The acid–base edit has 2 matches (the two water hydrogens are exchanged by an automorphism of L), |Aut L| = 2, and one occurrence. Gillespie's propensities h_r = Π C(n_i, α_ir) count occurrences, which is the convention physis's kinetics follow.

**Result 4 (derivations of chains).** From n units labelled 1..n, each with a free amine terminus and a free carboxyl terminus, the coupling rule joins a free carboxyl to a free amine of another chain. Its terminal graphs are the n! directed chains, each reached by exactly (n − 1)! derivations, n!(n − 1)! in all. They form one orbit of S_n.

*Proof.* Each unit gains at most one outgoing and one incoming link, and the condition forbids closing a cycle, so a terminal graph is a directed Hamiltonian path: a linear order of the units. A derivation of a given order adds its n − 1 links one at a time. Every order of adding them is valid, since a link of a path always joins two different chains, and each step has exactly one match. ∎

For n = 3 this is Figs. grammar-growth and grammar-derivation-tree made complete: 6 labelled backbones, 12 derivations, 2 per backbone, and 2 H₂O created. The figures show the two derivations of one backbone. Because a peptide bond is directed, joining units 1 and 2 can give Gly₁–Gly₂ or Gly₂–Gly₁.

### 2.4 Accounting

Eq. (stoich), AN = 0 and zᵀN = 0, is checked by physis's `balance()` for:
- the acid–base equilibrium;
- the condensation 3 Gly → Gly₃ + 2 H₂O (C₂H₅NO₂ ×3 = C₆H₁₁N₃O₄ + 2 H₂O);
- Frank's network of App. I.4. The left kernel of its stoichiometric matrix is spanned by (1, 1, 1, 2), so a + ℓ + r + 2d is conserved, exactly the idealized accounting the appendix assumes.

## 3. Dynamics and Selection

### 3.1 The Jump Generator

**Proposition 4.2** is correct: (𝓛_{gC} f)(gx) = 𝓛_C(f ∘ g)(x), by reindexing the sum over rule–match occurrences.

The engine builds the generator of Frank's network on the states reachable from (A, L, R, D) = (24, 1, 1, 0), 1872 states, and checks that it commutes exactly, in rational arithmetic, with the mirror L ↔ R. The law of the absorbing state is solved exactly by the embedded jump chain:
- from (8, 1, 1, 0) it is mirror-symmetric in rational arithmetic;
- from (24, 1, 1, 0) it is symmetric to rounding;
- every outcome except mutual extinction (probability 0.039) holds one enantiomer only.

This is the manuscript's distinction between symmetry of the law and symmetry of a sample (Remark, App. B.3), computed. The Selection view draws the bimodal law and Gillespie paths, each choosing a sign.

### 3.2 Stabilizers Along Flows

**Lemma 4.1** states that along a unique G_C-equivariant flow, H_{x₀} ⊆ H_{Φᵗx₀}. A stronger statement holds.

**Result 5.** If Φ is a flow (Φᵗ invertible with inverse Φ⁻ᵗ, as for every ODE with unique solutions both ways), then H_{Φᵗx₀} = H_{x₀} for every t where the flow is defined.

*Proof.* Let g ∈ H_{Φᵗx₀}, with gC = C. By equivariance, g x₀ = g Φ⁻ᵗ Φᵗ x₀ = Φ⁻ᵗ(g Φᵗ x₀) = Φ⁻ᵗ Φᵗ x₀ = x₀. With the lemma's inclusion, the two stabilizers are equal. ∎

So the stabilizer is constant along a trajectory, as the isotropy of an equivariant flow is (Golubitsky, Stewart and Schaeffer 1988, *Singularities and Groups in Bifurcation Theory II*, ch. XIII). A larger stabilizer can only be approached as t → ∞, or along non-invertible semiflows and maps.

In Frank's model integrated by calx's Dormand–Prince method, a racemic start stays exactly racemic (m ≡ 0 to the last bit). From ℓ₀ > r₀, the excess stays positive.

### 3.3 Frank's Network and the Normal Form

The equation ṁ = (kᵢs/2) m(1 − m²) of App. I.4 is correct. The engine checks it along the integrated solution to a relative 1.5 × 10⁻¹⁶. More holds.

**Result 6.** In ℓ̇ = kaℓ − kᵢℓr, ṙ = kar − kᵢℓr, ȧ = −ka(ℓ + r), with k, kᵢ, a₀ > 0 and ℓ₀ > r₀ ≥ 0:

1. δ = ℓ − r obeys δ̇ = kaδ, so δ(t) = δ₀ exp(k ∫₀ᵗ a) > 0;
2. a(t) → 0;
3. r(t) → 0;

so m(t) → 1. By symmetry, m → −1 when ℓ₀ < r₀.

*Proof.*

1. ℓ̇ − ṙ = ka(ℓ − r).
2. s = ℓ + r ≥ δ ≥ δ₀, so ȧ ≤ −kδ₀a and a ≤ a₀e^{−kδ₀t}.
3. ℓ = r + δ ≥ δ₀, so ṙ ≤ r(ka − kᵢδ₀). Since a decays exponentially, ka < kᵢδ₀/2 after some time, and r then decays exponentially.

Finally m = δ/(δ + 2r) → 1. ∎

The decay rate of the minor enantiomer is about kᵢδ_∞, so convergence is slow when the initial excess is small: from m₀ = 0.005, m(400) = 0.96; from m₀ = 0.1, m(400) = 1 to rounding.

With uncatalysed production A → L, A → R at rate k₀ (as in Frank 1953, *Biochim. Biophys. Acta* 11, 459), the equation gains a term: ṁ = (kᵢs/2) m(1 − m²) − (2k₀a/s) m.

The cubic normal form (§7.3, App. E.2) is correct: f′(0) = a and f′(±√(a/b)) = −2a.

### 3.4 Selection Rules and Point Groups

**Proposition 7.1** and its extension in App. E.4 are correct. They are Neumann's principle: a response transforming by ρ lies in the fixed subspace of ρ(H), whose dimension is the mean of the character,

  dim V^H = (1/|H|) Σ_{h∈H} χ_ρ(h),

with χ = tr Q for polar vectors and χ = det Q · tr Q for axial ones.

physis's `point_group()` finds the group of a molecule from the positions of its nuclei. Its check recovers all 26 groups from C1 to Ih from generic orbits, with:
- their orders and Schoenflies symbols;
- their numbers of conjugacy classes;
- the polar groups (C1, Cs, Cn, Cnv), the groups that fix an axial vector (C1, Ci, Cs, Cn, Cnh, S2n), and the chiral groups (C1, Cn, Dn, T, O, I).

morphe computes each group's character table exactly.

On idealized geometries (VSEPR angles, bond lengths as sums of covalent radii):

| Molecule | Point group | Fixed vectors | Consequence |
|---|---|---|---|
| (Z)-1,2-dichloroethene | C2v | 1 | a dipole is allowed |
| (E)-1,2-dichloroethene | C2h | 0 | no dipole |
| CHBrClF | C1 | 3 | chiral |
| hydrogen peroxide, skew | C2 | 1 | chiral and polar |
| benzene | D6h | 0 | — |
| methane | Td | 0 | — |
| SF₆ | Oh | 0 | — |

The (Z)/(E) result confirms App. I.1: one constitutional graph, a dipole allowed in one configuration and forbidden in the other.

## 4. Kernels, Readouts and Inference

### 4.1 Contrastive Models

**Proposition 3.1** (non-entailment) and **Proposition 8.1** (restriction-relative observational equivalence) are correct, with complete proofs in Apps. A and F. The two-state models of App. F agree in C₀ and differ by |q − p| in C*.

### 4.2 Readouts, Factorization and Contraction

All four are correct:
- **Proposition 5.2** (covariance of generated structures), by a change of variables. The engine checks it on a finite ensemble with an equivariant readout. As App. C.3 remarks, a constant readout is equivariant and erases the broken sectors.
- **Proposition 6.1** (structural sufficiency), with the final σ-algebra. Example 6.1's distance is δ.
- **Proposition 10.1**: the data-processing inequality for total variation, and its mixture bound. Both are checked exactly in rationals on 200 random finite cases each.
- **Proposition 11.1** (unrestricted readout).

## 5. Formalism II: The Tetrahedral Field Theory

### 5.1 The Field, the Measure and the Regulator

The interaction of Example 5.1, V[ψ] = Σ ψ(a,b,c) ψ(c,d,e) ψ(e,b,f) ψ(f,d,a), is Boulatov's tetrahedral vertex (Boulatov 1992, *Mod. Phys. Lett. A* 7, 1629), which should be cited. Its four factors are triangles, any two sharing one index, an edge of the tetrahedron. The values V[1] = 64 and V[(−1)^a] = 0 over Z₂ are correct with counting measure. §5.2 announces a normalized Haar measure, under which V[1] = 1; the measure should be stated.

**Result 7 (the regulator is necessary).** V is homogeneous of degree 4, even, and takes negative values, for example V = −368 at ψ = (−2, −2, −2, −2, 2, −2, 1, −2) over Z₂ in the order (000, 001, …, 111). So V is unbounded below and e^{−λV} is integrable for no λ ≠ 0. Z_Λ exists only with the regulator; the series in λ is formal.

*Proof.* V[tψ] = t⁴ V[ψ]. Along a ψ with V[ψ] < 0, −λV[tψ] → +∞ as t → ∞ if λ > 0. If λ < 0, use ψ = 1, where V = 64. ∎

### 5.2 The Symmetry Group of the Interaction

**Proposition C.1** (invariance under simultaneous translation) is correct for every finite group, abelian or not. The full symmetry is larger.

**Result 8.** Let (U_g ψ)(x₁, x₂, x₃) = ψ(g₁⁻¹x₁, g₂⁻¹x₂, g₃⁻¹x₃) for g ∈ G³. Then V ∘ U_g = V for all ψ exactly when g₁ = g₃. The group of slot translations preserving V is {(g, h, g)} ≅ G × G.

*Proof.* The strands of the vertex join slot 1 to slot 3 (the indices a, c, e, f) and slot 2 to slot 2 (b, d). If g₁ = g₃, substituting g₁⁻¹ or g₂⁻¹ in every index is a bijection of the sums.

Conversely, for the field δ_y supported on one point y, V[δ_y] = 1 if y₁ = y₃ and 0 otherwise: all four factors must equal y, which forces a = c = e = f = y₁ = y₃. Take y = (1, 1, 1). Then V[δ_y] = 1, while V[U_g δ_y] = V[δ_{gy}] = [g₁ = g₃]. ∎

The same holds for right translations. V is also invariant under reversing the slots, ψ(a,b,c) ↦ ψ(c,b,a). The engine checks invariance under (g, h, g) for Z₂, Z₃ and S₃, and its failure elsewhere. Over Z₂ the invariant translations are exactly 000, 010, 101 and 111. The manuscript's U_p is the diagonal (111), one element of this group.

### 5.3 Global and Gauge Symmetry for Abelian Groups

§5.2 says the diagonal right-invariance of group field theories is a gauge redundancy, "not to be conflated with a chemically observable symmetry that may be broken". §5.3 then breaks the simultaneous translation U_p of Z₂.

**Result 9.** Let G be abelian, P = |G|⁻¹ Σ_h R_h the gauge projection with (R_hψ)(x) = ψ(x₁h, x₂h, x₃h), and η a linear functional with η ∘ U_p = χ(p) η for a nontrivial character χ of G (for Z₂ and χ(p) = −1: an odd functional). Then η ∘ P = 0.

*Proof.* For abelian G, R_h = U_{h⁻¹}, so

  |G| η[Pψ] = Σ_h η[U_h ψ] = Σ_h η[U_{ph} ψ] = χ(p) Σ_h η[U_h ψ].

Since χ(p) ≠ 1, the sum vanishes. ∎

So the symmetry §5.3 breaks is the redundancy §5.2 warns about. On the gauge-invariant fields of a group field theory (Pψ = ψ), every odd order parameter, η included, vanishes identically; the engine checks this over Z₂.

**Amendment.** Either declare that the model has no gauge projection, so that U_p is global and §5.2's caution does not apply to it; or take G non-abelian, where diagonal left translations commute with the right projection and remain global.

### 5.4 Order Parameters

Over Z₂ a field is eight numbers, expanded in the characters χ_s(x) = (−1)^{s·x}, s ∈ Z₂³. The simultaneous translation multiplies χ_s by (−1)^{|s|}. So the odd part of the field is four-dimensional, spanned by s ∈ {100, 010, 001, 111}.

**Result 10.** V[χ_s] = 64 if s₁ = s₃ and 0 otherwise.

*Proof.* In Π_k χ_s(factor k), each index contributes (−1) to the power of its value times the sum of the bits of s at its two slots:
- a, c, e and f sit in slots 1 and 3, giving s₁ + s₃;
- b and d sit in slot 2 twice, giving 2s₂ ≡ 0.

So V[χ_s] = Σ (−1)^{(s₁+s₃)(a+c+e+f)} · 4, which is 2⁶ if s₁ = s₃ and 0 otherwise. ∎

More generally, expanding ψ = Σ ψ̂(s) χ_s gives V[ψ] = 64 · V[ψ̂]: the tetrahedral vertex is self-dual under the Fourier transform of Z₂³.

Consequences for §5.3:

1. **The order parameter lies in a flat direction.** η = 2⁻³ Σ (−1)^a ψ is the coefficient of χ_100, and V vanishes identically along it: V[tχ_100] = 0 for all t. The manuscript's test value V[(−1)^a] = 0 is this fact. So the quartic term of V_eff(η) = αη² + βη⁴ cannot come from the interaction; along that line the action is ½μ·8t² plus the regulator.

2. **A correct derivation.** Along the odd character χ_010 = (−1)^b, V[tχ_010] = 64t⁴. With S = ½μ⟨ψ, ψ⟩ + λV, the action restricted to ψ = η(−1)^b is

   S = 4μη² + 64λη⁴:

   the Landau form of Eq. (effective-potential) with α = 4μ and β = 64λ, derived rather than postulated. The double well appears for μ < 0, λ > 0, with minima ±√(−μ/32λ). The Field view draws both lines.

3. **η is incomplete.** Proposition 5.1's "every configuration with η ≠ 0 has trivial stabilizer" is true, but the converse fails: ψ = (−1)^b is odd and has η = 0. A field is fixed by p exactly when its odd part ψ − U_pψ vanishes. The order parameter should be that projection, or η should be said to detect breaking along one of four odd directions.

Proposition 5.1 itself (the minima ±√(−α/2β), with V(η±) = −α²/4β and V″(η±) = −4α) is correct.

### 5.5 The Wick Expansion

**Result 11.** Let the values ψ(x), x ∈ G³, be independent standard normals: the Gaussian measure of the kinetic term ½⟨ψ, ψ⟩, whose propagator identifies the slots of two factors one by one. Then

  ⟨Vⁿ⟩ = Σ_{π} |G|^{F(π)},

the sum over the (4n − 1)!! perfect matchings π of the 4n field factors (Wick pairings), F(π) the number of closed strands (faces) of π.

*Proof.* By Isserlis's theorem, E[Π ψ(x_i)] = Σ_π Π_{(i,j)∈π} δ_{x_i x_j}. Under π, the propagators identify index values slot by slot, and the strands of the vertex copies close into F(π) cycles. An assignment of the indices satisfies all the deltas exactly when it is constant on each cycle, which gives |G|^{F(π)} assignments. ∎

The group Aut(v) ≀ S_n acts on the pairings and preserves F. Here Aut(v) is the factor permutations of the vertex that carry strands to strands slot by slot: the cyclic group of order 4. A class Γ therefore holds 4ⁿ n!/|Aut Γ| pairings, and

  Z/Z₀ = Σ_n ((−λ)ⁿ/n!) ⟨Vⁿ⟩ = Σ_Γ (−4λ)^{n(Γ)} |G|^{F(Γ)} / |Aut Γ|

as a formal series. This is Eq. (complex-expansion) made exact. The sum runs over stranded Feynman graphs Γ, with |Aut Γ| of the stranded graph and the vertex factor 4. Different stranded graphs may have isomorphic 2-complexes, so the sum is not over 2-complexes.

physis checks the theorem against the Gaussian moment computed by brute force:

| Order | G = 1 | Z₂ | Z₃ | S₃ |
|---|---|---|---|---|
| ⟨V⟩ | 3 | 24 | 99 | 1368 |
| ⟨V²⟩ | 105 | 1792 | 17 145 | — |
| ⟨V³⟩ | 10 395 | 261 120 | — | — |

It also checks the class sizes, with 2, 10 and 54 classes at orders 1, 2 and 3.

The classes carry topology. At order 1, the one-vertex graph with two faces has H₁ = Z/4, the homology of the one-tetrahedron triangulation of L(4,1). At order 2 the classes and the homology of their 2-complexes (morphe) are:

| Faces | Pairings | \|Aut Γ\| | Connected | H₀, H₁, H₂ | BF amplitude over Z₂ |
|---|---|---|---|---|---|
| 8 | 1 | 32 | no | Z², 0, Z⁴ | 16 |
| 6 | 4 | 8 | no | Z², Z/4, Z² | 8 |
| 6 | 4 | 8 | yes | Z, 0, Z³ | 8 |
| 5 | 8 | 4 | yes | Z, 0, Z² | 4 |
| 4 | 4 | 8 | no | Z², Z/4 ⊕ Z/4, 0 | 4 |
| 4 | 4 | 8 | yes | Z, Z/2, Z | 4 |
| 3 | 32 | 1 | yes | Z, Z/4, 0 | 2 |
| 3 | 16 | 2 | yes | Z, Z/8, 0 | 2 |
| 3 | 16 | 2 | yes | Z, Z, Z | 2 |
| 3 | 16 | 2 | yes | Z, Z, Z | 2 |

### 5.6 The Complexes of Case Study II

**The six-cycle and its filled complex** (App. L.2) are correct: Betti numbers (1, 1) and (1, 0, 0), with ∂₁∂₂ = 0.

**The two-vertex complex** (Proposition 9.2, App. L.3) is correct:
- (V, E, F) = (2, 4, 6), χ = 4;
- the boundary matrices as printed, ∂₁∂₂ = 0;
- homology Z, 0, Z³, torsion-free.

It is the dual 2-skeleton of the 3-sphere cut into two tetrahedra. The two tetrahedra are glued face to face by the identity, and a tetrahedron glued to itself that way is the double of a 3-ball, S³. Its six faces are the six edges of the triangulation. H₂ = Z³ = Z^{4−1}, from the triangulation's four vertices (the dual 3-cells missing from the 2-skeleton) less one relation. Its BF amplitude with gauge projection is |G|^{F−E+V−1} = 8 over Z₂ and 216 over S₃.

Two amendments concern how the complex is said to arise.

1. **"Pairing corresponding field slots" picks one of 24 gluings.** Of the 24 ways to join two vertices factor to factor:
   - 4 give this complex (the four automorphisms of the vertex);
   - 4 give four faces, with H₁ = Z/2 (homology consistent with RP³);
   - 16 give three faces, with H₁ = Z and H₂ = Z (consistent with S² × S¹).

   Whether those gluings are manifolds, with spherical vertex links, is not checked here. The order-λ² term holds all of them, and 72 more pairings with self-loops (the table of §5.5).

2. **The complex is not generated by a field configuration.** Eq. (case-dipole-generated) writes 𝔉 = (Z₂, 3, ψ, V_tet) ⇝ k_dip: "one chooses the field data …, evaluates the quartic interaction, and glues". Feynman graphs arise from the Gaussian integral over all fields, by Wick contraction of ⟨V²⟩₀, not from a chosen ψ. The test fields ψ₀ and ψ₁ play no part in the gluing. The dipole's weight is |G|^F = 2⁶ = 64 without the gauge projection; that this equals V[1] = 2⁶ is a coincidence of two different counts.

## 6. Scientific Framing and Language

**Conjecture 7.1 seems already instantiated.** As an existence claim under the manuscript's own operational definitions (Definition 7.1 with an odd observable, Proposition 7.1), known experiments appear to satisfy it:
- stirred crystallization of sodium chlorate gives nearly homochiral batches whose optical rotation, an odd observable, the achiral solution lacks (Kondepudi, Kaufman and Singh 1990, *Science* 250, 975);
- the Soai reaction amplifies enantiomeric excess autocatalytically (Soai, Shibata, Morioka and Choji 1995, *Nature* 378, 767).

The conjecture should cite such instances and be sharpened to what they do not settle, for example a function that is not itself the odd observable.

**Notation.** Several symbols carry several meanings:
- **p**: the reflection, the productions p±, the ensemble p_C, a protocol p ∈ 𝒫, and a probability (App. F); "pp₊ = p₋" is hard to read.
- **ν**: the state lift, the stoichiometric increment, and the endpoint law.
- **K**: the DPO interface, the kinetic kernel K_Λ, the response kernel K_{C,J}, and the Markov kernels K₀, K*.
- **ℛ**: the matched events in Definition 4.3, but the rule family in Eq. (gram-generator) and App. H.
- **p_C (§5.4) and π_C (the glossary, §10)**: one ensemble under two names; likewise ℋ_𝔾 and ℋ_gr.
- **D**: a domain, a matrix, a product, a determinant and a set.

**Two sentences.**
- §4.5: "it is represented separately from the grammar derivation by a grammar's derivation sequence" does not parse. Proposed: "it is represented separately from the grammar's derivation sequence."
- §2.2: "The parser is an explanatory mathematical representation." A parser recognizes; the grammar or generator is meant.

## 7. What the Revision Should Change

**Errors.**
1. Definition 4.2 cannot change a kept atom's charge: adopt DPO with relabelling (§2.1).
2. For abelian G, §5.3 breaks the gauge symmetry of §5.2: declare no projection, or take G non-abelian (§5.3).
3. §5.3's η lies along a flat direction of V: use (−1)^b, which gives S = 4μη² + 64λη⁴ (§5.4).
4. Eq. (case-dipole-generated): the complex comes from the Wick expansion, not from ψ (§5.6).

**Amendments.**
5. Give η's incompleteness, or use the odd projection (§5.4).
6. Make the regulator part of Z_Λ, and call the series formal (§5.1).
7. Write Eq. (complex-expansion) over stranded graphs, with the factor 4 = |Aut v| (§5.5).
8. Say which of the 24 gluings is meant, and that the λ² term holds all of them (§5.6).
9. Note that τ is relative to the ligands' order; name X₊ = (R), X₋ = (S) (§2.2).
10. State the measure, counting or Haar (§5.1).
11. Cite instances of Conjecture 7.1 and sharpen it (§6).
12. Fix the notation and the two sentences (§6).

**Stronger statements available.**
13. Lemma 4.1: stabilizers are constant along invertible flows (§3.2).
14. App. I.4: the sign of the excess is invariant and m → ±1 (§3.3).
15. App. C: the symmetry of V is G × G (§5.2).

**Also to add.**
- Name closure under G among Proposition 4.1's hypotheses (§2.3).
- Write L_term by τ-labelled graphs (§2.2).
- Cite Boulatov 1992, Frank 1953, and Habel and Plump 2002.

## 8. What Was Kept, and Where It Lives

The manuscript's philosophy is kept whole: a symmetry-first, explicitly provisional modelling, with formal derivation, physical realization and ontological attribution kept apart. Its two formalisms are executed as it defines them, amended only where §7 says.

General mathematics went to the foundations:

| Foundation | What it gained |
|---|---|
| morphe (decision 0007) | graph rewriting with relabelling, matches, derivation spaces and equivariance |
| physis (decision 0011) | point groups and the responses they allow; real stranded field theories on finite groups and their Wick expansion; reaction networks with mass action and their exact Markov chains, in its new `chemistry` module |

The plug-in (`../engine/`) holds only what names molecules, provisionally, until a chemistry domain holds it:
- nine elements;
- idealized geometries;
- the CIP rule for a single-atom centre;
- the manuscript's examples and the checks of its claims.

## References

- Boulatov, D. V. (1992). A model of three-dimensional lattice gravity. *Mod. Phys. Lett. A* 7, 1629–1646.
- Cotton, F. A. (1990). *Chemical Applications of Group Theory*, 3rd ed. Wiley.
- Cordero, B. et al. (2008). Covalent radii revisited. *Dalton Trans.*, 2832–2838.
- Ehrig, H., Ehrig, K., Prange, U. and Taentzer, G. (2006). *Fundamentals of Algebraic Graph Transformation*. Springer.
- Frank, F. C. (1953). On spontaneous asymmetric synthesis. *Biochim. Biophys. Acta* 11, 459–463.
- Gillespie, D. T. (1977). Exact stochastic simulation of coupled chemical reactions. *J. Phys. Chem.* 81, 2340–2361.
- Golubitsky, M., Stewart, I. and Schaeffer, D. G. (1988). *Singularities and Groups in Bifurcation Theory*, vol. II. Springer.
- Habel, A. and Plump, D. (2002). Relabelling in graph transformation. *ICGT 2002*, LNCS 2505, 135–147.
- IUPAC (2013). *Nomenclature of Organic Chemistry: IUPAC Recommendations and Preferred Names 2013*, P-92.
- Isserlis, L. (1918). On a formula for the product-moment coefficient of any order of a normal frequency distribution in any number of variables. *Biometrika* 12, 134–139.
- Kondepudi, D. K., Kaufman, R. J. and Singh, N. (1990). Chiral symmetry breaking in sodium chlorate crystallization. *Science* 250, 975–976.
- Nye, J. F. (1957). *Physical Properties of Crystals*. Oxford.
- Soai, K., Shibata, T., Morioka, H. and Choji, K. (1995). Asymmetric autocatalysis and amplification of enantiomeric excess of a chiral molecule. *Nature* 378, 767–768.
