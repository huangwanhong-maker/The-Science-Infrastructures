# Group Grammars for Genomes: the GE-CSG Framework Made Exact

*Status: draft, 2026-10-24.*

This note serves learning and research. It takes the draft of the user, W. Huang, *Towards a Group Context-Sensitive Grammar for Genome Analysis* (2026; "the draft"), and states its framework as definitions with proved results, in the form in which the infrastructure implements and checks it.

The code:

- **morphe** holds the general machinery:
  - `algebra/action.hpp`: actions of finite groups on sets and words;
  - `algebra/characters.hpp` and `algebra/group_algebra.hpp`: characters, group algebras and their Wedderburn data;
  - `formal/grammar.hpp`, `formal/earley.hpp`, `formal/inside_outside.hpp`, `formal/equivariant.hpp` and `formal/breaking.hpp`: grammars, their parsers, group-invariant grammars and breakings.
- **The genomics domain** holds the nucleotides:
  - `grammar/nucleotides.hpp`, `strands.hpp`, `structures.hpp` and `code_symmetry.hpp`;
  - `rna/`: stochastic grammars of RNA structure and the equivariant EM (§2.6).

Every proposition below that can be checked numerically is checked, in `morphe/formal/tests/check_grammars.cpp`, `genomics/grammar/tests/check_grammar.cpp` and the `rna` module's checks. The numbers quoted are theirs.

The draft's thesis is that a finite group G is the primitive of a grammar of genomes. The alphabet, the rules and the parsing all derive from G, its group algebra ℂ[G] and the latter's Wedderburn decomposition, and concrete symbols arise by a symmetry-breaking map. Much of this survives exactly; some of it needs restating before it can be proved; a few statements are not true as written. §8 lists the amendments; §9 says where the draft's open questions stand.

## 1. Groups Acting on Grammars

### 1.1 Actions on Alphabets and Words

**Definition 1.1.** An *action* of a finite group G on a finite alphabet Σ is a homomorphism G → Sym(Σ), written x ↦ g·x. A *sign* is a homomorphism ε: G → {±1}. The *action on words* is
- g·(x₁…xₙ) = (g·x₁)…(g·xₙ) when ε(g) = +1,
- g·(x₁…xₙ) = (g·xₙ)…(g·x₁) when ε(g) = −1.

Reversal commutes with every letter permutation, so g ↦ (letters, ε(g)) is a homomorphism into Sym(Σ) × ℤ/2. This is an action of G on Σ* by automorphisms and anti-automorphisms of the free monoid. The strand symmetries of DNA are of this kind: complement acts by letters, reversal by sign.

### 1.2 Actions on Grammars and Invariant Languages

**Definition 1.2.** Let Γ = (V, Σ, R, S) be a grammar.
- G *acts on Γ* when it acts on Σ (with a sign) and on V, fixing S. An element carries a rule α → β to g·α → g·β, the forms mapped symbol by symbol and reversed when ε(g) = −1.
- Γ is *G-invariant* (the draft's "G-equivariant", Definition 3.3) when g·R = R for every g ∈ G.

**Proposition 1.3.** If Γ is G-invariant, then S ⇒* u implies S ⇒* g·u. Hence L(Γ) is G-invariant. If A ∈ V is fixed by G, then A ⇒* u if and only if A ⇒* g·u.

*Proof.* Let u = u₀ ⇒ u₁ ⇒ … be a derivation, each step replacing an occurrence of α in uᵢ by β for a rule α → β.
- When ε(g) = +1, g·uᵢ ⇒ g·uᵢ₊₁ by the rule g·α → g·β ∈ R.
- When ε(g) = −1, g·uᵢ contains g·α (the reversed image) where uᵢ contained α, and replacing it by g·β gives g·uᵢ₊₁.

The image of a derivation is therefore a derivation. For the last claim, apply this to g and to g⁻¹. ∎

### 1.3 The Number of Rule Orbits

**Proposition 1.4** (Burnside's lemma; Rotman 1995, *An Introduction to the Theory of Groups*, Theorem 3.22). The rules of a G-invariant grammar fall into |R/G| = (1/|G|) Σ_g |Fix_R(g)| orbits. This number is at least |R|/|G|, with equality exactly when no element other than 1 fixes a rule.

*Proof.* The orbit count is Burnside's lemma. Fix_R(1) = R and every other term is non-negative, so the sum is at least |R|. It equals |R| exactly when every other term vanishes. ∎

Storing one representative per orbit compresses R by a factor |R|/|R/G| ≤ |G|. Checked on the closures of random grammars:
- under S₃ on three letters, 250 rules fall into 83 orbits;
- under the Klein group on four, 259 rules fall into 70.

The draft states the inequality the other way (§8, item 1).

## 2. Parsing Up to Symmetry

### 2.1 Matching Terminals by Their Orbits

Let π: Σ → Σ/G be the orbit map, and Γ/G the image of Γ under it (terminals renamed by their orbits; `formal::orbit_grammar`). For context-free Γ, L(Γ/G) = π(L(Γ)), since the context-free languages are closed under letter-to-letter homomorphisms with the image grammar as witness (Hopcroft and Ullman 1979, Theorem 6.2).

**Proposition 2.1.** Earley's parser whose scan rule matches a rule's terminal a against an input letter x by x ∈ G·a alone (the draft's §5.2.2) accepts exactly π⁻¹(L(Γ/G)) = π⁻¹(π(L(Γ))).

This set is the closure of L(Γ) under the group Gⁿ acting position by position, independent elements at each position, while Γ is invariant only under the diagonal G. It contains L(Γ), in general strictly.

*Proof.* With the scan rule so weakened, the parser is Earley's parser for Γ/G run on π(w). ∎

*Example.* S → aa | bb is invariant under a ↔ b and has one rule orbit. The orbit-matching parser accepts ab; the language does not contain it. On the random invariant grammars of the checks, it accepted 312 words outside the language under S₃ and 264 under the Klein group.

### 2.2 The Equivariant Earley Parser

A parser that uses only one rule per orbit must bind the group element of each rule once. Assume G fixes every nonterminal and no element reverses words. An *equivariant item* [r̄ → α • β, i, C] carries a set C ⊆ G:
- prediction creates [r̄ → • γ, j, G] for each orbit representative r̄;
- scanning a terminal a against x replaces C by {g ∈ C : g·a = x};
- completion leaves C unchanged.

The children's elements are independent of the parent's, since nonterminals are fixed and R contains g·r for every g separately.

**Lemma 2.2.** For a given item, C is either empty or a left coset g₀K of K = ⋂_t Stab(a_t), the intersection over the terminals a_t that r̄ has scanned before the dot.

*Proof.* {g : g·a = x} is empty or a left coset of Stab(a), and the intersection of left cosets of subgroups is empty or a left coset of their intersection. ∎

**Proposition 2.3.** The equivariant item [r̄ → α • β, i, C] in the set S_j stands for the plain Earley items [g·r̄ → α • β, i] for g ∈ C. This correspondence is a bijection between the plain chart and the disjoint union, over the equivariant chart, of these sets. Hence:
- the equivariant parser accepts w exactly when Earley's parser of the full rule set does;
- the equivariant chart is never larger than the plain one.

*Proof.* By induction on the order in which items are made, each rule of the parser is the image of the corresponding plain rule. Two items with the same (r̄, dot, i, j) and different sets have, by Lemma 2.2, disjoint cosets of K. Since Stab(r̄) ⊆ K, distinct cosets give disjoint sets of rules {g·r̄}. So no plain item is represented twice. ∎

**What the equivariant parser saves.**
- **Rules and predictions.** It stores |R/G| rules and makes one prediction per orbit.
- **Not the worst case.** Its chart is the plain chart taken up to the action, so the worst-case time stays O(n³ |R|), not O(n³ |R|/|G|) as the draft's Proposition 5.1 states (§8, item 6).
- **Measured.** The chart had 209,157 items against 389,304 on the S₃ grammars, and 84,176 against 198,416 on the Klein ones.
- **Hairpins.** On hairpin grammars under complement (§5) it made 132,064 items against 209,306.

### 2.3 Sharing Spans Across Orbits

By Proposition 1.3, a nonterminal A fixed by G derives a span u exactly when it derives g·u. The inside algorithm may therefore decide the original nonterminals once per *class* of spans and reuse the answer for every span of the class (`formal::span_recognition`). With the trivial group the class is the span's content (repeated substrings); the group widens it to their images.

This is exact, and it is where the draft's "share parse-tree caches across symmetry-equivalent subsequences" (§1.2) holds as stated.
- On the checks' random grammars, sharing by orbit reused 25,128 spans against 11,232 by content alone (S₃), and 13,440 against 5,088 (the Klein group).
- On hairpins under the strand group it reused 4,884 of 15,566 spans.

### 2.4 Context-Sensitive Grammars

A monotone grammar derives w only through forms no longer than w, so membership is decided by searching them. This takes time exponential in the worst case, and the problem is PSPACE-complete in general (Hopcroft and Ullman 1979, §9.3; Karp 1972).

**Proposition 2.4.** Let G fix S and Γ be G-invariant. Searching the canonical forms (least images) instead of the forms visits the orbits of the derivable forms. It decides the same question, visiting between |forms|/|G| and |forms| classes.

*Proof.* By Proposition 1.3 the set of forms derivable from S is a union of orbits, and the successors of g·u are the images of u's. The canonical search therefore explores this set's orbits. The bounds are Proposition 1.4's, applied to the derivable set. ∎

The factor is at most |G|, so the complexity class does not change (§8, item 10). Measured on invariant monotone grammars under the Klein group: 61,248 classes against 241,248 forms.

### 2.5 Graduated Breaking, Measured

The coarse-to-fine theorem of `formal/breaking.hpp` makes the draft's graduated breaking (§6.3) a pruning that loses nothing: an image grammar's inside–outside pass keeps every span of every complete fine parse.

Measured on ten of *Escherichia coli*'s tRNA genes, with the grammar of a Watson–Crick hairpin (stem ≥ 4, loop 3–8) anywhere:
- parsed first over IUPAC's R/Y alphabet (its image under the breaking of ⟨τ⟩), it keeps 37.7% of the spans open;
- parsed first over W/S, it keeps 40.2%;
- in both cases the fine inside count over the kept spans equals the unrestricted count exactly.

### 2.6 Weighted Grammars, and the Equivariant EM

A weighted grammar is G-invariant when w(g·r) = w(r). Its inside values then obey I(A, u) = I(A, g·u) for fixed A, by Proposition 1.3 applied to derivations with their weights.

For a stochastic grammar with probabilities tied along orbits, θ_r = θ_O for r ∈ O, maximum likelihood by EM must respect the tie. Write c_r for a rule's expected count (the E-step: inside–outside, `formal::expected_counts`), C_O = Σ_{r∈O} c_r, and A for the shared left side of an orbit O, which is shared since nonterminals are fixed.

**Proposition 2.6.** Under the constraints Σ_{O ⊆ R_A} |O| θ_O = 1, the expected complete log-likelihood Σ_r c_r log θ_r is maximised by

  θ_O = C_O / (|O| · Σ_{O′ ⊆ R_A} C_O′).

So the tied EM is an EM algorithm for the quotient model and never decreases the likelihood (Dempster, Laird and Rubin 1977, J. R. Statist. Soc. B 39, 1–38, Theorem 1).

*Proof.* The objective is Σ_O C_O log θ_O. The Lagrangian condition for each A gives C_O/θ_O = λ_A |O|, so θ_O = C_O/(λ_A |O|). The constraint fixes λ_A = Σ_{O′ ⊆ R_A} C_O′. The objective is concave in log θ on the constraint set, so this point is the maximum. ∎

This answers the draft's Q1 as far as convergence goes: the orbit-tied inside–outside algorithm is an EM algorithm. The `rna` module runs it on *Escherichia coli*'s tRNAs, checks that the likelihood never decreases, and compares tied and untied models.

## 3. The Group of the Nucleotides

### 3.1 The Bases as a Torsor of the Klein Group

**Proposition 3.1.** The two involutions below commute, and with their product σ = κτ (A ↔ C, G ↔ T) they form the Klein four-group V₄ = {1, κ, τ, σ} ≅ C₂ × C₂, acting freely and transitively on {A, C, G, T}:
- Watson–Crick complementarity κ (A ↔ T, C ↔ G);
- transition τ (A ↔ G, C ↔ T: the exchange within the purines and within the pyrimidines).

The bases are a V₄-torsor. With the domain's codes A = 0, C = 1, G = 2, T = 3, the action is exclusive or: κ = ⊕3, τ = ⊕2, σ = ⊕1. The two-bit encoding is therefore the torsor's trivialisation at the base point A.

*Proof.* By computation. The checks verify the action axioms and that each stabiliser is trivial. ∎

### 3.2 IUPAC's Codes as Cosets

**Proposition 3.2.** The cosets of the subgroups of V₄ are exactly the four bases, IUPAC's six two-letter codes and N (Cornish-Bowden 1985, Nucleic Acids Research 13, 3021–3030):
- ⟨τ⟩: R = {A, G}, Y = {C, T};
- ⟨κ⟩: W = {A, T}, S = {C, G};
- ⟨σ⟩: M = {A, C}, K = {G, T}.

The three-letter codes B, D, H and V (the complements of single bases) are cosets of no subgroup.

*Proof.* V₄ has five subgroups: 1, the three of order 2, and V₄. Their cosets are the sets listed. A set of three elements cannot be a coset, since its size does not divide 4. ∎

### 3.3 The Classifications as Characters

**Proposition 3.3.** V₄ has four irreducible characters, all real and ±1-valued, three of them nontrivial, each determined by its kernel, a subgroup of order 2. Through the torsor at a base point σ₀ (χ(g·σ₀) := χ(g)), each nontrivial character is a function on the bases, constant on the cosets of its kernel. The three are the classical binary classifications:
- purine–pyrimidine (kernel ⟨τ⟩);
- weak–strong, two or three hydrogen bonds in a pair (kernel ⟨κ⟩);
- amino–keto (kernel ⟨σ⟩).

They form a group, χ_RY χ_WS = χ_MK, so two classifications name a base.

Moving the base point to g₀·σ₀ multiplies each character by χ(g₀). The partition is canonical; which class the character calls +1 is the one choice a breaking makes. In the two-bit codes, purine–pyrimidine is the low bit, amino–keto the high bit, and weak–strong their parity.

*Proof.* The characters of C₂ × C₂ are the products of those of the factors (Serre 1977, *Linear Representations of Finite Groups*, §3.2). The rest is computation, which the checks carry out. ∎

### 3.4 Idempotents, Letters and the Breaking Chain

**Proposition 3.4.** The minimal idempotents of ℂ[V₄] are the four e_χ = ¼ Σ_g χ(g) g. They are the characters, not the letters.

Under the torsor the letters are the elements g ↦ g·σ₀ (the basis of ℂ[V₄]), and the Fourier transform carries this basis to the idempotents. "Terminal symbols arise as minimal idempotents" (the draft's abstract and §2.3) is therefore, for the nucleotides:
- the classifications arise as idempotents;
- the letters arise from the trivialisation of the torsor;
- the two are Fourier dual (§8, item 4).

*Proof.* ℂ[V₄] is commutative, so its Wedderburn blocks all have dimension 1, and their unit elements are the e_χ (Serre 1977, §2.6, Theorem 8). ∎

A chain of subgroups V₄ ⊇ H₁ ⊇ 1 gives the graded breaking of the draft's §7.2. The alphabets run from N, through one of R/Y, W/S or M/K, to the bases (`grammar::quotient`). The residual symmetry acting on Σ/H is V₄/H, and on the bases it is V₄ itself. Complementarity is the element with χ_RY(κ) = −1, χ_WS(κ) = +1 and χ_MK(κ) = −1:
- it exchanges purines and pyrimidines;
- it keeps the number of hydrogen bonds;
- it exchanges amino and keto.

### 3.5 Codons

Codons are a torsor of V₄³ = (C₂)⁶, acting by exclusive or on their six bits (§6).

## 4. Strands

The group of order 8, V₄ × ⟨ρ⟩, acts on words: ρ, reversal, has sign −1. Its subgroup {1, κ, ρ, κρ}, the strand group, contains the reverse complement κρ, which carries a word to the other strand read in its own direction.

By Burnside's lemma with the fixed-word counts of Definition 1.1, the canonical k-mers number:
- under ⟨κρ⟩: (4^k + 4^{k/2})/2 for even k (the reverse-complement palindromes are fixed), and 4^k/2 for odd k;
- under the strand group: (4^k + 4^⌈k/2⌉ + [k even] 4^{k/2})/4.

The checks verify these against enumeration up to k = 6, and against the closed forms up to k = 8.

**Which symmetries a genome keeps.** The total variation distance TV(f, g·f) between a genome's k-mer frequencies f (one strand) and their image under g is 0 exactly when the statistics are invariant under g. For *Escherichia coli* K-12 and 4-mers:

| Element | TV |
|---|---|
| reverse complement κρ | 0.0054 |
| reversal ρ | 0.167 |
| complement κ | 0.168 |
| σρ | 0.195 |
| transition τ | 0.195 |
| τρ | 0.240 |
| σ | 0.240 |

Of the eight elements, only the identity and the reverse complement nearly keep the statistics. This is Chargaff's second parity rule (Rudner, Karkas and Chargaff 1968, PNAS 60, 921–922, for the bases), found here as a group: the empirical symmetry group of the genome's statistics is {1, κρ} ≅ C₂.

The same fact reads through the group's characters (§7). The 4-mer frequency vector splits into eight orthogonal channels, one per character of the group of order 8 (`grammar::kmer_channels`):
- the channels odd under the reverse complement hold 5.0 × 10⁻⁵ of its energy;
- those odd under the complement alone hold 4.6 × 10⁻², some 900 times more.

This is the draft's question Q2, group discovery, answered for genomic k-mer statistics by measurement. As the draft's own Q6 anticipates, it is an empirical fact, not a structural one. The complement alone and reversal alone, which the draft's G_DNA lists as separate factors, are *not* symmetries of the statistics; only their product is.

## 5. Hairpins and Their Symmetry

The hairpin language H(s, l₁, l₂) is the set of words u v u′ with:
- |u| = |u′| ≥ s;
- the pairs (u_i, u′_{|u|+1−i}) all allowed;
- l₁ ≤ |v| ≤ l₂.

The grammar of `grammar::hairpin_grammar` emits stem and loop from the outside in (its rules are in the header), so that every rule is symmetric under reversal. Its language is checked against the definition on every word up to length 8.

**Proposition 5.1.**
- With Watson–Crick pairs only, the grammar is invariant under the strand group: (κa, κb) is a Watson–Crick pair when (a, b) is, and the reversal of a pair rule is the reversed pair's rule.
- With the wobble pairs G·T and T·G, it is invariant under ⟨ρ⟩ but not under κ, since κ(G·T) = C·A pairs in neither.

*Proof.* By inspection of the rules; `formal::invariant` checks both. ∎

The chemistry of the G·U pair thus breaks the complement symmetry of the hairpin grammar, and the invariance test records which part of the symmetry survives. This is a small, exact instance of the draft's "breaking is not total" (§6.1).

**Hairpins in real tRNAs.** Of *Escherichia coli*'s 86 tRNA genes, 98.8% contain a stem of at least 5 Watson–Crick pairs closing a loop of 3–9 bases, recognised by Earley's parser on the flanked grammar. Only 40.7% of random windows of the same lengths do. The anticodon and acceptor stems of the cloverleaf are what the grammar finds.

## 6. The Genetic Code on (C₂)⁶

A genetic code is a function aa on the codon torsor of V₄³. Its symmetries are the m with aa(c ⊕ m) = aa(c) for every c; they form a subgroup, the stabiliser of aa. The invariance profile I(m) = #{c : aa(c ⊕ m) = aa(c)} measures how far each element is a symmetry.

**Proposition 6.1** (from NCBI's tables, `sequence/genetic_code.hpp`).
- The standard code (table 1) has trivial stabiliser. The third-position transition τ₃ keeps the amino acid of 60 codons of 64, failing at exactly two pairs: ATA/ATG (Ile/Met) and TGA/TGG (stop/Trp).
- The vertebrate mitochondrial code (table 2), which reads ATA as Met and TGA as Trp, has stabiliser {1, τ₃}.

*Proof.* By evaluation of the tables; the checks do it. ∎

The draft's "universal grammar, different breakings" (§9.3, Q4) has here a small exact instance:
- the two codes are two functions on the same group;
- the code of vertebrate mitochondria keeps one more symmetry than the standard code;
- the wobble of the third position is invariance under a subgroup.

I(m) is the autocorrelation of the amino acids' indicator functions on the group, so its Walsh–Hadamard transform is their power spectrum: P(u) = Σ_a f̂_a(u)² = Σ_m (−1)^{u·m} I(m) (Terras 1999, *Fourier Analysis on Finite Groups and Applications*, ch. 2). The checks verify this identity exactly, in integers, for both codes, with morphe's transform.

**Proposition 6.2.** A code is invariant under m exactly when its spectrum vanishes on every character u with u·m odd.

*Proof.* f_a(c ⊕ m) = f_a(c) for every a and c if and only if f̂_a(u) = (−1)^{u·m} f̂_a(u) for every u, that is, f̂_a(u) = 0 when u·m is odd. ∎

The vertebrate mitochondrial code's spectrum vanishes on every character odd under τ₃, as Propositions 6.1 and 6.2 together require; this is checked.

Of the standard code's non-constant spectral energy:
- 73.4% lies in characters of the first two bases;
- 14.0% lies in characters of the first base alone.

This is the wobble of the third position, in Fourier terms.

## 7. The Wedderburn Decomposition, and What It Parallelises

**Theorem 7.1** (Wedderburn; Serre 1977, §2.6 and §6.2). ℂ[G] ≅ ⊕_{χ∈Ĝ} M_{χ(1)}(ℂ), with Σ_χ χ(1)² = |G|. The central primitive idempotents e_χ = (χ(1)/|G|) Σ_g χ(g⁻¹) g are the units of the blocks.

morphe computes the character tables exactly, by the Dixon–Schneider method with values in cyclotomic fields, together with the e_χ and the isotypic components of any element (`algebra/characters.hpp`, `group_algebra.hpp`). Over ℚ, which is what exact arithmetic sees, ℚ[C_n] ≅ ⊕_{d|n} ℚ(ζ_d): the blocks are cyclotomic fields, which split further over ℂ.

**Proposition 7.2.** Let Γ be G-invariant with fixed nonterminals and invariant weights. Then its weight function on Σⁿ, f(w) = Σ over derivations of w, lies in the trivial isotypic component of the permutation module ℂ[Σⁿ]: f · e_triv = f, and f · e_χ = 0 for χ ≠ 1.

*Proof.* f is constant on G-orbits by Proposition 1.3 with weights, and the trivial isotypic component of a permutation module is the orbit-constant functions. ∎

So the Wedderburn blocks do not split the recognition or weighting of a G-invariant language into |Ĝ| independent problems: everything lies in one block. They parallelise what is genuinely linear and equivariant.
- **Equivariant linear maps.** By Schur's lemma, an equivariant linear map of a G-module is block-diagonal in its isotypic decomposition (Serre 1977, §2.2). These are the draft's §8 encoders and layers, where block parallelism is exact.
- **Words of length n.** ℂ[G]^{⊗n} decomposes into |Ĝ|ⁿ blocks, one per tuple of irreducible representations, not |Ĝ| (§8, item 7).

**The channels, realised.** This is the idea's direction made exact, on a linear object: a genome's k-mer frequency vector f.
- **The projection.** The central idempotents project f onto the isotypic components of the permutation representation, P_χ f = (χ(1)/|G|) Σ_g χ(g)‾ g·f. These components are orthogonal and sum to f.
- **What a channel measures.** The energy ‖P_χ f‖² of each channel says how much of the statistics transforms by χ.
- **On *E. coli*,** under the group of order 8 (§4): the eight channels' energies sum to ‖f‖² within 8 × 10⁻¹⁶, and the channels odd under the reverse complement are nearly empty (`grammar/channels.hpp`, with characters from morphe's exact table).

**The algebra does not determine the group.** As an algebra, ℂ[G] determines only the multiset of degrees {χ(1)}. For example, ℂ[C₄] ≅ ℂ[C₂ × C₂] ≅ ℂ⁴, and ℂ[D₄] ≅ ℂ[Q₈] ≅ ℂ⁴ ⊕ M₂(ℂ). morphe's character tables of these groups check it. G is recovered as the group-like elements {x : Δx = x ⊗ x} of the *Hopf* algebra ℂ[G], which needs the coproduct Δ(g) = g ⊗ g. The draft's §4.2 says this; its Remark 2.1 and §9.2, which speak of the algebra alone, need the qualification (§8, item 3).

## 8. Amendments to the Draft

The draft is the user's, and these are proposals for its next version, each with the section it concerns and the reason.

1. **§3.2.1 and §9.3, the orbit count.** By Burnside's lemma, |R₁/G| ≥ |R₁|/|G|, not ≤. Equality holds when the action on the rules is free. The compression factor is at most |G|.
2. **Definition 2.2, matrix units.** The E^ρ_ij with i ≠ j are matrix units, nilpotent (E_ij² = 0), not idempotents. The minimal (primitive) idempotents are the E^ρ_ii and their conjugates. The relation e_{ρ,ij} e_{σ,kl} = δ_ρσ δ_jk e_{ρ,il} is that of matrix units, and is right.
3. **Remark 2.1, §4.2 and §9.2, G and ℂ[G].** As algebras they do not determine each other (ℂ[C₄] ≅ ℂ[C₂ × C₂]). As Hopf algebras they do, by the group-like elements, as §4.2 says.
4. **Example 2.1, §2.3 and §7.2, purines and pyrimidines.** C₂ acting by complementarity has orbits {A, T} and {C, G}; its idempotents e± split functions on the bases into complement-even and complement-odd parts. The purine–pyrimidine classes are not these. They are the cosets of ⟨τ⟩, the kernel of a character of the Klein group generated by complementarity and transition (Proposition 3.3).
   - **Derivation.** The classification is then derived from two chemical facts: base pairing, and that transitions keep the ring class.
   - **Letters and idempotents.** Letters correspond to group elements through a trivialisation of the torsor; the minimal idempotents are the characters, that is, the classifications (Proposition 3.4).
   - **The hydrogen bonds.** The draft's B₂ step (A(2), G(3), C(3), T(2)) is the weak–strong character.
5. **§5.2.2, the scan rule.** Matching a terminal by its orbit alone accepts π⁻¹(L(Γ/G)), which strictly contains L(Γ) in general (Proposition 2.1). The parser must bind each rule's group element, as a coset, and narrow it at each scan (§2.2).
6. **Proposition 5.1, complexity.** The equivariant chart is the plain chart up to the action (Proposition 2.3), so the worst case stays O(n³ |R|). What is compressed is the rule store, by |R|/|R/G| ≤ |G|, and the predictions; the savings on the chart are measured, not guaranteed.
7. **§3.2.1, §5.2.4 and §9.3, block parallelism.**
   - **Recognition.** A G-invariant language's weight function lies in the trivial block (Proposition 7.2), so the Wedderburn blocks do not give independent parsing problems for it. "No rule in R^ρ₁ can produce or consume elements of R^{ρ′}₁" needs rules that are linear maps of ℂ[G]-modules, which rewriting rules on words are not.
   - **Words.** Over words, the blocks of ℂ[G]^{⊗n} number |Ĝ|ⁿ.
   - **Where it is exact.** Block parallelism is exact for equivariant linear maps (Schur's lemma), the draft's §8.
8. **§6.2, the vacuum manifold.** The breakings of a G-torsor at its base points form a torsor themselves. Two of them differ by a right translation, an automorphism of the G-set G, not by an outer automorphism of G (`formal::torsor_breaking`, checked).
9. **§7.1, G_DNA.** Complement acts on letters and reversal on words by anti-automorphism, so ℤ₂ × ℤ₂ acts on Σ* as in Definition 1.1. ℤ_n (tandem shifts) and S_k (exon permutations), however, act on positions or on segments, not on letters. The product acts on segmented words, and the Wedderburn decomposition of ℂ[G_DNA] does not act on the nucleotide alphabet.
   - **The letters' group.** The group of the letters is the Klein group V₄ (§3).
   - **The statistics' group.** The group a genome's statistics keep is {1, κρ} (§4): of complement and reversal, the statistics keep only their product, not either alone.
10. **§1.2, the complexity wall.** Orbit compression changes the cost of context-sensitive membership by a factor of at most |G| (Proposition 2.4). The problem stays PSPACE-complete; the wall is lowered, not restructured.
11. **§3.3, two phases.** Making it exact, phase I is a grammar over an alphabet Γ (for a torsor, the group's elements), and phase II a letter-to-letter map β. Then L = β(L(Γ)), and β(Γ) is a context-free grammar for it (`formal::image`, checked). "x ∈ ℂ[G]*" becomes a weighted language over Γ, a formal power series, when superpositions are meant.

## 9. The Draft's Open Questions

- **Q1, grammar induction.** Orbit-tied inside–outside is an EM algorithm (Proposition 2.6); it is run on real tRNAs in the `rna` module. Convergence to a global maximum is not claimed; EM gives monotonicity and convergence to a stationary point under regularity conditions (Wu 1983, Ann. Statist. 11, 95–103).
- **Q2, group discovery.**
  - The symmetries of a genome's k-mer statistics are measured element by element (§4): {1, κρ} for *E. coli*.
  - Tied and untied models of RNA pairing are compared by likelihood in the `rna` module.
  - Searching larger groups, for instance subgroups of the hyperoctahedral group acting on k-mers, is open.
- **Q3, treewidth.** Open.
- **Q4, cross-species universality.** Genetic codes are functions on one group, with different stabilisers (§6). Whether organisms differ by breakings in a stronger sense is open.
- **Q5, non-abelian groups.** morphe now computes the characters, Wedderburn data and Fourier transforms of non-abelian groups (S_n by Young's seminormal form). An exon-permutation grammar with S_k acting on segments is the natural first case; open.
- **Q6, the origin of G.** Open, and partly answered by §4: the symmetry of the statistics is empirical.

## 10. References

- Cornish-Bowden, A. (1985). Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984. *Nucleic Acids Research* 13, 3021–3030.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *J. R. Statist. Soc. B* 39, 1–38.
- Earley, J. (1970). An efficient context-free parsing algorithm. *Comm. ACM* 13, 94–102.
- Hopcroft, J. E. and Ullman, J. D. (1979). *Introduction to Automata Theory, Languages, and Computation*. Addison-Wesley.
- Huang, W. (2026). Towards a group context-sensitive grammar for genome analysis. Draft, 15 March 2026.
- Rotman, J. J. (1995). *An Introduction to the Theory of Groups*, 4th ed. Springer.
- Rudner, R., Karkas, J. D. and Chargaff, E. (1968). *PNAS* 60, 921–922 (the second parity rule, for the bases of the separated strands of *B. subtilis* DNA).
- Serre, J.-P. (1977). *Linear Representations of Finite Groups*. Springer.
- Terras, A. (1999). *Fourier Analysis on Finite Groups and Applications*. Cambridge.
- Wu, C. F. J. (1983). On the convergence properties of the EM algorithm. *Ann. Statist.* 11, 95–103.
