# RNA Secondary Structure by Grammars

*Status: draft (2026-10-25). Serves learning, use and research.* Secondary structures and their number; the one recursion that maximises additive weights over them, and Nussinov's maximal pairing as its simplest case; stochastic context-free grammars whose derivations are the structures, each derived once; folding by the best parse, pair probabilities by inside and outside, and the maximum expected accuracy structure; training from sequences alone by expectation maximisation, with the probabilities tied along the rule orbits of a group (the equivariant EM of the user's question Q1), the tied M-step proved exact; and what *Escherichia coli*'s tRNAs say when the grammars are trained on them, including whether base pairing is symmetric under reversal (the user's question Q2, asked of RNA). The code is the module [`rna`](../rna/): [`structure.hpp`](../rna/include/genomics/rna/structure.hpp), [`pairing.hpp`](../rna/include/genomics/rna/pairing.hpp), [`scfg.hpp`](../rna/include/genomics/rna/scfg.hpp), [`training.hpp`](../rna/include/genomics/rna/training.hpp), [`symmetry.hpp`](../rna/include/genomics/rna/symmetry.hpp). The grammars are morphe's (`formal/grammar.hpp`), and so are their inside, outside, expected counts, best parse, derivations and rule orbits (`formal/inside_outside.hpp`, `formal/equivariant.hpp`); the symmetries of the bases and strands are the [group grammars](group-grammars.md)'.

## 1. Secondary Structures

**Definition 1.1 (secondary structure).** A *secondary structure* on the positions [n] = {0, …, n − 1} is a set y of pairs (i, j), i < j, such that

- every position lies in at most one pair (y is a partial matching), and
- no two pairs cross: there are no (i, j), (k, l) ∈ y with i < k < j < l.

It has *minimum loop* m when j − i − 1 ≥ m for every pair. Its *pair table* gives each position its partner, or −1 (`Structure::partner`). A sequence x ∈ Σⁿ, Σ = {A, C, G, U}, with a structure y is a *folded sequence* (x, y); the pair (i, j) is *canonical* when (x_i, x_j) is a Watson–Crick pair (A·U, U·A, G·C, C·G) or a wobble pair (G·U, U·G). Which pairs are canonical is chemistry, stated once (`canonical`); the grammars of §3 take it as a choice of which pairs have rules, not as an axiom.

*Remark.* Non-crossing is what makes structures context-free: the pairs nest like brackets, and pseudoknots, which cross, are outside every grammar here.

**Definition 1.2 (dot-bracket).** The dot-bracket word of y writes '(' at i and ')' at j for each pair (i, j), and '.' at every unpaired position.

**Proposition 1.3.** The dot-bracket word is a bijection from the secondary structures on [n] onto the words of length n over {(, ), .} whose brackets are balanced.

*Proof.* Because pairs nest, reading left to right every ')' closes the most recently opened pair still open, so the brackets are balanced. Conversely, given a balanced word, a non-crossing matching of its brackets must pair each ')' with the nearest unmatched '(' to its left: any other choice would enclose an opening bracket whose partner lies outside, a crossing. So the matching exists and is unique. ∎

`parse_dot_bracket` is this matching by a stack; `dot_bracket` the inverse.

**Proposition 1.4 (the number of structures).** Let S_m(n) be the number of secondary structures on [n] with minimum loop m, every pair allowed. Then S_m(0) = 1 and

S_m(n) = S_m(n − 1) + Σ_{k=0}^{n−m−2} S_m(k) · S_m(n − k − 2).

*Proof.* Either n − 1 is unpaired, leaving a structure on [n − 1], or it is paired with some k. No pair crosses (k, n − 1), so the other pairs lie in [0, k) or in (k, n − 1), independently; the second interval has n − k − 2 ≥ m positions. ∎

**Proposition 1.5.** S_0(n) = Σ_k C(n, 2k) · Cat_k, with Cat_k = C(2k, k)/(k + 1) the Catalan number.

*Proof.* Choose the 2k paired positions; by Proposition 1.3 the non-crossing perfect matchings of 2k points in a line are the balanced bracket words of length 2k, of which there are Cat_k. ∎

These are the Motzkin numbers; `check_genomics_rna_structure` holds `structure_count` to the enumeration for m ≤ 3, n ≤ 14, and to Proposition 1.5 for m = 0, in morphe's exact integers.

## 2. One Recursion for Additive Weights

Many objectives over structures are additive: a weight u_i for each unpaired position and w(i, j) for each pair.

**Proposition 2.1.** For weights u_i (i ∈ [n]) and w(i, j) (or "not allowed"), let M(i, i) = 0 and, over half-open intervals [i, j),

M(i, j) = max( u_i + M(i + 1, j),  max_{k ∈ [i + m + 1, j), w(i, k) allowed} w(i, k) + M(i + 1, k) + M(k + 1, j) ).

Then M(0, n) = max_y ( Σ_{i unpaired} u_i + Σ_{(i, j) ∈ y} w(i, j) ) over the structures with minimum loop m and allowed pairs.

*Proof.* Every structure on [i, j) falls into exactly one case by the state of position i: unpaired, leaving a structure on [i + 1, j); or paired with exactly one k, and then, since no pair crosses (i, k), the remaining pairs lie in [i + 1, k) or [k + 1, j) independently. The objective is additive over this decomposition, so the maximum over the structures of [i, j) is the maximum over the cases of the sums of the sub-maxima, by induction on j − i. ∎

The traceback (`best_structure`) recomputes the cases and prefers the unpaired one, then the least k, on ties; time O(n³), space O(n²).

*Remark (one sum in several semirings).* Because the decomposition is unique, the same recursion with (max, +) replaced by (+, ×) sums Π u_i Π w(i, j) over the structures, a partition function, and with (+, ×) on the integers counts them (Proposition 1.4 is its case u = w = 1). A structure grammar (§3) is this recursion written as rules: its derivations decompose structures uniquely, which is what lets morphe's inside–outside, generic in the semiring, compute every quantity of §4 and §5.

### 2.1 Nussinov's Maximal Pairing

Nussinov, Pieczenik, Griggs and Kleitman (1978) and Nussinov and Jacobson (1980) fold a sequence by maximising its number of canonical pairs: u_i = 0, w(i, j) = 1 when (x_i, x_j) is canonical and not allowed otherwise, with a minimum loop (`nussinov`). The check holds its count and its traceback to the best of every structure of 300 random sequences of up to 14 bases, minimum loops 0 to 3, with and without G·U.

*Caveat.* Counting pairs ignores everything that makes a structure stable (the stacking of helices, the cost of loops), so Nussinov's structures are rarely the real ones; the method is here as the simplest exact case of the recursion, and as the chemistry of canonical pairs made explicit.

## 3. Structure Grammars

### 3.1 The Grammars

**Definition 3.1 (the grammars).** Let P ⊆ Σ × Σ be a set of pairs. *Knudsen and Hein's grammar* (Knudsen and Hein 1999; the grammar G6 of Dowell and Eddy 2004) on P has nonterminals S (the start), L and F and the rules

- S → L S | L,
- L → a for every a ∈ Σ,
- L → a F b and F → a F b for every (a, b) ∈ P,
- F → L S.

Its *reversible* presentation replaces the list rules by

- S → L | L L | L S L,  F → L L | L S L.

Either is *without lonely pairs* when L → a F b is replaced by L → a H b and H → a F b for every (a, b) ∈ P, with a new nonterminal H. In the code, `structure_grammar(shape, lonely_pairs, pairs)`; `knudsen_hein()` is the published grammar with P = Σ × Σ (39 rules, 66 in morphe's binary form), and the trained grammar of §6 is the reversible one without lonely pairs on the six canonical pairs (27 rules, 42 binary). Each emission is its own rule, so the probability of the rule L → a F b is the joint probability of opening a pair and of its bases being (a, b).

**Definition 3.2 (the structure of a derivation).** A derivation tree of x in one of these grammars has the structure y whose pairs are the spans [i, j] of its pair-rule nodes (L → a F b, L → a H b, H → a F b, F → a F b): a pair rule emits its first terminal at i and its last at j.

### 3.2 Each Structure Derived Once

**Theorem 3.3.** For each grammar of Definition 3.1, the map from derivation trees to folded sequences (x, y) is a bijection onto the set D of (x, y) with |x| ≥ 1 and

1. every pair encloses at least two positions (j − i − 1 ≥ 2);
2. every pair has (x_i, x_j) ∈ P;
3. without lonely pairs: every pair (i, j) is stacked on a neighbour, (i + 1, j − 1) ∈ y or (i − 1, j + 1) ∈ y.

*Proof.* Call a *component* of an interval a maximal unit of its top level: an unpaired position, or a pair together with everything inside it. The top level of any interval is a sequence of components.

*Lists.* S derives exactly the non-empty sequences of L's, each in one way. In Knudsen and Hein's shape S ⇒* L^k only by k − 1 uses of S → L S and one of S → L. In the reversible shape, by induction on k: L¹ only by S → L, L² only by S → L L (S → L S L yields at least three), and L^k, k ≥ 3, only by S → L S L with S ⇒* L^{k−2}.

*Components.* L derives one component: an unpaired base (L → a), or a pair (i, j) with (x_i, x_j) ∈ P whose interior F derives (L → a F b); without lonely pairs, L → a H b and H → a′ F b′ force the pair (i + 1, j − 1) inside it.

*Interiors.* F derives, over an interval I, either a pair spanning I whose interior F derives (F → a F b), or a list of at least two components (Knudsen and Hein: F → L S; reversible: F → L L | L S L, two components or two plus a non-empty list). The two cases exclude each other: in the first the ends of I pair with each other, so I is one component; in the second it has at least two. Every F derives at least two positions (two components, or a pair around at least two), which is condition 1, and the pair rules carry only pairs of P, which is condition 2. Without lonely pairs, a pair spanning an F-interval is stacked on the pair around it, and every pair that begins a component has its inner neighbour by the H rule: condition 3. Conversely, given (x, y) ∈ D, the decomposition of each interval into components, and of each interior into "one spanning pair" or "two or more components", is determined by y, and gives the rules of the derivation one by one; condition 3 provides the inner pair the H rule needs. By induction on the length of the interval, every (x, y) ∈ D has exactly one derivation (`derivation` builds it), and every derivation yields an element of D. ∎

**Corollary 3.4.** With every rule weighted 1 in morphe's counting semiring, the inside value of x is the number of structures y with (x, y) ∈ D; for P = Σ × Σ with lonely pairs allowed it is S_2(|x|).

`check_genomics_rna_scfg` holds this for five grammars on words of 1 to 12 bases, against the enumeration filtered by `derivable`, and against Proposition 1.4.

*Remark.* Unambiguity is what makes the probabilities mean what they say: P(x) = Σ_y P(x, y) sums over structures, not over derivations of one structure counted twice, and the best derivation is the best structure. The distinction between ambiguous and unambiguous structure grammars is the one Dowell and Eddy (2004) draw among lightweight grammars; here it is proved for the grammars used, not assumed.

### 3.3 Knudsen and Hein's Grammar inside the Reversible One

**Proposition 3.5.** Let θ be probabilities for Knudsen and Hein's grammar, p = θ(S → L S), q = θ(F → L S). Give the reversible grammar (same P, same choice on lonely pairs) the probabilities θ′(S → L) = 1 − p, θ′(S → L L) = p(1 − p), θ′(S → L S L) = p², θ′(F → L L) = q(1 − p), θ′(F → L S L) = qp, and θ's emissions. Then θ′ sums to 1 on each left side, and P_θ′(x, y) = P_θ(x, y) for every (x, y).

*Proof.* By Theorem 3.3 both grammars derive (x, y) once, with the same emission and pair rules; they differ only in the rules that build lists. A top-level list of k components costs p^{k−1}(1 − p) in Knudsen and Hein's grammar. In the reversible one it costs 1 − p for k = 1, p(1 − p) for k = 2, and p² times the cost of the middle k − 2 for k ≥ 3, which is p^{k−1}(1 − p) by induction. An interior list of k ≥ 2 components costs q · p^{k−2}(1 − p) in the first (F → L S, then a list of k − 1); in the second q(1 − p) for k = 2 and qp · p^{k−3}(1 − p) = q p^{k−2}(1 − p) for k ≥ 3. The sums: (1 − p) + p(1 − p) + p² = 1, and q(1 − p) + qp = q. ∎

So Knudsen and Hein's family is a sub-model of the reversible one, which has two more free parameters (the reversible lists' lengths have a free parity). `as_reversible` is the map; the check holds the log-likelihoods of random sequences, and the probabilities of every structure of short ones, equal to 10⁻¹² relative.

### 3.4 Symmetries

Let W = V₄ × ⟨ρ⟩ be the group of order 8 of [group grammars](group-grammars.md) §4: the Klein group of the bases (complement κ: A ↔ U, C ↔ G; transition τ; their product σ) acting on letters, and reversal ρ. An element g acts on a rule X → α of a grammar by morphe's rule (`GrammarAction`, nonterminals fixed): the letters of α are permuted and, when g reverses (sign −1), α is reversed. On a pair rule X → a Y b this is X → (g a) Y (g b), or X → (g b) Y (g a) when g reverses; write g·(a, b) for that pair.

**Proposition 3.6.** (a) Knudsen and Hein's rule set is invariant under g exactly when g does not reverse and g·P = P. (b) The reversible rule set, with or without lonely pairs, is invariant under g exactly when g·P = P.

*Proof.* Rules without terminals are fixed by every letter permutation, and L → a is carried to L → g a, which is a rule. (a) If g reverses, S → L S goes to S → S L, which is no rule. Otherwise the pair rules go to pair rules exactly when g·P = P. (b) The list rules S → L | L L | L S L and F → L L | L S L are palindromes of nonterminals and are fixed by reversal; the pair rules are carried to pair rules exactly when g·P = P. ∎

On the canonical pairs: ρ keeps them (A·U ↔ U·A, G·C ↔ C·G, G·U ↔ U·G); κ carries G·U to C·A and κρ carries it to A·C, so neither keeps them, and the closure of the canonical pairs under κ, under κρ and under the strand group {1, κ, ρ, κρ} is the canonical pairs with A·C and C·A (`closure`). The complement κ maps the Watson–Crick pairs to themselves, reversed (κ(A, U) = (U, A)), which is why it coincides with ρ on them and departs from it only at G·U.

## 4. Folding

**Viterbi.** By Theorem 3.3 the most probable derivation is the most probable structure: `viterbi` is morphe's `best_parse` over (max, +) on log θ, its tree read back as pairs.

### 4.1 Pair Probabilities

**Proposition 4.1.** For i < j,

P(i·j | x) = Σ_{pair rules r = X → x_i Y x_j} O(X, i, j + 1) · θ_r · I(Y, i + 1, j) / P(x),

with I and O the inside and outside values of morphe's charts over half-open spans.

*Proof.* By Theorem 3.3, (i, j) ∈ y exactly when the derivation of (x, y) has a pair-rule node spanning [i, j], and it has at most one, since a position is in one pair. The derivations with node r over [i, j + 1) weigh, in sum, the derivations around the span (the outside value of X there) times θ_r times the derivations of Y over the interior (its terminals weigh 1): the inside–outside decomposition (Lari and Young 1990; Goodman 1999), which morphe's `outside` computes in any commutative semiring. Summing over the pair rules that can emit (x_i, x_j) and dividing by P(x) gives the posterior. ∎

This is the grammar's analogue of McCaskill's (1990) base-pair probabilities, with rule probabilities in place of Boltzmann weights and P(x) in place of the partition function (`pair_probabilities`). The check holds it to the enumerated marginals to 10⁻¹² on sixty words, for five grammars.

### 4.2 Maximum Expected Accuracy

**Definition 4.2.** With p_ij = P(i·j | x) and q_i = 1 − Σ_j p_ij, the *expected accuracy* of a structure ŷ at weight γ > 0 is A_γ(ŷ) = Σ_{i unpaired in ŷ} q_i + 2γ Σ_{(i, j) ∈ ŷ} p_ij (Do, Woods and Batzoglou 2006, who weigh the pairs by γ to trade specificity for sensitivity).

**Proposition 4.3.** A_1(ŷ) is the expected number of positions at which ŷ agrees with a structure Y drawn from P(· | x), agreement at i meaning that both leave i unpaired or both pair it with the same partner.

*Proof.* By linearity of expectation: a position unpaired in ŷ agrees with probability q_i; a pair (i, j) of ŷ makes i agree with probability p_ij and j with probability p_ij. ∎

The MEA structure maximises A_γ; by Proposition 2.1 with u_i = q_i and w(i, j) = 2γ p_ij on the pairs with p_ij > 0, it is one recursion (`maximum_expected_accuracy`; `mea_fold` uses minimum loop 2, as the grammars do).

**Proposition 4.4.** Let γ₁ < γ₂ and ŷ₁, ŷ₂ maximise A_{γ₁} and A_{γ₂}. With B(ŷ) = Σ_{(i, j) ∈ ŷ} 2 p_ij (the expected number of correct paired positions) and U(ŷ) = Σ_{i unpaired} q_i, B(ŷ₂) ≥ B(ŷ₁) and U(ŷ₂) ≤ U(ŷ₁).

*Proof.* Optimality gives U₁ + γ₁ B₁ ≥ U₂ + γ₁ B₂ and U₂ + γ₂ B₂ ≥ U₁ + γ₂ B₁. Adding, (γ₂ − γ₁)(B₂ − B₁) ≥ 0; then the first inequality gives U₁ − U₂ ≥ γ₁(B₂ − B₁) ≥ 0. ∎

The check holds the MEA optimum to the best of every structure (200 instances, six γ) and B to Proposition 4.4.

## 5. Training from Sequences Alone

### 5.1 Expectation Maximisation

Let x⁽¹⁾, …, x⁽ᴺ⁾ be sequences, ℓ(θ) = Σ_k log P_θ(x⁽ᵏ⁾), and n_r(x, y) the number of uses of rule r in the derivation of (x, y), so that log P_θ(x, y) = Σ_r n_r log θ_r. For the current θ, the *expected counts* are C_r(θ) = Σ_k E_θ[n_r | x⁽ᵏ⁾], which morphe's `expected_counts` computes by inside and outside (Lari and Young 1990), and

Q(θ′ | θ) = Σ_k E_θ[log P_θ′(x⁽ᵏ⁾, Y) | x⁽ᵏ⁾] = Σ_r C_r(θ) log θ′_r.

**Lemma 5.1 (Dempster, Laird and Rubin 1977).** ℓ(θ′) − ℓ(θ) ≥ Q(θ′ | θ) − Q(θ | θ).

*Proof.* For one sequence, log P_θ′(x) − log P_θ(x) = log Σ_y P_θ(y | x) · P_θ′(x, y)/P_θ(x, y) ≥ Σ_y P_θ(y | x) log (P_θ′(x, y)/P_θ(x, y)) by Jensen's inequality (log is concave), and the right side is the difference of the Q terms. Sum over the sequences. ∎

So any θ′ with Q(θ′ | θ) ≥ Q(θ | θ) is at least as likely as θ. The untied M-step θ_r = C_r / C_A, with C_A = Σ_{r ∈ R_A} C_r over the rules R_A of the left side A, maximises Q by Gibbs' inequality.

### 5.2 The Equivariant EM

Let a finite group G act on the grammar fixing every nonterminal, its rule set invariant (Proposition 3.6), and let R_A/G be the orbits of the rules of A. The *tied model* is Θ_G = {θ : θ constant on every orbit, Σ_{r ∈ R_A} θ_r = 1 for every A}.

**Theorem 5.2 (the tied M-step).** For counts C_r ≥ 0 with C_A > 0, the maximum of Q(θ) = Σ_r C_r log θ_r over Θ_G is attained at

θ_r = C_O / (|O| · C_A)  for r ∈ O ∈ R_A/G,  where C_O = Σ_{r ∈ O} C_r,

uniquely on every orbit with C_O > 0 (and with θ_O = 0 on the others).

*Proof.* For θ ∈ Θ_G write θ_O for its value on O and φ_O = |O| θ_O, the probability of the orbit. The map θ ↦ (φ_O) is a bijection of Θ_G onto the product over A of the simplices on R_A/G, since Σ_{O ∈ R_A/G} φ_O = Σ_{r ∈ R_A} θ_r = 1. In these coordinates

Q(θ) = Σ_A Σ_{O ∈ R_A/G} C_O log(φ_O / |O|) = Σ_A Σ_O C_O log φ_O − Σ_A Σ_O C_O log |O|,

the last sum constant. For each A, (C_O / C_A)_O is a probability vector, and Gibbs' inequality gives Σ_O C_O log φ_O ≤ Σ_O C_O log(C_O / C_A), with equality exactly when φ_O = C_O / C_A on every O with C_O > 0 (the remaining φ_O then sum to zero). Hence θ_O = φ_O / |O| = C_O / (|O| C_A). ∎

**Corollary 5.3.** The tied EM, θᵗ⁺¹ = argmax_{Θ_G} Q(· | θᵗ) started in Θ_G, never decreases the likelihood, and every iterate is tied.

*Proof.* θᵗ ∈ Θ_G, so Q(θᵗ⁺¹ | θᵗ) ≥ Q(θᵗ | θᵗ); apply Lemma 5.1. The iterates are tied by Theorem 5.2. ∎

The tied EM is thus the ordinary EM of the model whose parameters are the orbit probabilities φ_O, in which P_φ(x, y) = Π_O (φ_O/|O|)^{n_O(x, y)}. `train` projects its start onto Θ_G (orbit averages, `project`) and iterates `expectation` and `maximisation`; `check_genomics_rna_training` holds the tied M-step to its stationarity conditions C_O/(|O| θ_O) = C_A, to every random tied θ and every tied θ near it (four groups), and the iteration to monotonicity.

**Proposition 5.4 (invariance).** If θ ∈ Θ_G then P_θ(g·x) = P_θ(x) for every g ∈ G and every sequence x, and P_θ(g·x, g·y) = P_θ(x, y), where g·y is y mirrored when g reverses.

*Proof.* Let g act on derivation trees by replacing each node's rule r by g·r and, when g reverses, reversing the order of every node's children. The image is a derivation tree (g·r is a rule, and its right side is the children in their new order), its yield is g·x (letters permuted, and the word reversed when every node's children are), and the map is a bijection with inverse g⁻¹. It preserves weights, Π θ_{g·r} = Π θ_r, since θ is constant on orbits. Sum over derivations; the second statement follows from Theorem 3.3. ∎

So a tied model is a G-invariant distribution of folded sequences. The converse fails in general (an untied θ may give an invariant distribution), so the test below is a test of the tie, the parametric form of the symmetry.

### 5.3 The Test of a Tie

The untied model has Σ_A (|R_A| − 1) free parameters and the tied one Σ_A (|R_A/G| − 1), so the tie removes d = Σ_A (|R_A| − |R_A/G|), where |R_A/G| is Burnside's count (morphe's `rule_orbits` checks it). With ℓ̂ the maximised log-likelihoods on the same sequences, Λ = 2(ℓ̂_untied − ℓ̂_tied) is asymptotically χ²_d under the tie (Wilks 1938), and `likelihood_ratio` gives its tail by calx's incomplete gamma function. On the canonical grammar without lonely pairs, the reversal tie pairs A·U with U·A, G·C with C·G and G·U with U·G on each of L, H and F: d = 9.

*Caveat.* Wilks's theorem assumes maximised likelihoods and parameters in the interior. EM finds local maxima, and its fit depends on its start (§6.2); the untied fit is started from the tied one, so that Λ ≥ 0, and is compared with the untied fit from the uniform start, the larger taken. Some fitted probabilities lie at or near zero, on the boundary, where the χ² law is only approximate. The check calibrates the test on sequences simulated from the grammar: with a reversal-symmetric truth Λ averaged 7.05 over six replicates of 40 sequences (8.4, 6.3, 7.7, 5.1, 4.1, 10.7; d = 9, standard score of the mean −1.1), and with an asymmetric one Λ = 194.6, p = 4.5 × 10⁻³⁷.

## 6. Escherichia coli's tRNAs

### 6.1 The Data

The RefSeq record NC_000913.3 in the data cache has 86 tRNA genes (features of key `tRNA`, extracted by `genomics/sequence`'s GenBank reader, T read as U), of which 48 sequences are distinct; every distinct one ends in CCA (measured, not assumed). The distinct sequences, in the order of the genome, are dealt alternately into a training half (24 tRNAs, 1,902 bases) and a held-out half (24, 1,906). Isoacceptors resemble each other, so the halves are not independent families: held-out likelihoods measure prediction of related sequences, not of a new kind of RNA.

No reference structure is typed anywhere. What is checked is what can be measured: the *acceptor stem*, the seven pairs (k, n − 5 − k), k = 0, …, 6, which join the first seven bases to the seven before the discriminator base and the CCA, and shuffles of each tRNA (three each, Fisher–Yates on calx's streams), which keep its composition and lose its order.

### 6.2 Structure Learnt without Structures

The reversible grammar without lonely pairs on the six canonical pairs, trained by EM from a start that weighs every canonical pair alike (`starting_parameters`, nothing known of which pairs are common), converges in 38 iterations from ℓ = −2855.2 to −2668.1 on the training half. On the held-out half (`check_genomics_rna_trna`):

| Measure | Held-out tRNAs | Their shuffles |
|---|---|---|
| Viterbi structures with ≥ 6 of the 7 acceptor-stem pairs | 18 of 24 (10 with all 7) | 0 of 72 |
| MEA structures (γ = 1) with ≥ 6 of the 7 pairs | 17 of 24 | 0 of 72 |
| mean posterior probability of an acceptor-stem pair | 0.634 | 0.010 |
| more probable than the shuffle | 72 of 72 comparisons | |
| mean log-odds against the shuffle | 6.9 nats | |

The grammar's held-out log-likelihood is −2660.7 nats (2.01 bits a base), while base composition alone, estimated from the training half, gives −2610.0: as a model of sequences the grammar is worse than composition, since it spends probability on structure, but it tells the tRNAs from their shuffles every time, and its structures carry their acceptor stems.

Knudsen and Hein's grammar as published (all sixteen pairs, lonely pairs allowed), trained the same way from the uniform start, rises from −2814.6 to −2666.4, predicts the held-out half worse (−2698.1), and finds the acceptor stem in none of the 24 Viterbi structures (mean pair probability 0.000).

*Remark.* In development runs its Viterbi structures were full of lonely and non-canonical pairs closing loops of two bases: with sixteen pair rules and lonely pairs, EM explains local correlations of the sequence by short "pairs" rather than by helices. Neither restriction alone suffices: from the same uniform start, the reversible grammar on the canonical pairs with lonely pairs allowed gave the acceptor stem (six pairs or more) in 4 of the 24 held-out Viterbi structures, and on all sixteen pairs without lonely pairs in none; the two restrictions together give 18 (development runs, not repeated by the checks).

*Remark (local maxima).* The untied canonical fit started from the reversal-tied optimum (§6.3) reaches ℓ = −2650.2 on the training half, 17.9 nats above the fit from the uniform start. EM's answer depends on where it starts, and the uniform start is not the best one.

### 6.3 The Data's Verdict on the Symmetries of Base Pairing

Each symmetry is tested on the grammar whose pair set is the closure of the canonical pairs under it (Proposition 3.6): the canonical pairs for ρ, and the canonical pairs with A·C and C·A for κ and κρ.

| Tie | Pairs | d | Λ (training) | p | Held-out ℓ, tied | Held-out ℓ, untied |
|---|---|---|---|---|---|---|
| ρ (reversal) | canonical | 9 | 34.4 | 7.7 × 10⁻⁵ | −2666.1 | −2652.9 |
| κ (complement) | canonical, A·C, C·A | 14 | 43.3 | 7.7 × 10⁻⁵ | −2674.4 | −2658.1 |
| κρ (reverse complement) | canonical, A·C, C·A | 8 | 18.7 | 0.017 | −2659.0 | −2657.8 |

**Base pairing in these tRNAs is not reversal-symmetric.** The test rejects the reversal tie, and the untied model predicts the held-out half better by 13.2 nats. The asymmetry is in which base of a G·C pair comes first, and it is strongest at the ends of helices. In the untied fit the pairs that open a helix (L, the outer pair, its 5′ base first) are G·C 0.548 against C·G 0.167; the second pairs (H) are C·G 0.631 against G·C 0.170; the pairs further in (F) are G·C 0.483 against C·G 0.217. A·U and U·A are nearly balanced (0.061 and 0.067 on L, 0.148 and 0.146 on F); G·U and U·G are lopsided too (0.156 and 0.000 on L, 0.000 and 0.065 on H). The complement tie is rejected as strongly (it forces G·U to equal C·A). The reverse-complement tie, which fixes every Watson–Crick pair and ties G·U with A·C and the unpaired A with U and C with G, has p = 0.017, rejected at the 5% level but not at 1%, and its held-out likelihood is within 1.2 nats of the untied one.

*Caveat.* These are 24 tRNAs of one genome, a family whose helices are few and similar; the verdict is about them, under this grammar, at EM's maxima. It says that the statistics of these helices distinguish a pair from its reverse, not why.

## 7. Limits

- **Length.** Probabilities are morphe's `Arithmetic<double>`; they underflow beyond about 350 bases, and the code refuses rather than return a wrong number. A scaled or logarithmic inside–outside would belong in morphe.
- **Cost.** O(n³ |R_binary|) time per sequence for every computation, O(n² |V_binary|) space; tRNAs take tens of milliseconds, a 1,500-base ribosomal RNA would take minutes and underflow.
- **The model.** No stacking energies, no loop lengths, no pseudoknots; minimum hairpin loop two. The structures are those of a lightweight grammar, and the acceptor stem is the one feature checked.
- **EM.** Local maxima; the verdicts of §6.3 are likelihood-ratio tests at the maxima found, with the boundary caveat of §5.3.

## References

- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), Maximum likelihood from incomplete data via the EM algorithm, *Journal of the Royal Statistical Society B* 39, 1–38.
- Do, C. B., Woods, D. A. and Batzoglou, S. (2006), CONTRAfold: RNA secondary structure prediction without physics-based models, *Bioinformatics* 22, e90–e98.
- Dowell, R. D. and Eddy, S. R. (2004), Evaluation of several lightweight stochastic context-free grammars for RNA secondary structure prediction, *BMC Bioinformatics* 5, 71.
- Durbin, R., Eddy, S. R., Krogh, A. and Mitchison, G. (1998), *Biological Sequence Analysis*, Cambridge University Press, ch. 9–10.
- Goodman, J. (1999), Semiring parsing, *Computational Linguistics* 25, 573–605.
- Knudsen, B. and Hein, J. (1999), RNA secondary structure prediction using stochastic context-free grammars and evolutionary history, *Bioinformatics* 15, 446–454.
- Lari, K. and Young, S. J. (1990), The estimation of stochastic context-free grammars using the Inside-Outside algorithm, *Computer Speech and Language* 4, 35–56.
- McCaskill, J. S. (1990), The equilibrium partition function and base pair binding probabilities for RNA secondary structure, *Biopolymers* 29, 1105–1119.
- Nussinov, R., Pieczenik, G., Griggs, J. R. and Kleitman, D. J. (1978), Algorithms for loop matchings, *SIAM Journal on Applied Mathematics* 35, 68–82.
- Nussinov, R. and Jacobson, A. B. (1980), Fast algorithm for predicting the secondary structure of single-stranded RNA, *PNAS* 77, 6309–6313.
- Wilks, S. S. (1938), The large-sample distribution of the likelihood ratio for testing composite hypotheses, *Annals of Mathematical Statistics* 9, 60–62.
