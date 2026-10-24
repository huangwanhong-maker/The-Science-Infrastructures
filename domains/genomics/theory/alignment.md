# Alignment as a Sum over Paths

*Status: draft. Serves learning and use.* The pairwise alignment of sequences as a sum over the paths of a graph in a semiring: the edit graph and its modes; the semirings that give scores, costs, counts, probabilities and parametric families; the classical algorithms as schedules of one sum; the polytope algebra and parametric alignment as tropical geometry; the statistics of local scores; and the pair hidden Markov model. The definition is morphe's (`foundations/morphe/formal/`: `semiring.hpp`, `alignment.hpp`, `polytope.hpp`); the biology and the fast schedules are the module [`align`](../align/) of the genomics domain. The design is the idea note [The Genome Infrastructure](../../../development/ideas/2026-10-24-the-genome-infrastructure.md), §1.1.

## 1. The Edit Graph

Let Σ be a finite alphabet and x = x₁ ⋯ x_m, y = y₁ ⋯ y_n words over it; in the code x is the *query* and y the *target*.

**Definition 1.1 (edit graph).** The edit graph G(x, y) is the directed graph on the vertices V = {0, …, m} × {0, …, n} with three kinds of edge into a vertex (i, j):

- M, from (i − 1, j − 1), labelled by the pair (x_i, y_j);
- I, from (i − 1, j), labelled by (x_i, −): a letter of x against a gap;
- D, from (i, j − 1), labelled by (−, y_j).

G is acyclic, and the row-major order of V is a topological order.

**Definition 1.2 (alignment).** An alignment of x and y in a mode is a directed path of G that begins at a *source* and ends at a *sink* of the mode:

| Mode | Sources | Sinks |
|---|---|---|
| global (Needleman and Wunsch 1970) | (0, 0) | (m, n) |
| local (Smith and Waterman 1981) | every vertex | every vertex |
| semiglobal (x end to end, y's ends free) | the first row {(0, j)} | the last row {(m, j)} |
| overlap (free end gaps on both) | the first row and the first column | the last row and the last column |

A path is determined by its start and its word over {M, I, D}, read along it; the words of global paths are those with |w|_M + |w|_I = m and |w|_M + |w|_D = n.

**Proposition 1.3 (Delannoy).** The number of global alignments of words of lengths m and n is the Delannoy number

D(m, n) = Σ_k C(m, k) C(n, k) 2^k.

*Proof.* A global path with d steps M has m − d steps I and n − d steps D, so the words number Σ_d (m + n − d)! / (d! (m − d)! (n − d)!) = Σ_d C(m + n − d, m) C(m, d), the lattice paths from (0, 0) to (m, n) with steps (1, 0), (0, 1) and (1, 1). That this equals Σ_k C(m, k) C(n, k) 2^k is a classical identity of the Delannoy numbers (Banderier and Schwer 2005); both satisfy D(m, n) = D(m − 1, n) + D(m, n − 1) + D(m − 1, n − 1) with D(m, 0) = D(0, n) = 1, the recurrence of the last step. ∎

*Remark.* D(n, n) = 1, 3, 13, 63, 321, … (OEIS A001850) grows like (3 + 2√2)^n: exhaustive search over alignments is hopeless beyond toy lengths, and the dynamic programme of §2 is what makes alignment possible. morphe's check `alignment` counts D(60, 60), a number of 45 digits, exactly.

## 2. Semirings and the Sum over Paths

**Definition 2.1 (semiring).** A semiring (K, ⊕, ⊗, 0, 1) is a set with two associative operations, ⊕ commutative with identity 0, ⊗ with identity 1, ⊗ distributing over ⊕ on both sides, and 0 ⊗ a = a ⊗ 0 = 0 (Golan 1999, ch. 1). It is *idempotent* when a ⊕ a = a and *selective* when a ⊕ b ∈ {a, b}.

The semirings of alignment (morphe's `semiring.hpp`):

| Semiring | Carrier, ⊕, ⊗ | Selective | What the sum over paths is |
|---|---|---|---|
| Boolean | {0, 1}, ∨, ∧ | yes | whether an alignment exists |
| tropical (max, +) | ℤ ∪ {−∞}, max, + | yes | the optimal score |
| tropical (min, +) | ℤ ∪ {+∞}, min, + | yes | the least cost; the edit distance with unit costs |
| counting | ℕ, +, × | no | the number of alignments |
| probability | ℝ≥0, +, × | no | the probability of the pair HMM (§8) |
| tropical with multiplicities | (ℤ ∪ {−∞}) × ℕ | no | the optimal score and the number of optimal alignments |
| polytope algebra | lattice polygons, conv ∪, Minkowski + | no (idempotent) | every optimal alignment for every parameter (§6) |

**Definition 2.2 (weighted edit graph, alignment value).** A scheme of weights w assigns to every edge e of G an element w(e) ∈ K; the weight of a path π = e₁ ⋯ e_k is w(π) = w(e₁) ⊗ ⋯ ⊗ w(e_k) (the empty path weighs 1). The *alignment value* of x and y in a mode is

A(x, y) = ⊕ { w(π) : π a path from a source to a sink }.

The linear scheme weighs an M edge by a pair weight s(x_i, y_j) and a gap edge by a gap weight of its letter. In the tropical semiring with s the substitution score and every gap letter −e this is the score of alignment with linear gaps; with (min, +), unit mismatch costs and unit gaps, the edit distance.

**Proposition 2.3 (the dynamic programme).** Define, in row-major order,

V(i, j) = [(i, j) is a source] ⊕ V(i−1, j−1) ⊗ w(M_{ij}) ⊕ V(i−1, j) ⊗ w(I_{ij}) ⊕ V(i, j−1) ⊗ w(D_{ij}),

where [source] is 1 at a source and 0 elsewhere and terms with a negative index are absent. Then V(i, j) is the sum of the weights of the paths from a source to (i, j), and A(x, y) = ⊕_{sinks v} V(v).

*Proof.* By induction along the topological order. A path to (i, j) is the empty path (when (i, j) is a source) or a path to a predecessor u followed by the edge (u, (i, j)). Grouping the paths by their last edge and using right distributivity, ⊕_π w(π) ⊗ w(e) = (⊕_π w(π)) ⊗ w(e) = V(u) ⊗ w(e) (Mohri 2002, Theorem 1, for acyclic graphs). ∎

The cost is O(mn) operations of K, and two rows of values suffice for the sum. This is `morphe::formal::path_sum`; the check `alignment` compares it, in every mode and both gap models, with the enumeration of every path of short random words (`enumerate_paths`), in five semirings.

*Remark (local alignment).* In the local mode every vertex is a source, so the recurrence ⊕-adds the semiring's one at every vertex. In (max, +), whose one is 0, this is Smith and Waterman's V(i, j) = max(0, …): a local alignment may begin anywhere because beginning costs nothing. The same sentence in the counting semiring counts all local alignments, the empty ones at each vertex included.

**Proposition 2.4 (homomorphisms commute with the sum).** Let h : K → L be a homomorphism of semirings. Then h(A_w(x, y)) = A_{h∘w}(x, y).

*Proof.* A_w is a polynomial expression in the edge weights built from ⊕ and ⊗ (a sum of products), and h preserves both operations and both identities. ∎

Three consequences carry the design. Forgetting the multiplicities, (s, c) ↦ s, is a homomorphism from the tropical semiring with multiplicities to (max, +), so the count of optimal alignments is attached to the very score Gotoh computes. The map x ↦ −x is an isomorphism (max, +) → (min, +): maximising scores and minimising penalties are one problem. Evaluating a polygon's support function at w is a homomorphism from the polytope algebra to (max, +) (§6): one polygon carries the tropical optimum at every parameter. The map (s, c) ↦ c, by contrast, preserves ⊗ but not ⊕ ((1, 1) ⊕ (0, 1) = (1, 1)), and the check `semiring` exhibits the failure: the number of optimal alignments is not a function of the number of alignments.

**Proposition 2.5 (traceback).** In a selective semiring every sum equals one of its terms, so an optimal path is found by following back from a sink whose value equals A the first term of each recurrence that equals its value, until a source whose value is its start term.

*Proof.* Selectivity gives, by induction on the number of terms, that a finite ⊕-sum equals one of its terms. Each step of the walk moves to a predecessor whose value, multiplied by the edge weight, is the current value, so the product of the weights along the walk is A. ∎

This is `optimal_path`, and in the fast algorithms the directions they keep are a record of which term won. In the tropical semiring with multiplicities, which is not selective, there is no traceback of its own: its projection to (max, +) has one.

## 3. Affine Gaps and Three States

A gap of length L scores −(o + L e), with o the open and e the extend (the convention of BLAST and minimap2; `types.hpp`). This is not a sum of edge weights on G: the weight of an I edge depends on whether the edge before it was I.

**Definition 3.1 (the three-state scheme).** Let A be the automaton with states M, I, D, a begin and an end, a weight T[r][s] on each transition r → s and an emission weight on each state: the pair weight s(x_i, y_j) on M, a letter weight on I and D. The weighted graph is the product G × A: its vertices are pairs ((i, j), s), an edge of G of kind s entering state s. *Affine gaps* are the transitions T[r][I] = −o for r ≠ I, T[I][I] = 0 (and the same for D), all others 0, with emission −e on I and D, in (max, +).

**Proposition 3.2.** Every word over {M, I, D} is the label of exactly one path of A, so the paths of G × A from a source to a sink are in bijection with the alignments of the mode, and the weight of a path is the affine score of its alignment.

*Proof.* The automaton is deterministic on the alphabet {M, I, D} with every transition present, so a word has one run. On a run, the open is paid exactly at the first letter of each maximal block of I or of D, and the extension at every letter. ∎

So the counts, the enumeration and the semantics of §2 hold unchanged, and in the tropical semiring the recurrences of Proposition 2.3 on the product become Gotoh's (1982):

- I(i, j) = max(H(i − 1, j) − o − e, I(i − 1, j) − e),
- D(i, j) = max(H(i, j − 1) − o − e, D(i, j − 1) − e),
- H(i, j) = max(start(i, j), H(i − 1, j − 1) + s(x_i, y_j), I(i, j), D(i, j)),

with H the best over the three states. Gotoh's paper took the minimum of I and D over H only, which Altschul and Erickson (1986) showed can lose the optimum when it is traced back; with H including I, and o ≥ 0, a block is never split, because closing and reopening a gap costs o more than extending it. The code refuses o < 0 for this reason (`gotoh.hpp`).

*Remark.* The same three-state product with other transitions is the pair HMM (§8), whose X and Y may not follow each other: T[I][D] = T[D][I] = 0, the zero of the semiring. Generic code over the product serves both (`morphe::formal::three_state`).

## 4. Algorithms as Schedules

Every algorithm of the module computes the tropical value of §3 (or a sub-problem of it) in an order and a representation of its own. This is calx's separation of *what* from *how* (calx architecture §1.2) in the terms of alignment: the definition is one sum; a schedule chooses which vertices to visit, in what order, what to keep, and on what hardware.

### 4.1 Gotoh Row by Row

`gotoh.hpp` visits the vertices row by row, keeping two rows of H and I (D runs along the row in a register) and one byte of directions per vertex: the source of H in two bits, and whether I and D extended. The traceback reads the directions back from the optimal sink. Ties are broken by a fixed order (start, pair, I, D; open before extend; the first optimal sink in row-major order), so every machine returns the same alignment; reading back from the end, preferring the pair puts gaps as far left as they can go, which is the normal form of indels in VCF. `gotoh_score` and `gotoh_end` drop the directions and run in O(n) memory; `gotoh_field` keeps every value for drawing.

### 4.2 Linear Space

The directions take mn bytes. Hirschberg (1975) observed that the optimum through a middle row can be found from two passes that keep one row each, and Myers and Miller (1988) carried the idea to affine gaps.

**Proposition 4.1 (the middle row).** Let i* be a row strictly between 0 and m. For each column j let F(j) be the best score of a path from the start to (i*, j) and F_I(j) of one that ends with an I step; let B(j) and B_I(j) be the same for paths from (i*, j) to the end, beginning with an I step for B_I. Then the optimum is

max_j max(F(j) + B(j), F_I(j) + B_I(j) + o).

*Proof.* An optimal path π meets row i*; let (i*, j) be the last vertex of π on it, so that π leaves it by an M or an I step. Unless π also enters it by an I step, cutting it there splits it into a path counted by F(j) and one counted by B(j), whose scores add up to π's. If it does, a vertical block runs through (i*, j), and the cut divides it into two, each of which paid an open in F_I and B_I, and the + o restores the single open of π. Conversely each term is the score of some path, or less (a split horizontal block pays its open twice), so the maximum is at most the optimum. ∎

The recursion splits the problem at the maximising (i*, j); in the second case it removes the two vertical steps through (i*, j), which carry the open, and recurses with the open of a vertical gap at the shared boundary set to zero, Myers and Miller's tb and te. The time is O(mn) (each level reads half the area of the level above, so about twice Gotoh's) and the memory O(m + n). `hirschberg.hpp` also aligns in the other modes in linear space: the end of an optimal path from a linear-memory pass, its start from a pass backwards from that end, and the path between them by the global recursion on the two substrings, since every path between those two vertices is a path of the mode.

### 4.3 The Band

Restricting the recurrences to the vertices with j − i ∈ [low, high] is the sum over the sub-graph of paths inside the band: never more than the optimum, and equal to it when some optimal path lies inside (Chao, Pearson and Miller 1992).

**Proposition 4.2 (a certificate).** In global alignment, let s_max be the largest substitution score (or 0), d = n − m, and suppose the band contains the diagonals 0 and d. A path that leaves the band has at least g_out = min(2(high + 1) − d, d − 2(low − 1)) gap columns, and a path with g gap columns scores at most U(g) = s_max (m + n − g)/2 − o − g e. If max(U(g_out), U(m + n)) is at most the banded score S, the banded alignment is optimal.

*Proof.* A path reaches diagonal k > max(0, d) only after at least k D steps and returns to d only after at least k − d I steps, so it has at least 2k − d gap columns; similarly below. Its pairs number (m + n − g)/2, each scoring at most s_max, and it opens at least one gap. U is affine in g, so its maximum over g ≥ g_out is at an end of the range. Every path outside the band then scores at most S. ∎

`banded_global` doubles the band's half-width from a small start until the certificate holds, the iteration of Fickett (1984) and Ukkonen (1985): the cost is that of the last band, O(m w), and the answer is always optimal.

### 4.4 Pruned Frontiers

The X-drop extension (Zhang et al. 2000; Altschul et al. 1997) computes the recurrences from an anchor outwards and discards every vertex whose score falls more than X below the best seen; the frontier of live vertices follows the alignment as a band of variable width. It is not exact: a path that dips by more than X and recovers is lost, which is the point (two sequences diverge after the homologous region, and nothing should be computed there). With X = ∞ it is the best alignment anchored at the anchor with a free end; the check verifies this by brute force over all prefixes. minimap2's Z-drop (Li 2018) adds a stop when the best cell of the current front falls more than Z + e |Δi − Δj| below the best seen: the e |Δi − Δj| forgives the extension of a long gap, so the rule separates an alignment at a block of mismatches (an inversion, a divergent insertion) while letting it cross a clean indel.

### 4.5 Schedules on Lanes

The kernels of the second builder are further schedules of the same sums (their headers and checks give the details): Myers' bit vectors (1999) hold a column of the unit-cost (min, +) table as differences in {−1, 0, +1}, two bits per vertex, and advance 64 vertices per word operation; Farrar's striped kernel (2007) computes the local (max, +) recurrences on vector lanes in an order that makes the vertical dependence rare; the wavefront algorithm (Marco-Sola et al. 2021) visits vertices in order of increasing penalty along the diagonals, so that its work is proportional to the penalty of the optimum rather than to mn. Each is checked against Gotoh, and so, through the checks above, against the definition.

## 5. Counting Optimal Alignments

The tropical semiring with multiplicities has elements (s, c), a score and the number of paths attaining it, with (s, c) ⊕ (t, d) the pair of the larger score, adding the counts on a tie, and (s, c) ⊗ (t, d) = (s + t, cd). It is the image of the Laurent polynomials ℕ[t, t⁻¹] under the leading term, which is a homomorphism because ℕ has neither negatives nor zero divisors: an alignment of score s is the monomial t^s, the sum over alignments is the generating polynomial of the scores, and its leading term is the optimum with its multiplicity (Pachter and Sturmfels 2005, §2.1).

*Example.* For x = A^k and y = A^{k−1} with linear gaps, an optimal global alignment has k − 1 matches and one gap, which may stand in any of k places: the count is k (`reference_count`), and Gotoh's tie rule chooses the leftmost (check `gotoh`).

## 6. The Polytope Algebra and Parametric Alignment

The dependence of the optimal alignment on the scores is an old problem (Waterman, Eggert and Lander 1992; Gusfield, Balasubramanian and Naor 1994); Pachter and Sturmfels (2004) solved it as a sum over paths in one more semiring.

**Definition 6.1 (polytope algebra).** Let P₂ be the set of convex lattice polygons in ℤ², the empty set, points and segments included, with P ⊕ Q = conv(P ∪ Q), P ⊗ Q = P + Q = {p + q : p ∈ P, q ∈ Q}, 0 = ∅ and 1 = {(0, 0)}.

**Proposition 6.2.** P₂ is a commutative, idempotent semiring.

*Proof.* Convex hull of a union and Minkowski sum are associative and commutative, with the identities given, and ∅ annihilates the sum. Distributivity: conv(P ∪ Q) + R = conv((P + R) ∪ (Q + R)), since both are the convex hull of {p + r} ∪ {q + r}, the Minkowski sum of convex sets commuting with convex hulls. Idempotence: conv(P ∪ P) = P. ∎

**Definition 6.3 (the alignment polygon).** Weigh each edge of G by a point: a pair of equal letters by (0, 0), a mismatch by (1, 0), a gap column by (0, 1). Then the weight of an alignment is the point (its mismatches, its gap columns), and the alignment value in P₂ is the alignment polygon P(x, y), the convex hull of the points of all alignments. Its computation by Proposition 2.3 is *polytope propagation* (`alignment_polygon`, `parametric`).

**Proposition 6.4 (evaluation is tropical).** For w ∈ ℤ², the support function h_P(w) = max_{p ∈ P} ⟨w, p⟩ is a homomorphism P₂ → (max, +), and h_{P(x, y)}(w) is the optimal score of the tropical alignment whose edge weights are ⟨w, ·⟩ of the points.

*Proof.* h_{conv(P ∪ Q)} = max(h_P, h_Q) because a linear function attains its maximum over a convex hull at a point of the generating set; h_{P+Q} = h_P + h_Q because the maximum of ⟨w, p + q⟩ separates; h_∅ = −∞ and h_{{0}} = 0. Proposition 2.4 does the rest. ∎

In global alignment the points lose nothing: 2 (matches) + 2 (mismatches) + (gap columns) = m + n, so the triple (matches, mismatches, gaps) is determined by (x, g), and the score with match 1, mismatch −α and gap −β is

score = (m + n)/2 − (1 + α) x − (1/2 + β) g.

**Proposition 6.5 (the normal fan is the parametric decomposition).** An alignment with point v is optimal at (α, β) if and only if v minimises ⟨c, ·⟩ over P(x, y) for c = (1 + α, 1/2 + β), that is, −c lies in the outer normal cone N(v); for c in the interior of a vertex's cone, the optimal alignments are exactly those with that point. The regions of the (α, β) plane are therefore cones with common apex (−1, −1/2), one for each vertex.

*Proof.* The score is an affine function of the point with linear part −c, and a linear function's minimisers over a polygon are a face, a single vertex exactly when the direction is interior to that vertex's normal cone. ∎

So the polygon with its normal fan is the whole parametric family. A witness for each vertex comes from the tropical traceback at an interior direction of its cone, where the maximiser is unique (`alignment_polygon`); the check `parametric` finds the same vertex with the domain's own Gotoh at that direction. In the language of tropical geometry (Pachter and Sturmfels 2004, §2): the alignment value is a tropical polynomial in the parameters, the maximum of the linear functions ⟨c, v⟩ over the points; P(x, y) is its Newton polygon; the parameters at which the optimum changes form the tropical curve of that polynomial, which is the codimension-one skeleton of the normal fan.

**Proposition 6.6 (few vertices).** A convex lattice polygon contained in an N × N box has O(N^{2/3}) vertices (Andrews 1963), so the alignment polygon has O((m + n)^{2/3}) vertices, and polytope propagation costs O(mn V log V) with V of that order.

*Sketch.* The edges of a convex lattice polygon have distinct directions, each a multiple of a primitive vector, and their lengths sum to the perimeter, at most 4N. The number of primitive vectors of length at most r grows as (6/π) r², so V edges need total length of order V^{3/2}, and V^{3/2} ≲ N gives V ≲ N^{2/3}. The sharp constant is 12 (4π²)^{−1/3} ≈ 3.52 (Acketa and Žunić 1995). ∎

*Remark.* This is why parametric alignment is practical: the number of genuinely different optimal alignments, as the scores vary over the whole plane, is a few dozen for sequences of a few hundred letters (the check reports 12 for a pair of 120 and 132 bases), while the number of alignments is astronomical (Proposition 1.3).

## 7. The Statistics of Local Scores

The size of the best local score of two unrelated sequences, and so the significance of an observed one, was found by Karlin and Altschul (1990) for ungapped alignment.

Let the letters of the two sequences be independent with frequencies p and q, and let P(s) = Σ_{a,b} p_a q_b [s_ab = s] be the distribution of the score of one aligned pair; suppose E[s] < 0 (otherwise long random alignments score well and local alignment degenerates into global) and P(s) > 0 for some s > 0.

**Proposition 7.1 (λ).** The equation Σ_s P(s) e^{λs} = 1 has exactly one positive root.

*Proof.* φ(λ) = Σ_s P(s) e^{λs} is a strictly convex function with φ(0) = 1 and φ′(0) = E[s] < 0, so it dips below 1 for small λ > 0, and grows without bound because some s > 0 has positive probability. A strictly convex function meets the level 1 at most twice, once at 0. ∎

The root has a meaning: the numbers p_a q_b e^{λ s_ab} are positive and sum to 1, a distribution of aligned pairs (the *target frequencies*), and λ s_ab = log(p_ab / (p_a q_b)) is the log-odds score of an aligned pair against chance (Altschul 1991). A substitution matrix, scaled by its λ, is the log-odds matrix of some target distribution, and its relative entropy H = Σ p_ab log(p_ab / p_a q_b) = λ Σ_s s P(s) e^{λs} is the information per aligned pair, in nats. BLOSUM62's file states 0.6979 bits relative to its own background; under Robinson and Robinson's frequencies H = 0.4012 nats.

**Theorem 7.2 (Karlin and Altschul 1990).** For sequences of lengths m and n, the number of distinct local alignments with score at least S is asymptotically Poisson with mean E = K m n e^{−λS}, so the best score M satisfies P(M ≥ S) ≈ 1 − exp(−K m n e^{−λS}), an extreme-value law of location ln(K m n)/λ and scale 1/λ. The constant K is

K = λ δ e^{−2σ} / (H (1 − e^{−λδ})),   σ = Σ_{k ≥ 1} (1/k) [E(e^{λ S_k}; S_k < 0) + P(S_k ≥ 0)],

with δ the greatest common divisor of the scores of positive probability and S_k the sum of k independent scores (Karlin and Altschul 1990, Appendix; Karlin and Dembo 1992).

`statistics.hpp` finds λ by Brent's method, H by its sum, and K by the series, each P_k by convolution, as NCBI BLAST sums it, with the closed forms BLAST uses when the reduced scores reach only −1 or only +1 (the check confirms the series against them). The bit score S′ = (λS − ln K)/ln 2 makes the law universal, E = m n 2^{−S′}. With BLOSUM62 and Robinson and Robinson's frequencies the code gives λ = 0.3176, K = 0.1337, H = 0.4012 against NCBI BLAST's 0.3176, 0.134, 0.4012, and agrees with BLAST's table for the seven other matrices (BLOSUM80 in the half-bit form BLAST builds in, which differs from NCBI's third-bit file); a simulation of 3000 random protein pairs reproduces the tail and the mean of the extreme-value law once the lengths are corrected for the edge, m′ = m − ln(K m n)/H (Altschul and Gish 1996).

*Caveat.* For gapped alignment no such theory is proved; λ and K are estimated by simulation, and the module takes BLAST's table of them (BLOSUM62 with open 11 and extend 1: λ = 0.267, K = 0.041).

## 8. The Pair Hidden Markov Model

A probabilistic model of two related sequences turns the scores of alignment into probabilities (Durbin et al. 1998, ch. 4).

**Definition 8.1 (pair HMM).** States M, X, Y, a begin that behaves as M, and an end. M emits an aligned pair (a, b) with probability p_ab, X emits a letter a of x against a gap with probability q_a, Y a letter b of y with q_b. Transitions: M → M with 1 − 2δ − τ, M → X and M → Y with δ, X → X and Y → Y with ε, X → M and Y → M with 1 − ε − τ, every state to the end with τ; X and Y do not follow each other.

It is the three-state scheme of §3 with these transitions and emissions, so by Propositions 2.3 and 2.4:

- the probability P(x, y) = Σ_π P(x, y, π) is the alignment value in the probability semiring (+, ×): the *forward* algorithm;
- the most probable alignment is the value in (max, +) of the log probabilities: the *Viterbi* algorithm; and the traceback of Proposition 2.5 is Viterbi's.

The two are the same programme in two semirings, related by the tropical limit: T log Σ e^{a/T} → max a as T → 0, so that the log-sum-exp semiring at temperature T degenerates to (max, +). With a random model for unrelated sequences, the Viterbi log-odds is an affine-gap alignment score, with s(a, b) = log(p_ab / q_a q_b) up to a constant, an open and an extend determined by δ, ε and τ (Durbin et al. 1998, §4.1): the scores of §3 are the logarithms of a model.

**Proposition 8.2 (posterior of an aligned pair).** With f^M(i, j) the forward value of the paths to (i, j) ending in M and b^M(i, j) the backward value of what follows from (i, j) in M,

P(x_i ◇ y_j | x, y) = f^M(i, j) b^M(i, j) / P(x, y).

*Proof.* Every alignment that aligns x_i with y_j passes through the vertex (i, j) in state M, and is the concatenation of a path counted by f^M(i, j) and a continuation counted by b^M(i, j); the Markov property makes the probability factor. ∎

For a fixed i the events "x_i is aligned with y_j", j = 1 … n, are disjoint, so the posteriors along a row sum to the probability that x_i is aligned at all, at most 1; their total over (i, j) is the expected number of aligned pairs, which the check computes independently by enumerating every alignment. The *maximum expected accuracy* alignment (Holmes and Durbin 1998; ProbCons, Do et al. 2005) maximises the sum of the posteriors of its aligned pairs, by the dynamic programme of Proposition 2.3 in (max, +) with the posteriors as pair weights and zero gap weights; it aligns what is probably aligned, rather than the single most probable path, and it is the step from which consistency-based multiple alignment proceeds. In the code all of it runs in log space, log(e^a + e^b) = max(a, b) + log1p(e^{−|a − b|}), so that P(x, y) of two kilobase sequences, far below the smallest double, is still exact to rounding in its logarithm.

## 9. Where the Code Is

| Section | Code |
|---|---|
| §1, §2 the edit graph, the modes, the sum over paths, traceback, enumeration | `morphe/formal/alignment.hpp` (`path_sum`, `optimal_path`, `enumerate_paths`, `delannoy`) |
| §2, §5 the semirings and their homomorphisms | `morphe/formal/semiring.hpp` |
| §3 affine gaps, three states | `morphe/formal/alignment.hpp` (`affine`, `three_state`); `align/reference.hpp` |
| §4.1 Gotoh | `align/gotoh.hpp` |
| §4.2 linear space | `align/hirschberg.hpp` |
| §4.3 the band | `align/banded.hpp` |
| §4.4 X-drop and Z-drop | `align/xdrop.hpp` |
| §4.5 lanes, bit vectors, wavefronts | `align/bitparallel.hpp`, `align/wavefront.hpp`, and the second builder's further headers |
| §6 polytope algebra, parametric alignment | `morphe/formal/polytope.hpp`; `align/parametric.hpp` |
| §7 Karlin–Altschul | `align/statistics.hpp`, `align/scoring.hpp` |
| §8 the pair HMM | `align/pairhmm.hpp` |

## References

- Acketa, D. M. and Žunić, J. D. (1995), On the maximal number of edges of convex digital polygons included into an m × m grid, *J. Combin. Theory A* 69, 358–368.
- Altschul, S. F. (1991), Amino acid substitution matrices from an information theoretic perspective, *J. Mol. Biol.* 219, 555–565.
- Altschul, S. F. and Erickson, B. W. (1986), *Bull. Math. Biol.* 48, 603–616.
- Altschul, S. F. and Gish, W. (1996), Local alignment statistics, *Methods Enzymol.* 266, 460–480.
- Andrews, G. E. (1963), A lower bound for the volume of strictly convex bodies with many boundary lattice points, *Trans. Amer. Math. Soc.* 106, 270–279.
- Banderier, C. and Schwer, S. (2005), Why Delannoy numbers?, *J. Statistical Planning and Inference* 135, 40–54.
- Durbin, R., Eddy, S. R., Krogh, A. and Mitchison, G. (1998), *Biological Sequence Analysis*, Cambridge University Press.
- Farrar, M. (2007), Striped Smith–Waterman speeds database searches six times over other SIMD implementations, *Bioinformatics* 23, 156–161.
- Golan, J. S. (1999), *Semirings and their Applications*, Kluwer.
- Gotoh, O. (1982), *J. Mol. Biol.* 162, 705–708.
- Gusfield, D., Balasubramanian, K. and Naor, D. (1994), *Algorithmica* 12, 312–326.
- Hirschberg, D. S. (1975), A linear space algorithm for computing maximal common subsequences, *CACM* 18, 341–343.
- Karlin, S. and Altschul, S. F. (1990), *PNAS* 87, 2264–2268; Karlin, S. and Dembo, A. (1992), *Adv. Appl. Prob.* 24, 113–140.
- Marco-Sola, S., Moure, J. C., Moreto, M. and Espinosa, A. (2021), Fast gap-affine pairwise alignment using the wavefront algorithm, *Bioinformatics* 37, 456–463.
- Mohri, M. (2002), Semiring frameworks and algorithms for shortest-distance problems, *J. Automata, Languages and Combinatorics* 7, 321–350.
- Myers, E. W. and Miller, W. (1988), *CABIOS* 4, 11–17; Myers, G. (1999), A fast bit-vector algorithm for approximate string matching based on dynamic programming, *J. ACM* 46, 395–415.
- Pachter, L. and Sturmfels, B. (2004), *PNAS* 101, 16138–16143; (2005), *Algebraic Statistics for Computational Biology*, Cambridge University Press, ch. 2.
- Smith, T. F. and Waterman, M. S. (1981), *J. Mol. Biol.* 147, 195–197; Needleman, S. B. and Wunsch, C. D. (1970), *J. Mol. Biol.* 48, 443–453.
- Waterman, M. S., Eggert, M. and Lander, E. (1992), *PNAS* 89, 6090–6093.
- Zhang, Z., Schwartz, S., Wagner, L. and Miller, W. (2000), *J. Comput. Biol.* 7, 203–214; Li, H. (2018), *Bioinformatics* 34, 3094–3100.
