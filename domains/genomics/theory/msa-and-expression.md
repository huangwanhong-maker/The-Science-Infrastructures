# Multiple Alignment and Expression

*Status: draft. Serves learning and use.* Two problems of many sequences at once. The first is the multiple alignment of a family: its sum-of-pairs objective and why it is hard, progressive alignment as a greedy order on a tree, profiles, refinement, UPGMA, and the centre star with its bound. The second is the expression of a transcriptome: pseudoalignment as a meet in a lattice of subsets, the abundances as the maximum of a likelihood on a simplex reached by the EM algorithm, the negative binomial as a gamma–Poisson mixture with its generalised linear model, and the Wald test with the false discovery rate. The code is the modules [`msa`](../msa/) and [`expression`](../expression/); pairwise alignment, on which the first rests, is [alignment.md](alignment.md); the design is the idea note [The Genome Infrastructure](../../../development/ideas/2026-10-24-the-genome-infrastructure.md), §4.

## 1. The Sum-of-Pairs Objective

Let Σ be an alphabet, − ∉ Σ the gap, and x₁, …, x_k words over Σ.

**Definition 1.1 (multiple alignment).** A multiple alignment of x₁, …, x_k is a k × L matrix A over Σ ∪ {−} with no column of gaps only, whose i-th row with its gaps deleted is x_i. Write A_i for the i-th row.

**Definition 1.2 (projection).** The projection of A on rows i and j is the 2 × L′ matrix left when the columns in which both A_i and A_j are gaps are deleted; it is a pairwise alignment of x_i and x_j (a path of their edit graph, [alignment.md](alignment.md) §1).

**Definition 1.3 (sum of pairs).** For a pairwise score S (a substitution matrix s and affine gaps, a maximal run of L gaps in one row scoring −(o + Le)), the sum-of-pairs score of A is

SP(A) = Σ_{i<j} S(π_{ij}(A)),

π_{ij} the projection (Carrillo and Lipman 1988). Charging each run of the projection, rather than each run of a row, is Altschul's *natural* gap cost (1989).

*Remark.* SP counts every pair of rows, so a family of many close relatives weighs heavily against a few distant ones; sequence weights (§2.4) correct for the redundancy.

**Proposition 1.4 (hardness).** Finding an alignment of maximum SP score is NP-hard, already for a fixed alphabet and a metric cost (Wang and Jiang 1994, by reduction from the shortest common supersequence problem; Just 2001 extends it to a wide class of scoring schemes).

**Proposition 1.5 (the exact dynamic programme).** The alignments of k words of length n are the paths of a k-dimensional lattice of (n + 1)^k vertices, each with 2^k − 1 predecessors; the SP optimum with linear gaps is a dynamic programme over it in O(2^k n^k) time.

*Proof.* A column is a nonempty subset of the rows advancing by one; the score of a column with linear gaps is a sum over pairs and depends only on the column, so the edit graph's argument of [alignment.md](alignment.md) §2 applies in k dimensions. ∎

The exponent is why exact SP alignment stops at a handful of sequences, and why every method below is a heuristic with a reason.

## 2. Progressive Alignment

### 2.1 A Greedy Order on a Tree

**Definition 2.1 (progressive alignment).** Let T be a rooted binary tree whose leaves are the sequences. Progressive alignment assigns to each leaf its sequence as a one-row alignment and to each internal node v the alignment of its children's alignments, A(v) = merge(A(a), A(b), P), P an optimal path of the profile alignment of A(a) and A(b) (§2.2); the result is A(root) (Feng and Doolittle 1987; Higgins and Sharp 1988).

The merge never changes what a child aligned: it inserts columns of gaps into each side. This is the rule "once a gap, always a gap".

**Proposition 2.2.** For every node v, the alignment A(root) restricted to the rows below v, with its columns of gaps only deleted, is A(v).

*Proof.* By induction on the distance from v to the root. At v it is A(v) by definition. If it holds at the parent u of v, A(u) = merge(A(v), A(w), P), and merge places A(v)'s columns in order among columns that are gaps in all of A(v)'s rows; deleting those columns gives A(v) back. ∎

*Remark.* The tree fixes the order of a greedy algorithm: the closest pairs are aligned first, when the signal is clearest, and their alignment is never revised. A wrong early decision is never undone, which is what refinement (§2.3) addresses.

### 2.2 Profiles and Their Alignment

**Definition 2.3 (profile).** The profile of an alignment X of rows with weights w_r > 0, W = Σ w_r, holds for each column i the weighted frequencies f_i(a) = Σ_{r : X_{ri} = a} w_r / W of the letters, and the occupancy o_i = Σ_a f_i(a), the weighted fraction of rows with a residue there.

**Proposition 2.4 (column scores are expected pair scores).** Let X and Y be alignments with profiles f and g. If column i of X and column j of Y are placed in one column of the merge, the sum over the pairs (r, q) of a row of X and a row of Y, each weighted w_r w′_q, of the substitution scores between residues is W W′ · s(i, j), where

s(i, j) = Σ_a Σ_b f_i(a) g_j(b) s(a, b).

*Proof.* Σ_{r,q} w_r w′_q s(X_{ri}, Y_{qj}) over the pairs where both are residues is Σ_{a,b} s(a, b) (Σ_{r: X_{ri}=a} w_r)(Σ_{q: Y_{qj}=b} w′_q) by distributivity. ∎

So the profile alignment's pair score is the sum-of-pairs objective restricted to the pairs across the two sides, per pair of rows. Its gap costs cannot be made exact so cheaply: whether a new gap in Y opposite column i of X opens a gap in the projection of (r, q) depends on whether q had a residue in the column before, and on the whole history of both rows; the exact problem of aligning two alignments under SP with affine gaps is NP-complete (Kececioglu and Starrett 2004). The code (`progressive.hpp`) charges a column of x against a new gap in y, inserted after y's column j:

- the extension e · o_i: each residue of x's column pays an extension against every row of y;
- half the opening, (o/2) · o_i · o′_j, at the run's first column, and half, (o/2) · o_i · o′_{j+1}, at its last: the pairs whose y-row had residues on both sides of the insertion are those whose gap is new (o′ = 1 beyond either end), MUSCLE's gap-open and gap-close.

**Proposition 2.5.** For two single sequences, the profile alignment is Gotoh's global alignment: the same optimal score and an optimal path.

*Proof.* All weights are 1 and all occupancies 1, so s(i, j) = s(x_i, y_j), a column against a gap costs e, and a run costs o/2 + o/2 + Le = o + Le: the recurrences are Gotoh's ([alignment.md](alignment.md) §3). ∎

The check `progressive` holds the code to it on 300 pairs.

### 2.3 Refinement

**Definition 2.6 (restricted partitioning).** An edge of the tree splits the rows into the leaves S below it and the rest S̄. Re-aligning the edge replaces A by merge(A|_S, A|_{S̄}, P), A|_S the sub-alignment of S with its columns of gaps only deleted and P an optimal path of their profile alignment; the result is kept if its SP score is higher (Edgar 2004, stage 3).

**Proposition 2.7.** Re-aligning an edge changes only the scores of the pairs across the split: SP(A′) − SP(A) = Σ_{r∈S, q∈S̄} (S(π_{rq}(A′)) − S(π_{rq}(A))).

*Proof.* For r, r′ ∈ S, the projection π_{rr′}(A′) deletes the columns where both are gaps; A′ restricted to S is A|_S with columns of gaps inserted, and A|_S is A restricted to S with columns of gaps deleted, so π_{rr′}(A′) = π_{rr′}(A). The same holds within S̄. ∎

So the refinement computes only |S| · |S̄| pair scores per trial (`scores.hpp`, `sp_across`), and it is an ascent: the SP score never falls, and since it rises strictly at each accepted trial and the alignments are finitely many, it stops; here after a few rounds or a round without change.

### 2.4 The Guide Tree

**Definition 2.8 (UPGMA).** Starting from singletons, merge the two clusters of least distance, the distance to a union being the size-weighted mean, d(a ∪ b, k) = (|a| d(a, k) + |b| d(b, k)) / (|a| + |b|), the mean distance between their members (Sokal and Michener 1958). The node of a merge at distance δ has height δ/2.

**Proposition 2.9 (reducibility).** Average linkage satisfies d(a ∪ b, k) ≥ min(d(a, k), d(b, k)).

*Proof.* A convex combination of two numbers is at least their minimum. ∎

Reducibility is what makes the nearest-neighbour chain correct (Bruynooghe 1977; Murtagh 1983; Müllner 2011): if a and b are each other's nearest neighbours, merging them makes no other cluster's nearest neighbour nearer, so pairs of reciprocal nearest neighbours may be merged in any order, and following a chain of nearest neighbours until it closes finds one in amortised O(n) steps. The chain gives the greedy algorithm's merges in another order; sorted by distance they are its dendrogram (`guide_tree.hpp`, held to the greedy algorithm on 400 inputs).

**Proposition 2.10.** If d is an ultrametric, d(x, z) ≤ max(d(x, y), d(y, z)), realised by a rooted tree with leaves at height 0 (d twice the height of the least common ancestor), UPGMA returns that tree.

*Proof.* The two leaves of least distance have an ancestor of least height whose clusters are all at that distance; every member of one child is at the same distance 2h from every member of the other, so the mean is 2h and merging preserves the ultrametric on the quotient. Induction on the number of clusters. ∎

The distances that build the first tree are MUSCLE's k-mer distance: F(x, y) the number of k-mers the two share (with multiplicity) over the number in the shorter, d = 1 − F, computed by merging sorted lists in O(L) (Edgar 2004). The second tree is built from the first alignment's identities, corrected for multiple substitutions: under Jukes and Cantor's model (1969) a site differs after distance d with probability p = (3/4)(1 − e^{−4d/3}), so d = −(3/4) ln(1 − 4p/3); for protein Kimura's empirical d = −ln(1 − p − p²/5) (1983).

**Definition 2.11 (CLUSTAL W weights).** In a rooted tree with branch lengths, the weight of a leaf is Σ over the branches on its path to the root of the branch's length divided by the number of leaves below it (Thompson, Higgins and Gibson 1994).

*Remark.* A branch shared by m leaves is divided among them, so ten copies of one sequence share one weight, and an isolated sequence keeps its whole path.

## 3. The Centre Star

**Definition 3.1.** Let D(x, y) be the cost of an optimal pairwise alignment. The centre c minimises Σ_j D(c, x_j). The centre-star alignment aligns each x_j optimally to c and merges the k − 1 alignments by giving each gap of c (before each residue and at the end) as many columns as the largest number of residues any x_j places there (Gusfield 1993).

**Proposition 3.2.** The centre-star alignment induces each pairwise alignment of x_j to c exactly.

*Proof.* Its columns are c's residues and, in each gap of c, columns holding the residues the alignments place there, left-justified. Deleting the columns where both x_j and c are gaps leaves x_j's residues against c's residues and against gaps in the order of its alignment to c. ∎

**Proposition 3.3 (Gusfield's bound).** Let the cost of a pairwise alignment be a sum over its columns of δ(a, b) for a, b ∈ Σ ∪ {−}, with δ a metric on Σ ∪ {−} and δ(−, −) = 0 (edit distance; linear gaps). Then the SP cost of the centre-star alignment of k sequences is at most (2 − 2/k) times the optimal SP cost.

*Proof.* Write C for the centre-star alignment and D_C(i, j) for the cost of its projection on i and j. Summing δ over the columns of C, with δ(−, −) = 0 so that columns of two gaps add nothing, the triangle inequality column by column gives D_C(i, j) ≤ D_C(i, c) + D_C(c, j), and by Proposition 3.2 D_C(i, c) = D(i, c). Hence

SP(C) = Σ_{i<j} D_C(i, j) ≤ Σ_{i<j} (D(i, c) + D(c, j)) = (k − 1) Σ_j D(c, j),

since each index lies in k − 1 pairs (D(c, c) = 0). For any alignment A, D_A(i, j) ≥ D(i, j), so

OPT ≥ Σ_{i<j} D(i, j) = ½ Σ_i Σ_j D(i, j) ≥ ½ Σ_i Σ_j D(c, j) = (k/2) Σ_j D(c, j),

the last because c minimises the inner sum. The ratio is at most (k − 1)/(k/2) = 2 − 2/k. ∎

*Remark.* The proof needs a cost that decomposes over columns; affine gaps do not, and neither do similarity scores, so with them the bound is not claimed. The check `progressive` verifies it under edit distance against the lower bound Σ_{i<j} D(i, j), on 30 families: the largest ratio of cost to bound was 1.42, against 2 − 2/k ≥ 1.33 … 1.83.

## 4. Pseudoalignment as a Meet

Let T = {t₁, …, t_n} be the transcripts, each a word over {A, C, G, T}, and k fixed. For a word w of length k let w̄ be its reverse complement.

**Definition 4.1 (the class of a k-mer).** κ(w) = {t ∈ T : w or w̄ is a factor of t}. Since κ(w) = κ(w̄), it is a function of the canonical k-mer min(w, w̄).

**Definition 4.2 (the class of a read).** For a read r with k-mers K(r), its class is

C(r) = ⋂ { κ(w) : w ∈ K(r), κ(w) ≠ ∅ },

and r is pseudoaligned when some κ(w) is nonempty and C(r) ≠ ∅ (Bray et al. 2016). For a pair, K is the union of both mates' k-mers.

The sets κ(w) generate, under nonempty intersection, a ∩-semilattice L of subsets of T, and C(r) is the meet in L of the generators its k-mers name. The index stores the generators, one per distinct class; reads meet the others, which `pseudoalign` adds to its table in the order of the reads.

**Proposition 4.3 (the source is compatible).** If r is an error-free fragment of t, on either strand, then t ∈ C(r).

*Proof.* Every k-mer of r is a factor of t or of t̄, so t ∈ κ(w) for every w ∈ K(r), and every κ(w) is nonempty. ∎

*Caveat.* With errors the inclusion can fail: a k-mer changed by an error is usually in no transcript (κ = ∅, skipped), but if it is in another transcript, the intersection excludes the source. This is the price of compatibility being exact, and why kallisto skips k-mers absent from the index rather than failing the read.

**Proposition 4.4 (the incremental index).** Process the transcripts in order t₁, …, t_n, giving a canonical k-mer met first in t_i the class {t_i} and a k-mer of class C met again in t_i the class C ∪ {t_i}. After t_i, each k-mer w seen so far has the class κ(w) ∩ {t₁, …, t_i}. Every k-mer of class C met in t_i receives the same class C ∪ {t_i}, so a memo of the classes made during t_i builds each class once.

*Proof.* Induction on i: a k-mer met in t_i either was unseen (its class restricted to t₁ … t_i is {t_i}) or had κ(w) ∩ {t₁, …, t_{i−1}} = C, and t_i is greater than every member of C, so the union is C with t_i appended, a sorted list. ∎

The check `pseudoalign` holds the index to an independent construction (all (k-mer, transcript) pairs sorted) and every read's class to the intersection computed by brute force, and on a small transcriptome to Definition 4.2 itself, searching each k-mer in every transcript and its reverse complement.

## 5. Abundances as a Maximum on a Simplex

### 5.1 The Likelihood

A sample has molecules of transcript t in fraction ρ_t. A molecule of length l holds Σ_{f ≤ l} p(f)(l − f + 1) = P(F ≤ l)(l − μ_l + 1) fragments in expectation, p the law of fragment lengths and μ_l its mean truncated to [1, l]; kallisto's effective length is l̃ = l − μ_l + 1 (Pachter 2011; Bray et al. 2016). A fragment comes from t with probability α_t = ρ_t l̃_t / Σ_u ρ_u l̃_u and from a given position of it with probability 1/l̃_t, so a fragment of class c has probability Σ_{t∈c} α_t / l̃_t, and the n_c fragments of each class give

L(α) = Π_c (Σ_{t∈c} α_t / l̃_t)^{n_c},   α ∈ Δ = {α ≥ 0, Σ_t α_t = 1}.

**Proposition 5.1.** log L is concave on Δ.

*Proof.* Each factor is the logarithm of a linear form, concave, and n_c ≥ 0. ∎

### 5.2 The EM Algorithm

Let the latent transcript of each fragment be unobserved. Given α, a fragment of class c is from t ∈ c with posterior probability γ_{ct} = (α_t / l̃_t) / Σ_{u∈c} α_u / l̃_u. The EM map is

E(α)_t = (1/N) Σ_{c ∋ t} n_c γ_{ct},   N = Σ_c n_c.

**Proposition 5.2 (monotonicity).** L(E(α)) ≥ L(α), with equality only at a fixed point.

*Proof.* For any α′ ∈ Δ, by Jensen's inequality applied to the concave logarithm with the weights γ_{ct},

log Σ_{t∈c} α′_t / l̃_t − log Σ_{t∈c} α_t / l̃_t = log Σ_{t∈c} γ_{ct} (α′_t / α_t) ≥ Σ_{t∈c} γ_{ct} log(α′_t / α_t).

Summing with the weights n_c, log L(α′) − log L(α) ≥ Σ_t (Σ_{c∋t} n_c γ_{ct}) log(α′_t / α_t) = N Σ_t E(α)_t log(α′_t / α_t). At α′ = E(α) the right side is N times the Kullback–Leibler divergence of E(α) from α, which is nonnegative and zero only when E(α) = α. ∎

This is Dempster, Laird and Rubin's Theorem 1 (1977) in the case at hand, the EM as a minorise–maximise step; the check `quantify` records log L over 811 rounds and finds it never decreasing.

**Proposition 5.3 (fixed points).** α ∈ Δ is a fixed point of E if and only if, for every t with α_t > 0, ∂ log L / ∂α_t = Σ_{c∋t} n_c (1/l̃_t) / Σ_{u∈c} α_u / l̃_u = N. A fixed point at which also ∂ log L / ∂α_t ≤ N for every t with α_t = 0 is a maximum of L on Δ.

*Proof.* E(α)_t = α_t · (∂ log L / ∂α_t) / N, so E(α) = α exactly when the derivative is N on the support. The second set of conditions, with the first, are the Karush–Kuhn–Tucker conditions of maximising the concave log L over Δ (the multiplier of Σ α_t = 1 is N, since Σ_t α_t ∂_t log L = N by homogeneity), and for a concave problem they are sufficient. ∎

kallisto stops the iteration when every transcript of estimated count above 0.01 changes by less than one per cent, after at least 50 rounds (`quantify.hpp`).

### 5.3 Identifiability

L depends on α only through the class probabilities π_c = Σ_{t∈c} α_t / l̃_t, a linear map M: α ↦ π. Where M is not injective on Δ, the maximum is a face, not a point.

*Remark (an isoform no k-mer identifies).* Let a gene have exons e₁, e₂, e₃ and isoforms t₁ = e₁e₂e₃, t₂ = e₁e₂, t₃ = e₂e₃, t₄ = e₁e₃. Every k-mer of t₁ lies in another isoform: those of e₁ in t₂ and t₄, of e₂ in t₂ and t₃, of e₃ in t₃ and t₄, of the junction e₁e₂ in t₂, of e₂e₃ in t₃ (unless e₂ is shorter than k). So t₁ has no class of its own, and its reads can be redistributed among t₂, t₃, t₄ with little or no change in L. The simulated transcriptome of the check makes such genes; the EM leaves some of them at 0 with their reads given to their sub-isoforms. The check therefore holds the estimate to the truth where it is identifiable, on transcripts with a tenth of their reads in classes of their own (mean relative error 0.038, the bootstrap's ±2 deviations covering the truth for 226 of 237) and on the genes' totals (0.003), and reports the rest (mean 0.31, median 0.074).

## 6. The Negative Binomial Model

### 6.1 A Gamma–Poisson Mixture

**Proposition 6.1.** Let λ ~ Γ(shape r, scale μ/r) and, given λ, K ~ Poisson(λ). Then

P(K = y) = Γ(y + r) / (Γ(r) y!) · (r/(r + μ))^r (μ/(r + μ))^y,

the negative binomial NB(μ, α) with α = 1/r, of mean μ and variance μ + αμ².

*Proof.* P(K = y) = ∫₀^∞ e^{−λ} λ^y / y! · λ^{r−1} e^{−λr/μ} (r/μ)^r / Γ(r) dλ = (r/μ)^r / (Γ(r) y!) · Γ(y + r) / (1 + r/μ)^{y+r}, which rearranges to the stated form. By the laws of total expectation and variance, E K = E λ = μ and Var K = E λ + Var λ = μ + μ²/r. ∎

The gamma is the biological variation of a gene's expression between replicates; the Poisson, the sampling of fragments. For a gene i in sample j, K_ij ~ NB(μ_ij, α_i), μ_ij = s_j q_ij, s_j the size factor, and log q_ij = Σ_r x_jr β_ir (Love, Huber and Anders 2014).

### 6.2 Size Factors

If most genes do not change, the ratio of a gene's count in sample j to its geometric mean over the samples is s_j over the geometric mean of the s, up to noise, for each such gene; the median over genes is robust to the minority that change (Anders and Huber 2010, eq. 5):

ŝ_j = median_i k_ij / (Π_v k_iv)^{1/m}.

On exactly scaled counts it returns s_j over the geometric mean of the s (check `differential`, to 10⁻¹²).

### 6.3 Iteratively Reweighted Least Squares

The log-likelihood of one gene, ℓ(β) = Σ_j log P(y_j; μ_j), with log μ_j = log s_j + x_jᵀβ, has score U = Σ_j x_j (y_j − μ_j)/(1 + αμ_j) and expected information I = Xᵀ W X, W = diag(μ_j/(1 + αμ_j)).

**Proposition 6.2.** The Fisher scoring step β + I⁻¹U is the weighted least squares solution (XᵀWX)⁻¹XᵀWz with the working response z_j = x_jᵀβ + (y_j − μ_j)/μ_j.

*Proof.* XᵀWz = XᵀWXβ + Σ_j x_j (μ_j/(1 + αμ_j)) (y_j − μ_j)/μ_j = Iβ + U. ∎

At the fit the score vanishes; the check holds Σ_j x_j(y_j − μ_j)/(1 + αμ_j) to zero, and at equal size factors in a two-group design the closed form β₀ = log ȳ₀, β₁ = log(ȳ₁/ȳ₀), which solves U = 0 there.

### 6.4 The Dispersion

Maximum likelihood estimates α with the β fitted at each α, and is biased down because β's fit absorbs part of the variation. The Cox–Reid adjusted profile likelihood (Cox and Reid 1987; McCarthy, Chen and Smyth 2012) subtracts half the log-determinant of β's information,

ℓ_CR(α) = ℓ(α, β̂(α)) − ½ log det(XᵀW(α)X),

the adjustment for a nuisance parameter orthogonal to α (the NB GLM's information is block-diagonal between β and α). The code maximises it over log α at the fitted means, by a grid and Brent's method on its derivative (digamma terms, and ½ α tr((XᵀWX)⁻¹ Xᵀ diag(μ²/(1 + αμ)²) X) from the adjustment). At α = 0.1 with 5 + 5 samples the median estimate of 400 genes is 0.092, against 0.072 by maximum likelihood.

With few samples each gene's estimate is poor, but many genes estimate how dispersion falls with the mean. DESeq2 fits the trend α_tr(μ̄) = a₀ + a₁/μ̄ (Poisson-like genes at high counts keep only the biological a₀) and takes the maximum a posteriori under a normal prior on log α centred on the trend, whose variance is the spread of the log estimates about the trend less the sampling variance of a log estimate, approximated by ψ′((m − p)/2): for a scaled χ²_ν variable, Var log(χ²_ν) = ψ′(ν/2) (Love et al. 2014, Methods). Genes far above the trend keep their own estimate, since shrinking them would make their tests liberal.

## 7. Testing Many Genes

### 7.1 The Wald Test

Under H₀: β_r = 0, the maximum likelihood estimate is asymptotically normal with the variance [(XᵀWX)⁻¹]_rr, so z = β̂_r / SE(β̂_r) is approximately N(0, 1), and the two-sided p-value is P(|Z| ≥ |z|) = P(χ²₁ ≥ z²) = Q(½, z²/2), the regularised upper incomplete gamma function. The check finds the p-values of 4000 null genes at a Kolmogorov–Smirnov distance 0.017 from uniform (the 1 % critical value for 4000 is 0.026), with 6.4 % below 0.05: with estimated dispersions and ten samples the normal approximation is slightly liberal.

### 7.2 The False Discovery Rate

Of m tests, let V be the number of true null hypotheses rejected and R the number of rejections; the false discovery rate is E[V / max(R, 1)].

**Theorem 7.1 (Benjamini and Hochberg 1995, Theorem 1).** With p_(1) ≤ … ≤ p_(m) the ordered p-values of independent tests, rejecting the hypotheses of p_(1), …, p_(k̂), where k̂ = max{i : p_(i) ≤ i q / m}, controls the false discovery rate at (m₀/m) q ≤ q, m₀ the number of true nulls.

**Proposition 7.2.** Let q_(i) = min_{k ≥ i} min(1, m p_(k)/k). Then q_(i) ≤ q exactly when i ≤ k̂: the adjusted p-values reject what the step-up procedure rejects, at every q at once.

*Proof.* If i ≤ k̂, then q_(i) ≤ m p_(k̂)/k̂ ≤ q. If i > k̂, every k ≥ i has p_(k) > k q/m, so every term exceeds q and so does their minimum. ∎

On Benjamini and Hochberg's own example of fifteen p-values (from Neuhaus et al. 1992), k̂ = 4 at q = 0.05: p_(4) = 0.0095 ≤ 4 · 0.05/15 = 0.0133, while p_(5) = 0.0201 > 0.0167 and no later p-value meets its bound; q_(4) = 15 · 0.0095/4 = 0.035625 and q_(5) = 15 · 0.0201/5 = 0.0603. The check reproduces both.

## 8. Where the Code Is

| Section | Code |
|---|---|
| §1 alignments, projections, the SP score | `msa/msa.hpp`, `msa/scores.hpp` (`pair_score`, `sp_score`) |
| §2.1–2.2 progressive alignment, profiles, their alignment | `msa/progressive.hpp` (`make_profile`, `align_profiles`, `merge`, `align_along`, `progressive`) |
| §2.3 refinement | `msa/progressive.hpp` (`refine`), `msa/scores.hpp` (`sp_across`) |
| §2.4 k-mer and corrected distances, UPGMA, weights | `msa/guide_tree.hpp` |
| §3 the centre star | `msa/centre_star.hpp`; the accuracy against a truth, `msa/scores.hpp` (`accuracy`) |
| §4 pseudoalignment | `expression/pseudoalign.hpp` (`ClassTable`, `Index`, `pseudoalign`) |
| §5 effective lengths, the likelihood, the EM, the bootstrap | `expression/transcriptome.hpp` (`effective_lengths`), `expression/quantify.hpp` |
| §6 the negative binomial, size factors, IRLS, dispersions | `expression/differential.hpp` |
| §7 the Wald test, Benjamini–Hochberg | `expression/differential.hpp` (`differential`, `benjamini_hochberg`) |

## References

- Altschul, S. F. (1989), Gap costs for multiple sequence alignment, *J. Theor. Biol.* 138, 297–309.
- Anders, S. and Huber, W. (2010), Differential expression analysis for sequence count data, *Genome Biology* 11, R106.
- Benjamini, Y. and Hochberg, Y. (1995), Controlling the false discovery rate, *J. R. Statist. Soc. B* 57, 289–300.
- Bray, N. L., Pimentel, H., Melsted, P. and Pachter, L. (2016), Near-optimal probabilistic RNA-seq quantification, *Nature Biotechnology* 34, 525–527.
- Bruynooghe, M. (1977), Méthodes nouvelles en classification automatique de données taxinomiques nombreuses, *Statistique et Analyse des Données* 2(3), 24–42.
- Carrillo, H. and Lipman, D. (1988), The multiple sequence alignment problem in biology, *SIAM J. Appl. Math.* 48, 1073–1082.
- Cox, D. R. and Reid, N. (1987), Parameter orthogonality and approximate conditional inference, *J. R. Statist. Soc. B* 49, 1–39.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), Maximum likelihood from incomplete data via the EM algorithm, *J. R. Statist. Soc. B* 39, 1–38.
- Edgar, R. C. (2004), MUSCLE: multiple sequence alignment with high accuracy and high throughput, *Nucleic Acids Research* 32, 1792–1797; and *BMC Bioinformatics* 5, 113.
- Feng, D.-F. and Doolittle, R. F. (1987), *J. Mol. Evol.* 25, 351–360; Higgins, D. G. and Sharp, P. M. (1988), *Gene* 73, 237–244.
- Gusfield, D. (1993), Efficient methods for multiple sequence alignment with guaranteed error bounds, *Bull. Math. Biol.* 55, 141–154; (1997), *Algorithms on Strings, Trees, and Sequences*, §14.6.
- Jukes, T. H. and Cantor, C. R. (1969), in *Mammalian Protein Metabolism* III, 21–132; Kimura, M. (1983), *The Neutral Theory of Molecular Evolution*, Cambridge University Press.
- Just, W. (2001), Computational complexity of multiple sequence alignment with SP-score, *J. Comput. Biol.* 8, 615–623.
- Kececioglu, J. and Starrett, D. (2004), Aligning alignments exactly, *Proceedings of RECOMB 2004*.
- Love, M. I., Huber, W. and Anders, S. (2014), Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2, *Genome Biology* 15, 550.
- McCarthy, D. J., Chen, Y. and Smyth, G. K. (2012), *Nucleic Acids Research* 40, 4288–4297.
- Müllner, D. (2011), Modern hierarchical, agglomerative clustering algorithms, arXiv:1109.2378; Murtagh, F. (1983), *The Computer Journal* 26, 354–359.
- Pachter, L. (2011), Models for transcript quantification from RNA-Seq, arXiv:1104.3889.
- Sokal, R. R. and Michener, C. D. (1958), *University of Kansas Science Bulletin* 38, 1409–1438.
- Thompson, J. D., Higgins, D. G. and Gibson, T. J. (1994), CLUSTAL W, *Nucleic Acids Research* 22, 4673–4680; Thompson, J. D., Plewniak, F. and Poch, O. (1999), *Nucleic Acids Research* 27, 2682–2690.
- Wang, L. and Jiang, T. (1994), On the complexity of multiple sequence alignment, *J. Comput. Biol.* 1, 337–348.
