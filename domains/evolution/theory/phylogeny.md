# Phylogeny: Trees, Metrics and Markov Chains on Trees

*Status: draft. Serves learning and research.* The mathematics behind the evolution domain's `phylogeny` module: trees as systems of splits and as metrics, the reconstruction of a tree from its metric, continuous-time Markov chains on a tree and the likelihood they give a set of aligned sequences, Felsenstein's pruning as a sum–product over a semiring of which parsimony is the tropical case, rates across sites, the search of tree space, the bootstrap, and the molecular clock. The code is [`phylogeny/`](../phylogeny/); each section says which header implements it.

## 1. Trees and Their Splits

Let X be a finite set of taxa, |X| = n ≥ 3.

**Definition 1.1 (X-tree).** An *unrooted X-tree* is a tree T together with a bijection from X onto its leaves, whose internal vertices have degree at least three. It is *binary* when every internal vertex has degree three. A *rooted* X-tree has in addition a distinguished vertex ρ, the root, of degree at least two, every other internal vertex of degree at least three.

**Definition 1.2 (splits).** Removing an edge e of an unrooted X-tree T disconnects it into two subtrees, and so partitions X into two non-empty blocks A | X∖A, the *split* of e. The set of splits of T is Σ(T). A split is *trivial* when one block is a single taxon; every X-tree has the n trivial splits, and a binary one has n − 3 others.

**Definition 1.3 (compatibility).** Two splits A | A′ and B | B′ are *compatible* when at least one of A ∩ B, A ∩ B′, A′ ∩ B, A′ ∩ B′ is empty.

**Theorem 1.4 (the splits-equivalence theorem; Buneman 1971; Semple and Steel 2003, §3.1).** A set Σ of X-splits containing the trivial ones is Σ(T) for some X-tree T if and only if its splits are pairwise compatible, and then T is unique up to isomorphism.

*Proof (sketch).* Two splits of one tree are compatible, since the four blocks lie on the two sides of two edges and the side of one edge away from the other contains a whole block of the other. Conversely, fix x₀ ∈ X and represent each split by its block not containing x₀, a *cluster*. Compatibility of A | A′ and B | B′ with x₀ ∈ A′ ∩ B′ means that A ∩ B = ∅, or A ⊆ B, or B ⊆ A: the clusters form a laminar family, a hierarchy, and a hierarchy is the set of clusters of a unique tree rooted at x₀, each cluster the leaves below a vertex. ∎

*Remark (the code).* `splits.hpp` holds a split as the side without the first taxon, a bit set over the taxa; `compatible` is the test of containment or disjointness in the proof, and `tree_from_splits` is the proof's construction: clusters sorted by size, each attached under the smallest cluster containing it, the first taxon next to the root.

**Definition 1.5 (Robinson–Foulds distance).** For X-trees T₁, T₂, d_RF(T₁, T₂) = |Σ(T₁) △ Σ(T₂)|, counting non-trivial splits (Robinson and Foulds 1981).

It is a metric on X-trees, being the symmetric-difference metric on their sets of splits, and by Theorem 1.4 it vanishes only between isomorphic trees. For binary trees it is at most 2(n − 3); `normalised_robinson_foulds` divides by that.

**Proposition 1.6 (majority-rule consensus).** Let T₁, …, T_m be X-trees. The splits that occur in more than m/2 of them are pairwise compatible, and so (Theorem 1.4) are the splits of a tree, the majority-rule consensus tree.

*Proof.* If two splits each occur in more than m/2 of the trees, some tree has both, and splits of one tree are compatible. ∎

Margush and McMorris (1981) characterised this consensus axiomatically: it is the median of the profile for the Robinson–Foulds metric among trees whose splits are majority splits. Below one half, the splits are no longer guaranteed compatible; `majority_consensus` with a lower threshold adds them greedily by frequency, keeping each one compatible with those already kept.

## 2. Tree Metrics and Their Reconstruction

**Definition 2.1.** A metric d on X is a *tree metric* when there is an X-tree T with positive edge lengths ℓ such that d(x, y) is the length of the path from x to y. It is an *ultrametric* when moreover T can be rooted so that every leaf lies at the same distance from the root.

**Theorem 2.2 (the four-point condition; Buneman 1974; Semple and Steel 2003, ch. 7).** A metric d on X is a tree metric if and only if for every x, y, z, w ∈ X, of the three sums

d(x, y) + d(z, w),  d(x, z) + d(y, w),  d(x, w) + d(y, z),

the largest two are equal. It is an ultrametric if and only if for every x, y, z the largest two of d(x, y), d(x, z), d(y, z) are equal (the three-point condition).

The four-point condition reads the quartet tree off the sums: if xy | zw is the split of the quartet, the paths from x to z, x to w, y to z and y to w all cross the central edge of length c, so the two crossing sums exceed d(x, y) + d(z, w) by 2c.

**Proposition 2.3 (UPGMA on an ultrametric).** If d is an ultrametric, UPGMA returns its tree, with its heights.

*Proof.* The closest pair is a cherry of the ultrametric tree, of height half their distance. By the three-point condition every other taxon is equidistant from the two, so the mean of their distances (average linkage, or any weighting) is the distance from the cherry's cluster, and the reduced matrix is the ultrametric of the tree with the cherry contracted to a leaf. Induction on n. ∎

On a matrix that is not ultrametric, UPGMA still returns an ultrametric tree (each join's height exceeds its children's, since the merged distances are means of distances no smaller than the joined one), which is a good guide tree but a biased phylogeny when rates differ among lineages. Code: `distance_trees.hpp`, `upgma`.

**Neighbour joining.** With r clusters and R_i = Σ_k d(i, k), let

Q(i, j) = (r − 2) d(i, j) − R_i − R_j.

**Theorem 2.4 (Saitou and Nei 1987; Studier and Keppler 1988).** If d is a tree metric, a pair minimising Q is a cherry of its tree; joining it at the new vertex u with d(i, u) = d(i, j)/2 + (R_i − R_j)/(2(r − 2)) and d(u, k) = (d(i, k) + d(j, k) − d(i, j))/2 leaves the tree metric of the tree with that cherry replaced by u. So neighbour joining recovers a tree from its tree metric exactly.

*Remark (what Q measures).* Saitou and Nei introduced the criterion as the sum S_ij of branch lengths of the tree in which i and j form a cherry and every other taxon hangs from one central vertex, estimated by least squares; in Studier and Keppler's form, S_ij = Q(i, j)/(2(r − 2)) + Σ_{k<l} d(k, l)/(r − 2), so minimising Q minimises S. Gascuel and Steel (2006, Mol. Biol. Evol. 23, 1997–2000) showed that neighbour joining is the greedy algorithm for the *balanced* minimum-evolution criterion, which explains its consistency, and Bryant (2005, "On the uniqueness of the selection criterion in neighbor-joining", Journal of Classification 22, 3–15) that Q is, up to positive affine change, the only criterion of its kind that is consistent.

**Theorem 2.5 (Atteson 1999).** If T has positive edge lengths with shortest edge ℓ_min, and the matrix D satisfies |D(x, y) − d_T(x, y)| < ℓ_min/2 for all x, y, then neighbour joining returns the topology of T. The radius 1/2 is the largest any method can have.

This is why NJ works on estimated distances: they need only be close to a tree metric, entry by entry. The checks perturb an additive matrix by 0.45 ℓ_min and recover its topology.

BIONJ (Gascuel 1997) keeps NJ's choice of pair and its branch lengths and replaces the reduction by the convex combination d(u, k) = λ(d(i, k) − d(i, u)) + (1 − λ)(d(j, k) − d(j, u)), with λ minimising the variance of the new estimate under a model in which the variance of a distance is proportional to it and covariances to the shared path; on a tree metric every λ gives the same reduction, so BIONJ is exact too.

**Distances as inverted expectations.** Under a Markov model (§3), the expected proportion of differing sites between two sequences separated by a path of length t is a function p(t); a distance estimator is its inverse applied to the observed proportion. For JC69, P_xx(t) = 1/4 + (3/4)e^{−4t/3}, so p(t) = (3/4)(1 − e^{−4t/3}) and

d = −(3/4) ln(1 − 4p/3).

If the rate r of a site is Γ(α, α)-distributed (mean 1), the expectation of e^{−4rt/3} is the moment generating function (1 + 4t/(3α))^{−α}, and inverting gives the Γ form d = (3/4) α [(1 − 4p/3)^{−1/α} − 1]: each −ln x becomes α(x^{−1/α} − 1). The same substitution in K80 and TN93 gives their Γ forms (Jin and Nei 1990; Tamura and Nei 1993, eq. 15). Code: `distance.hpp`.

*Remark (LogDet).* For any Markov chain on a tree, the joint frequency matrix of two leaves is J = Π_x P_{path}, a product of the matrices along the path, and det is multiplicative: −ln det J − (½)(ln det Π_x + ln det Π_y) is additive along paths whatever the rates, even without stationarity, which makes the LogDet distance consistent for tree building where other distances fail (Steel 1994, Appl. Math. Lett. 7, 19–24; Lockhart et al. 1994). Under a stationary model det P(t) = e^{tr(Q) t}, so the distance (with the normalisation −¼) is −tr(Q) t/4: t itself only when −tr Q = 4, as under JC69. The checks verify this identity under TN93.

## 3. Continuous-Time Markov Chains on a Tree

**Definition 3.1.** A *rate matrix* on the states {A, C, G, T} is a real 4 × 4 matrix Q with Q_ij ≥ 0 for i ≠ j and rows summing to zero. Its transition matrices P(t) = e^{Qt} form a semigroup, P(s + t) = P(s)P(t), of stochastic matrices. A distribution π is *stationary* when πQ = 0, and the chain is *reversible* when π_i Q_ij = π_j Q_ji, that is, when ΠQ is symmetric (Π = diag π).

**Proposition 3.2 (spectral form).** If Q is reversible with π > 0, then S = Π^{1/2} Q Π^{−1/2} is symmetric, Q has real eigenvalues 0 = λ₄ > λ₃ ≥ λ₂ ≥ λ₁ (for an irreducible chain), and with S = U Λ Uᵀ,

P(t) = Π^{−1/2} U e^{Λt} Uᵀ Π^{1/2},  dP/dt = Q P(t),  d²P/dt² = Q² P(t).

*Proof.* (Π^{1/2} Q Π^{−1/2})_ij = π_i^{1/2} Q_ij π_j^{−1/2}, and reversibility makes this symmetric. Q is similar to S, so their exponentials are similar. The eigenvalue 0 has the eigenvector Π^{1/2}𝟙 = √π of S; the others are negative since Q generates a contraction on the complement of the stationary direction. ∎

`models.hpp` computes U by calx's symmetric eigensolver once per model and then P(t) for every t at the cost of four exponentials; it sets the eigenvalue 0 and its vector √π exactly, so that P(t) → 𝟙π (every row π) as t → ∞ to rounding. A reversible rate matrix is Q_ij = s_ij π_j with symmetric exchangeabilities s (GTR, Tavaré 1986); JC69, K80, F81, HKY85 and TN93 are its special cases. Q is scaled so that −Σ π_i Q_ii = 1, the expected number of substitutions per unit time at stationarity: branch lengths are then expected substitutions per site.

**Definition 3.3 (the tree model).** Let T be a rooted tree with edge lengths t_e. A state at the root drawn from π, and along each edge u → v a state at v drawn from row x_u of P(t_e), give a joint distribution of states at all vertices,

Pr(x) = π(x_ρ) Π_{edges u→v} P_{x_u x_v}(t_e),

and its marginal at the leaves is the probability of a column of the alignment. Sites are independent given the tree.

**Proposition 3.4 (the pulley principle; Felsenstein 1981).** For a reversible chain started at stationarity, the distribution of the leaves depends only on the unrooted tree and its edge lengths, not on where the root is placed.

*Proof.* Moving the root from ρ to a child c across an edge of length t replaces π(x_ρ) P_{x_ρ x_c}(t) by π(x_c) P_{x_c x_ρ}(t), which is equal by reversibility (π_i P_ij(t) = π_j P_ji(t) follows from π_i Q_ij = π_j Q_ji by the power series). Placing it inside an edge splits P(t) into P(s)P(t − s), the semigroup property. ∎

Hence only the unrooted tree with its lengths is identifiable from data, and the two edges at a root of degree two are seen only through their sum. The likelihood code works on unrooted trees and fits each branch as though the root lay on it (§4).

## 4. Pruning as a Sum–Product

Fix a column. The likelihood is a sum over the states of the internal vertices,

L = Σ_{x_internal} π(x_ρ) Π_{u→v} P_{x_u x_v}(t_e),

with the leaves' states fixed (a missing state summed over, contributing 1). There are 4^{n−1} terms for n − 1 internal vertices of a binary rooted tree.

**Proposition 4.1 (Felsenstein's pruning algorithm, 1981).** Define, from the leaves up, L_v(x) = 1[x is v's state] at a leaf and

L_v(x) = Π_{c child of v} Σ_y P_xy(t_c) L_c(y)

at an internal vertex. Then L = Σ_x π(x) L_ρ(x), computed in O(n · 16) operations.

*Proof.* The distributive law: the sum over the states of the subtree below v factors through v's state, since the subtrees of its children share no vertex. Induction from the leaves. ∎

*Remark (the semiring).* Nothing in the proof uses more than a commutative semiring (R, ⊕, ⊗): the same recursion computes ⊕ over labellings of ⊗ along edges of edge weights w_e(x, y). In (ℝ≥0, +, ×) with w = P it is the likelihood; in (ℝ ∪ {∞}, min, +) with w = c(x, y) a cost it is Sankoff's parsimony (§7); with the costs −ln P_xy(t) and (min, +) it is the most probable joint reconstruction of ancestral states (Pupko et al. 2000, Mol. Biol. Evol. 17, 890–896). This is the sum–product (belief propagation) algorithm on a tree, exact because a tree has no cycles, and the same structure as the genomics domain's alignments over semirings (the design note, §1).

**Underflow.** A column's likelihood shrinks geometrically with the number of leaves (below 10⁻³⁰⁸, the smallest double, beyond about 500 leaves). `likelihood.hpp` multiplies a partial by 2²⁵⁶ whenever its largest entry falls below 2⁻²⁵⁶ and counts the rescalings, which is to compute in the semiring of pairs (mantissa, exponent), a representation of (ℝ>0, ×). The checks evaluate a tree of 700 leaves whose single columns have log-likelihoods near −930.

**Messages from both sides.** For an edge v below p, let D_v = L_v be the partial from v's subtree and W_v the conditional likelihood at p of everything outside it (the product of p's other children's messages and of p's own outer message through its edge). By Proposition 3.4,

L(t_v) = Σ_x π_x W_v(x) Σ_y P_xy(t_v) D_v(y),

the likelihood as a function of one branch length with the rest fixed. With the spectral form, P_xy(t) = Σ_k A_xk e^{λ_k t} B_ky, so

L(t) = Σ_c w_c Σ_k e^{λ_k r_c t} G_ck,  G_ck = (Σ_x π_x W_v(x) A_xk)(Σ_y B_ky D_v(y)),

a sum of exponentials whose derivatives in t are immediate. The code forms G once per branch and runs Newton's method on d ln L/dt inside a bracket of its change of sign; a sweep visits the branches in preorder, recomputing W for each child as it descends and D for each vertex as it returns, so every branch is fitted on the current tree and the likelihood never decreases. Coordinate ascent of this kind converges to a stationary point of the likelihood in the branch lengths, which for a fixed topology is typically, though not always, the global maximum (Steel 1994, "The maximum likelihood point for a phylogenetic tree is not unique", Syst. Biol. 43, 560–564).

**Proposition 4.2 (two sequences).** For two sequences under JC69 with k differences among L sites, the maximum-likelihood length of the path between them is the JC69 distance −(3/4) ln(1 − 4k/(3L)), when k/L < 3/4.

*Proof.* The log-likelihood is (L − k) ln(1/4 + (3/4)e^{−4t/3}) + k ln(1/4 − (1/4)e^{−4t/3}) + const, a function of p(t) = (3/4)(1 − e^{−4t/3}) through the binomial log-likelihood (L − k) ln(1 − p) + k ln(p/3); it is maximised at p = k/L, and t ↦ p(t) is a bijection of (0, ∞) onto (0, 3/4). ∎

The checks verify both the closed form of this likelihood and that the fitted length equals the JC69 distance.

## 5. Rates Across Sites

Sites evolve at different rates: a site's rate r scales every branch, P(rt). Yang (1994) took r ~ Γ(α, α) (mean 1, variance 1/α) and replaced the integral over r by K equiprobable categories, each represented by its mean.

**Proposition 5.1 (the category means).** Let b₀ = 0 < b₁ < … < b_K = ∞ be the quantiles of Γ(α, α) at k/K. The mean rate of category k is

r_k = K [P(α + 1, α b_k) − P(α + 1, α b_{k−1})],

where P is the regularised lower incomplete gamma function.

*Proof.* With g(x; a, β) = β^a x^{a−1} e^{−βx}/Γ(a), x g(x; α, α) = g(x; α + 1, α) (since Γ(α + 1) = αΓ(α)), so ∫_{b_{k−1}}^{b_k} x g(x; α, α) dx = P(α + 1, αb_k) − P(α + 1, αb_{k−1}), and each category has probability 1/K. ∎

For α = 1 the rate is exponential, the quantiles are ln(4/3), ln 2, ln 4, and the four means are exactly 1 − 3 ln(4/3), 1 + 3 ln(4/3) − 2 ln 2, 1 and 1 + ln 4; the checks hold the code to these and to Yang's published values for α = 0.5 (0.0334, 0.2519, 0.8203, 2.8944). The likelihood of a column is the mixture Σ_c w_c L(column | rates r_c). A proportion p of invariant sites is one more category, of rate 0, whose likelihood is Σ_x π_x [x agrees with every observed state], with the other rates divided by 1 − p to keep the mean 1. Code: `models.hpp`, `discrete_gamma` and `gamma_rates`.

## 6. Searching the Space of Trees

**Proposition 6.1.** There are (2n − 5)!! = 1 · 3 · 5 ⋯ (2n − 5) unrooted binary X-trees for |X| = n ≥ 3.

*Proof.* A binary tree on k leaves has 2k − 3 edges; inserting leaf k + 1 on any of them gives a binary tree on k + 1 leaves, and every such tree arises from exactly one tree on k leaves (delete the leaf and suppress its neighbour) and one edge. So the count is Π_{k=3}^{n−1} (2k − 3). ∎

The proof is also an algorithm: inserting each taxon on a uniformly chosen edge gives the uniform distribution on topologies (`random_unrooted_tree`; the checks test it by χ² on the 15 trees of five leaves), and recursing over every edge enumerates the trees (the parsimony checks score all 10,395 of eight taxa). Since finding the maximum-parsimony tree (Foulds and Graham 1982) and the maximum-likelihood tree (Roch 2006, IEEE/ACM TCBB 3, 92–94) is NP-hard, the code searches locally.

**Nearest-neighbour interchanges.** An internal edge separates four subtrees, (A, B | C, D); the other two binary trees on them, (A, C | B, D) and (A, D | B, C), are its NNI neighbours, and the NNI graph on binary topologies is connected. `likelihood.hpp` scores each neighbour with the five branches around the edge refitted from the four subtrees' partials, makes the gaining interchanges that touch disjoint parts of the tree, refits, and climbs to a local maximum (Guindon and Gascuel 2003). The quartet's fit is a lower bound on the refitted tree's likelihood, so the best single interchange always gains when it is predicted to.

**Branch and bound.** Adding a taxon to a tree never lowers its parsimony score: restricting an optimal labelling of the larger tree to the smaller (suppressing the new leaf's neighbour) leaves at most the same number of changes. So a partial tree whose score reaches the best complete score found cannot lead to a better one, and the exact search (Hendy and Penny 1982) prunes it.

## 7. Parsimony: the Tropical Case

In the semiring (ℝ ∪ {∞}, min, +), Proposition 4.1 becomes Sankoff's algorithm: with a cost c(x, y) of a change,

S_v(x) = Σ_{c child of v} min_y [c(x, y) + S_c(y)],  score = min_x S_ρ(x),

and an optimal labelling is read off by tracing the minimising y's back from the root. With unit costs (c = 1 − δ), the functions S_v take only two values that matter, the minimum m_v and m_v + 1, and S_v is determined by m_v and the set F_v = argmin S_v: this is Fitch's algorithm, F_v = F_a ∩ F_b at no cost if non-empty, F_a ∪ F_b at cost one otherwise (Fitch 1971), and for vertices of higher degree, the states in the most children's sets at a cost of the number of children whose sets lack them (Hartigan 1973). The checks compare all three with enumeration of every internal labelling.

*Remark (parsimony as a limit of likelihood).* With weights e^{−c(x,y)/ε}, the (+, ×) likelihood becomes, as ε → 0, −ε ln L → min over labellings of the total cost: Maslov dequantisation turns the sum–product into the min–sum. Parsimony is the likelihood's tropical shadow, and it inherits the shadow's weakness: it ignores branch lengths, and is inconsistent when two long branches are separated by a short one (Felsenstein 1978, Syst. Zool. 27, 401–410).

## 8. The Bootstrap

The columns of an alignment are taken as a sample from an unknown distribution over columns; a bootstrap replicate is a sample of L columns drawn with replacement from the L observed, and a tree is built from each replicate by the method under study (Felsenstein 1985). Since every method here depends on the alignment only through its pattern counts, a replicate is the vector of multinomial(L; n_q/L) counts of the patterns, and `bootstrap.hpp` draws it from the stream (seed, BOOTSTRAP, b), builds each replicate's tree as a task of a calx task graph, and counts in a final task, in order, the fraction of replicate trees that contain each split of the reference tree. The result does not depend on the number of threads.

A bootstrap support is the probability that the method recovers a split on data like these, not the probability that the split is true; Hillis and Bull (1993) found it conservative above 70 % under their conditions, and Efron, Halloran and Holmes (1996, PNAS 93, 13429–13434) explained the bias geometrically, as the curvature of the boundaries between regions of tree space.

## 9. The Molecular Clock and Dated Tips

Under a strict clock of rate u, a tip sampled at time t_i lies at distance d_i = u(t_i − t_r) from the root, sampled at t_r (Rambaut et al. 2016, eq. 1). For sequences sampled over time (viruses, ancient DNA), the regression of root-to-tip distance on sampling date estimates u as its slope and t_r as the date where it meets the axis.

The root must first be placed. Put it on the edge above v at distance x from v: a tip below v has distance a_i + x, any other a_i − x, where a_i is its distance from v. The regression's sums are then polynomials in x, S_td(x) = c₀ + c₁x (linear) and S_dd(x) = v₀ + v₁x + v₂x² (quadratic), so that the residual sum of squares S_dd − S_td²/S_tt is a quadratic in x, minimised at its vertex or an end; and R²(x) = S_td²/(S_tt S_dd), whose derivative vanishes where S_td = 0 or where

2S_td′ S_dd − S_td S_dd′ = 0 ⟺ (2c₁v₀ − c₀v₁) + (c₁v₁ − 2c₀v₂) x = 0,

a linear equation, its quadratic terms cancelling. Every branch is tried and the best root kept, each in closed form (`dating.hpp`). The points of the regression are not independent, since tips share branches, so its residuals measure clock-likeness for exploration, not a likelihood.

## 10. References

- Buneman, P. (1971), The recovery of trees from measures of dissimilarity, in *Mathematics in the Archaeological and Historical Sciences*, 387–395; (1974), A note on the metric properties of trees, *J. Combin. Theory B* 17, 48–50.
- Semple, C. and Steel, M. (2003), *Phylogenetics*, Oxford University Press.
- Felsenstein, J. (1981), *J. Mol. Evol.* 17, 368–376; (1985), *Evolution* 39, 783–791; (2004), *Inferring Phylogenies*, Sinauer.
- Saitou, N. and Nei, M. (1987), *Mol. Biol. Evol.* 4, 406–425; Studier, J. A. and Keppler, K. J. (1988), *Mol. Biol. Evol.* 5, 729–731; Atteson, K. (1999), *Algorithmica* 25, 251–278; Gascuel, O. (1997), *Mol. Biol. Evol.* 14, 685–695; Gascuel, O. and Steel, M. (2006), *Mol. Biol. Evol.* 23, 1997–2000.
- Yang, Z. (1994), *J. Mol. Evol.* 39, 306–314; Yang, Z. (2014), *Molecular Evolution: A Statistical Approach*, Oxford University Press.
- Tavaré, S. (1986), *Lectures on Mathematics in the Life Sciences* 17, 57–86; Lockhart, P. J., Steel, M. A., Hendy, M. D. and Penny, D. (1994), *Mol. Biol. Evol.* 11, 605–612.
- Fitch, W. M. (1971), *Syst. Zool.* 20, 406–416; Hartigan, J. A. (1973), *Biometrics* 29, 53–65; Sankoff, D. (1975), *SIAM J. Appl. Math.* 28, 35–42; Hendy, M. D. and Penny, D. (1982), *Math. Biosci.* 59, 277–290; Foulds, L. R. and Graham, R. L. (1982), *Adv. Appl. Math.* 3, 43–49.
- Guindon, S. and Gascuel, O. (2003), *Syst. Biol.* 52, 696–704.
- Rambaut, A., Lam, T. T., Max Carvalho, L. and Pybus, O. G. (2016), *Virus Evolution* 2, vew007.
