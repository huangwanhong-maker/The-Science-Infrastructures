# Population Genetics

*Status: draft. Serves learning and use.* The models of genetic variation in a population that the module [`population/`](../population/) implements, and the reasoning behind its simulators and estimators: the Wright–Fisher model and its diffusion limit, the coalescent as its genealogy, mutation and the estimators of θ, Tajima's D, recombination and the ancestral recombination graph, linkage disequilibrium, Hardy–Weinberg equilibrium, F_ST and the island model, principal components, and admixture. The module's README lists the checks that hold each statement below to numbers; §12 says where the code implements it.

**Conventions.** A population has N diploid individuals, so 2N gene copies. Coalescent time is measured in units of 2N generations. θ = 4Nμ and ρ = 4Nr, with μ and r per generation for the region studied. Write a_n = Σ_{i=1}^{n−1} 1/i and b_n = Σ_{i=1}^{n−1} 1/i². A site is biallelic, with alleles A₁ and A₂, or ancestral and derived.

## 1. The Wright–Fisher Model

**Definition 1.1 (the Wright–Fisher chain).** Let X_t ∈ {0, 1, …, 2N} be the number of copies of A₁ in generation t, and p_t = X_t/(2N). Given p_t = p, the next generation draws its 2N genes independently from a gamete pool in which A₁ has frequency ψ(p):

X_{t+1} ~ Binomial(2N, ψ(p)).

Under neutrality ψ(p) = p. The module's life cycle composes three maps (Fisher 1930; Wright 1931):

- **selection,** with fitnesses 1 + s, 1 + hs, 1 for A₁A₁, A₁A₂, A₂A₂: ψ_s(p) = (p²(1 + s) + pq(1 + hs))/w̄;
- **mutation,** A₁ → A₂ at rate μ and A₂ → A₁ at rate ν: ψ_μ(p) = p(1 − μ) + qν;
- **migration** from the pool of all demes at rate m (the island model): ψ_m(p) = (1 − m)p + m p̄.

**Proposition 1.2 (neutral fixation).** Under neutrality, (p_t) is a martingale, absorbed at 0 or 1 in finite time almost surely, and P(fixation of A₁ | p₀) = p₀.

*Proof.* E[p_{t+1} | p_t] = p_t, since the binomial's mean is 2N p_t. From any interior state the chain is absorbed in one step with positive probability, so absorption happens in finite time almost surely. By optional stopping, p₀ = E[p_τ] = P(p_τ = 1). ∎

**Proposition 1.3 (the decay of heterozygosity).** Under neutrality, E[2p_{t+1}q_{t+1} | p_t] = (1 − 1/(2N)) · 2p_t q_t. Hence E[H_t] = (1 − 1/(2N))^t H₀ exactly, with H = 2pq.

*Proof.* For X ~ Binomial(2N, p),

E[X(2N − X)] = 2N · 2Np − (2Np(1 − p) + 4N²p²) = 2N(2N − 1) p(1 − p).

Divide by (2N)². ∎

*Remark.* Proposition 1.3 is an identity of the chain, not an approximation. The check `wright_fisher` therefore holds the simulation to it within sampling error alone.

**Proposition 1.4 (identity in the island model; Wright 1931).** Let F be the probability that two genes drawn from one deme are identical by descent, migrants being unrelated to residents (an infinite pool). Then

F′ = (1 − m)² (1/(2N) + (1 − 1/(2N)) F),

and its fixed point is F* = (1 − m)² / (2N − (2N − 1)(1 − m)²) ≈ 1/(1 + 4Nm) for small m.

*Proof.* Both genes descend from residents with probability (1 − m)². Given that, they copy the same parental gene with probability 1/(2N), and otherwise two distinct genes of the previous generation, which are identical with probability F. Solving F* = (1 − m)²(1/(2N) + (1 − 1/(2N))F*) gives F*. With (1 − m)² ≈ 1 − 2m, the denominator is about 1 + 2m(2N − 1) ≈ 1 + 4Nm. ∎

At equilibrium F* is F_ST = 1 − H_S/H_T among the demes (§9). The check `wright_fisher` finds 0.501 with 4Nm = 1 among a hundred islands.

## 2. The Diffusion Limit and the Probability of Fixation

Let N → ∞ with Ns, Nμ and Nm fixed, and measure time in units of 2N generations. The chain's frequency then converges to a diffusion on [0, 1] (Kimura 1962; Ewens 2004, ch. 4). Per generation it is characterised by the mean and variance of its change,

M(x) = E[Δp | p = x] = s x(1 − x)(h + (1 − 2h)x) + O(s²),  V(x) = Var[Δp | p = x] = x(1 − x)/(2N).

The mean comes from Δp = pq(p(w₁₁ − w₁₂) + q(w₁₂ − w₂₂))/w̄ with the fitnesses above.

**Proposition 2.1 (the probability of fixation; Kimura 1962).** The probability u(p) that A₁ fixes from frequency p is

u(p) = ∫₀ᵖ G(x) dx / ∫₀¹ G(x) dx,  G(x) = exp(−∫₀ˣ 2M/V) = exp(−2Ns(2hx + (1 − 2h)x²)).

For genic selection, h = ½, write s_K = s/2 for the advantage of the allele itself (fitnesses ≈ 1 + 2s_K, 1 + s_K, 1). Then

u(p) = (1 − e^{−4N s_K p}) / (1 − e^{−4N s_K}).

*Proof.* u(p_t) is a martingale: u(p) = E[u(p′) | p]. Expanding to second order in Δp gives the backward equation ½V u″ + M u′ = 0, with u(0) = 0 and u(1) = 1. Its solution is the normalised integral of the scale density G. At h = ½, G(x) = e^{−4N s_K x}. ∎

*Remark.* A new mutant has p = 1/(2N). For 4Ns_K ≫ 1 and s_K ≪ 1, u ≈ 1 − e^{−2s_K} ≈ 2s_K, which is Haldane's 2s. Without selection the neutral mean time to absorption is t̄(p) = −4N(p ln p + q ln q) generations (Kimura and Ohta 1969). Both are checked against the chain.

*Caveat.* The diffusion approximates the chain to first order in s and 1/N. At N = 100 and 4Ns_K between −1 and 4 the simulated fixation probabilities agree with u(p) within their sampling error, about one percent. Strong selection (Ns of order N) leaves the diffusion's domain.

## 3. The Coalescent

Sample n genes now and follow their ancestry backwards. In the neutral Wright–Fisher population two genes share a parent with probability 1/(2N) per generation. Among k genes some pair does so with probability C(k, 2)/(2N), and two pairs at once with probability O(N⁻²). In units of 2N generations this becomes a continuous-time process.

**Definition 3.1 (Kingman's n-coalescent).** The n-coalescent is the Markov process (Π_t)_{t≥0} on the equivalence relations of [n] = {1, …, n}, starting from the discrete relation, in which each pair of classes merges at rate 1, independently of the others.

**Proposition 3.2.** Let T_k be the time during which there are k classes.

1. The T_k are independent, with T_k ~ Exp(C(k, 2)).
2. The sequence of merges is independent of the times, and at each step every pair is equally likely.
3. E[T_MRCA] = 2(1 − 1/n) and Var[T_MRCA] = Σ_{k=2}^{n} 4/(k(k − 1))².
4. The total branch length L = Σ_k k T_k has E[L] = 2a_n and Var[L] = 4b_n.

*Proof.* Parts 1 and 2 follow from the competition of C(k, 2) independent exponential clocks. For part 3, Σ 2/(k(k − 1)) telescopes to 2(1 − 1/n). For part 4, k T_k ~ Exp((k − 1)/2), whose mean is 2/(k − 1) and whose variance is 4/(k − 1)². ∎

**Theorem 3.3 (Kingman 1982).** The genealogy of n genes sampled from the neutral Wright–Fisher population, with time in units of 2N generations, converges in law as N → ∞ to the n-coalescent. The same holds for any exchangeable reproduction law with finite offspring variance σ², with time in units of 2N/σ² generations.

*Remark (consistency and exchangeability).* The restriction of the n-coalescent to [m] ⊂ [n] is the m-coalescent. The family is projective, with an exchangeable limit on the partitions of ℕ. This is why a simulation may follow only the sample. It is also why a subsample's genealogy, and so its site frequency spectrum, can be read off the sample's, which is the projection of §4.

**Proposition 3.4 (time rescaling; Griffiths and Tavaré 1994).** Let the population size at time t in the past be N₀λ(t), with λ positive and piecewise continuous. Then the genealogy is the n-coalescent run on the clock

Λ(t) = ∫₀ᵗ du/λ(u).

That is, the merge times are Λ⁻¹(τ_i), where the τ_i are those of the constant-size process.

*Proof.* In a smaller population pairs meet faster, at rate 1/λ(t) per pair. With k lineages, P(no merge in (t, t′]) = exp(−C(k, 2)(Λ(t′) − Λ(t))). This is the constant-size waiting law on the clock Λ. ∎

The module's demography is piecewise exponential, λ(t) = λ_e e^{−α_e(t − t_e)} on epoch e, so Λ and Λ⁻¹ are closed forms on each epoch. For two genes, P(T₂ > t) = e^{−Λ(t)}, which the check holds the simulation to under growth, a bottleneck and three epochs. Under strong exponential growth Λ rises steeply, every merge happens at about the same time, and the tree is nearly a star (Slatkin and Hudson 1991). External branches are then long and singletons in excess (§5).

## 4. Mutation and the Estimators of θ

**Definition 4.1 (the infinite-sites model).** Given the genealogy, mutations form a Poisson process of rate θ/2 per unit of branch length (each lineage mutates at 2Nμ = θ/2 per 2N generations), each at a new site. A mutation on the branch above node v is carried by exactly the leaves below v. The **site frequency spectrum** ξ_i, for i = 1, …, n − 1, is the number of sites at which i of the n sequences carry the derived allele.

**Proposition 4.2 (Watterson 1975).** The number of segregating sites S = Σ_i ξ_i has

E[S] = θa_n,  Var[S] = θa_n + θ²b_n.

*Proof.* S | L ~ Poisson(θL/2). By the laws of total expectation and total variance, E[S] = (θ/2)E[L] and Var[S] = (θ/2)E[L] + (θ/2)²Var[L]. Then apply Proposition 3.2. ∎

**Proposition 4.3 (π; Tajima 1983).** The mean number of differences between two sequences, π = Σ_i 2i(n − i)/(n(n − 1)) · ξ_i, has E[π] = θ and Var[π] = b₁θ + b₂θ², where b₁ = (n + 1)/(3(n − 1)) and b₂ = 2(n² + n + 3)/(9n(n − 1)).

*Proof of the mean.* π is the average over the C(n, 2) pairs of their differences. A pair's genealogy is the 2-coalescent: two branches of length T₂ ~ Exp(1), carrying Poisson(θT₂) differences, whose mean is θ. The variance is Tajima's. ∎

**Proposition 4.4 (the expected spectrum; Fu 1995).** E[ξ_i] = θ/i for i = 1, …, n − 1.

*Reason.* E[ξ_i] = (θ/2) E[L_i], where L_i is the length of the branches above exactly i leaves. While there are k lineages, a given one of them subtends i leaves with the probability of a Pólya urn, C(n − i − 1, k − 2)/C(n − 1, k − 1). Summing k · E[T_k] times that probability over k gives E[L_i] = 2/i (Fu 1995). ∎

*Remark (the hyperplane of unbiased estimators).* A linear estimator θ̂_w = Σ_i w_i ξ_i is unbiased exactly when Σ_i w_i/i = 1: an affine hyperplane in the space of weights. The estimators of the module are points of it:

| Estimator | w_i |
|---|---|
| θ_W (Watterson) | 1/a_n |
| π (Tajima) | 2i(n − i)/(n(n − 1)) |
| θ_H (Fay and Wu) | 2i²/(n(n − 1)) |

Under the neutral, constant-size coalescent the difference of any two has mean zero. A test of neutrality is such a difference, normalised by its standard deviation. It is sensitive where the two weightings differ: rare alleles for θ_W against π, high-frequency derived alleles for π against θ_H (Achaz 2009).

*Remark (missing data and projection).* A site typed in N_j ≥ n sequences, k of them derived, contributes to the spectrum of a random subsample of n the hypergeometric law of k in n draws without replacement from N_j. By exchangeability this is exactly the expected spectrum of a subsample (Marth et al. 2004). The module projects every site to the smallest sample size before computing D and H; where nothing is missing the projection is the identity.

## 5. Tajima's D

Let d = π − S/a_n. Both terms estimate θ, so E[d] = 0, and Tajima (1989) shows that Var[d] = c₁θ + c₂θ², where

c₁ = b₁ − 1/a₁,  c₂ = b₂ − (n + 2)/(a₁n) + a₂/a₁²,

with a₁ = a_n and a₂ = b_n. Since E[S] = a₁θ and E[S(S − 1)] = Var[S] + E[S]² − E[S] = (a₁² + a₂)θ², the statistic e₁S + e₂S(S − 1), with e₁ = c₁/a₁ and e₂ = c₂/(a₁² + a₂), is an unbiased estimator of Var[d]. Tajima's statistic is

D = d / √(e₁S + e₂S(S − 1)).

**Proposition 5.1 (Tajima's bounds).** For fixed n, as S → ∞,

- D → D_min = (2/n − 1/a₁)/√e₂ when every site is a singleton;
- D → D_max = (n/(2(n − 1)) − 1/a₁)/√e₂ when every site is at frequency n/2 (n even).

*Proof.* A singleton contributes 2/n to π, and a site at n/2 contributes n/(2(n − 1)). The S² term dominates the denominator. ∎

*Remark (reading D).* A negative D means an excess of rare variants. This is what population growth, a recent selective sweep or purifying selection produce: star-like genealogies (§3) put mutations on external branches. A positive D means an excess of intermediate frequencies, from balancing selection, population structure or a recent contraction. Fay and Wu's H = π − θ_H is negative when derived alleles at high frequency are in excess, the signature of hitchhiking near a sweep. Both are tests of the whole model (neutrality, constant size, no structure), and none says alone which assumption failed.

*Caveat.* D is a ratio of dependent quantities and is not exactly standardised. Under neutrality its mean is slightly negative and its variance below 1. The check finds −0.07 and 0.85 for n = 10 and θ = 5.

## 6. Recombination and the Ancestral Recombination Graph

With recombination, different positions of a chromosome have different genealogies. Backwards in time each lineage carries **ancestral material**: the set A ⊆ [0, 1) of positions at which it is ancestral to some sampled gene, a finite union of intervals. Two events change the lineages:

- **recombination** at a breakpoint x splits a lineage into two, carrying A ∩ [0, x) and A ∩ [x, 1). Each lineage recombines at rate (ρ/2) times the length of its chromosome. In Hudson's algorithm only breakpoints inside (min A, max A) are simulated, since outside it one part carries nothing the sample can see;
- **common ancestry** merges two lineages into one carrying A ∪ B. On A ∩ B the two coalesce in the marginal genealogies, and a position carried by only one remaining lineage has found its most recent common ancestor and is dropped.

The process, until no material is left, is the ancestral recombination graph restricted to the sample (Hudson 1983; Griffiths and Marjoram 1997).

**Definition 6.1 (marginal tree).** The marginal tree at x ∈ [0, 1) is the genealogy traced by the lineages that carry x.

**Proposition 6.2 (Hudson 1983).** Every marginal tree is distributed as the n-coalescent.

*Proof.* Let K_x(t) be the number of lineages carrying x. A recombination gives x to exactly one of the two parts, so it leaves K_x unchanged. Each pair of lineages meets at rate 1, in particular each pair carrying x. So K_x is the block-counting process of the n-coalescent, and by symmetry each merge joins a uniform pair. ∎

The check holds E[T_MRCA] and E[L] at a fixed position to Kingman's values at ρ = 20.

**Proposition 6.3 (two loci; Griffiths 1981).** For two sequences and two loci separated by ρ, measured as 4Nr between them,

Cov[T_A, T_B] = (ρ + 18)/(ρ² + 13ρ + 18),

while Var[T_A] = Var[T_B] = 1.

*Reason.* Backwards, the pair of sequences moves among three configurations:

- two lineages each carrying both loci;
- three lineages, one carrying both loci and two carrying one each;
- four lineages each carrying one locus.

Coalescence and recombination move it between them. E[T_A T_B] solves the linear equations of first-step analysis on this chain (Griffiths 1981; Durrett 2008, ch. 3). ∎

The covariance falls from 1 at ρ = 0 towards 0 as ρ grows. This is the genealogical root of the decay of linkage disequilibrium (§7).

*Remark (the tree sequence).* The module records the graph as nodes with times and edges (left, right, parent, child): node `parent` is the parent of `child` on [left, right). Adjacent marginal trees share most of their edges, so the record costs little more than the differences between consecutive trees. The trees are read by sweeping the chromosome, removing the edges that end at x and inserting those that begin there (Kelleher, Etheridge and McVean 2016, Algorithm T). This is the succinct tree sequence of msprime. Mutations fall on edges, at positions in the edge's interval, so each belongs to exactly one marginal tree.

## 7. Linkage Disequilibrium

Let x_AB, x_Ab, x_aB, x_ab be the frequencies of the four haplotypes of two biallelic sites, and p_A, p_B the frequencies of A and B. Then

D = x_AB − p_A p_B = x_AB x_ab − x_Ab x_aB,

the covariance of the indicators of A and B on a random haplotype. The two forms agree because x_AB(1 − x_AB − x_Ab − x_aB) = x_AB x_ab. The other measures normalise D:

- **r²** = D²/(p_A q_A p_B q_B) is the squared correlation (Hill and Robertson 1968). r² n is approximately χ²₁ under independence, and r² bounds the power of one site to tag another.
- **D′** = D/D_max divides D by its extreme value at the given allele frequencies, which is min(p_A q_B, q_A p_B) for D > 0 and min(p_A p_B, q_A q_B) for D < 0. |D′| = 1 exactly when at most three haplotypes occur (Lewontin 1964).

In an infinite population D decays as (1 − r)^t. In a finite one drift regenerates it, and the ratio of expectations σ²_d = E[D²]/E[p_A q_A p_B q_B] tends to (10 + C)/(22 + 13C + C²), with C = 4Nr (Ohta and Kimura 1971). This decline with distance is what the check observes among the coalescent's sites.

**Proposition 7.1 (unphased genotypes).** Let a 3 × 3 table count the dosages at two sites. Then:

1. The likelihood of the haplotype frequencies depends on the data only through the table.
2. Only the double heterozygotes, d in number, have an ambiguous phase: AB/ab or Ab/aB.
3. The allele frequencies are observed exactly, so the likelihood is a function of x_AB alone on [max(0, p_A + p_B − 1), min(p_A, p_B)].

The EM algorithm iterates

x_AB ← (c_AB + d φ)/(2n),  φ = x_AB x_ab / (x_AB x_ab + x_Ab x_aB),

where c_AB counts the AB haplotypes known from unambiguous genotypes and n is the number of individuals. Each step raises the likelihood, and its fixed points satisfy a cubic equation in x_AB (Hill 1974).

*Proof.* Each individual's allele count at each site is observed, so the update preserves p_A and p_B. φ is the posterior probability of the phase AB/ab. Monotonicity is that of every EM algorithm (Dempster, Laird and Rubin 1977). Clearing the quadratic denominator of φ in the fixed-point equation gives the cubic. ∎

A cubic has at most two local maxima on an interval, so the module starts the EM at linkage equilibrium and near both ends, and keeps the best. The check confirms that the result is the maximum over a fine grid of x_AB.

## 8. Hardy–Weinberg Equilibrium and the Exact Test

Random union of gametes gives the genotype frequencies p², 2pq and q². Conditioned on the allele counts, the test asks whether the genotypes are a random pairing of the alleles.

**Proposition 8.1 (Levene 1949).** Given N individuals carrying n_A copies of allele A and n_B = 2N − n_A of B, under Hardy–Weinberg equilibrium

P(n_AB | N, n_A) = 2^{n_AB} N! / (n_AA! n_AB! n_BB!) · n_A! n_B! / (2N)!.

*Proof.* Conditionally the 2N alleles are arranged uniformly among the 2N slots, two per individual, which gives (2N)!/(n_A! n_B!) arrangements. Those with given genotype counts number N!/(n_AA! n_AB! n_BB!), times 2^{n_AB} for the order within each heterozygote. ∎

The ratio of consecutive terms gives the recurrences of Wigginton, Cutler and Abecasis (2005), used from the mode outwards.

**Proposition 8.2 (the exact test keeps its size).** Let P_HWE(k) = Σ_{j : P(j) ≤ P(k)} P(j), and let the test reject when P_HWE ≤ α. Then the probability of rejection under equilibrium is at most α.

*Proof.* The rejection region R is a down-set for the order by probability: if k ∈ R and P(j) ≤ P(k), then P_HWE(j) ≤ P_HWE(k) ≤ α. Let k* be an outcome of R of largest probability. Then R ⊆ {j : P(j) ≤ P(k*)}, so P(R) ≤ P_HWE(k*) ≤ α. ∎

The χ² test has no such guarantee. In Wigginton's Table 1 (N = 100, n_A = 21) it rejects with probability 0.070 at nominal 0.05, while the exact test does so with probability 0.010. The check reproduces both figures.

## 9. F_ST, the Island Model and Balding–Nichols

Wright's fixation index compares the variance of allele frequencies among populations with its maximum:

F_ST = Var(p_k) / (p̄(1 − p̄)) = 1 − H_S/H_T,

with H_S the mean heterozygosity within populations and H_T that of the pooled population. Hudson et al. (1992) write it through the probabilities that two genes differ, within populations (H_w) and between them (H_b): F_ST = 1 − H_w/H_b.

**Proposition 9.1 (Hudson's estimator is unbiased term by term).** For two populations with sample frequencies p̂₁ and p̂₂ from n₁ and n₂ independent genes,

N = (p̂₁ − p̂₂)² − p̂₁q̂₁/(n₁ − 1) − p̂₂q̂₂/(n₂ − 1)

has expectation (p₁ − p₂)², and D = p̂₁q̂₂ + p̂₂q̂₁ has expectation p₁q₂ + p₂q₁.

*Proof.* E[(p̂₁ − p̂₂)²] = (p₁ − p₂)² + p₁q₁/n₁ + p₂q₂/n₂, and E[p̂q̂] = pq(1 − 1/n), so E[p̂q̂/(n − 1)] = pq/n. The samples are independent, so E[p̂₁q̂₂] = p₁q₂. ∎

So Σ_j N_j / Σ_j D_j is a ratio of unbiased sums, and is consistent as sites accumulate. A mean of per-site ratios is not: rare variants, with D_j small, dominate it (Bhatia et al. 2013).

Weir and Cockerham's θ̂ comes from an analysis of variance of the allele indicators, nested as gametes within individuals within populations. Its components a, b, c are unbiased for the variances between populations, between individuals within populations, and between the gametes of an individual. θ̂ = Σa/Σ(a + b + c) estimates the coancestry of genes within a population relative to the whole.

**Proposition 9.2 (Balding–Nichols).** In the island model with migrant pool frequency p, the stationary law of a deme's frequency is

Beta(4Nm p, 4Nm q) = Beta(p(1 − F)/F, q(1 − F)/F),  F = 1/(1 + 4Nm),

with mean p and variance F pq (Wright 1931; Balding and Nichols 1995).

*Proof.* The stationary density of the diffusion with drift m(p − x) and variance x(1 − x)/(2N) is proportional to x^{4Nmp − 1}(1 − x)^{4Nmq − 1}. With a = p(1 − F)/F and b = q(1 − F)/F, a + b = (1 − F)/F, and the Beta variance ab/((a + b)²(a + b + 1)) equals pq F. ∎

The Balding–Nichols model takes this law as a model of divergence from an ancestral frequency, with F as its parameter, and the module's simulator draws populations from it. Two populations each diverged by F from the same ancestor have E[(f₁ − f₂)²] = 2Fpq and E[f₁q₂ + f₂q₁] = 2pq, so both estimators estimate F. The check finds both within 10⁻⁴ of F = 0.1.

## 10. Principal Components as the Spectrum of a Kinship Matrix

Standardise the genotypes, M_ij = (g_ij − μ_j)/√(p_j(1 − p_j)) with p_j = μ_j/2, and form the n × n matrix X = MMᵀ/m over m sites. Entry X_ik is the average over sites of the standardised product of the dosages of samples i and k, an estimate of their genetic relatedness: X is a kinship matrix. Each column of M is centred, so 1ᵀM = 0 and X1 = 0, and the all-ones vector is a structural null vector.

**Without structure,** the columns are independent, and X is approximately a Wishart matrix with m degrees of freedom on the (n − 1)-dimensional complement of 1, divided by m. Its eigenvalues fill the Marčenko–Pastur law, and the largest fluctuates at the law's upper edge on the scale of the Tracy–Widom law TW₁ (Johnstone 2001).

**With K populations,** the expected kinship adds a matrix of rank K − 1 to the noise: a block structure, since samples of one population are more related than samples of two. Patterson et al. (2006) show that for divergence τ ≪ 1 and enough samples the covariance has K − 1 eigenvalues that grow with the sample, while the rest stay at 1 + O(τ). This is a spiked Wishart model. A spike is visible only above a threshold that falls as m grows, the phase change of Baik, Ben Arous and Péché (2005), which Patterson discusses.

**The test.** Patterson normalises the largest eigenvalue to the edge and spread of the noise, and estimates the effective number of independent markers n′ by moments (their eq. 10). Linkage between markers makes n′ smaller than m. The statistic x is then compared with TW₁, and the k-th eigenvalue is tested after removing the first k − 1. The module computes x and n′. Its check finds statistics above 70 for the two axes of three populations, and a null statistic distributed about as TW₁: a mean of −1.03 against −1.21, and 9 of 200 samples above the 95% point.

*Remark.* Patterson's n′ = (m + 1)(Σλ)²/((m − 1)Σλ² − (Σλ)²), with m the number of individuals, differs from the exact Wishart moment match by one part in m. The module uses it as printed.

## 11. Admixture as a Mixture Model

**Definition 11.1 (the admixture model).** Each individual i has ancestry proportions q_i on the simplex Δ^{K−1}, and each population k an allele frequency f_kj at each site j. Each of the two alleles of i at j comes from population k with probability q_ik, independently, and is the alternate allele with probability f_kj. Hence

g_ij ~ Binomial(2, π_ij),  π_ij = Σ_k q_ik f_kj,

and the log-likelihood is ℓ(q, f) = Σ_ij [g_ij log π_ij + (2 − g_ij) log(1 − π_ij)] on (Δ^{K−1})ⁿ × [0, 1]^{Km} (Pritchard, Stephens and Donnelly 2000; Tang et al. 2005; Alexander et al. 2009).

**Proposition 11.1 (the EM).** Let the latent datum be the population of origin of each allele. Then:

- the E step gives an alternate allele the posterior q_ik f_kj / π_ij of coming from k, and a reference allele the posterior q_ik(1 − f_kj)/(1 − π_ij);
- the M step sets q_ik to the expected fraction of i's alleles from k, and f_kj to the expected fraction of alternate alleles among those from k.

Each step raises ℓ.

*Proof.* The complete-data likelihood is multinomial in the origins and Bernoulli in the alleles given their origins. Its maximisers are the expected proportions named. Monotonicity follows from Dempster, Laird and Rubin (1977). ∎

*Remark (symmetry).* The symmetric group S_K acts on the parameters by permuting populations, (q, f) ↦ (qσ, σ⁻¹f), and ℓ is invariant. Maximum-likelihood estimates are therefore orbits, and an estimate is compared with the truth up to this action: the check minimises the error over the six permutations of S₃. The likelihood also has local maxima, so the module offers random restarts.

*Remark (acceleration).* EM converges linearly, with a rate set by the fraction of missing information, which is large when populations are close. SQUAREM (Varadhan and Roland 2008) extrapolates along the secant of two EM steps, θ′ = θ₀ − 2αr + α²v. It projects back onto the parameter space, takes one EM step to stabilise, and falls back to the plain step if the likelihood would fall. So it keeps EM's monotonicity and global convergence while taking larger steps. ADMIXTURE uses quasi-Newton steps for the same purpose.

## 12. Where the Code Implements It

| Section | Header | Functions |
|---|---|---|
| 1, 2 | `wright_fisher.hpp` | `wright_fisher_trajectory`, `wright_fisher_replicates`, `kimura_fixation`, `fixation_probability`, `heterozygosity_decay`, `island_fst`, `island_fst_exact` |
| 3 | `coalescent.hpp` | `Demography` (Λ, Λ⁻¹), `coalescent`, `expected_tmrca`, `expected_tree_length` |
| 4, 5 | `diversity.hpp` | `theta_watterson`, `nucleotide_diversity`, `theta_h`, `site_frequency_spectrum` (projection), `tajima_constants`, `tajima_d`, `fay_wu_h` |
| 6 | `coalescent.hpp` | `coalescent` (Hudson's algorithm), `for_each_tree` (Algorithm T), `haplotypes`, `marginal_trees`, `tmrca_at` |
| 7 | `linkage.hpp` | `linkage_phased`, `PackedHaplotypes`, `linkage_unphased`, `ld_decay` |
| 8 | `hardy_weinberg.hpp` | `heterozygote_law`, `hwe_exact`, `hwe_chi_square` |
| 9 | `diversity.hpp`, `structure.hpp` | `fst_hudson`, `fst_weir_cockerham`, `window_statistics`, `simulate_structure` |
| 10 | `structure.hpp` | `kinship`, `pca`, `tracy_widom_statistic` |
| 11 | `structure.hpp` | `admixture`, `admixture_log_likelihood` |
| [haplotypes.md](haplotypes.md) §5 | `sweep.hpp` | `sweep`: a selective sweep forward in time, its one-generation law, fixation and hitchhiking |

## 13. References

- Achaz, G. (2009), Frequency spectrum neutrality tests: one for all and all for one, *Genetics* 183, 249–258.
- Baik, J., Ben Arous, G. and Péché, S. (2005), Phase transition of the largest eigenvalue for nonnull complex sample covariance matrices, *Annals of Probability* 33, 1643–1697.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), Maximum likelihood from incomplete data via the EM algorithm, *J. R. Statist. Soc. B* 39, 1–38.
- Durrett, R. (2008), *Probability Models for DNA Sequence Evolution*, 2nd ed., Springer.
- Ewens, W. J. (2004), *Mathematical Population Genetics I*, 2nd ed., Springer.
- Fisher, R. A. (1930), *The Genetical Theory of Natural Selection*, Oxford.
- Griffiths, R. C. and Marjoram, P. (1997), An ancestral recombination graph, in *Progress in Population Genetics and Human Evolution*, IMA Vol. 87, 257–270.
- Johnstone, I. M. (2001), On the distribution of the largest eigenvalue in principal components analysis, *Annals of Statistics* 29, 295–327.
- Levene, H. (1949), On a matching problem arising in genetics, *Annals of Mathematical Statistics* 20, 91–94.
- Marth, G. T., Czabarka, E., Murvai, J. and Sherry, S. T. (2004), The allele frequency spectrum in genome-wide human variation data, *Genetics* 166, 351–372.
- Ohta, T. and Kimura, M. (1971), Linkage disequilibrium between two segregating nucleotide sites under the steady flux of mutations in a finite population, *Genetics* 68, 571–580.
- Pritchard, J. K., Stephens, M. and Donnelly, P. (2000), Inference of population structure using multilocus genotype data, *Genetics* 155, 945–959.
- The papers of the module's README (§7): Alexander et al. 2009; Balding and Nichols 1995; Bhatia et al. 2013; Excoffier and Slatkin 1995; Fay and Wu 2000; Fu 1995; Griffiths 1981; Griffiths and Tavaré 1994; Hill 1974; Hill and Robertson 1968; Hudson 1983, 1990, 2002; Hudson, Slatkin and Maddison 1992; Kelleher, Etheridge and McVean 2016; Kimura 1962; Kimura and Ohta 1969; Kingman 1982; Lewontin 1964; Nei and Li 1979; Patterson, Price and Reich 2006; Slatkin and Hudson 1991; Tajima 1983, 1989; Tang et al. 2005; Varadhan and Roland 2008; Watterson 1975; Weir and Cockerham 1984; Wigginton, Cutler and Abecasis 2005; Wright 1931.
