# Haplotypes: Copying Models, Phasing, Imputation and Haplotype Scans

*Status: draft. Serves learning and use.* The models behind the module [`haplotypes/`](../haplotypes/) and the sweep simulator of [`population/`](../population/):

- Li and Stephens's copying model of a haplotype on a reference panel, and the algebra that makes its recursions linear in the panel;
- imputation as the model's posterior predictive;
- the pair model of a diploid and phasing as a joint maximum;
- the extended haplotype homozygosity and the integrated haplotype score;
- the selective sweep that the scans are checked against.

The module's README lists the checks that hold each statement to numbers; §6 says where the code implements it.

**Conventions.**

- A panel is K haplotypes h_1, …, h_K at N biallelic sites, with alleles 0 and 1 (0 ancestral where polarity matters) and genetic positions g_1 ≤ … ≤ g_N.
- A target is x = (x_1, …, x_N) with x_j ∈ {0, 1, −1}, −1 missing; a diploid's genotype is g_j ∈ {0, 1, 2, −1}.
- [K] = {1, …, K}, u is the uniform distribution on [K], J the K × K matrix of ones and P = J/K.
- Writing α for a row vector over [K], αP = (Σ_k α_k) u.

## 1. The Copying Model

**Definition 1.1 (Li and Stephens 2003).** Let ρ ≥ 0 and ε ∈ [0, ½]. The copying chain on a panel is the Markov chain Z_1, …, Z_N on [K] with Z_1 ~ u and, for j ≥ 2, with d_j = g_j − g_{j−1},

T_j = P(Z_j = · | Z_{j−1} = ·) = (1 − r_j) I + r_j P,  r_j = 1 − e^{−ρ d_j / K}.

Given the path, the alleles are independent, x_j = h_{Z_j}(j) with probability 1 − ε and 1 − h_{Z_j}(j) otherwise. A missing x_j contributes the factor 1. Li and Stephens write the miscopying through a mutation parameter θ: with probability θ/(K + θ) a mutation makes the allele uniform on {0, 1}. So ε = θ/(2(K + θ)), and the value they suggest is θ̃ = (Σ_{m=1}^{K−1} 1/m)⁻¹.

*Remark (what the model approximates).* The conditional distribution of a further haplotype given K sampled ones, under the coalescent with recombination, has no tractable form. Li and Stephens's chain keeps its two qualitative features: the new haplotype is locally close to one of the others, and which one changes along the chromosome at the rate of recombination. The scale ρ/K reflects that with more haplotypes the nearest one is closer in time, so fewer recombinations separate it.

**Proposition 1.2 (the jumps form a semigroup).** Let T(t) = e^{−t} I + (1 − e^{−t}) P for t ≥ 0. Then T(t) = exp(t(P − I)), T(s)T(t) = T(s + t), and T_j = T(ρ d_j / K). So the copying chain is the continuous-time Markov chain with generator (ρ/K)(P − I) on the map, observed at the sites.

*Proof.* P is idempotent, so for the commuting pair −I and P, exp(t(P − I)) = e^{−t} exp(tP) = e^{−t}(I + (e^t − 1)P) = T(t). The semigroup law follows. ∎

**Corollary 1.3 (missing sites are harmless).** Inserting a site at which the target is missing, or deleting one, leaves the joint law of the other sites' states and alleles unchanged.

*Proof.* At a missing site the emission is 1, so summing out its state composes the two adjacent transitions: T(a)T(b) = T(a + b), the transition across the joined interval. ∎

This is why imputation (§2) may run the chain over every panel site, typed or not.

**Proposition 1.4 (the O(K) step).** Let α be a probability vector and β any vector over [K]. Then αT_j = (1 − r_j)α + r_j u and T_jβ = (1 − r_j)β + r_j (mean β) 𝟙.

*Proof.* αP = (Σ_k α_k) u = u, and Pβ = (mean β) 𝟙. ∎

The scaled forward variables α̂_j sum to one (Rabiner 1989, §V.A), so the forward step is α_j(l) = e_j(l)((1 − r_j) α̂_{j−1}(l) + r_j/K), with e_j(l) the emission. The backward step uses the second identity on e_j β̂_j. The posterior is γ_j(k) ∝ α̂_j(k) β̂_j(k). Each costs O(K) per site, so O(NK) in all.

**Proposition 1.5 (the Viterbi step).** For any δ ∈ [0, ∞)^K, max_k δ(k) T_j(k, l) = max((1 − r_j + r_j/K) δ(l), (r_j/K) max_k δ(k)).

*Proof.* T_j(l, l) = 1 − r_j + r_j/K and T_j(k, l) = r_j/K for k ≠ l, so the left side is max((1 − r_j + r_j/K)δ(l), (r_j/K) max_{k≠l} δ(k)). Enlarging the second maximum to all k adds only (r_j/K)δ(l). Since r_j ≤ 1, this is at most (1 − r_j + r_j/K)δ(l), already in the first term. ∎

The same argument shows that the back-pointer is exact. When the jump term wins strictly, its argmax is not l, so the transition it records is the one whose weight it used.

**Proposition 1.6 (calibration under the model).** If x is drawn from the model, E[γ_j(Z_j)] = E[Σ_k γ_j(k)²] for every j.

*Proof.* Given x, Z_j has law γ_j, so E[γ_j(Z_j) | x] = Σ_k γ_j(k)². Take expectations. ∎

The check draws mosaics from the model, hides 30% of their alleles, and holds the two sides to each other within sampling error. This tests the posterior itself, not only the likelihood.

## 2. Imputation

**Proposition 2.1 (the posterior predictive).** Let the target be typed on a set O of sites, and let j ∉ O. Then

P(x_j = 1 | x_O) = Σ_k γ_j(k) P(1 | h_k(j)),  P(1 | h) = 1 − ε if h = 1, ε otherwise,

where γ_j is the posterior of the chain run with every site outside O missing.

*Proof.* Given Z_j, x_j is independent of x_O and of the other states. So P(x_j = 1 | x_O) = Σ_k P(Z_j = k | x_O) P(x_j = 1 | Z_j = k). By Corollary 1.3, P(Z_j = · | x_O) is the posterior of the chain with the untyped sites missing. ∎

**Corollary 2.2.** Under the model the dosage d_j = P(x_j = 1 | x_O) is unbiased: E[x_j − d_j] = 0. For a diploid, d_j = E[g_j | g_O] from the pair posterior (§3), and E[g_j − d_j] = 0.

*Proof.* d_j = E[x_j | x_O]; apply the tower property. ∎

Against the truth of a real population, imputation is only as good as the panel's relatedness to the target. The check measures it on the coalescent's samples:

- the concordance of the best guess;
- the squared error of the dosage;
- r², the squared correlation of dosage and truth across targets at a site, the usual measure of an imputed site's quality.

Each is compared with what the panel's allele frequency alone gives.

## 3. The Pair Model and Phasing

**Definition 3.1 (the pair model).** A diploid's two haplotypes x and y copy the panel by independent chains with the same jumps and ε. The hidden state is the ordered pair (A_j, B_j) ∈ [K]², with transition T_j ⊗ T_j. The observation is g_j = x_j + y_j, so

P(g_j | a, b) = Σ_{x+y=g_j} P(x | h_a(j)) P(y | h_b(j)).

A model of unphased genotypes as pairs of panel haplotypes is how the first version of IMPUTE imputed them (Marchini et al. 2007).

**Proposition 3.2 (the O(K²) step).** Identify a vector α over [K]² with the K × K matrix α(a, b). Then α(T ⊗ T) = TᵀαT = TαT, and with T = (1 − r)I + rP,

TαT = (1 − r)² α + (1 − r) r (Pα + αP) + r² PαP.

Here (αP)(a, b) is the row sum Σ_{b′} α(a, b′)/K, (Pα)(a, b) the column sum Σ_{a′} α(a′, b)/K, and PαP the total over K².

*Proof.* (α(T ⊗ T))(c, d) = Σ_{a,b} α(a, b) T(a, c) T(b, d) = (TᵀαT)(c, d), and T is symmetric. Expand the product. ∎

So a step costs O(K²): one pass for the row and column sums, one for the update. The dense recursion over the K² × K² matrix, O(K⁴) per site, is the reference.

**Proposition 3.3 (the maximum step).** With s = 1 − r + r/K and q = r/K,

max_{a,b} δ(a, b) T(a, c) T(b, d) = max{ s² δ(c, d), s q max_b δ(c, b), s q max_a δ(a, d), q² max_{a,b} δ(a, b) }.

*Proof.* The weight T(a, c)T(b, d) is s² when (a, b) = (c, d), sq when exactly one coordinate agrees, and q² otherwise. Each relaxed maximum on the right includes some pairs at a weight below their own. Since q ≤ s, every such term is bounded by the term in which the pair has its own weight. ∎

**Proposition 3.4 (exchangeability of the copies).** The swap of the two copies, σ: (A, B, x, y) ↦ (B, A, y, x), preserves the joint law of the states, the haplotypes and the genotypes. Hence the posterior satisfies γ_j(a, b) = γ_j(b, a). At a heterozygous site the posterior probability that x carries allele 1 is exactly ½, whatever the data.

*Proof.* T ⊗ T commutes with the swap of coordinates, the initial law is uniform on [K]², the two copies miscopy alike, and g = x + y = y + x is unchanged. So the joint law is σ-invariant, and so is the posterior given g. In particular P(x_j = 1 | g) = P(y_j = 1 | g). At a heterozygous site x_j + y_j = 1, so the two probabilities sum to 1 and each is ½. ∎

Phase is therefore not a marginal quantity: it lives in the correlation between sites, which a path carries.

**Definition 3.5 (phasing).** The phase of a diploid is the joint maximum

(Â, B̂, x̂, ŷ) = argmax P(A, B, x, y)  over x_j + y_j = g_j at the typed sites.

The alleles maximise per site given the states, so this is the Viterbi recursion of Proposition 3.3 with the emission replaced by its largest term, max_{x+y=g} P(x | h_a) P(y | h_b). At a missing genotype the alleles copied (x = h_a, y = h_b) multiply every pair's emission by the same (1 − ε)², so they are chosen without changing the path. Ties at a heterozygous site whose copied haplotypes carry the same allele give allele 1 to x.

**Definition 3.6 (switch errors).** Among the sites heterozygous in both the truth and the inference, a switch is a pair of consecutive such sites whose relative phase differs. The switch error rate is switches/(sites − 1). Random phase gives ½.

*Remark (larger panels).* The pair model costs O(NK²), which limits it to panels of a few hundred haplotypes. For thousands, the copying is restricted to the panel haplotypes that share long matches with the target near each site. The positional Burrows–Wheeler transform (Durbin 2014) finds the longest such matches in time linear in the size of the panel's N × K matrix. It is being written in morphe's `formal`; restricting the copying chain to its matches is future work, and nothing here depends on it.

## 4. Extended Haplotype Homozygosity and the Integrated Haplotype Score

**Definition 4.1 (EHH; Sabeti et al. 2002).** Let C_a be the haplotypes carrying allele a at the core site c, with n_a = |C_a| ≥ 2. For a site j, let π_j be the partition of C_a by the alleles at the sites between c and j inclusive. Then

EHH_a(j) = Σ_{B ∈ π_j} C(|B|, 2) / C(n_a, 2),

the probability that two distinct carriers drawn at random are identical over the whole interval. A carrier missing an allele in the interval is put in a block of its own.

**Proposition 4.2.** EHH_a(c) = 1. Going away from the core on either side, the partitions π_j form a descending chain in the lattice of partitions of C_a, and EHH_a is non-increasing.

*Proof.* At the core every carrier has allele a, so π_c is the single block C_a. Extending the interval by a site refines each block by its allele there. The number of pairs within blocks can only fall under refinement. ∎

The walk computes π_j from π_{j−1} by one pass over the carriers still in blocks of two or more, since singletons never contribute again. The check holds it to the definition by pairs, exactly.

**Definition 4.3 (iHH).** iHH_a is the integral of EHH_a over the genetic map on both sides of the core, by the trapezoid rule over the sites walked:

iHH_a = Σ_sides Σ_{i ≥ 1} |g_{j_i} − g_{j_{i−1}}| (EHH_a(j_{i−1}) + EHH_a(j_i))/2,  j_0 = c.

Each side's walk stops at the first site where EHH_a has fallen to the cut-off or below, that interval included. The cut-off is 0.05 (Voight et al. 2006). A side that reaches the end of the data first is *truncated*: its integral is a lower bound, and a scan leaves such cores out unless asked to keep them.

**Definition 4.4 (iHS; Voight, Kudaravalli, Wen and Pritchard 2006).** The unstandardised score of a core site is U = ln(iHH_A / iHH_D), with A the ancestral and D the derived allele. Let the derived frequencies be binned. For the sites in a bin b, with m_b their mean of U and s_b their standard deviation (divisor the bin's count n_b),

iHS = (U − m_b)/s_b.

**Proposition 4.5 (standardisation).** Within each bin with n_b ≥ 2 and s_b > 0, the iHS of the bin's own sites have mean 0 and variance 1 (divisor n_b).

*Proof.* The map U ↦ (U − m_b)/s_b is affine. It sends the bin's mean to 0 and scales its variance s_b² to 1. ∎

*Remark (reading the sign).* The derived allele is the younger, so under neutrality it sits on longer haplotypes the rarer it is. Binning by frequency removes that dependence. What is left is the departure of one site from sites of the same frequency. A strongly negative iHS means the derived allele's haplotypes are unusually long for its frequency: it rose faster than drift would carry it, the signature of an incomplete sweep. A complete sweep leaves no polymorphism at the selected site, and iHS cannot see it there.

## 5. A Sweep Forward in Time

The scans need a truth with selection in it. `population/sweep.hpp` simulates one forward in time, in the model below. Write M = 2N.

**Definition 5.1 (the sweep).**

- *Start.* The population is M chromosomes [0, 1). Its standing variation is a coalescent sample of all M, with θ = 4Nμ and ρ = 4Nr for the whole chromosome. At generation 0 a beneficial allele arises at one position on one chromosome drawn uniformly.
- *A generation.* Each chromosome of the next generation is formed from two of the present one, drawn independently with probability proportional to fitness (1 + s with the allele, 1 without). It copies them alternately between k ~ Poisson(ρ/(4N)) crossovers at uniform positions. No new mutation arises.
- *Stop.* An attempt that loses the allele starts again from the same standing variation. The first that reaches the stop frequency is sampled.

**Proposition 5.2 (one generation).** Given the frequency p of the beneficial allele, the next generation's count is Binomial(M, ψ(p)) with ψ(p) = p(1 + s)/(1 + sp). The deterministic recursion p ↦ ψ(p) multiplies the odds by 1 + s, so it takes 2 ln(M − 1)/ln(1 + s) generations from 1/M to 1 − 1/M.

*Proof.* The selected site of a child comes from one of its two parents (which one depends on the crossovers). Each parent is a carrier with probability p(1 + s)/(p(1 + s) + 1 − p) = ψ(p), independently of everything else in that child, and the children are independent. Then ψ(p)/(1 − ψ(p)) = (1 + s) p/(1 − p), and the odds go from 1/(M − 1) to M − 1. ∎

**Corollary 5.3 (fixation).** In the diffusion limit the drift is ψ(x) − x = s x(1 − x)/(1 + sx) ≈ s x(1 − x), and the variance is x(1 − x)/M. By Kimura's formula (theory/population.md, Proposition 2.1, with s_K = s) the probability of fixing from one copy is

u(1/M) = (1 − e^{−2s})/(1 − e^{−4Ns}).

The reciprocal of the mean number of attempts estimates it, and the check holds the two together.

**Proposition 5.4 (drift on neutral sites).** Let a neutral site's allele in each child come from a parent drawn independently of that allele, as for a site unlinked to the selected one. Then the expected pairwise diversity of the population falls by the factor (1 − 1/M) each generation, as in theory/population.md, Proposition 1.3.

*Proof.* The child's allele is then a draw from the parental frequency, and the children are independent. The count is binomial, and Proposition 1.3 of the population note applies. ∎

Near the selected site the neutral alleles are not drawn independently of fitness. They travel with the beneficial allele until recombination separates them, and the sweep removes diversity there (Maynard Smith and Haigh 1974): this is hitchhiking. The check measures both regions against the same standing variation:

- near the site, π falls to about 14% of drift's share (Proposition 5.4);
- far from it, π is drift's share within a few per cent.

*Remark (why iHS sees a partial sweep).* The carriers of the beneficial allele descend from one founder within the G generations of the sweep. Their haplotypes differ only where recombination with non-carriers broke them, so EHH_D decays over map distances of order 1/(rate of crossovers × time), against the far shorter scale of neutral haplotypes of the same frequency. The check scores partial sweeps to 0.7 against a neutral background of the same θ and ρ. 30 of 39 selected sites have iHS < −2, against 4.6% of neutral sites with |iHS| > 2.

*Caveat (no new mutation).* Without mutation during the sweep, drift lowers the neutral diversity everywhere by (1 − 1/M)^G, about 20% for G = 200 and M = 1000. The checks compare against that factor rather than against equilibrium. The structured coalescent of a sweep (Kaplan, Hudson and Langley 1989) would avoid the cost of forward simulation, and new mutation would remove the loss; neither is written.

## 6. Where the Code Implements It

| Section | Header | Functions |
|---|---|---|
| 1 | `haplotypes/copying.hpp` | `copying_chain`, `li_stephens_theta`, `miscopying`, `forward`, `backward`, `log_likelihood`, `posterior`, `viterbi`, `path_log_probability`, `forward_dense`, `sample_mosaic` |
| 2 | `haplotypes/imputation.hpp` | `impute_haplotype`, `impute_genotypes`, `dosage_r2` |
| 3 | `haplotypes/phasing.hpp` | `pair_log_likelihood`, `pair_log_likelihood_dense`, `for_each_pair_posterior`, `pair_posterior`, `phase`, `phased_log_probability`, `switch_errors` |
| 4 | `haplotypes/ehh.hpp` | `ehh`, `ehh_at`, `ihh`, `ihs_scan`, `fit_ihs_standardisation`, `standardise` |
| 5 | `population/sweep.hpp` | `sweep` |

## 7. References

- Durbin, R. (2014), Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (PBWT), *Bioinformatics* 30, 1266–1272.
- Kaplan, N. L., Hudson, R. R. and Langley, C. H. (1989), The "hitchhiking effect" revisited, *Genetics* 123, 887–899.
- Kimura, M. (1962), On the probability of fixation of mutant genes in a population, *Genetics* 47, 713–719.
- Li, N. and Stephens, M. (2003), Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data, *Genetics* 165, 2213–2233.
- Marchini, J., Howie, B., Myers, S., McVean, G. and Donnelly, P. (2007), A new multipoint method for genome-wide association studies by imputation of genotypes, *Nature Genetics* 39, 906–913.
- Maynard Smith, J. and Haigh, J. (1974), The hitch-hiking effect of a favourable gene, *Genetical Research* 23, 23–35.
- Rabiner, L. R. (1989), A tutorial on hidden Markov models and selected applications in speech recognition, *Proceedings of the IEEE* 77, 257–286.
- Sabeti, P. C. et al. (2002), Detecting recent positive selection in the human genome from haplotype structure, *Nature* 419, 832–837.
- Voight, B. F., Kudaravalli, S., Wen, X. and Pritchard, J. K. (2006), A map of recent positive selection in the human genome, *PLoS Biology* 4, e72.
