# Structural Variants, Copy Number and Loss of Heterozygosity

*Status: draft. Serves learning and use.* How a genome that departs from its reference by more than a few bases is written, what such a departure leaves in the reads of a sample, and how it is read back from them: a rearranged genome as a set of adjacencies of oriented breakends; the signatures each rearrangement leaves in pairs, split reads and depth; the probabilities of breakpoints and the likelihood of genotypes; copy number as a step function and its segmentation; the alleles beside the depth, which tell a loss of heterozygosity from a deletion; and how calls are compared with a truth. The code is the module [`structural`](../structural/) of the genomics domain, with its simulation in [`simulate`](../simulate/) (`structural.hpp`, `pairs.hpp`); it stands on the mapper and the caller of [mapping and variants](mapping-and-variants.md), on the hidden Markov chains of [`hmm`](../hmm/) and on tyche's change points.

## 1. A Rearranged Genome as Adjacencies of Oriented Breakends

### 1.1 Breakends

Let a reference be a finite family of words r₁, …, r_m (its contigs) over {A, C, G, T, N}, and write |r_c| for the length of contig c. Positions are 0-based.

**Definition 1.1 (breakend).** A *breakend* is a triple (c, p, o) with 0 ≤ p < |r_c| and o ∈ {+, −}: the side of base p of contig c that faces a junction. (c, p, +) is the junction after base p, the reference to its left kept; (c, p, −) the junction before it, the reference to its right kept.

The reference itself is the set of its *reference adjacencies* {(c, p, +), (c, p + 1, −)} for 0 ≤ p < |r_c| − 1: each base is joined to the next.

**Definition 1.2 (oriented segment, chromosome).** An *oriented segment* is an interval [a, b) of a contig read forward, r_c[a, b), or backward and complemented, \overline{r_c[a, b)}. A *chromosome* of a sample is a word s₁ s₂ ⋯ s_k of oriented segments and novel words; a *sample* is a finite multiset of chromosomes.

A segment read forward is entered at (c, a, −) and left at (c, b − 1, +); read backward it is entered at (c, b − 1, +) and left at (c, a, −). Two consecutive segments s_i s_{i+1} are therefore joined by the adjacency {exit(s_i), entry(s_{i+1})}, with the novel bases between them if a novel word separates them.

**Definition 1.3 (junction).** A *junction* of a sample is an adjacency of one of its chromosomes that is not a reference adjacency. A rearrangement is described by its junctions together with the copies of the segments between them.

This is the language of the VCF specification's breakends ("an arbitrary rearrangement event can be summarized as a set of novel adjacencies. Each adjacency ties together 2 breakends", VCFv4.2 §5.4), and of the breakpoint graph of genome rearrangements (Bafna and Pevzner 1996; the graphs of JaBbA and of LINX, Hadi et al. 2020, Cameron et al. 2021). The code is `structural/breakend.hpp` (`Breakend`, `Adjacency`) and, for the simulation, `simulate/structural.hpp` (`Piece`, `Chromosome`, `Junction`).

### 1.2 Copy Number and Adjacencies Balance

Let the copy number c(x) of a reference base x be the number of times it occurs in the sample's chromosomes (counting segments read in either direction). Let m(J) be the multiplicity of an adjacency J among the chromosomes.

**Proposition 1.4 (balance at a breakend).** For every breakend e = (c, p, o) that is not at the end of a chromosome of the sample, the copy number of base p equals the sum of the multiplicities of the adjacencies at e:

c(p) = Σ_{J ∋ e} m(J).

*Proof.* Every occurrence of base p in a chromosome is inside an oriented segment. If p is not an end of that segment, the occurrence lies inside a reference adjacency at e (the one joining p to p + 1 when o = +, to p − 1 when o = −). If it is an end on the side o, the occurrence is the exit or the entry of the segment at e, and it lies in exactly one adjacency of the chromosome, a junction or a novel word's. Conversely, each adjacency at e in a chromosome contains exactly one occurrence of p on its side o. ∎

The balance is what ties depth to junctions: a heterozygous deletion of [a, b) gives the flank base a − 1 two copies, one in its reference adjacency and one in the junction {(a − 1, +), (b, −)}, and the deleted base a one copy, its reference adjacency on the other haplotype. The check `simulate/structural` counts the copies of every base of a rearranged genome by brute force over the lifted bases and compares them with the step function of copies the simulator reports.

### 1.3 Microhomology and the Canonical Form of a Junction

A junction is read by a read that crosses it. The bases the read shows are the derivative chromosome's, and several placements of the junction can give the same derivative.

Index the derivative's bases by k ∈ ℤ, the junction between k = −1 and k = 0, read in the direction that leaves the reference at a and enters it at b. The side left at a = (c, p, +) is λ(k) = r_c[p + 1 + k]; at a = (c, p, −) it is the complement of r_c[p − 1 − k]. The side entered at b = (c′, q, −) is ρ(k) = r_{c′}[q + k]; at b = (c′, q, +) it is the complement of r_{c′}[q − k]. The derivative is λ(k) for k < 0 and ρ(k) for k ≥ 0, and both formulas make sense for every k.

**Definition 1.5 (shift).** The *shift* of the junction by m ∈ ℤ moves it m bases later in the derivative: (p, +) ↦ (p + m, +) and (p, −) ↦ (p − m, −) on the leaving side, (q, −) ↦ (q + m, −) and (q, +) ↦ (q − m, +) on the entering side.

**Proposition 1.6 (the homology interval).** The shifts that leave the derivative unchanged are the integers of an interval [−ℓ, ρ], where ρ is the largest number with λ(k) = ρ(k) for 0 ≤ k < ρ and ℓ the largest with λ(−k) = ρ(−k) for 1 ≤ k ≤ ℓ.

*Proof.* Shifting by m replaces the derivative's bases at offsets between 0 and m − 1 (for m > 0) by λ(k) instead of ρ(k), and leaves the others; so the derivative is unchanged exactly when λ and ρ agree on those offsets. The condition for m implies the one for every m′ between 0 and m, so the set of good shifts is an interval containing 0. ∎

The bases λ(−ℓ), …, λ(ρ − 1) are the junction's *microhomology*; most junctions made by non-homologous end joining or by microhomology-mediated repair carry a few such bases, and they make the base-pair position of a breakpoint ambiguous by exactly this interval. The code chooses one member: the adjacency written with a ≤ b, shifted to the end of the interval that puts a leftmost on the reference (`canonical`). For a deletion and a tandem duplication this is VCF's left alignment, and Tan, Abecasis and Kang's lattice argument (mapping-and-variants.md, §7.1) is Proposition 1.6 in the special case of one contig and one orientation. The checks enumerate, for 900 random junctions with planted repeats (deletions, tandem duplications, translocations), every placement that makes the same derivative, and find the interval of Proposition 1.6 and the same canonical form for all of them.

### 1.4 Signatures and the Events of Two Junctions

**Definition 1.7 (signature).** The signature of an adjacency {a, b} with a ≤ b is: *translocation* when a and b lie on different contigs; otherwise *deletion* when (o_a, o_b) = (+, −), *duplication* when (−, +), and an *inversion's* junction when (+, +) or (−, −).

| Rearrangement of [a, b) | Words | Junctions |
|---|---|---|
| deletion | ⋯ [·, a) [b, ·) ⋯ | (a − 1, +) ~ (b, −) |
| tandem duplication, k extra copies | ⋯ [·, b) [a, b)^k [b, ·) ⋯ | (b − 1, +) ~ (a, −), k times |
| inversion | ⋯ [·, a) \overline{[a, b)} [b, ·) ⋯ | (a − 1, +) ~ (b − 1, +) and (a, −) ~ (b, −) |
| novel insertion at a | ⋯ [·, a) u [a, ·) ⋯ | (a − 1, +) ~ (a, −) through u |
| dispersed copy of [a, b) at p | ⋯ [·, p) [a, b) [p, ·) ⋯ | (p − 1, +) ~ (a, −) and (b − 1, +) ~ (p, −); reversed copy (p − 1, +) ~ (b − 1, +) and (a, −) ~ (p, −) |
| reciprocal translocation at p_A, p_B | A[0, p_A) B[p_B, ·) and B[0, p_B) A[p_A, ·) | (A: p_A − 1, +) ~ (B: p_B, −) and (A: p_A, −) ~ (B: p_B − 1, +) |
| copy-neutral LOH of [a, b) | ⋯ [·, a)_h [a, b)_{h′} [b, ·)_h ⋯ | none |

The events of one junction are named by its signature alone. An event of two junctions is recognised by how its four breakends pair.

**Proposition 1.8 (two junctions).** Let two junctions u, v be written u = {u_a, u_b}, v = {v_a, v_b}, and suppose that {u_a, v_a} and {u_b, v_b} each consist of one + and one − breakend of one contig (call them the two *sides*). Call a side *adjacent* when its − breakend is the base after its + breakend, and *spanning* when its − breakend is at a and its + breakend at b − 1 > a. Then, in the table above:

1. both sides adjacent on one contig, with u and v of signatures (+, +) and (−, −), is an inversion of the segment between the sides;
2. both sides adjacent on two contigs is a balanced reciprocal translocation (or, with (+, +) and (−, −) junctions, its inverted form);
3. one side adjacent at p and one spanning [a, b) is a copy of [a, b) inserted at p, reversed when the junctions are (+, +) and (−, −).

*Proof.* Read the junctions of the table. An inversion's junctions have sides {(a − 1, +), (a, −)} and {(b − 1, +), (b, −)}, both adjacent; a translocation's have sides {(p_A − 1, +), (p_A, −)} on A and {(p_B − 1, +), (p_B, −)} on B; a dispersed copy's have the side {(p − 1, +), (p, −)} at its destination, adjacent, and {(a, −), (b − 1, +)} at its source, spanning. Conversely an adjacent side is a cut of the reference whose two ends are both re-joined, and a spanning side is a segment whose two ends are both joined elsewhere; joining the ends of a cut to the ends of a segment, in the two ways the orientations allow, is a copy of the segment into the cut, and joining the ends of two cuts to each other is an exchange of arms (on one contig, of the segment between them: an inversion). ∎

Adjacency is asked within a tolerance, because the two junctions of a real event rarely meet exactly (target-site duplications, small deletions at the breakpoints) and because breakpoints from pairs alone are imprecise. This is the reasoning that GRIDSS2 and LINX apply to chains of any number of breakends (Cameron et al. 2021); the module applies it to pairs (`structural/call.hpp`, `interpret`).

## 2. What Each Variant Leaves in the Reads

### 2.1 Pairs

A library's fragments have lengths F drawn from a law with survival function S(x) = P(F ≥ x); both ends are read inwards, so that a fragment of an unrearranged region maps forward at its left end and reverse at its right (FR orientation), at a distance F.

**Proposition 2.1 (the template across a junction).** Let a fragment of length F cross a junction {a, b} with x bases on the side left at a and F − x on the side entered at b. Its mates map near a and near b, the mate on a's side mapped forward when a is + and reverse when a is −, and the mate on b's side reverse when b is − and forward when b is +. If a = (c, p, +) and b = (c, q, −) with q > p (a deletion of q − p − 1 bases), the pair maps in FR orientation with template length F + (q − p − 1).

*Proof.* The derivative read in the direction leaving at a is the fragment's strand; the mate read from its start reads the bases λ(−x), …, which are r_c[p + 1 − x, …) read forward when a is +, so that it maps forward, and their reverse complement when a is −. The other mate reads the other strand from the fragment's other end, which is ρ(F − x − 1), …, read back towards the junction: on the reference this is reverse when b is − and forward when b is +. For a deletion, the mates' outer ends are p + 1 − x and q + F − x on contig c. ∎

So a pair is *discordant* by orientation (everted RF pairs for a tandem duplication, FF and RR pairs for an inversion's two junctions, mates on two contigs for a translocation) or by distance (an FR pair whose template exceeds the library's range for a deletion). A deletion is seen in pairs when its length exceeds the spread of F; with DELLY's cut-off, median + 9 MAD (Rausch et al. 2012), and a library of 450 ± 45 bases, from about 280 bases.

A pair with one mate in a novel insertion longer than a read has that mate unmapped: it is *one-end anchored*, its mapped mate pointing to the insertion.

### 2.2 Split Reads and Clips

A read that crosses a junction with enough bases on both sides maps in two parts, the mapper's primary and a supplementary alignment (mapping-and-variants.md, §2: a second primary chain). Its parts, in the read's order, leave the reference at the end of the first and enter it at the start of the second; an overlap of the two on the read is microhomology, a gap between them novel bases. A read with too few bases on one side for a second chain (fewer than about 40 with minimizers of k = 21, w = 11) is soft-clipped instead: an exact breakend whose other side is the clipped bases. A gap of 50 or more bases inside one alignment is a deletion or an insertion seen whole.

**Proposition 2.2 (the expected support of a junction).** In a sample of N pairs of reads of length L over chromosomes of total length G, from fragments of length F, a junction on one chromosome is crossed by mates with at least m bases on each side at the rate 2N(L − 2m + 1)/G, and is spanned by pairs with both mates wholly on one side each at the rate N E[(F − 2L + 1)⁺]/G.

*Proof.* A fragment starts at a uniform position, and each of its two mates starts at one of L − 2m + 1 positions that put m bases on each side of the junction; a pair spans it when the junction lies between the inner ends of its mates, F − 2L + 1 positions of the fragment's start. ∎

At 30× with L = 150, F = 450 ± 45 and m = 20 these are about 22 mates and 15 pairs for each copy of a junction (G counts both haplotypes); a heterozygous junction has about 37 pieces of evidence in expectation, fewer where its reads are not split by the mapper, which is why a threshold of 4 (LUMPY's in its simulations) is far from the noise and far below the signal.

### 2.3 Depth

**Proposition 2.3 (the depth of copies).** If fragments start at reference position x at the rate λ(g(x)) per copy, g(x) the GC content of the fragment that would start there (Benjamini and Speed 2012), the count of fragments starting in a bin B of the sample is Poisson with mean (c/2) Σ_{x ∈ B} λ(g(x)) when every base of B has c copies, relative to a baseline of two.

*Proof.* Each copy of each position contributes fragments independently at its rate (Poisson sampling of fragments), and the c copies of a base add their rates. ∎

The exposure e_B = Σ λ(g(x)) is estimated from the sample itself: λ_g is the number of fragments whose start has GC stratum g over the number of mappable positions of that stratum (`depth.hpp`). Under the model the dispersion D = Σ_B (n_B − e_B)²/e_B of the bins at the baseline has expectation the number of bins, which the checks verify on draws from the model and measure on simulated reads.

## 3. Probabilities of Breakpoints and Genotypes

### 3.1 Breakpoints as Products of Probability Vectors

LUMPY represents a breakpoint by two intervals of positions with a probability vector on each (Layer et al. 2014). A discordant pair with a forward mate starting at s and a reverse mate ending at e′ points to a + breakend at i and a − breakend at j with the fragment length (i − s + 1) + (e′ − j); its probability, given the pair, is proportional to P(F = (i − s + 1) + (e′ − j)).

**Proposition 3.1 (the marginal of a side).** With the other mate of length ℓ′ wholly on its side, the marginal probability of the breakend at i is proportional to Σ_j P(F = i − s + 1 + e′ − j) over j ≤ e′ − ℓ′, which is S(i − s + 1 + ℓ′).

*Proof.* The sum over j of the probabilities of the lengths i − s + 1 + e′ − j for e′ − j ≥ ℓ′ is the probability that F ≥ i − s + 1 + ℓ′. ∎

LUMPY takes P(i) ∝ P(F > i − s), which omits ℓ′; the code adds it (`evidence.hpp`, a *Change*). The evidence of one junction, if its pieces are independent given the junction and the prior on positions is flat, has the likelihood ∏_e P_e(i) on each side, which is LUMPY's default *product* rule; its sum rule is the mixture of the pieces, robust to a piece of another junction. The code combines by the product over the evidence that contains the maximum of the sum, repeatedly (`cluster.hpp`), and reports the maximum and the interval of 95% (VCF's CIPOS and CIEND). A split read is a junction exactly; when a breakpoint has one, its position is the most frequent exact junction of its split reads.

### 3.2 Genotypes

SVTyper's genotype model (Chiang et al. 2015) is binomial: with R reads supporting the reference and A the alternate, L(g) = C(R + A, A) p_g^A (1 − p_g)^R, p_g the fraction of alternate observations genotype g expects.

**Proposition 3.2 (the fraction a junction expects).** Let an adjacency be made on h of a diploid's two haplotypes, k times on each, and let the reference adjacency at its breakends survive on r copies. If every copy of a junction is crossed by reads at one rate, the expected fraction of alternate observations is hk/(hk + r).

*Proof.* The alternate observations come from the hk copies of the junction, the reference observations from the r copies of the reference adjacency, at one rate. ∎

A deletion, an inversion's or a translocation's junction has r = 2 − h, so it expects 1/2 heterozygous and 1 homozygous (SVTyper's 1/2 and 0.9); a tandem duplication of k extra copies keeps its reference adjacencies on both haplotypes, r = 2, and expects k/(k + 2) heterozygous and k/(k + 1) homozygous, where SVTyper takes 0.2 and 1/3. The code uses Proposition 3.2 (`genotype.hpp`) and checks bayes_gt against SVTyper's own code with SVTyper's fractions. The reference support of a junction is the mean over its two breakends of the reads and proper pairs that cross the reference there; adding the two breakends, as SVTyper does, counts the reference allele twice against the junction's one.

## 4. Copy Number as a Step Function

The copy number of a contig is a step function c: [0, |r|) → ℕ, constant between consecutive breakends (Proposition 1.4). Its observations are the bins' counts n_B with means (c_B/2) e_B (Proposition 2.3).

**Definition 4.1 (segmentation).** A segmentation of the bins 0, …, n − 1 is a set of change points 0 < τ₁ < ⋯ < τ_m < n; its penalised cost is Σ_k C(τ_k, τ_{k+1}) + βm, with C(a, b) twice the negative maximised Poisson log-likelihood of bins [a, b) with a common rate, 2(Y − Y log(Y/E)) for their sums Y = Σ n and E = Σ e, up to terms free of the segmentation.

PELT finds the exact minimiser in expected linear time (Killick, Fearnhead and Eckley 2012); tyche's `changepoint.hpp` implements it with the Poisson cost with exposures, and the module takes the BIC penalty β = 2 log n.

**Proposition 4.2 (the integer copy number of a segment).** On a segment with sums Y and E, ℓ(c) = Y log μ_c − μ_c E with μ_c = c/2 is a concave function of μ with its maximum at μ̂ = Y/E; so the best integer copy number c ≥ 1 is ⌊2Y/E⌋ or ⌈2Y/E⌉.

*Proof.* ℓ″(μ) = −Y/μ² ≤ 0. A concave function of one variable is increasing up to its maximum and decreasing after it, so over the integers its maximum is at a neighbour of the real maximum. ∎

c = 0 is compared apart, with μ₀ a small floor (the rate of mismapped fragments). Neighbouring segments of one integer are merged: the segmentation of overdispersed data cuts levels the integers cannot tell apart, and the integer is the quantity of interest. A contig nearly all of whose bins are at one copy number other than two is aneuploid.

## 5. The Alleles Beside the Depth

### 5.1 B-allele Fractions

At a germline heterozygous site in a segment with c copies of which m carry one allele (the minor copies), the fraction of reads showing the B allele is m/c or 1 − m/c, as the alternate lies on the minor or the major haplotype.

| State | c | m | Heterozygous sites | Depth |
|---|---|---|---|---|
| normal | 2 | 1 | BAF 1/2 | 1 |
| copy-neutral LOH | 2 | 0 | none (BAF 0 or 1) | 1 |
| hemizygous deletion | 1 | 0 | none | 1/2 |
| homozygous deletion | 0 | — | none (no reads) | 0 |
| single-copy gain (trisomy) | 3 | 1 | BAF 1/3 and 2/3 | 3/2 |
| copy-neutral LOH with gain | 3 | 0 | none | 3/2 |

Neither depth nor B-allele fractions alone tell the second row from the third: both lose every heterozygous site, and a copy-neutral LOH keeps the depth. Together they do (`allelic.hpp`).

### 5.2 Runs of Homozygosity

BCFtools/RoH is a hidden Markov model over the variant sites of a sample (Narasimhan et al. 2016): two states, AZ (autozygous, one haplotype's alleles twice) and HW, emissions from the genotype likelihoods L(RR), L(RA), L(AA) of each site and an allele frequency f,

P(D | AZ) = (1 − f) L(RR) + f L(AA),   P(D | HW) = (1 − f)² L(RR) + 2f(1 − f) L(RA) + f² L(AA),

and transitions P(AZ | HW) = t_A, P(HW | AZ) = t_H per site (6.7 · 10⁻⁸ and 5 · 10⁻⁹ by default). The code builds it as a chain with emissions given step by step (`hmm/chain.hpp`), decodes it by Viterbi and scores each run by its posterior.

**Proposition 5.1 (the shortest run).** Among heterozygous sites with clean likelihoods (L(RA) = 1 at a heterozygous site, L(AA) = 1 at a homozygous one, the others negligible), the Viterbi path enters AZ for a run of k homozygous-alternate sites exactly when k log(1/f) > −log t_A − log t_H, up to terms of order k max(t_A, t_H).

*Proof.* At a homozygous-alternate site P(D | AZ)/P(D | HW) = f/f² = 1/f; at a heterozygous one the ratio is negligible, so a path in AZ over any heterozygous site is never best. The only competitors are the path in HW throughout and the path that enters AZ before the run and leaves after it; the second gains k log(1/f) in emissions and pays log t_A + log t_H for its two transitions against log(1 − t_A)² for staying, and log(1 − t_H)^{k−1} against log(1 − t_A)^{k−1} inside the run. ∎

With f = 0.4 and the defaults the threshold is 38.9: a run of 39 homozygous sites is called and one of 38 is not, which the checks verify. At a heterozygosity of 10⁻³ per base this is a resolution of a few tens of kilobases.

## 6. Evaluation

Calls are compared with the truth as Truvari compares them (English et al. 2022): one type, starts and ends within 500 bases, size similarity at least 0.7, matched one to one in decreasing score, and precision, recall and genotype concordance by type and by size, the strata of GIAB's benchmark (Zook et al. 2020). Breakpoints are compared after both are made canonical (§1.3), so microhomology costs nothing. Copy number is compared bin by bin, and loss of heterozygosity base by base.

## 7. The State of the Art and What Is Written

| Method | Its idea | Here |
|---|---|---|
| DELLY (Rausch et al. 2012) | discordant pairs by orientation, a library by median and MAD, split-read refinement | the library and its cut-offs (median ± 5 and 9 MAD), the classes of pairs; split reads from the mapper's supplementary parts and clip piles instead of a consensus alignment |
| LUMPY (Layer et al. 2014) | breakpoints as probability vectors, combined by product | written, with the other mate's length in the pair's probability |
| Manta (Chen et al. 2016), GRIDSS (Cameron et al. 2017, 2021) | local assembly of breakend regions; single breakends; complex events | the two-junction events of Proposition 1.8; insertions by a two-sided consensus of clips; assembly planned |
| SVTyper (Chiang et al. 2015) | binomial genotypes from split and spanning reads | written, with the fractions of Proposition 3.2 and the reference averaged over the breakends |
| CNVnator (Abyzov et al. 2011) | read depth in bins, GC-corrected, segmented | written, with Benjamini and Speed's fragment-level GC model and PELT (tyche) instead of mean shift |
| CBS (Olshen et al. 2004), PELT (Killick et al. 2012) | segmentation of copy-number signals | tyche's, used here |
| BCFtools/RoH (Narasimhan et al. 2016) | an HMM of autozygosity over variant sites | written on hmm/chain.hpp, without population allele frequencies |
| Sniffles2 (Smolka et al. 2024) | long-read signatures, repeat-aware clustering | gaps of ≥ 50 bases inside alignments are evidence; long-read calling planned |
| Truvari (English et al. 2022), GIAB (Zook et al. 2020) | matching and stratified benchmarks | written for the simulated truth |

## 8. Where the Code Is

| Section | Code |
|---|---|
| §1 breakends, adjacencies, microhomology, canonical forms | `structural/breakend.hpp`; `simulate/structural.hpp` (the words of pieces, junctions, copy number) |
| §1.4 the events of two junctions | `structural/call.hpp` (`interpret`) |
| §2 pairs, splits, clips, depth | `simulate/pairs.hpp`; `structural/pairs.hpp`, `evidence.hpp`, `depth.hpp` |
| §3 breakpoint probabilities, genotypes | `structural/evidence.hpp`, `cluster.hpp`, `genotype.hpp` |
| §4 segmentation and copy number | `structural/copy_number.hpp`; tyche `time_series/changepoint.hpp` |
| §5 B-allele fractions, runs of homozygosity | `structural/allelic.hpp`; `hmm/chain.hpp`; the caller `variant/caller.hpp` |
| §6 evaluation | `structural/evaluate.hpp` |
| the calls said in VCF, HGVS and ISCN | `structural/vcf.hpp`, `structural/describe.hpp`; `nomenclature/describe.hpp`, `karyotype/describe.hpp` |

## References

- Abyzov, A., Urban, A. E., Snyder, M. and Gerstein, M. (2011), CNVnator: an approach to discover, genotype, and characterize typical and atypical CNVs from family and population genome sequencing, *Genome Research* 21, 974–984.
- Bafna, V. and Pevzner, P. A. (1996), Genome rearrangements and sorting by reversals, *SIAM Journal on Computing* 25, 272–289.
- Benjamini, Y. and Speed, T. P. (2012), Summarizing and correcting the GC content bias in high-throughput sequencing, *Nucleic Acids Research* 40, e72.
- Cameron, D. L., Schröder, J., Penington, J. S. et al. (2017), GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly, *Genome Research* 27, 2050–2060.
- Cameron, D. L., Baber, J., Shale, C. et al. (2021), GRIDSS2: comprehensive characterisation of somatic structural variation using single breakend variants and structural variant phasing, *Genome Biology* 22, 202.
- Chen, X., Schulz-Trieglaff, O., Shaw, R. et al. (2016), Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications, *Bioinformatics* 32, 1220–1222.
- Chiang, C., Layer, R. M., Faust, G. G. et al. (2015), SpeedSeq: ultra-fast personal genome analysis and interpretation, *Nature Methods* 12, 966–968.
- English, A. C., Menon, V. K., Gibbs, R. A., Metcalf, G. A. and Sedlazeck, F. J. (2022), Truvari: refined structural variant comparison preserves allelic diversity, *Genome Biology* 23, 271.
- Hadi, K., Yao, X., Behr, J. M. et al. (2020), Distinct classes of complex structural variation uncovered across thousands of cancer genome graphs, *Cell* 183, 197–210.
- Killick, R., Fearnhead, P. and Eckley, I. A. (2012), Optimal detection of changepoints with a linear computational cost, *Journal of the American Statistical Association* 107, 1590–1598.
- Layer, R. M., Chiang, C., Quinlan, A. R. and Hall, I. M. (2014), LUMPY: a probabilistic framework for structural variant discovery, *Genome Biology* 15, R84.
- Narasimhan, V., Danecek, P., Scally, A., Xue, Y., Tyler-Smith, C. and Durbin, R. (2016), BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data, *Bioinformatics* 32, 1749–1751.
- Olshen, A. B., Venkatraman, E. S., Lucito, R. and Wigler, M. (2004), Circular binary segmentation for the analysis of array-based DNA copy number data, *Biostatistics* 5, 557–572.
- Rausch, T., Zichner, T., Schlattl, A., Stütz, A. M., Benes, V. and Korbel, J. O. (2012), DELLY: structural variant discovery by integrated paired-end and split-read analysis, *Bioinformatics* 28, i333–i339.
- Smolka, M., Paulin, L. F., Grochowski, C. M. et al. (2024), Detection of mosaic and population-level structural variants with Sniffles2, *Nature Biotechnology* 42, 1571–1580.
- Tan, A., Abecasis, G. R. and Kang, H. M. (2015), Unified representation of genetic variants, *Bioinformatics* 31, 2202–2204.
- The SAM/BAM and VCF Format Specification Working Groups, *The Variant Call Format Specification, VCFv4.2*, samtools/hts-specs.
- Yoon, S., Xuan, Z., Makarov, V., Ye, K. and Sebat, J. (2009), Sensitive and accurate detection of copy number variants using read depth of coverage, *Genome Research* 19, 1586–1592.
- Zook, J. M., Hansen, N. F., Olson, N. D. et al. (2020), A robust benchmark for detection of germline large deletions and insertions, *Nature Biotechnology* 38, 1347–1355.
