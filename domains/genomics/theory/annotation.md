# The Annotation of Genomes

*Status: draft. Serves learning and use.* What can be read off a genome's sequence about its parts: features and their locations on a circle; motifs as distributions, with their information and the exact law of their scores; the discovery of a motif as inference of a hidden alignment, by Gibbs sampling and by EM; the statistics of codons that find genes; and exact repeats as nodes of the suffix tree. The code is the module [`annotation`](../annotation/): [`features.hpp`](../annotation/include/genomics/annotation/features.hpp), [`motifs.hpp`](../annotation/include/genomics/annotation/motifs.hpp), [`genes.hpp`](../annotation/include/genomics/annotation/genes.hpp), [`repeats.hpp`](../annotation/include/genomics/annotation/repeats.hpp); hidden Markov models are [hmm.md](hmm.md), suffix arrays [indexes.md](indexes.md).

Throughout, a genome is a word x = x₀ ⋯ x_{n−1} over {A, C, G, T} (codes 0–3, N = 4 for an unknown base), coordinates are 0-based and half-open, and rc(x) is the reverse complement, rc(x)_i = 3 − x_{n−1−i}.

## 1. Features and Their Locations

**Definition 1.1 (location).** A location is a word of the context-free language

location ::= span | complement(location) | join(location, …) | order(location, …),

whose spans are intervals of the genome (INSDC Feature Table, §3.4). Its *value* is a list of parts, each an interval with a strand: complement reverses the list and the strands of what it encloses; join and order concatenate. The value of a CDS is read in order: its bases are the concatenation of its parts, a reversed part read as its reverse complement (`sequence::extract`).

**Definition 1.2 (crossing the origin).** On a circular genome of length n, list the parts of a feature on one strand in the forward strand's order (as read for the forward strand, reversed for the reverse strand). The feature *crosses the origin* when some part begins before the part listed before it. Its extent is then the arc [start, n) ∪ [0, stop), start the beginning of the first part and stop the end of the last.

*Remark.* A feature cannot say "from 2901 to 203" with one span, since a span is an interval of the line; join(2901..3000,1..203) is how GenBank writes the arc; phiX174's genes A, A* and B are written so (A is join(3981..5386,1..136) in NC_001422.1), the six features the check finds across its origin. On the circle the arc is a single interval of ℤ/nℤ, which is what the circular view draws; the linear view draws its two segments (`TrackFeature::segments`). A GFF3 file writes the same arc as one interval with end − start beyond the landmark's length (the specification's bacteriophage f1 example, gene II from 6006 to 7238 on a genome of 6407), and `annotation::segments` cuts it at n.

**Definition 1.3 (the part-of graph of GFF3).** The features of a GFF3 file are the classes of lines sharing an ID (a line without one is its own class). Parent=p on a line of feature f gives an edge f → p. The resulting graph must be acyclic (GFF3 1.26, "Parent (part_of) Relationships"); it is not in general a forest, since an exon may be part of several transcripts (exon00004 of the canonical gene is part of all three of EDEN's mRNAs).

`read_gff3` builds the graph, reports every Parent that names no ID, and finds cycles by a depth-first search in three colours: an edge into a vertex still on the stack closes a cycle. O(V + E).

**Lanes.** For drawing, the features of a kind are placed in lanes so that no two in a lane overlap: a greedy colouring of the interval graph in order of start, which on intervals of a line uses the least possible number of lanes (an interval graph is perfect, and the greedy colouring in order of left ends is optimal; Golumbic 1980, ch. 8). Arcs across the origin take a lane each before the rest, which keeps the colouring valid on the circle though not always optimal there.

## 2. Motifs as Distributions

**Definition 2.1 (PWM).** A motif of width W is a product distribution p = (p_i)_{i<W} on words of length W, p_i a distribution on the bases. Against a background q, the *weight* of base b at position i is w_{i,b} = log₂(p_{i,b}/q_b), and the *score* of a word is S(u) = Σ_i w_{i,u_i} = log₂ P(u | motif)/P(u | background): a log-likelihood ratio, in bits (Stormo 2000).

**Proposition 2.2 (information is the mean score).** The mean score of words drawn from the motif is the relative entropy

E_p[S] = Σ_i D(p_i ‖ q) = Σ_i Σ_b p_{i,b} log₂(p_{i,b}/q_b) ≥ 0,

with equality only when p = q. For a uniform background, D(p_i ‖ q) = 2 − H(p_i): a column of one base carries 2 bits, two bases equally 1 bit.

*Proof.* Linearity of expectation over positions, and Gibbs' inequality. ∎

**Estimation.** From n aligned sites with counts n_{i,b}, p_{i,b} = (n_{i,b} + α q_b)/(n + α), the posterior mean under a Dirichlet prior α q: the pseudocounts keep a base never seen at a position from scoring −∞.

**Definition 2.3 (Schneider's R_sequence).** With observed frequencies f_{i,b} and H_i = −Σ f log₂ f, R_sequence(i) = H_g − (H_i + e(n)), where H_g is the background's entropy and e(n) = H_g − E[H_n] is the bias of the entropy of n bases drawn from the background (Schneider, Stormo, Gold and Ehrenfeucht 1986). A sequence logo draws at position i the letters with heights f_{i,b} R_sequence(i) (Schneider and Stephens 1990).

**Proposition 2.4.** For the uniform background, e(1) = 2 and e(2) = 5/4; and e(n) = 3/(2n ln 2) + O(n⁻²) (Miller 1955).

*Proof.* One base has entropy 0, so E[H₁] = 0. Two bases are equal with probability 1/4 (entropy 0) and differ with probability 3/4 (entropy 1 bit), so E[H₂] = 3/4. The asymptotic form is the second-order expansion of E[H_n] for a multinomial (Miller 1955; Basharin 1959). ∎

`small_sample_correction` sums E[H_n] exactly over the C(n + 3, 3) compositions of n for n ≤ 150, and the check holds its value at 150 to Miller's formula within 0.6%.

## 3. The Exact Law of Scores

**Definition 3.1 (the lattice).** Fix a granularity ε > 0, and round each weight to an integer W_{i,b} = round(w_{i,b}/ε). The rounded score of u is S′(u) = ε Σ_i W_{i,u_i}.

**Proposition 3.2 (Staden 1989).** Under the background, the generating function of the rounded score is the product

F(z) = Π_{i<W} (Σ_b q_b z^{W_{i,b}}),

so P(S′ = kε) is the coefficient of z^k; the dynamic programme P_i(s) = Σ_b q_b P_{i−1}(s − W_{i,b}) multiplies the factors in turn, in O(W R) for R lattice points.

*Proof.* The positions of a background word are independent, so the generating function of the sum is the product of the generating functions of the terms. ∎

**Proposition 3.3 (the rounding bound; after Touzet and Varré 2007).** For every word u, |S(u) − S′(u)| ≤ Wε/2, and hence for every t,

P(S′ ≥ t + Wε/2) ≤ P(S ≥ t) ≤ P(S′ ≥ t − Wε/2).

*Proof.* Each weight moves by at most ε/2 in rounding. If S(u) ≥ t then S′(u) ≥ t − Wε/2, which gives the upper bound; if S′(u) ≥ t + Wε/2 then S(u) ≥ t, the lower. ∎

The p-value of a word's own rounded score is exact for the rounded matrix, and converges to the true one as ε → 0. The check enumerates all 4⁶ and 4⁸ words of random PWMs: the lattice law equals the enumeration point by point (to 10⁻¹⁸), and every real tail lies within the bounds of 3.3.

## 4. Discovering a Motif

Given n sequences, each containing (or not) an occurrence of an unknown motif of width W, the unknowns are the motif p and the *alignment* a = (a_1, …, a_n) of site positions. Two strategies infer them.

### 4.1 The Gibbs Sampler

**Proposition 4.1 (the predictive update).** Under a Dirichlet prior B·q on each column, with the other n − 1 sites fixed and counts c_{i,b} among them, the probability that sequence z's site is at a, integrating out p, is proportional to

Π_i q̂_{i, x_{a+i}} / q_{x_{a+i}}, q̂_{i,b} = (c_{i,b} + B q_b)/(n − 1 + B).

*Proof.* The posterior of p_i given the other sites is Dirichlet(c_i + B q), whose predictive probability of a base b is q̂_{i,b}; the positions are independent; and the background of the rest of sequence z, Π q over the bases outside the window, is the same for every a up to the factor Π q over the window, which divides. ∎

So the Gibbs sampler of Lawrence et al. (1993) is a *collapsed* Gibbs sampler on the alignment alone: it removes one site, draws it again from Proposition 4.1, and repeats. Its stationary distribution is the posterior of the alignment. The phase shifts move every site by δ at once, δ drawn in proportion to the likelihood ratio of the shifted alignment, a Metropolis-within-Gibbs move between alignments that the one-site moves connect only through improbable states (a motif found two bases off). The objective tracked, and returned at its maximum, is

F = Σ_i Σ_b c_{i,b} log₂(q̂_{i,b}/q_b),

the log-likelihood ratio of the alignment (Lawrence et al. 1993). Restarts are independent streams of `calx::random::Stream`, run in parallel and compared in order, so the result does not depend on the threads.

### 4.2 Expectation Maximisation

**Definition 4.2 (OOPS and ZOOPS; Bailey and Elkan 1994).** In the *one occurrence per sequence* model, sequence j has its site at a uniformly random window of its m_j; in the *zero or one* model, it has one with probability γ, then uniformly placed. The rest is background. With R_jw = Π_i p_{i,x}/q_x the odds of window w,

log P(x_j | OOPS) − log P(x_j | background) = log((1/m_j) Σ_w R_jw), log P(x_j | ZOOPS) − log P(x_j | background) = log((1 − γ) + (γ/m_j) Σ_w R_jw).

EM (hmm.md, Theorem 4.3) with the site indicators as hidden data: the E-step gives the posterior Z_jw of each window, the M-step re-estimates p from the expected counts Σ Z_jw [x_{w+i} = b] and γ as the mean of Σ_w Z_jw. With pseudocounts, the log posterior never decreases (hmm.md, Proposition 4.4), which the check verifies.

*Remark (starting points).* EM climbs to a local maximum. MEME starts from each distinct word of the data, made a matrix with probability 0.5 for the word's base and 1/6 for the others, runs one iteration from each, and continues from the best (Bailey and Elkan 1995). Every start is independent, so they run in parallel.

**Proposition 4.3 (a background sequence is not evidence against a site).** If x_j is drawn from the background, E[Σ_w R_jw] = m_j.

*Proof.* E_q[R] = Σ_u q(u) p(u)/q(u) = 1 for each window. ∎

*Remark (why ZOOPS overcalls).* The mean of Σ_w R_jw is carried by rare high-scoring windows, so a typical background sequence has Σ R ≪ m_j and is evidence against a site; but the best window of each such sequence is a partial match that the maximum-likelihood motif can absorb. In the check, the ZOOPS fit takes 7 weak sites in sequences with none, at a log-likelihood ratio of 200.1 bits against the planted truth's 198.2: the excess calls are the likelihood's, not EM's failure. MEME chooses the number of sites by the significance of the motif instead.

## 5. The Statistics of Codons

**Definition 5.1 (periodic chains).** A three-periodic Markov chain of order k assigns each base of a coding sequence the probability P_c(x_t | x_{t−k} ⋯ x_{t−1}, t mod 3), the phase being the position in the codon; a homogeneous chain P_n(x_t | x_{t−k} ⋯ x_{t−1}) models non-coding sequence. The *log-likelihood ratio* of a reading frame from s to the stop at e is

Λ(s, e) = Σ_{t=s}^{e−4} log₂ P_c(x_t | ⋯, (t − s) mod 3) − log₂ P_n(x_t | ⋯).

**Proposition 5.2.** If the bases of [s, e) are drawn from the coding chain, E[Λ] is a sum of conditional relative entropies D(P_c(· | context, phase) ‖ P_n(· | context)) ≥ 0, so Λ grows linearly in the length of a true gene; on non-coding sequence it decreases linearly.

*Proof.* As Proposition 2.2, position by position, conditioning on the context. ∎

This is why the long ORFs carry the signal: a random frame of non-coding sequence meets one of the three stops of 64 codons every 21 codons on average, so an ORF of 300 codons is almost always a gene, and the first coding model can be learnt from them without annotation (Glimmer's long ORFs; GeneMarkS's iterations). Re-training on the predictions and re-predicting is *classification EM*: a hard E-step (each candidate a gene or not) followed by the maximum-likelihood chains of the classes; its predictions settle within a few rounds (on E. coli the number of genes is constant from the fourth).

**Estimation with back-off.** P_j(b | c) = (n(c, b) + α P_{j−1}(b | c′))/(n(c) + α), c′ the context c without its most distant base: the posterior mean under a Dirichlet prior centred on the chain one order lower (MacKay and Peto 1995). A context seen often is its own estimate; one seen rarely falls back toward shorter contexts, which is what Glimmer's interpolated Markov models do by other means (Salzberg et al. 1998).

**Starts.** Among the in-frame starts of an ORF, Λ(s, e) + log₂ P(start codon) is maximised. By Proposition 5.2 lengthening into non-coding sequence lowers Λ, so the maximum is near the true start; without a ribosome-binding-site model the finder gets 74% of E. coli's starts right, against 95% of its stops (the check).

**Proposition 5.3 (overlaps).** Let every interval be longer than 2h. Two intervals [b, e) and [b′, e′) overlap by at most 2h if and only if their shrunk forms [b + h, e − h) and [b′ + h, e′ − h) are disjoint.

*Proof.* With b ≤ b′, the overlap is min(e, e′) − b′. The shrunk intervals are disjoint iff min(e, e′) − h ≤ b′ + h, that is, iff the overlap is at most 2h. ∎

So the best set of genes with overlaps of at most `max_overlap` bases is a maximum-weight set of disjoint intervals, found exactly by weighted interval scheduling: sorted by right end, OPT(j) = max(OPT(j − 1), w_j + OPT(p(j))), p(j) the last interval ending before j begins (Kleinberg and Tardos 2006, §6.1). O(G log G).

## 6. Exact Repeats

**Definition 6.1 (Gusfield 1997, §7.12).** A *maximal pair* of x is a pair of occurrences i ≠ j of a factor α such that x_{i−1} ≠ x_{j−1} (or one of them is at the start) and x_{i+|α|} ≠ x_{j+|α|} (or one is at the end). α is a *maximal repeat* when it has a maximal pair, and a *supermaximal repeat* when, besides, it is a factor of no other maximal repeat.

**Theorem 6.2 (Gusfield, Theorem 7.12.1).** α is a maximal repeat if and only if it is the label of an internal node v of the suffix tree of x that is *left-diverse*: two leaves below v have different left letters (a suffix at the start counting as a letter of its own).

*Proof.* An internal node is right-maximal: two suffixes below it continue differently after α. If v is left-diverse, take leaves with left letters a ≠ a′; if they lie below different children of v, they form a maximal pair; otherwise a leaf below another child has a left letter differing from one of a, a′, and pairs with it. Conversely the two occurrences of a maximal pair continue differently, so α is a node, and their left letters differ, so it is left-diverse. ∎

**Proposition 6.3 (on the suffix array).** The internal nodes of the suffix tree are the lcp-intervals of the suffix array: intervals [i, j] with LCP[k] ≥ ℓ on (i, j], = ℓ somewhere there, and LCP[i], LCP[j + 1] < ℓ (Abouelhoda, Kurtz and Ohlebusch 2004). The node of [i, j] is left-diverse iff the Burrows–Wheeler letters x[SA[k] − 1], i ≤ k ≤ j, are not all equal; its children are all leaves iff LCP = ℓ on all of (i, j]; and it is a supermaximal repeat iff it is left-diverse, its children are all leaves and their left letters are pairwise distinct (Gusfield, §7.12).

`maximal_repeats` enumerates the lcp-intervals bottom-up with a stack in one pass over LCP, marking each interval that receives a child interval, and tests left-diversity in O(1) with the runs of the BWT precomputed: O(n) for the traversal after morphe's SA-IS and Kasai's LCP.

**Both strands.** On the word T = x #₀ rc(x) #₁, with #₀ and #₁ distinct letters occurring once, the map that sends an occurrence at p in x to the occurrence of the reverse complement at the mirrored place in rc(x) is an anti-automorphism of the factors of T: it carries α to rc(α), left letters to complemented right letters and back. So α is a (supermaximal) repeat of T iff rc(α) is, with the same places on the double strand; `maximal_repeats` reports the smaller of the two, its occurrences given on the forward strand with the strand they read on. A palindrome α = rc(α) occurring once in x is found twice in T, at one place on each strand: an inverted repeat of itself, the shape of a hairpin. The check verifies all of this against brute force over every factor and pair of occurrences.

## 7. Where the Code Is

| Concept | Code |
|---|---|
| Definitions 1.1–1.2, lanes | `feature_tracks`, `TrackFeature`, `Tracks` |
| Definition 1.3 | `read_gff3`, `Gff`, `GffEntity`, `segments` |
| Definitions 2.1, 2.3, Proposition 2.4 | `Pwm`, `pwm_from_counts`, `pwm_from_sites`, `information_content`, `small_sample_correction`, `logo` |
| Propositions 3.2–3.3 | `score_distribution`, `ScoreDistribution`, `scan` |
| Proposition 4.1 | `gibbs_motif` |
| Definition 4.2 | `meme_motif`, `motif_log_likelihood` |
| Definition 5.1, back-off | `PeriodicChain` |
| Proposition 5.3, the finder | `find_genes` |
| Theorem 6.2, Proposition 6.3 | `maximal_repeats`, `supermaximal_repeats`, `longest_repeat` |

## References

- Abouelhoda, M. I., Kurtz, S. and Ohlebusch, E. (2004), *Journal of Discrete Algorithms* 2, 53–86.
- Bailey, T. L. and Elkan, C. (1994), *Proceedings of ISMB* 2, 28–36; (1995), *Machine Learning* 21, 51–80.
- Basharin, G. P. (1959), On a statistical estimate for the entropy of a sequence of independent random variables, *Theory of Probability and Its Applications* 4, 333–336.
- Besemer, J., Lomsadze, A. and Borodovsky, M. (2001), *Nucleic Acids Research* 29, 2607–2618; Borodovsky, M. and McIninch, J. (1993), *Computers & Chemistry* 17, 123–133.
- Delcher, A. L., Harmon, D., Kasif, S., White, O. and Salzberg, S. L. (1999), *Nucleic Acids Research* 27, 4636–4641; Salzberg, S. L., Delcher, A. L., Kasif, S. and White, O. (1998), *Nucleic Acids Research* 26, 544–548.
- Golumbic, M. C. (1980), *Algorithmic Graph Theory and Perfect Graphs*, Academic Press.
- Gusfield, D. (1997), *Algorithms on Strings, Trees, and Sequences*, Cambridge University Press.
- Kleinberg, J. and Tardos, É. (2006), *Algorithm Design*, Addison-Wesley.
- Lawrence, C. E. et al. (1993), *Science* 262, 208–214.
- MacKay, D. J. C. and Peto, L. C. B. (1995), A hierarchical Dirichlet language model, *Natural Language Engineering* 1, 289–307.
- Miller, G. A. (1955), Note on the bias of information estimates, in *Information Theory in Psychology*, ed. H. Quastler, Free Press, 95–100.
- Schneider, T. D., Stormo, G. D., Gold, L. and Ehrenfeucht, A. (1986), *Journal of Molecular Biology* 188, 415–431; Schneider, T. D. and Stephens, R. M. (1990), *Nucleic Acids Research* 18, 6097–6100.
- Staden, R. (1989), *Computer Applications in the Biosciences* 5, 89–96; Touzet, H. and Varré, J.-S. (2007), *Algorithms for Molecular Biology* 2, 15.
- Stormo, G. D. (2000), *Bioinformatics* 16, 16–23.
- The Sequence Ontology, *Generic Feature Format Version 3*, version 1.26 (2020); the INSDC Feature Table Definition.
