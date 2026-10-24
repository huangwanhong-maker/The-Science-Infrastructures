# Mapping Reads and Calling Variants

*Status: draft. Serves learning and use.* How a read finds its place in a genome, how sure one can be of that place, and what the reads piled at a place say about the genome of the individual they came from. Seeds and their sampling by minimizers; chains as paths of a graph of anchors, scored by a dynamic programme; mapping quality as a posterior; the base-level alignment of a chain; pileups, genotype likelihoods and their posterior; the normal form of an indel as the canonical member of a class of equivalent representations; and the evaluation of calls. The code is the modules [`mapping`](../mapping/) and [`variant`](../variant/) of the genomics domain; the alignments they compute are those of [alignment.md](alignment.md), and their seeds come from the minimizer index of [indexes.md](indexes.md), §6.

## 1. Seeds

Let a genome be a finite family of words c₁, …, c_r (its contigs) over {A, C, G, T, N}, and a read x a word of length L over the same letters. Write ū for the reverse complement of a word u: the word read backwards with A ↔ T and C ↔ G exchanged. The genome is double-stranded, so x may come from a contig or from its reverse complement.

**Definition 1.1 (anchor).** For k ≤ 32, an *anchor* of x on the genome is a triple (i, j, s) of a read position, a contig position and a strand, with x[i, i + k) = c[j, j + k) when s = + and x[i, i + k) = c[j, j + k)‾ when s = −.

An anchor is the cheapest evidence of a place: k letters that agree exactly. A read has L − k + 1 k-mers, and looking up all of them is wasteful when consecutive ones overlap in k − 1 letters. The mapper looks up only its (w, k)-minimizers: of every w consecutive k-mers, the one of least hash (indexes.md, §6).

**Proposition 1.2 (the window guarantee).** If x and a contig share a word u of length w + k − 1 without N, they share an anchor inside u.

*Proof.* The w k-mers of u form a window in x and in the contig. The minimizer of a window is a function of its k-mers alone, so the two windows select the same k-mer at the same offset of u, and it is a minimizer of both. ∎

*Remark.* This is why the presets differ. A read with an error rate e per base has an error-free run of w + k − 1 bases at a given place with probability (1 − e)^{w+k−1}: 0.73 for short reads (e = 0.01, k = 21, w = 11), but only 0.14 for long reads (e = 0.08, k = 15, w = 10). Long reads therefore use shorter seeds; they make up for the low rate by their thousands of windows.

**Proposition 1.3 (the strand mapped).** With L the read's length, x[i, i + k) = c[j, j + k)‾ if and only if x̄[L − k − i, L − i) = c[j, j + k).

*Proof.* x̄[t] is the complement of x[L − 1 − t]. So x̄[L − k − i + m] is the complement of x[i + k − 1 − m], for m = 0, …, k − 1. This is letter m of x[i, i + k)‾. ∎

So every anchor on the reverse strand is a forward anchor of the read's reverse complement. A chain and its alignment are then always forward, of x or of x̄ (`chain.hpp`, `seeds_of`).

*Remark (seeds by SMEMs).* BWA-MEM seeds a read differently (Li 2013). A *maximal exact match* of x is a factor x[b, e) that occurs in the genome and cannot be extended to either side; a *super-maximal* one (SMEM) is contained in no other. There are at most L of them, and the FMD-index finds them all (indexes.md, §5.2). Each occurrence of an SMEM of at least 19 bases is a seed of the match's length. Such seeds are fewer and longer than minimizers: an error-free read is one seed. They are chained by the same programme (§2), which is why its coordinates are the seeds' ends (`MapperOptions::short_reads_smem`).

## 2. Chains and Their Dynamic Programme

Fix a strand and a contig, and write a seed as a = (x_a, y_a, w_a): its end x_a on the target and its end y_a on the strand of the read (one past its last base), and its span w_a (k for a minimizer). The letter x names target positions from here on, as in Li 2018, where an anchor is kept by its last base.

**Definition 2.1 (precedence, chain).** An anchor a *precedes* b, written a ≺ b, when

- 0 < y_b − y_a ≤ G and 0 < x_b − x_a ≤ G, and
- |l| ≤ B, where l = (y_b − y_a) − (x_b − x_a).

Here G is the maximum gap and B the bandwidth. A *chain* is a sequence a₁ ≺ a₂ ≺ ⋯ ≺ a_m. Its score is

S(a₁, …, a_m) = w_{a₁} + Σ_{t=2}^{m} [α(a_{t−1}, a_t) − β(a_{t−1}, a_t)],

with α(a, b) = min{y_b − y_a, x_b − x_a, w_b} and β(a, b) = γ(l), where γ(0) = 0 and

γ(l) = 0.01 · w̄ · |l| + 0.5 · log₂|l|   (l ≠ 0)

(Li 2018, §2.1.2, eqs. 1–2; w̄ is the mean span).

The two terms say what a step gains and costs:

- α(a, b) is the number of bases b adds to the chain: all of its w_b when the two anchors do not overlap, fewer when they do.
- β is the cost of the gap between them. l is the difference of the two distances: it is 0 when they lie on one diagonal, and an indel of |l| bases otherwise.

*Remark (ends, not starts).* For seeds of one span, the differences of their ends are those of their starts, and either may be used. For seeds of different spans only the ends make α the bases of b beyond a. Take a seed of 39 bases at 0 followed by one of 72 bases at 40: by starts, α = min{40, 72} = 40, as if the second seed overlapped the first in 32 bases; by ends, α = min{73, 72} = 72. The first version of the code used starts, and the same three SMEMs chained to different scores on the two strands of an exact repeat, one 147 and the other 114. A read in an exact repeat then got a high mapping quality, and the check of such reads caught it.

*Rationale.* The gap cost is concave: a linear part and a logarithmic one. Indel lengths have a heavy tail, so a long gap is cheaper per base than a short one. A chain can then cross a 300-base deletion that a linear cost would have cut.

**Proposition 2.2 (the dynamic programme).** Sort the anchors by (x, y), and let

f(i) = max{ w_i, max_{j ≺ i} [f(j) + α(j, i) − β(j, i)] }.

Then f(i) is the largest score of a chain ending at anchor i.

*Proof.* The relation ≺ is acyclic. If j ≺ i, then x_j < x_i, so j comes before i in the sorted order, which is therefore a topological order. A chain ending at i is either i alone (score w_i), or a chain ending at some j ≺ i extended by i (score S + α(j, i) − β(j, i)). Induction on the order. ∎

This is the sum over paths of alignment.md, §2, again. The (max, +) semiring is the same; the graph is now the DAG of anchors, with a weight on each vertex and on each edge. The code (`chain_dp`) restricts the predecessors of i to its h nearest ones (h = 50). This is a sub-DAG: it computes at most f, in O(N h) instead of O(N²) for N anchors. The check holds f equal to the best chain over every subset of up to twelve anchors.

**Proposition 2.3 (the integer gap cost).** Let

γ_ℤ(l) = ⌊(2 w̄ |l| + 100 ⌊log₂|l|⌋) / 200⌋.

Then γ(l) − 1.5 < γ_ℤ(l) ≤ γ(l) for l ≠ 0.

*Proof.* The quotient before the floor is 0.01 w̄ |l| + 0.5 ⌊log₂|l|⌋, which lies in (γ(l) − 0.5, γ(l)]. The floor takes off less than 1 more. ∎

The code uses γ_ℤ, so that scores are integers and are the same on every machine. The largest gap to γ measured over |l| ≤ 5000 is 1.485.

**Proposition 2.4 (backtracking).** Let p(t) be the predecessor attaining f(t). Follow p from an anchor i until an anchor k that is already in a chain, or none. The anchors passed form a chain whose score, counted from k, is f(i) − f(k) (with f(none) = 0).

*Proof.* f(t) = f(p(t)) + α(p(t), t) − β(p(t), t) along the walk, and the sum telescopes. ∎

Visiting the anchors in decreasing f and stopping at used anchors partitions them into disjoint chains (`chain`). The best chain comes first, and each later chain is the best use of what is left. A chain that overlaps a better one on the read is a rival placement of the same bases. It is *secondary* to it when the overlap is at least half the shorter of the two (`select_chains`). A chain overlapping no better one is another part of the read, such as the other side of a chimera: a *supplementary* primary.

*Remark (fragments).* The backtracking leaves behind the anchors that a better chain did not take, and they can form a short chain of the same placement: a parallel diagonal a few bases away, across an indel. Such a chain is not another place the read could come from, and counting it as a rival would lower the mapping quality of a read that is unique. The code drops a chain on the same contig and strand, overlapping a better one on the read, whose diagonal is within the bandwidth. This is a rule of the code that the paper does not state.

## 3. Mapping Quality as a Posterior

**Definition 3.1 (mapping quality).** The mapping quality of a placement is Q = −10 log₁₀ P(the placement is wrong | x), on the Phred scale (Li, Ruan and Durbin 2008).

**Proposition 3.2 (the posterior of placements).** Let a read have candidate placements u₁, …, u_n, equally likely a priori, with P(x | u_t) ∝ e^{λ S_t} for alignment scores S₁ > S₂ ≥ … ≥ S_n. Then

P(u₁ wrong | x) = Σ_{t ≥ 2} e^{λ(S_t − S₁)} / Σ_t e^{λ(S_t − S₁)},

and when one rival dominates the others, Q ≈ 10 λ log₁₀(e) · (S₁ − S₂).

*Proof.* Bayes' rule with a uniform prior. Dividing by e^{λS₁} gives the first formula. When the sums are dominated by t = 2, the probability is about e^{−λ(S₁ − S₂)}. ∎

*Remark.* The hypothesis says that a scoring scheme is a log-likelihood ratio at the scale λ. This is the content of Karlin and Altschul's theory (alignment.md, §7). Two equal placements give P(wrong) ≥ 1/2, so Q ≤ 3, and mappers report 0.

Li's mapping quality (2018, §2.1.4) is a heuristic of the same shape, read from the chain scores f₁ ≥ f₂ of the primary chain and its best rival:

mapQ = 40 · (1 − f₂/f₁) · min{1, m/10} · ln f₁,

for a primary of m anchors. It vanishes when f₂ = f₁ and grows with the margin. The factor ln f₁ makes a long, unique chain more certain than a short one, and the factor min{1, m/10} distrusts chains of few anchors. Once the chains are aligned, the alignment scores s₁ ≥ s₂ are finer evidence than the chain scores. The mapper takes the best alignment of a group of rival chains as the primary, and replaces f₂/f₁ by (f₂/f₁)(s₂/s₁), the product minimap2's implementation uses. Seeded by SMEMs, a chain has one to a few seeds, and the factor min{1, m/10} would cap the quality of every read. That path uses BWA-MEM's approximation instead (Li 2013), which is Proposition 3.2 with λ fixed by the match score a:

mapQ = 6.02 · (s₁ − s₂)/a · ι².

Here 6.02 = 10 log₁₀ 4, a matching base being four times likelier on the right placement than on a random one; s₂ is at least 19a (an SMEM's worth); and ι is the identity of the alignment.

Whether the result is calibrated is an empirical question, and the checks answer it with the curve of Li 2018, Fig. 1: the error rate of the mapped reads at each threshold of quality. On Escherichia coli, the short reads of MAPQ ≥ 20 are wrong 5 times in 10⁵ seeded by minimizers, and 2.5 times seeded by SMEMs.

## 4. Base-Level Alignment of a Chain

**Proposition 4.1 (the alignment of a chain).** Let (i₀, j₀), …, (i_m, j_m) be vertices of the edit graph G(x, y) with i_t ≤ i_{t+1} and j_t ≤ j_{t+1}, and let π_t be a path of the sub-graph G(x[i_t, i_{t+1}), y[j_t, j_{t+1})). Then the concatenation π₁ ⋯ π_m is a path of G(x, y) from (i₀, j₀) to (i_m, j_m). Its score is at least the sum of the scores of the pieces, with equality unless two gaps of one kind meet across a boundary.

*Proof.* The sub-graph on [i_t, i_{t+1}] × [j_t, j_{t+1}] is the induced sub-graph of G(x, y) on those vertices, so its paths are paths of G, and consecutive pieces share their end vertex. Each pair and each gap keeps its score. Where a run of gaps of one kind continues across a boundary, the concatenation opens it once instead of twice. ∎

The mapper uses this with the anchors of a chain. Each anchor is a run of k pairs. Between consecutive anchors there is a global alignment of the gap, by the wavefront algorithm or by bands certified optimal (alignment.md, §4.3–4.5). The ends are extended outwards (§4.2 below). The path is optimal among the paths through the anchors' vertices, not over the whole edit graph. The anchors' exact matches are trusted, and `rescore` is the arbiter of the final score, which the checks hold to SAM's AS.

### 4.2 The Ends and the Clipping Penalty

Beyond its first and last anchors, a read is extended by the X-drop rule with minimap2's Z-drop. This is a pruned frontier of Gotoh's recurrences that follows the alignment and stops when the score falls too far below its best (alignment.md, §4.4). What the extension does not reach is soft-clipped. The best cell of an extension is a *local* end: every vertex is a sink.

For short reads a second end is considered: the *glocal* end, whose sinks are the read's last row. Its score is compared with the local one. The glocal end is kept when it scores more than the local end less a clipping penalty p = 5 (BWA-MEM's rule; Li 2013). So a mismatch two bases from the end, which costs −8 against the +4 of the two matches beyond it, is aligned instead of clipped.

*Remark.* p is a prior on the ends. It is the log-odds, in score units, that a read's last bases are not genomic (an adapter, a chimera) against that they carry errors. A large p trusts the read to its end.

*Remark (the anchor forcing the path).* The glocal end is computed as a semiglobal alignment, the read end to end with the target's ends free, in a band. Q is the read's end followed by the anchor's k bases, and T the window of the target followed by the same k bases. A path that did not end on the anchor would lose k exact matches, 2k in score, so the optimum runs through the anchor. The code does not rely on this: it verifies that the path ends with k pairs on the anchor, and falls back to the X-drop end otherwise.

## 5. Pileups and Genotype Likelihoods

**Definition 5.1 (pileup).** Given alignments of reads to a reference, the *pileup* at position p is the multiset of the triples (b, q, s): the read's base b aligned to p by a pair, its Phred quality q and its strand s. The *indels anchored* at p are the insertions between p and p + 1 and the deletions of [p + 1, p + L].

The pileup keeps bases of quality at least 13 and reads of mapping quality at least 20 (§3: a read placed wrongly is a source of false variants).

The model of a diploid site has three assumptions (Li 2011, §2.1–2.2):

1. the reads are independent given the genotype;
2. each read comes from either of the two chromosomes with probability 1/2;
3. a base is read wrongly with probability ε = 10^{−q/10}, and the error lands on the other allele of the site (the biallelic simplification).

**Proposition 5.2 (Li's eq. 2).** At a site with reference allele R and alternate A, let a genotype have r copies of R. Then a read showing R has

P(R | r) = [r(1 − ε) + (2 − r) ε] / 2,   P(A | r) = [(2 − r)(1 − ε) + r ε] / 2,

and the likelihood of the reads D is L(r) = ∏_j P(b_j | r).

*Proof.* By (2) and (3), P(b | r) = Σ_{chromosome} (1/2) P(b | its allele), with P(b | a) = 1 − ε when b = a and ε otherwise. The product is (1). ∎

The same rule over n alleles gives the n(n + 1)/2 diploid genotypes, P(b | {a₁, a₂}) = [p(b | a₁) + p(b | a₂)]/2. It reduces to eq. 2 for n = 2, and the checks hold it to closed forms on small pileups.

*Remark.* Under (3) every error lands on the other allele. A less simplified model would spread an error over the three other bases, ε/3 each. The biallelic form is Li's, and it is the one checked here.

An indel is not a base, so it is counted rather than scored: of the reads that can show it, k carry it and the others do not, each misreporting with a fixed probability ε = 0.01. Which reads can show it is a question about words.

**Proposition 5.3 (the interval of ambiguity).** Let a deletion of L bases have the equivalent placements a ∈ [a₋, a₊]: ref[0, a) · ref[a + L, n) is the same word h for all of them. Then:

- a reference word ref[s, t) with t ≤ a₊ equals h[s, t);
- one with s ≥ a₋ + L equals h[s − L, t − L).

So a read that does not reach past a₊, or that begins at or after a₋ + L, reads the same on both haplotypes. An insertion with placements b ∈ [b₋, b₊] behaves the same way, with s ≥ b₋ and t ≤ b₊.

*Proof.* Take the placement a₊: h agrees with ref on [0, a₊). Take a₋: h[i] = ref[i + L] for i ≥ a₋. ∎

The caller therefore counts, as reference support, only the reads spanning the indel's whole interval: from its anchor, before the leftmost placement, to the base after the rightmost one. This is a sufficient condition for the read to tell the two apart. It raised the concordance of the insertions' genotypes in the check from 0.86 to 0.98: before it, reads ending inside a repeat counted for the reference.

## 6. The Posterior of a Genotype

**Proposition 6.1 (the infinite-sites prior; Fu 1995).** Under the standard neutral coalescent with infinite sites and θ = 4N_eμ per site, the expected number of sites per base at which the derived allele has i copies in a sample of n chromosomes is θ/i, for 1 ≤ i ≤ n − 1.

Taking the reference as ancestral gives GATK's prior for one diploid (DePristo et al. 2011):

- π(0/0) = 1 − 3θ/2, π(0/1) = θ, π(1/1) = θ/2;
- a genotype of two different alternate alleles, which needs two mutations at one site, has π(1/2) = θ², normalised.

The defaults are θ = 10⁻³ for SNVs and θ = 1.25 · 10⁻⁴ for indels.

**Definition 6.2.** The posterior is P(g | D) = L(g) π(g) / Σ_{g′} L(g′) π(g′). Write ĝ for the genotype of largest posterior. Its quality is GQ = −10 log₁₀ P(g ≠ ĝ | D), held to 99. The quality of the site is QUAL = −10 log₁₀ P(0/0 | D).

At a depth of a thousand reads L(g) underflows in linear space, so every quantity is a sum of logarithms. The constants are factored out before exponentiation, and P(g ≠ ĝ) is summed from the other genotypes, never computed as 1 − P(ĝ).

**Proposition 6.3 (Fisher's exact test; Fisher 1935).** Condition a 2 × 2 table [[a, b], [c, d]] on its margins r₁ = a + b, r₂ = c + d and c₁ = a + c. Under independence of rows and columns, a has the hypergeometric law

P(a) = C(r₁, a) C(r₂, c₁ − a) / C(r₁ + r₂, c₁).

*Proof.* Given the margins, the c₁ entries of the first column are a uniformly random c₁-subset of the r₁ + r₂ items. a counts those that fall in the first row. ∎

The two-sided p-value sums the probabilities of the tables no more probable than the observed one. For the reads by allele and strand it measures strand bias. FS = −10 log₁₀ p is GATK's FisherStrand, and sites above 60 (SNVs) or 200 (indels) are filtered. The check reproduces Fisher's tea-tasting experiment: three of four cups of each kind named rightly gives p = 34/70.

## 7. The Normal Form of an Indel

### 7.1 Representations and Their Classes

Let ref be a word of length n.

**Definition 7.1 (representation, equivalence).** A *representation* is a triple (p, R, A):

- p is a position;
- R = ref[p, p + |R|) is the reference allele;
- A is a word, the alternate allele.

It makes the haplotype h(p, R, A) = ref[0, p) · A · ref[p + |R|, n). Two representations are *equivalent* when they make the same haplotype.

Fix a haplotype h ≠ ref, and let δ = |h| − n. A representation of h is determined by p and r = |R|, since A is then the middle of h. Let P be the length of the longest common prefix of ref and h, and S that of their longest common suffix. Then (p, r) represents h if and only if

0 ≤ p ≤ P,   n − p − r ≤ S,   r ≥ 0,   r + δ ≥ 0,   p + r ≤ n.

The class of h is therefore the set of lattice points of a polygon in the (p, r) plane.

**Definition 7.2.** A representation is *parsimonious* when both alleles are non-empty and r is least among such representations of its class (so |R| + |A| = 2r + δ is least too). It is *left-aligned* when no equivalent parsimonious representation has a smaller p.

**Theorem 7.3 (the canonical form).** Each class has exactly one parsimonious, left-aligned representation. It is

r* = max{1, 1 − δ, n − S − P},   p* = max{0, n − S − r*}.

*Proof.* With both alleles non-empty, the constraints at a position p ≤ P are r ≥ max{1, 1 − δ, n − S − p}. The bound decreases in p, so its least value over 0 ≤ p ≤ P is r*, attained at p = P.

The positions attaining r* are those with n − S − p ≤ r* and p ≤ P, the interval [max{0, n − S − r*}, P]. It is non-empty, because n − S − r* ≤ P by the definition of r*. Its least point is p*. ∎

*Caveat.* The last constraint, p + r ≤ n, holds at (p*, r*) when h is not empty: p* + r* = n − S when p* > 0, and r* ≤ n because 1 − δ = n + 1 − |h| ≤ n.

*Remark.* In the language of alignment.md, §6: the class is a lattice polygon, and the canonical form is the point that the lexicographic order (r, then p) selects, a vertex of the polygon picked out by a linear functional. The equivalent representations of an indel in a repeat are the points of an edge of slope −1 (p + r constant on the side n − p − r = S). Normalising walks to that edge's end.

### 7.2 Algorithm 1

Tan, Abecasis and Kang's Algorithm 1 has two phases:

1. repeat until nothing changes: drop a last base common to all alleles; when an allele becomes empty, prepend the reference base before p to every allele and decrease p;
2. then drop a first base common to all alleles of length at least two.

Each step preserves the haplotype:

- a last base common to R and A is the base after the shortened R;
- a base prepended to both is ref[p − 1];
- a dropped first base is ref[p].

Phase 1 decreases p each time it extends, so it terminates; phase 2 shortens the alleles. Tan et al. (2015) prove that the result is parsimonious and left-aligned. By Theorem 7.3 it is therefore the canonical form.

The check confirms it by brute force. It takes 400 random indels in random tandem repeats and enumerates all 40 649 representations of their haplotypes. Every one of them normalises to one form, the shortest of least position.

*Remark.* At p = 0 nothing precedes the first base. When an allele becomes empty there, the code appends the base after the event instead, as VCF 4.2 pads an event at position 1.

### 7.3 The Gaps of an Alignment

A read that carries an indel, aligned to the reference, has a gap somewhere in the interval of ambiguity of Proposition 5.3. Different aligners, or different reads, may put it in different places.

**Definition 7.4 (shift).** Let a path of the edit graph have a run D^L of deletions, consuming target bases y[j, j + L), right after a pair (x_i, y_{j−1}). Its *shift* moves the run one column left: the pair becomes (x_i, y_{j+L−1}) and the run consumes y[j − 1, j + L − 1). The shift is allowed when y_{j−1} = y_{j+L−1}. The shift of a run of insertions is defined the same way, on the query's letters.

**Proposition 7.5.** An allowed shift preserves the path's score, and its lengths on both sequences.

*Proof.* The gap keeps its length, so it keeps its cost. The pair keeps its query letter x_i, and its target letter changes from y_{j−1} to y_{j+L−1}, which is the same letter. ∎

Each shift lowers the sum of the gaps' positions, so repeated shifting terminates (`left_align_gaps`). When a run meets another of the same kind, they merge and the score can only rise. For a read that agrees with the haplotype around a single deletion, the shift condition y_{j−1} = y_{j+L−1} is the condition for moving the deletion's placement left by one. So the left-aligned gap starts at a₋, the placement of the canonical form of Theorem 7.3. This is why every read carrying an indel reports it at one place in the pileup: 188 of the 196 truth indels in the check, the others being reads that end inside the repeat.

## 8. Evaluation

The two sets are compared allele by allele.

**Definition 8.1.** Normalise both sets by Theorem 7.3, and split multi-allelic calls into one call per alternate allele. A call is a *true positive* when a truth variant has the same contig, position, REF and ALT, a *false positive* otherwise. A truth variant matched by no call is a *false negative*. Then:

- precision = TP/(TP + FP);
- recall = TP/(TP + FN);
- F₁ = 2 · precision · recall / (precision + recall);
- the concordance of genotypes is the fraction of true positives whose number of alternate copies is the truth's.

*Remark.* This is the first level of the GA4GH benchmarking team's comparisons (Krusche et al. 2019). It does not compare haplotypes: a complex change written as an SNV and an indel in one set and as one block in the other would not match. For simulated diploids, whose variants are spaced apart, the two levels agree.

## 9. Where the Code Is

| Section | Code |
|---|---|
| §1 anchors, the strand mapped, SMEMs | `index/minimizer_index.hpp` (`anchors`), `index/fmd.hpp` (`smems`); `mapping/chain.hpp` (`seeds_of`), `mapping/mapper.hpp` (`seeds`) |
| §2 the chaining DP, the integer gap cost, backtracking, primaries, fragments | `mapping/chain.hpp` (`chain_dp`, `gap_cost`, `chain`, `select_chains`) |
| §3 mapping quality, its calibration | `mapping/chain.hpp` (`mapping_quality`), `mapping/mapper.hpp` (`map`), `mapping/evaluate.hpp` (`accuracy`) |
| §4 the alignment of a chain, the ends, clipping | `mapping/mapper.hpp` (`align_chain`); `align/wavefront.hpp`, `align/banded.hpp`, `align/xdrop.hpp`; SAM in `mapping/sam.hpp` |
| §5 pileups, eq. 2, the interval of ambiguity | `variant/pileup.hpp`; `variant/genotype.hpp` (`log10_likelihoods`); `variant/caller.hpp` (`indel_at`) |
| §6 prior, posterior, GQ, QUAL, Fisher's test | `variant/genotype.hpp` (`genotype_prior`, `posterior`, `fisher_exact`, `strand_bias`) |
| §7 the normal form, Algorithm 1, the gaps of alignments | `variant/normalize.hpp` (`normalize`, `apply`, `left_align_gaps`) |
| §8 evaluation | `variant/evaluate.hpp` (`compare`); VCF in `variant/vcf.hpp` |

## References

- Danecek, P., Auton, A., Abecasis, G. et al. (2011), The variant call format and VCFtools, *Bioinformatics* 27, 2156–2158.
- DePristo, M. A., Banks, E., Poplin, R. et al. (2011), A framework for variation discovery and genotyping using next-generation DNA sequencing data, *Nature Genetics* 43, 491–498.
- Fisher, R. A. (1935), *The Design of Experiments*, Oliver and Boyd, ch. 2.
- Fu, Y.-X. (1995), Statistical properties of segregating sites, *Theoretical Population Biology* 48, 172–197.
- Karlin, S. and Altschul, S. F. (1990), Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes, *PNAS* 87, 2264–2268.
- Krusche, P., Trigg, L., Boutros, P. C. et al. (2019), Best practices for benchmarking germline small-variant calls in human genomes, *Nature Biotechnology* 37, 555–560.
- Li, H. (2011), A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data, *Bioinformatics* 27, 2987–2993.
- Li, H. (2013), Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM, arXiv:1303.3997.
- Li, H. (2018), Minimap2: pairwise alignment for nucleotide sequences, *Bioinformatics* 34, 3094–3100.
- Li, H., Ruan, J. and Durbin, R. (2008), Mapping short DNA sequencing reads and calling variants using mapping quality scores, *Genome Research* 18, 1851–1858.
- Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M. and Yorke, J. A. (2004), Reducing storage requirements for biological sequence comparison, *Bioinformatics* 20, 3363–3369.
- Tan, A., Abecasis, G. R. and Kang, H. M. (2015), Unified representation of genetic variants, *Bioinformatics* 31, 2202–2204.
