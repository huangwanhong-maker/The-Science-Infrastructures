# Signals, Classification and Comparison of Genomes

*Status: draft. Serves learning and use.* This note covers four methods that read genomes whole rather than base by base: the Fourier spectrum of a DNA sequence and the period-3 signal of its genes; the classification of reads by exact k-mers and the lowest common ancestor; seeds chosen by their content (syncmers and strobemers); and the alignment of two genomes by their maximal unique matches. The code is in the genomics domain:

- [`signal/spectrum.hpp`](../signal/include/genomics/signal/spectrum.hpp), on calx's FFT ([`fft.hpp`](../../../foundations/calx/numerics/transforms/include/calx/transforms/fft.hpp));
- [`classify/taxonomy.hpp`](../classify/include/genomics/classify/taxonomy.hpp) and [`classify/kraken.hpp`](../classify/include/genomics/classify/kraken.hpp);
- [`kmers/strobemers.hpp`](../kmers/include/genomics/kmers/strobemers.hpp), beside the syncmers of [`sequence/kmer.hpp`](../sequence/include/genomics/sequence/kmer.hpp);
- [`comparative/mums.hpp`](../comparative/include/genomics/comparative/mums.hpp) and [`comparative/whole_genome.hpp`](../comparative/include/genomics/comparative/whole_genome.hpp), on morphe's suffix arrays ([`suffix_array.hpp`](../../../foundations/morphe/formal/include/morphe/formal/suffix_array.hpp)) and increasing subsequences ([`increasing.hpp`](../../../foundations/morphe/combinatorics/include/morphe/combinatorics/increasing.hpp)).

Nothing general is written twice: the transform is calx's, the suffix array, the LCP array and the heaviest increasing subsequence are morphe's (the last added to morphe's `combinatorics` for this purpose), and the open-addressing table of the classifier is the k-mer counter's.

Throughout, a DNA sequence is a word x = x₀ ⋯ x_{N−1} over {A, C, G, T, N}, with A = 0, C = 1, G = 2, T = 3 (the domain's codes), and N a base of unknown identity.

## 1. The Spectrum of a DNA Sequence

### 1.1 Indicator Sequences and the Fourier Transform of ℤ_N

**Definition 1.1.** For each base b ∈ {A, C, G, T}, the *indicator sequence* of x is u_b: ℤ_N → {0, 1}, u_b(n) = 1 if x_n = b and 0 otherwise (Voss 1992). A base N is 0 in all four. Let ν = Σ_b Σ_n u_b(n) be the number of bases other than N.

The positions are read in ℤ_N, the cyclic group of order N, so that the sequence becomes a function on a group.

**Definition 1.2.** The characters of ℤ_N are the homomorphisms χ_k: ℤ_N → ℂ^×, χ_k(n) = e^{2πi kn/N}, one for each k ∈ ℤ_N. The *Fourier transform* of f: ℤ_N → ℂ is f̂(k) = Σ_n f(n) \overline{χ_k(n)} = Σ_n f(n) e^{−2πi kn/N}. With U_b = û_b, the *spectrum* of x is

  S(k) = Σ_b |U_b(k)|²,  k ∈ ℤ_N.

The characters χ_k are the irreducible representations of ℤ_N, all of dimension one, and k ↦ χ_k is an isomorphism of ℤ_N with its dual group. The frequencies of the discrete Fourier transform are therefore exactly the irreducible representations of the cyclic group: the DFT *is* the Fourier transform of ℤ_N.

*Remark (the group Fourier transform).* For a finite group G the Fourier transform of f: G → ℂ takes each irreducible representation ρ to the matrix f̂(ρ) = Σ_g f(g) ρ(g), and Plancherel's theorem makes it an isomorphism of ℂ[G] with ⊕_ρ M_{dim ρ}(ℂ) (Serre 1977, §6.2; Diaconis 1988, ch. 2). The user's paper on group-equivariant grammars speaks of a group Fourier transform layer; the spectrum of a DNA sequence is its abelian case, for G = ℤ_N acting on positions by translation, where every ρ is a character and every block has dimension one. The non-abelian forms, and the exact transforms of finite abelian groups in cyclotomic fields, belong to morphe's `algebra` module (see also [group grammars](group-grammars.md), §§3 and 7); here the transform is calx's, in floating point, because a spectrum of a genome is a numerical object, except at N/3, where §1.3 makes it exact.

The code computes the four transforms by calx's FFT, radix 2 when N is a power of two and Bluestein's chirp otherwise, in O(N log N) (`indicator_transforms`, `spectrum`). The checks hold both to the direct sum over ℤ_N, with the phase kn reduced modulo N in integers, at N = 1024 and N = 999, to within 10⁻¹⁶ N.

### 1.2 Parseval and the Symmetries of the Spectrum

**Proposition 1.3 (Parseval).** Σ_k |f̂(k)|² = N Σ_n |f(n)|². Hence Σ_k S(k) = N ν, and the mean of the spectrum over the N frequencies is S̄ = ν.

*Proof.* Σ_k |f̂(k)|² = Σ_{n,m} f(n) \overline{f(m)} Σ_k e^{−2πi k(n−m)/N}, and the inner sum is N when n = m and 0 otherwise (orthogonality of the characters). For an indicator sequence |u_b(n)|² = u_b(n), and Σ_b Σ_n u_b(n) = ν. ∎

**Proposition 1.4.** (i) S(N − k) = S(k). (ii) The reverse complement y of x, y_n = x̄_{N−1−n} with Ā = T, C̄ = G, has the same spectrum.

*Proof.* (i) u_b is real, so U_b(N − k) = \overline{U_b(k)}. (ii) u^y_b(n) = u^x_{b̄}(N − 1 − n), so U^y_b(k) = Σ_m u^x_{b̄}(m) e^{−2πi k(N−1−m)/N} = e^{−2πi k(N−1)/N} \overline{U^x_{b̄}(k)}, and |U^y_b(k)| = |U^x_{b̄}(k)|. Summing over b, which b̄ permutes, gives S^y = S^x. ∎

So the spectrum, like a canonical k-mer, is a function of the double strand: a gene on the reverse strand shows the same signal as one on the forward strand.

**Proposition 1.5 (periodic sequences).** If x has period p and p divides N, then U_b(k) = 0 unless N/p divides k. In particular (ACG)^m, N = 3m, has S(N/3) = 3m² = N²/3 and no other non-zero value in (0, N/2].

*Proof.* Write n = r + jp with 0 ≤ r < p and 0 ≤ j < N/p. Then U_b(k) = Σ_r u_b(r) e^{−2πi kr/N} Σ_j e^{−2πi kj/(N/p)}, and the inner sum is N/p when (N/p) | k and 0 otherwise. For (ACG)^m each of A, C, G occupies one residue mod 3, so |U_b(m)| = m for three bases and 0 for T. ∎

The checks verify it on (ACG)^400 and on random words of periods 5, 7, 10 and 21 repeated 60 times.

### 1.3 The Period-3 Component Is a Transform on ℤ₃

Protein-coding DNA is read in codons, and the three positions of a codon use the bases unequally. Its spectrum has a peak at k = N/3 (Tiwari et al. 1997; Anastassiou 2001). When 3 divides N this component has an exact arithmetic of its own.

**Proposition 1.6.** Let 3 | N, ζ = e^{−2πi/3}, and let c_{b,r} be the number of positions n ≡ r (mod 3) with x_n = b. Then

  U_b(N/3) = c_{b,0} + c_{b,1} ζ + c_{b,2} ζ²,

an Eisenstein integer, and

  S(N/3) = Σ_b N(U_b(N/3)) = Σ_b (c_{b,0}² + c_{b,1}² + c_{b,2}² − c_{b,0}c_{b,1} − c_{b,1}c_{b,2} − c_{b,2}c_{b,0}) = ½ Σ_b Σ_{r<s} (c_{b,r} − c_{b,s})²,

where N(α) = α ᾱ is the norm of ℤ[ζ]. S(N/3) is an integer, unchanged by any relabelling of the three residues.

*Proof.* The character χ_{N/3}(n) = e^{2πi n/3} depends only on n mod 3: it is the composite of the quotient map π: ℤ_N → ℤ₃ with a character of ℤ₃. So U_b(N/3) = Σ_n u_b(n) ζ^{n mod 3} = Σ_r (π_* u_b)(r) ζ^r, the Fourier transform on ℤ₃ of the push-forward π_* u_b = (c_{b,0}, c_{b,1}, c_{b,2}). With 1 + ζ + ζ² = 0, c₀ + c₁ζ + c₂ζ² = (c₀ − c₂) + (c₁ − c₂)ζ, and N(a + bζ) = a² − ab + b² since ζ + ζ̄ = −1 and ζζ̄ = 1; expanding gives the quadratic form, which is symmetric in c₀, c₁, c₂ and equals half the sum of the squared differences. ∎

*Remark.* Since Σ_{r<s} (c_r − c_s)² = 3 Σ_r (c_r − c̄)² for three numbers of mean c̄, S(N/3) = (3/2) Σ_b Σ_r (c_{b,r} − c̄_b)²: the period-3 power is, up to the factor 3/2, the sum of squared deviations of the bases' counts across the codon positions, an analysis of variance of composition by codon position. It vanishes exactly when every base is spread evenly over the three positions.

So the period-3 signal needs no transform: twelve counts determine it, and a window sliding by one base changes two of them, so the profile along a genome of n bases costs O(n) in exact integers (`CodonCounts`, `period3_power`, `period3_profile`). The checks hold the exact value to S(N/3) from calx's FFT on 50 random stretches and the sliding profile to counts made afresh for every window.

### 1.4 The Strength of the Signal in Random and in Coding Sequence

The *relative strength* of the period-3 signal is P = S(N/3)/S̄ = S(N/3)/ν, the peak over the spectrum's mean (Tiwari et al. 1997 measure coding potential by this ratio; with the mean taken over all N frequencies, Proposition 1.3 makes it S(N/3)/ν exactly).

**Proposition 1.7.** Let 3 | N and let the bases be independent, base n being b with probability p_{b,n mod 3}. Then

  E S(N/3) = (N/3) Σ_{b,r} (p_{b,r} − p_{b,r}²) + (N²/9) Σ_b |Σ_r p_{b,r} ζ^r|².

In particular, when the composition does not depend on the position (p_{b,r} = p_b), E S(N/3) = N (1 − Σ_b p_b²), and E P = 1 − Σ_b p_b², which is 3/4 at equal frequencies.

*Proof.* E|Σ_n u_n ζ^n|² = Σ_n E u_n² + Σ_{n≠m} E u_n E u_m ζ^{n−m} = Σ_n p_n − Σ_n p_n² + |Σ_n p_n ζ^n|², using u_n² = u_n and independence. Each residue occurs N/3 times, so Σ_n p_n = (N/3) Σ_r p_r, Σ_n p_n² = (N/3) Σ_r p_r², Σ_n p_n ζ^n = (N/3) Σ_r p_r ζ^r. Sum over b. When p_{b,r} = p_b the last term vanishes because 1 + ζ + ζ² = 0, and Σ_b p_b = 1. ∎

The proposition separates a floor, of order N, from a signal of order N²: a position-dependent composition makes the expected strength grow linearly with the window. The checks verify the floor on 25 000 disjoint windows of random sequence at GC contents 0.3, 0.5 and 0.7, within four standard errors.

**What *Escherichia coli* shows.** On *E. coli* K-12 MG1655 (NC_000913.3), windows of W = 120 bases taken every third base were divided into those lying wholly inside annotated CDSs (pseudogenes excluded) and those touching no annotated gene of any kind (measured 2026-10-24):

| W | windows in CDSs | mean P | intergenic windows | mean P | ratio | AUC |
|---|---|---|---|---|---|---|
| 60 | 1 257 583 | 1.218 | 108 442 | 0.739 | 1.65 | 0.737 |
| 120 | 1 188 315 | 1.685 | 69 170 | 0.744 | 2.26 | 0.857 |
| 180 | 1 120 112 | 2.145 | 43 011 | 0.742 | 2.89 | 0.919 |
| 240 | 1 053 560 | 2.597 | 26 350 | 0.742 | 3.50 | 0.949 |
| 351 | 936 743 | 3.408 | 10 644 | 0.747 | 4.57 | 0.972 |

The AUC is the probability that a coding window outscores an intergenic one. Two facts stand out, and both are Proposition 1.7. The intergenic mean does not move with W and sits at 1 − Σ_b p_b² of the intergenic composition (0.741): between genes the sequence has no period 3, as for random sequence. The coding mean grows nearly linearly with W, by about 0.0075 a base, from an intercept near 0.79: under the model of Proposition 1.7 the signal term (1/9) Σ_b |Σ_r p_{b,r} ζ^r|² of the codon-position composition of *E. coli*'s genes would be near 0.0075, the slope. The checks assert, at W = 120, a coding mean more than twice the intergenic one, an AUC above 0.83, and the intergenic mean within 0.03 of 1 − Σ p_b²; at W = 351, a ratio above 4 and an AUC above 0.96.

*Caveat.* The signal measures coding potential, not genes. Short genes and short windows separate poorly (AUC 0.74 at 60 bases), a window across a gene's end is diluted, and the strength says nothing of the frame or the strand. A gene finder (`annotation/genes.hpp`) reads codons; the spectrum is a reference-free summary that needs no training.

## 2. Classification of Reads by k-mers and the Lowest Common Ancestor

### 2.1 The Taxonomy of Lineages

A GenBank record names its organism and its lineage, the taxa above it from the top down. The taxonomy built here is the trie of the records' lineages under one root, with each organism a leaf below its lineage (`Taxonomy::add`). No taxonomy is typed in: the fourteen cached genomes give a tree of 73 taxa, in which, for example, the human and Neanderthal mitochondria are two leaves of *Homo*.

**Definition 2.1.** On the taxa write a ≤ b when a lies on the path from the root to b (a is *above* b, or is b). The *clade* of a is {b : a ≤ b}.

**Proposition 2.2.** (T, ≤) is a meet-semilattice whose meet a ∧ b is the lowest common ancestor of a and b, and the root is its least element. Hence the LCA is associative, commutative and idempotent, and the root absorbs it.

*Proof.* The elements below b in ≤ (its ancestors) form a chain, the root path of b. The common lower bounds of a and b are the common part of two root paths, a chain with a greatest element, the deepest common ancestor. A meet in any poset is associative, commutative and idempotent. ∎

The LCA is computed by lifting the deeper taxon to the other's depth and then both together, O(depth) (`lca`). *Caveat.* A flat file's lineage has no ranks, so depth here is not a rank: Cercopithecinae sits between Cercopithecidae and *Macaca*, while Hominidae is followed directly by *Homo*; and the lineages of the primate mitochondria stop at the genus, so a species and a subspecies are siblings.

### 2.2 The Database

**Definition 2.3** (Wood and Salzberg 2014). Let each reference genome g carry a taxon t(g) and the set K(g) of its canonical k-mers. The database maps each k-mer x in ⋃_g K(g) to

  D(x) = ⋀ { t(g) : x ∈ K(g) },

the LCA of the taxa of the genomes that contain it.

**Proposition 2.4.** Building D one k-mer occurrence at a time, each meeting the value already stored, gives D whatever the order of the genomes and of their k-mers. D(x) is the deepest taxon whose clade contains every genome that has x.

*Proof.* The stored value after any sequence of insertions is the meet of the taxa inserted for x, which by Proposition 2.2 does not depend on their order or their repetitions. The second statement is the definition of the meet in a tree. ∎

The checks compare D with the brute-force LCA of the genomes containing each k-mer, on twelve sets of ten genomes sharing segments on both strands, and with the database built in the opposite order. On the fourteen cached genomes at k = 31 the database holds 4 772 316 canonical k-mers, 23 049 of them on an inner taxon and 2 958 at the root (shared between the bacterium or its phages and the eukaryotes or the other viruses), built in half a second.

### 2.3 The Heaviest Root-to-Leaf Path

**Definition 2.5** (ibid.). Let a read have h(t) of its k-mers with D(x) = t, and H = {t : h(t) > 0}. The *classification tree* is C = {a : a ≤ t for some t ∈ H}; its leaves are its maximal elements. A leaf ℓ has the weight W(ℓ) = Σ_{a ≤ ℓ} h(a). The read's label is the leaf of greatest weight, or, when several tie, the meet of those leaves. A read with H empty is unclassified.

**Proposition 2.6.** Define W(t) = Σ_{a ≤ t} h(a) for every t ∈ H. The maximum of W over the leaves of C equals its maximum over H, and every t ∈ H attaining it is a leaf of C.

*Proof.* The leaves of C lie in H, since a maximal element of C is below some t ∈ H and maximal, so equal to it. Let t ∈ H not be a leaf: some u ∈ C has t < u, and u lies below a leaf ℓ ∈ H, so t < ℓ and W(ℓ) ≥ W(t) + h(ℓ) > W(t), the hits being positive. So no non-leaf attains the maximum over H, and the maxima over H and over the leaves coincide. ∎

So the label is found among the hit taxa themselves, O(m² · depth) for m distinct taxa, without building C (`heaviest_path`). The checks compare it with the tree built explicitly, its leaves weighed one by one, on 3000 random sets of hits.

*Our thresholds.* A read can be required to have `min_hits` k-mers in the database (1 is Kraken's rule). With a confidence τ > 0, its label then climbs toward the root until the hits in its clade, Σ_{label ≤ a} h(a), reach τ times the read's k-mers; if even the root falls short, the read is unclassified.

### 2.4 What the Classification of a Read Can Say

**Proposition 2.7.** Let R be the meet of the taxa of the reference genomes that share at least one k-mer with a read. The read's label, if any, lies in the clade of R.

*Proof.* Each hit taxon D(x) is the meet of the taxa of some of those genomes, so R ≤ D(x). The label is a hit taxon or a meet of hit taxa, and a meet of elements of the clade of R lies in it, R being a lower bound of each and hence of their meet. ∎

So a read is never labelled outside the clade spanned by the genomes it shares k-mers with, and no more specifically than a single genome's taxon. Within that clade the method can be specific to a fault: a read of a genome absent from the database goes to whichever relatives share its k-mers, a relative's own leaf among them.

**What the fourteen genomes show** (measured 2026-10-24; 600 simulated reads of 150 bases a genome, with the simulator's errors of 0.1 to 1 per cent along the read):

- **Every genome in the database.** 0.9996 of the 8400 reads are labelled on their own lineage, 0.963 at their genome's leaf. Genomes with no close relative in the database (E. coli, lambda, phiX174, SARS-CoV-2, gorilla, gibbon, macaque, baboon) have 0.997 of their reads at their leaf; human and Neanderthal 0.79 each, the rest at *Homo*, because a 150-base read of one often carries no k-mer that tells it from the other. The labels are the same on one thread and on four.
- **Each genome left out in turn.** E. coli's, phiX174's and SARS-CoV-2's reads are all unclassified: nothing else in the database shares their 31-mers. Lambda's are unclassified (0.892) or labelled E. coli K-12, and nothing else: 3 863 positions of K-12 begin one of lambda's 31-mers, 1 024 of them inside genes whose products the K-12 record names as prophage proteins (BorD, NohA, RzpD, RzoR, StfQ, StfR, and the Rac and DLP12 prophages' LomR and RenD), 1 633 in their neighbours (the tail-fibre assembly proteins TfaD, TfaQ and TfaR, the DNA-packaging protein NohD, YbcV, YbcW and others) and 1 206 outside annotated CDSs: the sequence lambda shares lies in K-12's cryptic prophages and around them. A primate mitochondrion's reads are never labelled off the primates, as Proposition 2.7 requires (they share k-mers with primates only); of those labelled, 0.880 fall within the clade of the deepest ancestor still in the database, mostly on a relative's own leaf (a chimpanzee's reads on the bonobo, a human's on the Neanderthal), and the rest on more distant apes or monkeys, at k-mers that the closest relative has changed and a more distant one kept. The checks assert at least 0.99 unclassified for the three, lambda's at least 0.85 unclassified and otherwise on E. coli, no primate read off the primates, and at least 0.85 within the deepest shared clade.

*Caveat.* That a read of an absent organism goes only to a shared ancestor or nowhere is not a property of Kraken's rule: it goes to the leaf of whichever relative shares the most of its k-mers. Specificity beyond the database's taxa is never warranted, and the confidence threshold reduces such labels only when the relative's clade holds a small share of the read's k-mers.

## 3. Seeds Chosen by Their Content

### 3.1 Syncmers and Their Density

**Definition 3.1** (Edgar 2021). Fix s < k and an order on s-mers (here hash64 of the canonical s-mer). A k-mer with m = k − s + 1 s-mers is a *closed syncmer* when its least s-mer (the first, on ties) is its first or its last, and an *open syncmer* with offset t when it is the (t + 1)-th.

A syncmer is decided by the k-mer alone: whether a k-mer is chosen does not depend on its neighbours, so a k-mer conserved between two sequences is chosen in both or in neither, which minimizers do not guarantee. The code is `sequence/kmer.hpp`'s.

**Proposition 3.2.** (i) If the order on s-mers is a uniformly random total order independent of the sequence and the m s-mers of a k-mer are distinct, the least of them is at each of the m offsets with probability 1/m; closed syncmers then have density 2/m and open ones 1/m (Edgar 2021). (ii) For a fixed order, on a sequence of independent uniform bases, the expected density of closed syncmers is exactly the fraction of the 4^k k-mers that are closed syncmers.

*Proof.* (i) By symmetry of the random order over distinct elements. (ii) The k-mer at each position is uniform on the 4^k k-mers, and the density's expectation is, by linearity, the probability that a uniform k-mer is a closed syncmer. ∎

The check enumerates every k-mer for (k, s) = (9, 4), (11, 7) and (12, 8): the exact densities under hash64 are 0.33408, 0.40376 and 0.39952, beside Edgar's 1/3, 2/5 and 2/5, and a random sequence of two megabases agrees with the exact values within five standard deviations (the count of a (k − 1)-dependent sequence of indicators has variance at most (2k − 1) n p(1 − p)). The departures from 2/m come from the hash, a fixed order rather than a random one, and from s-mers repeated within a k-mer; they are below one per cent here.

### 3.2 Strobemers

**Definition 3.3** (after Sahlin 2021). Fix an order n ≥ 2, a strobe length ℓ and window bounds w_min ≤ w_max, with w = w_max − w_min + 1. The strobemer at position i has its first strobe the ℓ-mer at i and its j-th strobe (j ≥ 2) chosen among the ℓ-mers starting in W_j(i) = [i + w_min + (j − 2) w_max, i + (j − 1) w_max]:

- a *minstrobe* chooses the ℓ-mer of least hash h (leftmost on ties);
- a *randstrobe* chooses the ℓ-mer p minimising g(H_{j−1}, h(p)), where H₁ = h(strobe 1) and H_j = g(H_{j−1}, h(strobe j)); here g(H, h) = mix64(H ⊕ h), a bijection of 64-bit words, so that a randstrobe's later strobes depend on its earlier ones.

The seed's hash is H_n. Minstrobes cost O(L) for L bases, one sliding minimum serving every window; randstrobes O(L n w). *Remark.* The windows of order 2 are Sahlin's; for higher orders they are laid one after another at fixed offsets from i, and the function g is ours, where Sahlin's randstrobes use a different function of the same purpose.

**Proposition 3.4 (survival of an indel).** Take n = 2 and a one-base insertion or deletion at a position d with i + ℓ ≤ d < i + w_min, so between the first strobe and the window. Suppose the w + 1 ℓ-mers involved (the window's before and after the edit) are distinct, and the order on them is a uniformly random order independent of the sequence (for randstrobes, the order p ↦ g(H₁, h(p)), H₁ being unchanged). Then the strobemer at i of the edited sequence has the same strobes, as words, with probability (w − 1)/(w + 1) = 1 − 2/(w + 1). A k-mer spanning the strobemer contains d and is changed.

*Proof.* The first strobe is untouched. After a deletion the window holds the old window's ℓ-mers except its first, and the ℓ-mer after its end; after an insertion, all but its last, and the ℓ-mer before its start; the w − 1 others are common to both. The strobe is the same word exactly when the least of the w + 1 ℓ-mers lies among the w − 1 common ones: if it is the one lost, the old window chose it and the new one cannot; if it is the one gained, the new window chooses it and the old one could not. Under a random order the least is equally likely to be any of the w + 1. ∎

The checks verify the definitions against brute force for orders 2 to 4 on random sequences with bases N, construct a minstrobe that survives a deletion while the 85-base k-mer of its span is lost, and measure over 3000 random one-base indels (ℓ = 15, w = 31) a survival of 0.942 for minstrobes and 0.937 for randstrobes against 1 − 2/32 = 0.938, and of 0 for k-mers of the same span.

## 4. Whole-Genome Alignment by Maximal Unique Matches

### 4.1 Maximal Unique Matches

Let A and B be genomes, read on one strand of B (or on its reverse complement, for the other).

**Definition 4.1** (Delcher et al. 1999). A *unique match* is a word w, free of N, occurring exactly once in A, at p, and exactly once in B, at q. It is a *maximal unique match* (MUM) when no unique match strictly contains it.

**Proposition 4.2.** If a unique match u (at p, q) is a factor of a unique match v (at p′, q′), then p − p′ = q − q′: u is v's factor at one offset in both genomes. Hence a unique match is a MUM exactly when it is left- and right-maximal: the letters before its two occurrences differ (or one is missing or N), and so do the letters after.

*Proof.* Let u occur in v at offset o. Then u occurs in A at p′ + o and in B at q′ + o, and uniqueness gives p = p′ + o and q = q′ + o. So a unique match strictly contained in another extends, at one end at least, by equal letters in both genomes; conversely an extension by equal letters is again a unique match, since a word containing u occurs at most where u does. ∎

**Proposition 4.3 (MUMs in the suffix array).** Cut A and B into their runs free of N and join the runs into one word T, each followed by a separator of its own, distinct from every letter and from each other. Let SA and LCP be its suffix array and LCP array (LCP[i] the longest common prefix of the suffixes SA[i − 1] and SA[i]). The MUMs of length at least m are exactly the pairs of suffixes SA[i − 1], SA[i] such that

- LCP[i] ≥ m, LCP[i − 1] < LCP[i], and LCP[i + 1] < LCP[i] (where they exist);
- one suffix starts in a run of A and the other in a run of B;
- the letters before them differ, or one of them starts its run;

the MUM being their common prefix, of length LCP[i].

*Proof.* The suffixes that begin with a word w of length ℓ form an interval of SA, of as many suffixes as w has occurrences in the runs, which are its occurrences in A and B free of N, since no word of letters spans a separator. w is a unique match exactly when that interval is two suffixes, one from each genome: two neighbours i − 1, i with LCP[i] ≥ ℓ and LCP[i − 1], LCP[i + 1] < ℓ. Taking ℓ = LCP[i] makes w right-maximal, the next letters differing (or a separator, which matches nothing); conversely a right-maximal unique match has length exactly LCP[i]. Left-maximality is the third condition. By Proposition 4.2 these are the MUMs. ∎

One pass over SA finds them, in O(|A| + |B|) with SA-IS (`maximal_unique_matches`). The checks compare them with brute force (every pair of positions extended to its maximal exact match, kept when long enough, left-maximal and unique in both genomes) on 400 random pairs sharing segments on both strands, with repeats and bases N.

**Proposition 4.4.** No MUM's interval in A contains another MUM's interval in A, and likewise in B. Hence MUMs have distinct starts in A and in B, and a sequence of MUMs whose starts increase in A and in B has its ends increasing too.

*Proof.* Let [a₂, a₂ + ℓ₂) ⊆ [a₁, a₁ + ℓ₁). Then MUM₂'s word occurs in MUM₁'s at offset a₂ − a₁, so in B at b₁ + a₂ − a₁, and uniqueness in B gives b₂ = b₁ + a₂ − a₁: MUM₂ is a factor of MUM₁ at one offset in both genomes, a unique match strictly contained in another unless the two are equal, against maximality. The argument in B is the same. ∎

### 4.2 The Chain

**Definition 4.5.** A *chain* is a sequence of MUMs of one strand whose starts strictly increase in A and in B (B's reverse complement, for the reverse strand). Its *weight* is the sum of their lengths.

A chain is an alignment that keeps the order of both genomes. MUMmer orders the MUMs by their positions in one genome and takes an increasing subsequence of their positions in the other (Delcher et al. 1999); the heaviest chain is the heaviest increasing subsequence of the starts in B with the MUMs sorted by their starts in A, which morphe's `heaviest_chain` finds in O(m log m) by a Fenwick tree of prefix maxima over the ranks of the values (Jacobson and Vo 1992; Fredman 1975 for unit weights). The checks hold it to the quadratic recurrence over its MUMs on random related pairs with moved segments.

*Caveat.* Consecutive anchors of a chain may overlap, in a repeat or beside an indel in a homopolymer; the weight then counts the shared bases twice. *Rationale.* The mapper's chaining (`mapping/chain.hpp`, minimap2's) is a DP over a bounded number of predecessors within a gap cut-off and a band, made for the seeds of one read; between whole genomes the MUMs are few and their gaps as long as a rearrangement, and the chain wanted is the global one, exactly.

### 4.3 Gaps, Blocks and the Dot Plot

Between consecutive anchors lie length_a bases of A and length_b of B, negative when the anchors overlap. Equal lengths of one are a single-base substitution; equal greater lengths a stretch of substitutions; unequal lengths an indel of |length_b − length_a| bases, with or without substitutions beside it. A chain is cut into collinear *blocks* where a gap exceeds a bound in either genome. The dot plot draws every MUM as a segment of slope +1 (forward) or −1 (reversed), the chained ones marked.

One heaviest chain a strand recovers one collinear backbone a strand. An inversion is recovered whole: its MUMs are those of the reverse strand, collinear in the reverse complement's coordinates, and the forward chain jumps over it. Two inversions, a translocation or a duplication need local chains, joined afterwards.

### 4.4 What a Simulated Inversion and the Primate Mitochondria Show

**A known rearrangement.** A random genome of 200 kb, its segment [60 000, 90 000) inverted, then mutated with 232 SNVs and indels by `simulate` (whose liftover is the truth), aligned with MUMs of 20 bases or more (measured 2026-10-24): three blocks, +[0, 60 002), −[60 000, 90 000), +[89 998, 200 000), the inversion recovered to the base and the forward flanks overrunning it by two bases that match by chance. Every base of the chained anchors lies where the truth puts it, but for 31 of 199 760 within eight bases of an anchor's end: bases that match by chance across a breakpoint or beside an inserted base, and indels in repeats, whose placement is ambiguous, placed elsewhere than the simulator wrote them. Each of the 152 isolated SNVs is a single-base gap at its place, each of the 178 single-base gaps lies at an SNV, and each of the 38 isolated indels is a gap of its length near its place. The anchors cover 0.999 of both genomes.

**Human and chimpanzee mitochondria.** NC_012920.1 (16 569 bases) against NC_001643.1 (16 554), MUMs of 20 or more: 233 forward MUMs, none reversed, 224 of them chained over 15 976 human and 15 968 chimpanzee bases (0.964 of each), with offsets b − a between −586 and −576; the anchors cover 0.465 of each genome, and the gaps are 40 single-base substitutions, 173 stretches of substitutions and 10 indels. The nine MUMs left out of the chain lie on the same diagonal modulo the chimpanzee's length: the two records begin 582 bases apart on the circle, and a linear chain cannot pass the origin. With the chimpanzee's circle rotated by −582 the chain takes all 234 MUMs and spans 0.999 of both genomes. The checks assert a span of at least 0.95 (0.99 after rotation), every forward MUM on the diagonal modulo |B|, and no reversed MUM.

## 5. Notation and Code

| Symbol | Meaning | Code |
|---|---|---|
| u_b, U_b | indicator sequence of base b and its transform | `indicator`, `indicator_transforms` |
| S(k), S̄ | spectrum, its mean (= ν) | `spectrum` |
| c_{b,r} | count of b at positions ≡ r (mod 3) | `CodonCounts::c` |
| S(N/3), P | period-3 power, its strength S(N/3)/ν | `CodonCounts::power`, `Period3::strength` |
| ≤, ∧ | the order of the taxonomy (above), its meet (LCA) | `Taxonomy::above`, `Taxonomy::lca` |
| D(x) | the taxon of a k-mer in the database | `Database::lookup` |
| h, W | a read's hits per taxon, the weight of a root path | `heaviest_path`, `clade_hits` |
| W_j(i), ℓ, w | strobe windows, strobe length, window width | `StrobemerShape` |
| MUM (a, b, ℓ) | a maximal unique match | `Match` |
| chain, gap, block | heaviest increasing chain of MUMs and what lies between | `Chain`, `Gap`, `Block`, `chain_matches`, `blocks`, `align_genomes` |

## 6. References

- Anastassiou, D. (2001), Genomic signal processing, *IEEE Signal Processing Magazine* 18(4), 8–20.
- Delcher, A. L., Kasif, S., Fleischmann, R. D., Peterson, J., White, O. and Salzberg, S. L. (1999), Alignment of whole genomes, *Nucleic Acids Research* 27, 2369–2376.
- Diaconis, P. (1988), *Group Representations in Probability and Statistics*, Institute of Mathematical Statistics.
- Edgar, R. (2021), Syncmers are more sensitive than minimizers for selecting conserved k-mers in biological sequences, *PeerJ* 9, e10805.
- Fredman, M. L. (1975), On computing the length of longest increasing subsequences, *Discrete Mathematics* 11, 29–35.
- Jacobson, G. and Vo, K.-P. (1992), Heaviest increasing/common subsequence problems, *Combinatorial Pattern Matching (CPM 1992)*, LNCS 644.
- Sahlin, K. (2021), Effective sequence similarity detection with strobemers, *Genome Research* 31, 2080–2094.
- Serre, J.-P. (1977), *Linear Representations of Finite Groups*, Springer.
- Tiwari, S., Ramachandran, S., Bhattacharya, A., Bhattacharya, S. and Ramaswamy, R. (1997), Prediction of probable genes by Fourier analysis of genomic sequences, *CABIOS* 13, 263–270.
- Voss, R. F. (1992), Evolution of long-range fractal correlations and 1/f noise in DNA base sequences, *Physical Review Letters* 68, 3805–3808.
- Wood, D. E. and Salzberg, S. L. (2014), Kraken: ultrafast metagenomic sequence classification using exact alignments, *Genome Biology* 15, R46.
