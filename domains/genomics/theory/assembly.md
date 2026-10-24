# Assembly on de Bruijn Graphs

*Status: draft. Serves learning and use.* How a genome is read back from the short words that sequencing gives: the statistics of k-mers in a sample of reads, the de Bruijn graph of a double strand and its unitigs, the cleaning of what errors and heterozygous sites make, and assembly as an Eulerian path, whose ambiguity the BEST theorem measures. The code is the `kmers` and `assembly` modules of the genomics domain, with the graph theory and the exact determinant in morphe (`combinatorics/eulerian.hpp`, `algebra/determinant.hpp`); §7 says where each result is implemented.

Throughout, Σ = {A, C, G, T} with the complement A ↔ T, C ↔ G, a word x = x₁ … x_k has reverse complement x̄ = x̄_k … x̄₁, and Σ* is the free monoid on Σ. A genome is a word g of length G; a read is a subword of g or of ḡ, possibly with errors.

## 1. The k-mer Spectrum

**Definition 1.1 (composition, spectrum).** The *k-mer composition* of a word s is the multiset C_k(s) of its |s| − k + 1 subwords of length k. The *canonical form* of a k-mer x is min(x, x̄) in the lexicographic order. For a sample of reads, the *count* n(x) of a canonical k-mer is the number of its occurrences, in either orientation, in all reads; the *spectrum* is h(c) = #{x : n(x) = c}, c ≥ 1.

The canonical form is a section of the quotient of Σ^k by the involution x ↦ x̄: a k-mer and its reverse complement are one piece of a double strand, read from either side, and sequencing does not say which side it read.

**Proposition 1.2 (the coverage of a k-mer).** Let N reads of length L start at the points of a Poisson process of rate N/G on a genome of length G, on either strand, without errors. Then the number of reads that contain a given occurrence of a k-mer is Poisson distributed with mean

λ = N(L − k + 1)/G = c (L − k + 1)/L,

where c = NL/G is the base coverage. If each read independently contains the k-mer without error with probability q, the number of error-free copies is Poisson with mean qλ.

*Proof.* A read contains the occurrence at [p, p + k) exactly when it starts in [p + k − L, p], a window of L − k + 1 positions, so the count is Poisson of mean (N/G)(L − k + 1). Keeping each point of a Poisson process independently with probability q gives a Poisson process of rate q times the original (the thinning theorem; Kingman 1993, *Poisson Processes*, §5.1). ∎

The simulator's errors rise linearly along a read, so q is the mean over the L − k + 1 positions of a k-mer in a read of ∏ (1 − pᵢ)(1 − δ) over its bases, pᵢ the substitution rate at base i and δ the indel rate. The checks compare the fitted coverage with this qλ (kmers/tests/check_spectrum.cpp).

**Remark 1.3 (the shape of the spectrum).** A single-copy k-mer of a haploid genome is counted Poisson(λ) times; a k-mer of a repeat of j copies Poisson(jλ) times; and an error at a base of a read makes up to k k-mers that are, almost always, seen once. The spectrum is therefore a steep descent from c = 1, a first minimum, and a peak at λ with smaller peaks at 2λ, 3λ, …. In a diploid genome a window of k bases that holds a heterozygous site gives two k-mers, one from each haplotype, each at λ/2. The *error threshold* is the first minimum, the least c with h(c + 1) ≥ h(c); k-mers counted at least that often are *solid*.

**Proposition 1.4 (genome size).** Without errors, Σ_c c h(c) = N(L − k + 1), and N(L − k + 1)/λ = G, the number of positions of the genome, repeats counted with their multiplicity.

*Proof.* Each read has L − k + 1 k-mers, so the sum of the counts is N(L − k + 1) = λG. ∎

With errors, the k-mers below the threshold are discarded first: Ĝ = Σ_{c ≥ t} c h(c) / λ (Li and Waterman 2003). A heterozygous window contributes two k-mers at λ/2, which is λ in all, so Ĝ is the haploid length.

**The mixture model.** Above the threshold, h is fitted as the mixture Σ_j w_j Poisson(m_j λ), m ∈ {1/2, 1, 2, 3, 4}, by EM (Dempster, Laird and Rubin 1977): with responsibilities r_j(c) ∝ w_j P(c | m_j λ), the M-step is w_j = Σ_c h(c) r_j(c) / Σ_c h(c) and

λ = Σ_c c h(c) / Σ_c h(c) Σ_j m_j r_j(c),

which maximises Σ_c h(c) Σ_j r_j(c) log P(c | m_j λ) in λ (its derivative is Σ h r_j (c/λ − m_j)). This is a stated simplification of GenomeScope (Vurture et al. 2017), which fits negative binomials whose weights are functions of the heterozygosity and the duplication levels.

**Proposition 1.5 (heterozygosity).** Let heterozygous sites be independent with rate r per base. Then a window of k bases is heterozygous with probability f = 1 − (1 − r)^k. If W_j is the number of distinct k-mers in the component at m_j λ, then f is estimated by

f̂ = (W_{1/2}/2) / (W_{1/2}/2 + Σ_{j ≥ 1} j W_j), and r̂ = 1 − (1 − f̂)^{1/k}.

*Proof.* A heterozygous window gives two distinct k-mers at λ/2; a homozygous window of a j-copy repeat shares one k-mer at jλ with j − 1 others. Counting windows, W_{1/2}/2 are heterozygous and Σ j W_j homozygous. ∎

*Caveat.* A spectrum with peaks at x and 2x is a diploid genome with heterozygous k-mers at x or a haploid genome with duplications at 2x; the counts alone cannot tell them apart. The fit is made from λ₀ equal to the peak and to twice the peak, over one range, and the more likely kept, the first on a tie; a caller who knows the coverage gives it.

**Proposition 1.6 (Lander and Waterman 1988).** In the model of Proposition 1.2, let two reads that overlap by at least T = θL bases be joined. The expected number of islands (maximal joined runs of reads) is N e^{−cσ}, σ = 1 − θ.

*Proof.* A read is the last of its island when no read starts in the Lσ bases after its start; the number of starts there is Poisson of mean (N/G)Lσ = cσ, which is 0 with probability e^{−cσ}. Each island has exactly one last read. ∎

On the node-centric de Bruijn graph two error-free reads are joined when they overlap by k − 1 bases or more, since the last k-mer of one and the first of the other are then joined by an arc (§2), so T = k − 1. The checks count the contigs of error-free reads at 2×, 3× and 4× against N e^{−cσ} and against the exact expectation 1 + (N − 1)(1 − g/G′)^N of uniform starts on G′ positions with g = L − k + 1 (assembly/tests/check_assemble.cpp).

## 2. De Bruijn Graphs

### 2.1 Words and Line Graphs

**Definition 2.1.** The de Bruijn graph B(σ, n) has the words of length n over σ letters as vertices and an edge from x₁ … x_n to x₂ … x_{n+1} for each word x₁ … x_{n+1}.

**Proposition 2.2.** B(σ, n) is the line digraph of B(σ, n − 1): its vertices are the edges of B(σ, n − 1), and two are joined when the head of the first is the tail of the second. Every vertex has in- and out-degree σ, and A^n = J, the matrix of ones, for the adjacency matrix A.

*Proof.* An edge of B(σ, n − 1) is a word of length n, from its prefix to its suffix; two such words follow each other exactly when the suffix of the first is the prefix of the second, which is the edge relation of B(σ, n). A walk of length n from x to y spells x followed by y, and there is exactly one. ∎

The de Bruijn graph of a composition C_k(s) is the sub-multigraph of B(4, k − 1) with an edge, from prefix to suffix, for each element of C_k(s), counted with multiplicity (Pevzner 1989; Idury and Waterman 1995).

### 2.2 The Graph of a Double Strand

Reads come from both strands, so a k-mer and its reverse complement are one node.

**Definition 2.3 (bidirected de Bruijn graph).** Let K be a set of canonical k-mers. An *oriented node* is a pair (x, ε), x ∈ K, ε ∈ {+, −}, spelling x when ε = + and x̄ when ε = −; its *reverse* is (x, ε)⁻ = (x, −ε). There is an arc o → o′ when the last k − 1 letters spelled by o are the first k − 1 letters spelled by o′ (Medvedev, Georgiou, Myers and Brudno 2007; Compeau, Pevzner and Tesler 2011).

**Proposition 2.4.** (i) o → o′ is an arc if and only if o′⁻ → o⁻ is. (ii) If k is odd, no k-mer is its own reverse complement, so o ≠ o⁻ for every oriented node. (iii) Each oriented node has at most four successors.

*Proof.* (i) The reverse complement of the overlap of o and o′ is the overlap of o′⁻ and o⁻. (ii) If x = x̄ then x_{(k+1)/2} is its own complement, which no letter is. (iii) The successor is determined by the next letter. ∎

So the graph on the 2|K| oriented nodes with the involution o ↦ o⁻, which reverses arcs, is a *skew-symmetric graph* (Goldberg and Karzanov 1996); a pair of arcs {o → o′, o′⁻ → o⁻} is one *bidirected edge*, joining a side of each node. For even k the palindromes (ACGT is one) are their own reverse, and a walk can turn back on itself through them; odd k avoids this.

*Remark.* The graph is *node-centric*: every (k − 1)-overlap between k-mers of K is an arc, whether or not a read spans it. An edge-centric graph keeps only the overlaps seen in reads (the (k + 1)-mers); the node-centric one is a function of K alone, which is what makes the canonical form of §2.3 possible.

### 2.3 Unitigs

**Definition 2.5.** An arc o → o′ *extends* when it is the only arc out of o, the only arc into o′, and o and o′ are different nodes. A *unitig* is a maximal path o₁ → … → o_m of extending arcs, or a cycle of them (Chikhi, Limasset and Medvedev 2016, for the definition). It spells a word of m + k − 1 letters.

**Proposition 2.6.** (i) The extending arcs form disjoint paths and cycles on the oriented nodes, and the involution maps each to the reverse of another. (ii) The unitigs partition the nodes. (iii) An arc between two unitigs leaves the last node of one and enters the first node of the other (in suitable orientations).

*Proof.* (i) Each oriented node has at most one extending arc out and one in, so the extending arcs form paths and cycles; Proposition 2.4 (i) maps an extending arc o → o′ to the extending arc o′⁻ → o⁻. No path is its own reverse: its middle arc would join some o to o⁻, the same node. (ii) Each node lies on one path or cycle P and on its reverse P⁻, which are the same unitig read in the two directions. (iii) Let o → q be an arc from the last node o of a unitig to a node q of a unitig W. If q is not the first node of W in that orientation, q has a predecessor o′ in W along an extending arc, which is then the only arc into q; so o′ = o, the arc o → q extends, and o is not the last node of its unitig, a contradiction. ∎

The *compacted graph* has the unitigs as nodes and these arcs as its bidirected edges; a walk in the de Bruijn graph that begins and ends at unitig ends is a walk in the compacted graph, spelling the same word.

**Proposition 2.7 (canonical form).** Spell each unitig from the end whose first k-mer is the smaller integer, and number the unitigs by their least canonical k-mer. Then the compacted graph is a function of K. In particular, removing a set R of unitigs and compacting again gives the compacted graph of K minus the k-mers of R.

*Proof.* The de Bruijn graph is a function of K (node-centric), its unitigs a function of the graph, and the choices of spelling and numbering are functions of the unitigs. Removing whole unitigs removes their k-mers; the graph of the rest is the induced subgraph, whose extending arcs are found again by compaction, the maximal non-branching paths of unitigs being exactly the maximal non-branching paths of their nodes. ∎

**The algorithm.** Each node's two sides find their successors by four hash lookups each, in parallel over the nodes. The extending arcs then form paths, walked from every *head* (an oriented node with no extending arc into it) in parallel and kept from the smaller of their two heads, and cycles, walked from their least node. Then each unitig is spelled and the arcs between ends found. All of it is O(|K|) expected time, and the result is the same on any number of threads (assembly/debruijn.hpp). A unitig's coverage is the mean count of its k-mers.

## 3. Cleaning

Errors and heterozygous sites leave structures in the graph that a genome does not have. Where they come from is a statement about the k-mers that cover one base.

**Proposition 3.1 (the structures of one difference).** Let a read differ from the genome at one base, at position j of a read of length L (0-based), and let all its other k-mers be the genome's. The k-mers that cover the base are k-mers max(0, j − k + 1) to min(j, L − k): (i) if k − 1 ≤ j ≤ L − k, all k of them, and they form a path of k k-mers that leaves the genome's path and rejoins it: a *bubble*, whose arm spells 2k − 1 letters; (ii) otherwise fewer, and they form a path that leaves the genome's path and stops at the read's end: a *tip* of fewer than k k-mers, fewer than 2k − 1 letters.

*Proof.* The k-mers of the read that cover the base are consecutive and pairwise overlapping by k − 1; the one before the first and the one after the last are the genome's, and the first and last of the run overlap them. When the run reaches an end of the read, there is nothing after it. ∎

The same holds for a heterozygous SNV, with the read replaced by the haplotype that carries it: an arm of 2k − 1 letters on each haplotype, both covered λ/2.

**Proposition 3.2 (nearby sites).** Let two heterozygous SNVs lie d bases apart on a diploid genome, with no others near. (i) If d < k, they make one bubble, whose two arms are the haplotypes' and spell 2k − 1 + d letters. (ii) If d > k, they make two bubbles joined by the d − k shared k-mers between them. (iii) If d = k, they share no k-mer and one (k − 1)-mer: the two arms of the first site each join both arms of the second, a *crossed bubble*.

*Proof.* A k-mer starting at x covers neither site when x > p and x + k − 1 < p + d, which has d − k solutions; for d < k every k-mer between them covers a site, and since reads show only the two haplotypes, the k-mers covering both sites come in two versions, one per haplotype. For d = k the k-mer ending at p + k − 1 covers only the first site and the one beginning at p + 1 only the second, and they overlap in the k − 1 bases between the sites, which hold neither. ∎

The cleaning of clean.hpp is a set of heuristics built on these facts, in the tradition of Velvet (Zerbino and Birney 2008):

- *tips*: a unitig with one dead side, shorter than 2k, less covered than the best unitig leaving the same branching side, and *islands*, short unitigs with no neighbour at all;
- *bubbles*: arms that are single unitigs between the same two sides, of nearly equal length (Velvet's Tour Bus compares paths found by Dijkstra's algorithm; the simplified form pops single-unitig arms, and nested bubbles become simple as the rounds compact), and the crossed bubble of Proposition 3.2 (iii); the best covered arm is kept;
- *low coverage*: a short unitig covered less than a fifth of its neighbours.

Each step finds its removals on the graph as it stands, removes them all at once, and compacts; rounds repeat until nothing is removed. By Proposition 2.7 the result is the compacted graph of the k-mers that remain, which the checks verify. Popping a heterozygous bubble keeps one haplotype's arm: the contig that results is a *mosaic* of the haplotypes, not either of them.

## 4. Eulerian Paths

### 4.1 Strings with a Given Composition

**Proposition 4.1.** Let C = C_k(s) and let G be its de Bruijn graph, with m(x) parallel edges for each k-mer x of multiplicity m(x). Then the words t with C_k(t) = C are the spellings of the Eulerian trails of G, and

#{t : C_k(t) = C} = trails(G) / ∏_x m(x)!,

where trails(G) is the number of Eulerian trails counted as sequences of edges.

*Proof.* A word of length |s| with composition C visits, through its consecutive k-mers, every edge once along a walk: an Eulerian trail; conversely a trail spells a word with composition C. Parallel edges of G carry the same k-mer, since a prefix and a suffix of k − 1 letters that overlap in k − 2 determine it. The group ∏_x S_{m(x)}, permuting parallel copies, acts freely on trails (a trail uses each edge once, so a non-trivial permutation changes it), and its orbits are the fibres of the spelling map. ∎

### 4.2 Euler and Hierholzer

**Theorem 4.2 (Euler 1736; Hierholzer 1873; Good 1946 for digraphs).** A directed multigraph whose edges are weakly connected has an Eulerian circuit if and only if in(v) = out(v) at every vertex, and an Eulerian trail from s ≠ t if and only if out(s) = in(s) + 1, in(t) = out(t) + 1 and every other vertex is balanced.

*Proof (Hierholzer).* Necessity: a closed walk enters a vertex as often as it leaves it. Sufficiency, for the balanced case: walk from a vertex v along unused edges; at a vertex u ≠ v the walk has used one more edge in than out, so an unused edge out remains, and the walk can only stop at v. If unused edges remain, connectedness gives a vertex on the walk with an unused edge, and a closed walk from it on the unused edges (the graph of unused edges is still balanced) is spliced in. The trail case follows by adding an edge t → s. ∎

The iterative form keeps the current walk on a stack, emits edges as it backs up from vertices with no unused edge, and reverses the emission: O(n + m). Taking at each vertex the unused edge of least index makes the result a function of the graph (morphe/combinatorics/eulerian.hpp). It reconstructs s whenever s is the only word with its composition.

### 4.3 The Matrix-Tree Theorem

**Definition 4.3.** An *arborescence toward w* (an in-tree) is a set of edges containing, for each vertex v ≠ w, one edge out of v, such that following them from any vertex reaches w. The *out-degree Laplacian* is L = D_out − A, A_uv the number of edges u → v (loops cancel), and L_w is L with the row and column of w deleted.

**Theorem 4.4 (Tutte 1948; Borchardt 1860, Kirchhoff 1847 for undirected graphs).** The number of arborescences toward w is t_w(G) = det L_w.

*Proof (Chaiken 1982; Zeilberger 1985).* Expand det L_w = Σ_π sgn π ∏_v L_{v,π(v)} over permutations of V′ = V ∖ {w}, writing each diagonal entry as the sum, over the edges out of v that are not loops, of 1, and each off-diagonal entry as minus the sum over edges v → π(v). A term is a choice f of one edge out of each v ∈ V′, together with a set S of disjoint cycles of the functional graph of f (the non-trivial cycles of π, realised by edges); a cycle of length ℓ contributes sgn (−1)^{ℓ − 1} times (−1)^ℓ, that is −1. So

det L_w = Σ_f Σ_{S ⊆ cycles(f)} (−1)^{|S|},

and the inner sum is 0 unless f has no cycle. A function with no cycle, each vertex of V′ sent along one edge, reaches w from every vertex: it is an arborescence toward w. ∎

The determinant is computed exactly by Bareiss's fraction-free elimination (Bareiss 1968), whose intermediate entries are minors of the matrix by Sylvester's identity, so they are exact integers bounded by Hadamard's inequality (morphe/algebra/determinant.hpp).

### 4.4 The BEST Theorem

**Theorem 4.5 (van Aardenne-Ehrenfest and de Bruijn 1951; Smith and Tutte 1941).** For a connected Eulerian directed multigraph and any vertex w, the number of Eulerian circuits, as cyclic sequences of edges, is

ec(G) = t_w(G) ∏_v (out(v) − 1)!.

*Proof.* Fix an edge e out of w. Each cyclic circuit has exactly one rotation beginning with e, so ec(G) is the number of circuits beginning with e. For such a circuit and each v ≠ w let ℓ(v) be the last edge out of v. Following ℓ from v leads to a vertex v′ that is left again later than v was left for the last time, unless v′ = w, so the times of last exit increase along the path and it reaches w: the ℓ(v) form an arborescence T toward w. The circuit also gives, at each v ≠ w, the order in which its other out(v) − 1 edges are taken, and at w the order of its out(w) − 1 edges after e. Conversely, given T and these orders, walk from e, taking at each vertex the next unused edge in its order and the edge of T last. The walk ends at w, which it can leave no more; if an edge were unused, its tail v ≠ w would have its edge of T unused, so the head of that edge would have an unused edge out, and so on along T to w, which has none: a contradiction. So the walk is an Eulerian circuit beginning with e, whose last exits are T and whose orders are those given, and the two constructions are inverse. Counting the orders gives t_w ∏_{v ≠ w} (out(v) − 1)! (out(w) − 1)!. ∎

*Remark.* A vertex with one edge in and one out contributes 0! = 1 and no choice. morphe suppresses such vertices before the determinant, joining their two edges, so that a genome's graph of millions of vertices and a few hundred branchings costs a determinant of order a few hundred.

The number of Eulerian trails, as edge sequences, is m · ec(G) for an Eulerian graph (each circuit cut at each of its m edges, every rotation a different trail) and ec(G + (t → s)) for a graph with a trail from s to t (each circuit of the augmented graph passes the added edge once). With Proposition 4.1 this counts the words with a given composition exactly; the checks compare it with all 4⁹ words of length 9 and with the spellings of every trail of the textbook example TAATGCCATGGGATGTT, whose 3-mers spell two words, it and TAATGGGATGCCATGTT (Compeau and Pevzner 2015, ch. 3).

### 4.5 De Bruijn Sequences

**Corollary 4.6 (de Bruijn 1946; Flye Sainte-Marie 1894 for σ = 2).** The number of cyclic words of length σⁿ over σ letters in which every word of length n occurs exactly once is

(σ!)^{σ^{n−1}} / σⁿ, which is 2^{2^{n−1} − n} for σ = 2.

*Proof.* Such a word is an Eulerian circuit of B(σ, n − 1), read by the last letter of each edge. B(σ, m), m = n − 1, has σ^m vertices, each of out-degree σ. By Proposition 2.2, A^m = J has rank one, so the eigenvalues of A are σ (on the vector of ones) and 0 with multiplicity σ^m − 1, and those of L = σI − A are 0 once and σ with multiplicity σ^m − 1. For a connected Eulerian digraph every t_w is equal, and their sum is the product of the non-zero eigenvalues of L (the coefficient of the linear term of its characteristic polynomial; Stanley 1999, §5.6), so t_w = σ^{σ^m − 1}/σ^m. Then ec = σ^{σ^m − 1 − m} ((σ − 1)!)^{σ^m} = (σ!)^{σ^m}/σ^{m+1}. ∎

The checks compute these counts by Bareiss's determinant on the Laplacians of B(σ, n − 1) and compare them with the closed form and with enumeration: 1, 1, 2, 16, 2048 binary sequences of orders 1 to 5; 24 and 20 736 of order 2 over three and four letters; 189 321 481 108 517 289 984 DNA sequences of order 3.

## 5. Repeats and the Ambiguity of Assembly

**Proposition 5.1.** If s has no repeated (k − 1)-mer, then s is the only word with its k-mer composition.

*Proof.* Every vertex of the de Bruijn graph of C_k(s) is a (k − 1)-mer of s, occurring once, so the graph is a simple path and has one Eulerian trail. ∎

A repeat changes this only when a path can leave it by another exit. With s = A x B x C for one (k − 1)-mer x occurring twice, the loop x B x must be taken when the path first reaches x, or B is never used: one word. With two repeats whose occurrences interleave, s = A x B y C x D y E, the words A x D y C x B y E and s have one composition; so does a repeat of three occurrences, A x B x C x D and A x C x B x D; and a word that begins and ends with one (k − 1)-mer, x B x, shares its composition with its rotations. Ukkonen (1992) conjectured, and Pevzner (1995) proved, that two words share a composition exactly when one is obtained from the other by such transpositions and rotations, so that a word is determined by its composition when it has none of these.

The number of words with a composition, Proposition 4.1 and Theorem 4.5, is the measure of the ambiguity. It is a product over the branching vertices, (out(v) − 1)! each, times the arborescences among them: it grows factorially with the copies of a repeat and exponentially with the number of interleaved repeats. In a genome every repeat longer than k − 1 is such a vertex, and a contig must end there, since the unitigs are what every word with the composition shares: a unitig is a path that every Eulerian trail traverses whole. This is why the N50 of a short-read assembly is bounded by the spacing of the genome's repeats longer than k (assembly/tests/check_assemble.cpp), whatever the coverage, and why longer reads, read pairs or a larger k are needed to go further.

## 6. Approximate Membership

Two structures of kmers/filters.hpp hold k-mers in fixed memory at a stated error.

**Proposition 6.1 (Bloom 1970).** After n insertions into a Bloom filter of m bits with h hash positions per key, independent and uniform, a bit is 0 with probability (1 − 1/m)^{hn} ≈ e^{−hn/m}, and a key not inserted is reported present with probability about (1 − e^{−hn/m})^h, least at h = (m/n) ln 2.

The positions are drawn by double hashing, a(x) + i b(x) mod m, which has the same asymptotic rate (Kirsch and Mitzenmacher 2006); the checks find the measured rate within four binomial deviations of the formula for h = 1 to 12.

**Proposition 6.2 (Cormode and Muthukrishnan 2005, Theorem 1).** A count-min sketch of d = ⌈ln(1/δ)⌉ rows of w = ⌈e/ε⌉ counters, each row hashed by a pairwise independent function, estimates the count a(x) by â(x) = min_j over the rows, with a(x) ≤ â(x) and P[â(x) > a(x) + ε‖a‖₁] ≤ δ.

*Proof.* Each row's excess has expectation at most ‖a‖₁/w by pairwise independence, so by Markov's inequality it exceeds ε‖a‖₁ = e‖a‖₁/w with probability at most 1/e; the rows are independent. ∎

## 7. Where the Code Implements It

| Result | Code |
|---|---|
| Counting canonical k-mers by partitions (§1) | `kmers/counter.hpp`: `KmerCounter` |
| The spectrum, its threshold, the mixture, genome size and heterozygosity (Remark 1.3, Propositions 1.4, 1.5) | `kmers/spectrum.hpp`: `spectrum`, `fit_spectrum` |
| Lander–Waterman islands (Proposition 1.6) | `assembly/assemble.hpp`: `lander_waterman_contigs` |
| The bidirected graph, unitigs and the canonical form (§2.2–2.3) | `assembly/debruijn.hpp`: `build_graph`, `compact` |
| Tips, bubbles, crossed bubbles, low coverage (§3) | `assembly/clean.hpp`: `clean` |
| Strings with a composition (Proposition 4.1) | `assembly/eulerian.hpp`: `count_strings`, `reconstruct` |
| Hierholzer, matrix-tree, BEST (§4.2–4.4) | `morphe/combinatorics/eulerian.hpp`; the determinant in `morphe/algebra/determinant.hpp` |
| De Bruijn sequences (Corollary 4.6) | `assembly/eulerian.hpp`: `de_bruijn_sequences`, `de_bruijn_sequence` |
| Bloom filters and count-min sketches (§6) | `kmers/filters.hpp` |

## 8. References

- van Aardenne-Ehrenfest, T. and de Bruijn, N. G. (1951), Circuits and trees in oriented linear graphs, *Simon Stevin* 28, 203–217.
- Bareiss, E. H. (1968), Sylvester's identity and multistep integer-preserving Gaussian elimination, *Math. Comp.* 22, 565–578.
- Bloom, B. H. (1970), Space/time trade-offs in hash coding with allowable errors, *Comm. ACM* 13, 422–426.
- de Bruijn, N. G. (1946), A combinatorial problem, *Proc. Ned. Akad. Wet.* 49, 758–764.
- Chaiken, S. (1982), A combinatorial proof of the all minors matrix tree theorem, *SIAM J. Algebraic Discrete Methods* 3, 319–329.
- Chikhi, R., Limasset, A. and Medvedev, P. (2016), Compacting de Bruijn graphs from sequencing data quickly and in low memory, *Bioinformatics* 32, i201–i208.
- Compeau, P. E. C., Pevzner, P. A. and Tesler, G. (2011), How to apply de Bruijn graphs to genome assembly, *Nature Biotechnology* 29, 987–991.
- Compeau, P. and Pevzner, P. (2015), *Bioinformatics Algorithms*, Active Learning Publishers, ch. 3.
- Cormode, G. and Muthukrishnan, S. (2005), An improved data stream summary: the count-min sketch and its applications, *J. Algorithms* 55, 58–75.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), Maximum likelihood from incomplete data via the EM algorithm, *J. R. Stat. Soc. B* 39, 1–38.
- Goldberg, A. V. and Karzanov, A. V. (1996), Path problems in skew-symmetric graphs, *Combinatorica* 16, 353–382.
- Good, I. J. (1946), Normal recurring decimals, *J. London Math. Soc.* 21, 167–169.
- Hierholzer, C. (1873), Über die Möglichkeit, einen Linienzug ohne Wiederholung und ohne Unterbrechung zu umfahren, *Math. Ann.* 6, 30–32.
- Idury, R. M. and Waterman, M. S. (1995), A new algorithm for DNA sequence assembly, *J. Comput. Biol.* 2, 291–306.
- Kingman, J. F. C. (1993), *Poisson Processes*, Oxford University Press.
- Kirsch, A. and Mitzenmacher, M. (2006), Less hashing, same performance: building a better Bloom filter, *ESA*, LNCS 4168, 456–467.
- Lander, E. S. and Waterman, M. S. (1988), Genomic mapping by fingerprinting random clones: a mathematical analysis, *Genomics* 2, 231–239.
- Li, X. and Waterman, M. S. (2003), Estimating the repeat structure and length of DNA sequences using ℓ-tuples, *Genome Research* 13, 1916–1922.
- Medvedev, P., Georgiou, K., Myers, G. and Brudno, M. (2007), Computability of models for sequence assembly, *WABI*, LNBI 4645, 289–301.
- Pevzner, P. A. (1989), l-tuple DNA sequencing: computer analysis, *J. Biomol. Struct. Dyn.* 7, 63–73.
- Pevzner, P. A. (1995), DNA physical mapping and alternating Eulerian cycles in colored graphs, *Algorithmica* 13, 77–105.
- Smith, C. A. B. and Tutte, W. T. (1941), On unicursal paths in a network of degree 4, *Amer. Math. Monthly* 48, 233–237.
- Stanley, R. P. (1999), *Enumerative Combinatorics* 2, Cambridge University Press, §5.6.
- Tutte, W. T. (1948), The dissection of equilateral triangles into equilateral triangles, *Proc. Cambridge Philos. Soc.* 44, 463–482.
- Ukkonen, E. (1992), Approximate string-matching with q-grams and maximal matches, *Theoret. Comput. Sci.* 92, 191–211.
- Vurture, G. W. et al. (2017), GenomeScope: fast reference-free genome profiling from short reads, *Bioinformatics* 33, 2202–2204.
- Zeilberger, D. (1985), A combinatorial approach to matrix algebra, *Discrete Math.* 56, 61–72.
- Zerbino, D. R. and Birney, E. (2008), Velvet: algorithms for de novo short read assembly using de Bruijn graphs, *Genome Research* 18, 821–829.
