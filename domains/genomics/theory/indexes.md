# Indexes of Words and Genomes

*Status: draft. Serves learning and use.* This note covers the combinatorics of words that the indexes rest on, and the indexes themselves: suffix arrays, the Burrows–Wheeler transform read as a permutation, backward search, the FMD-index of a genome and its reverse complement, minimizers, and sketches read as estimators. The code is in two places:

- **Words and their indexes** name no object in nature, so they are morphe's `formal` module: [`words.hpp`](../../../foundations/morphe/formal/include/morphe/formal/words.hpp), [`suffix_array.hpp`](../../../foundations/morphe/formal/include/morphe/formal/suffix_array.hpp), [`bwt.hpp`](../../../foundations/morphe/formal/include/morphe/formal/bwt.hpp), [`succinct.hpp`](../../../foundations/morphe/formal/include/morphe/formal/succinct.hpp), [`fm_index.hpp`](../../../foundations/morphe/formal/include/morphe/formal/fm_index.hpp), [`aho_corasick.hpp`](../../../foundations/morphe/formal/include/morphe/formal/aho_corasick.hpp).
- **What knows the double strand** is the genomics domain's: [`index/fmd.hpp`](../index/include/genomics/index/fmd.hpp), [`index/minimizer_index.hpp`](../index/include/genomics/index/minimizer_index.hpp), [`sketch/minhash.hpp`](../sketch/include/genomics/sketch/minhash.hpp), [`sketch/hyperloglog.hpp`](../sketch/include/genomics/sketch/hyperloglog.hpp), [`sketch/distance.hpp`](../sketch/include/genomics/sketch/distance.hpp).

The contract is the design note [The Genome Infrastructure](../../../development/ideas/2026-10-24-the-genome-infrastructure.md), §§3.1 and 4.

Throughout, A is a finite totally ordered alphabet, here [0, σ) with the order of the integers. A* is the set of finite words over A, ε the empty word, |w| the length of w, and w[i] its letters from 0. Words are ordered lexicographically, a proper prefix being smaller.

## 1. Words

### 1.1 The Free Monoid and Its Morphisms

**Definition 1.1.** A* under concatenation is a monoid with unit ε. A *morphism* φ: A* → B* is a map with φ(uv) = φ(u)φ(v) and φ(ε) = ε.

**Proposition 1.2 (freeness).** A* is the free monoid on A: every map f: A → M into a monoid extends uniquely to a morphism A* → M, namely a₁ ⋯ a_n ↦ f(a₁) ⋯ f(a_n).

*Proof.* A morphism is determined by the images of the letters, since every word is a product of letters in exactly one way; and the formula is a morphism because concatenation is associative. ∎

So a substitution is given by one word for each letter, and `Morphism` holds exactly that. The checks verify φ(uv) = φ(u)φ(v) on random morphisms and words. Composition is applicative, (g ∘ f)(w) = g(f(w)), as morphe's glossary fixes it.

A morphism is *prolongable* on a when φ(a) = au with u ≠ ε. Then φⁿ(a) is a prefix of φⁿ⁺¹(a), and the words converge to an infinite fixed point φ^ω(a). Three fixed points recur in the checks:

- **The Thue–Morse word** t = μ^ω(0), μ: 0 ↦ 01, 1 ↦ 10. Its letter t_i is the parity of the number of ones in the binary expansion of i (Allouche and Shallit 2003, §1.6). It is *overlap-free*: it has no factor axaxa with a a letter (Thue 1912).
- **The Fibonacci word** f = φ^ω(0), φ: 0 ↦ 01, 1 ↦ 0. Its finite words f_m = φ^m(0) satisfy f_{m+1} = f_m f_{m−1}, with |f_m| = F_{m+2}.
- **The period-doubling word** d = δ^ω(0), δ: 0 ↦ 01, 1 ↦ 00, with d_i = ν₂(i + 1) mod 2.

### 1.2 Borders and Periods

**Definition 1.3.** A *border* of w is a word u ≠ w that is both a prefix and a suffix of w. An integer p ≥ 1 is a *period* of w when w[i] = w[i + p] for all 0 ≤ i < |w| − p.

**Proposition 1.4.** p is a period of w exactly when w has a border of length |w| − p.

*Proof.* Both say w[0, |w| − p) = w[p, |w|). ∎

**Proposition 1.5 (the failure function).** Let β(w) be the longest border of w. The borders of w are β(w), β(β(w)), …, ε.

*Proof.* A border of a border of w is a border of w. Conversely, a border u of w shorter than β(w) is a prefix of β(w), since both are prefixes of w, and likewise a suffix of β(w); so u is a border of β(w). Induction on length finishes the proof. ∎

The *border array* b[i] = |β(w[0..i])| is therefore computed in one left-to-right pass. To extend the border of w[0..i−1] by w[i], try it, then its longest border, and so on. Each fall shortens the current border, which grows by at most one per letter, so the pass makes fewer than 2n comparisons. This is the failure function of Knuth, Morris and Pratt (1977); run against a text instead of the word itself, it is KMP search (`border_array`, `find_all`).

The *Z-array*, Z[i] = lcp(w, w[i..]), carries the same information read from the left. Gusfield's algorithm keeps the rightmost interval [l, r) known to match a prefix, so that Z[i] starts from min(r − i, Z[i − l]) and only letters beyond r are compared (Gusfield 1997, §1.4). It also runs in O(n).

**Theorem 1.6 (Fine and Wilf 1965).** If w has periods p and q and |w| ≥ p + q − gcd(p, q), then gcd(p, q) is a period of w. The bound is sharp: for coprime p, q there are words of length p + q − 2 with periods p and q but not 1.

*Proof.* See Lothaire (1983), §1.3. The idea is Euclid's: if p < q, the factor w[q − p, |w|) has the periods p and q − p, and one recurses on (p, q − p). ∎

The checks verify the theorem on every binary word of length at most 14, and find a word attaining the bound for every coprime pair that fits. The extremal words are the central words of Sturmian theory, such as the prefixes of the Fibonacci word.

### 1.3 Primitive Words and Roots

**Definition 1.7.** A non-empty word w is *primitive* when w = u^e implies e = 1.

**Proposition 1.8 (Lyndon and Schützenberger 1962; Lothaire 1983, §1.3).**

1. Every non-empty word is u^e for a unique primitive u, its *primitive root*.
2. uv = vu if and only if u and v are powers of one word.
3. w is primitive if and only if w occurs in ww only at the positions 0 and |w|.

From (1) and Theorem 1.6, the primitive root of w has length p, the smallest period of w, when p divides |w|, and |w| otherwise. Indeed, if p | |w| then w = (w[0, p))^{|w|/p}, and a shorter root would give a shorter period. If p ∤ |w| and w = u^e with e ≥ 2, then |u| is a period of length at most |w|/2. So p + |u| ≤ |w|, Fine and Wilf make gcd(p, |u|) a period, and minimality gives p | |u| | |w|, a contradiction. This is `primitive_root_length`.

### 1.4 Lyndon Words and the Factorisation Theorem

**Definition 1.9.** A *Lyndon word* is a non-empty word strictly smaller than all of its proper rotations. Equivalently, it is primitive and the least in its conjugacy class (Lothaire 1983, §5.1).

**Theorem 1.10 (Chen, Fox and Lyndon 1958).** Every word w factors uniquely as w = ℓ₁ ℓ₂ ⋯ ℓ_m with each ℓ_i a Lyndon word and ℓ₁ ≥ ℓ₂ ≥ ⋯ ≥ ℓ_m. The first factor ℓ₁ is the longest prefix of w that is a Lyndon word.

*Proof.* Lothaire (1983), §5.1. ∎

*Remark (Duval's algorithm).* Duval (1983) reads the factorisation in one pass. It keeps the invariant that the scanned part w[i, j) has the form u^e u′, with u Lyndon and u′ a proper prefix of u. Let k index the letter of u that w[j] is compared with. Then:

- w[j] = w[k] continues the periodic pattern (k advances);
- w[j] > w[k] makes the whole of w[i, j] a Lyndon word (k returns to i);
- w[j] < w[k] ends the run, and the e copies of u are emitted as factors.

At most 2n comparisons are made, and the extra space is constant (`lyndon_factorisation`). The checks compare it, on every binary word up to length 12, with the greedy factorisation into longest Lyndon prefixes, computed by brute force.

**Proposition 1.11 (Witt's formula).** The number of Lyndon words of length n over k letters is

L_k(n) = (1/n) Σ_{d | n} μ(d) k^{n/d}.

*Proof.* By Proposition 1.8, each of the kⁿ words of length n is a power of a unique primitive word, whose length d divides n. The primitive words of length d fall into conjugacy classes of exactly d elements, since rotations of a primitive word are distinct, and each class holds exactly one Lyndon word. Hence kⁿ = Σ_{d|n} d L_k(d), and Möbius inversion gives the formula. ∎

The checks compare it with Duval's later algorithm, which enumerates the Lyndon words of length at most n in lexicographic order (Duval 1988): repeat the current word periodically to length n, strip trailing maximal letters, and increment the last letter. They also verify the identity kⁿ = Σ_{d|n} d L_k(d) exactly in morphe's integers, up to n = 60 (`for_each_lyndon_word`, `lyndon_count`).

### 1.5 Conjugacy

Two words are *conjugate* when they are xy and yx, that is, rotations of each other. Conjugacy is an equivalence relation, and its classes are the *necklaces*.

**Proposition 1.12.** u and v are conjugate exactly when they have equal length and equal least rotations.

*Proof.* The least rotation is a function of the class: rotations of rotations are rotations. ∎

Booth's algorithm (Booth 1980) finds the least rotation in O(n). It runs a failure function over ww, computed relative to a moving candidate start k, which jumps whenever a smaller letter shows that a later start is better. The code reduces its answer modulo the primitive root's length, since the rotations by r and r + |root| coincide (`least_rotation`, `are_conjugate`).

## 2. Suffix Arrays

### 2.1 Definitions and What They Count

**Definition 2.1.** The *suffix array* of w ∈ Aⁿ is the permutation SA of [0, n) with w[SA[0]..) < ⋯ < w[SA[n − 1]..). Its inverse is the *rank*. The *LCP array* is LCP[i] = lcp(w[SA[i − 1]..), w[SA[i]..)), with LCP[0] = 0.

**Proposition 2.2 (Kasai et al. 2001).** If the suffix at i has lcp h with its predecessor in SA, the suffix at i + 1 has lcp at least h − 1 with its predecessor. So LCP is computed in text order in O(n).

*Proof.* Let j be the predecessor of i, with lcp h ≥ 1. Then j + 1 precedes i + 1 in the order, and they share h − 1 letters. The predecessor of i + 1 lies between them in the order, so it shares at least as many. ∎

**Proposition 2.3 (distinct factors).** The number of distinct non-empty factors of w is n(n + 1)/2 − Σ_i LCP[i]. More finely, the number p(ℓ) of distinct factors of length ℓ is the number of i with LCP[i] < ℓ ≤ n − SA[i].

*Proof.* Each factor is a prefix of some suffix. A prefix of length ℓ of the suffix in row i is new, absent from earlier rows, exactly when ℓ > LCP[i]: the earlier rows share at most LCP[i] letters with it, by the sortedness. ∎

`factor_complexity` computes p(ℓ) for every ℓ by a difference array. Two theorems of the combinatorics of words become checks:

- **The Fibonacci word is Sturmian:** p(ℓ) = ℓ + 1 for every ℓ (Morse and Hedlund 1940), verified for ℓ ≤ 200 on a prefix of 10⁵ letters.
- **The Thue–Morse word** has Brlek's complexity (Brlek 1989; de Luca and Varricchio 1989). For n ≥ 3, write n − 1 = 2^r + q with 0 < q ≤ 2^r. Then p(n) = 6·2^{r−1} + 4q when q ≤ 2^{r−1}, and 8·2^{r−1} + 2q otherwise: 6, 10, 12, 16, 20, 22, 24, 28, ….

### 2.2 Induced Sorting

**Definition 2.4 (types).** A suffix w[i..) is *S-type* when it is smaller than w[i + 1..), and *L-type* otherwise. The last suffix is L-type, being larger than the empty one. Position i is *leftmost-S* (LMS) when it is S-type and i − 1 is L-type.

The type of i follows from the type of i + 1 by comparing w[i] with w[i + 1], since equal letters inherit the type. So one right-to-left pass computes the types.

**Lemma 2.5 (Ko and Aluru 2003).** Among the suffixes that begin with the same letter c, every L-type suffix is smaller than every S-type suffix.

*Proof.* Write an L-type suffix as c^p x… with x < c, or x the end of the word, and an S-type suffix as c^q y… with y > c. If p ≤ q, they differ at position p, where the first has x and the second c or y, both larger. If p > q, they differ at position q, where the first has c and the second y > c. ∎

So each *bucket*, the suffixes beginning with c, has its L-type suffixes first and its S-type suffixes last.

**Proposition 2.6 (induced sorting; Nong, Zhang and Chan 2011, §3).** Suppose the LMS suffixes are placed in their buckets in their correct relative order. A left-to-right scan that, on meeting w[i..), places w[i − 1..) at the next free slot of the L-part of its bucket when i − 1 is L-type, puts every L-type suffix in its place. A right-to-left scan that places w[i − 1..) at the last free slot of the S-part when i − 1 is S-type then puts every S-type suffix in its place.

*Proof.* An L-type suffix w[i − 1..) = c·w[i..) is larger than w[i..). It is therefore placed after w[i..) is met, and among suffixes beginning with c the order of the w[i..) is the order of the c·w[i..). The induction is on the order. The S-type scan is symmetric. ∎

Applied to LMS suffixes placed in an arbitrary order, the same two scans sort the LMS *substrings*: the factors from one LMS position to the next, inclusive (ibid., §3). The algorithm SA-IS then proceeds as follows:

1. Name the LMS substrings by rank.
2. Form the word of names, of length m ≤ ⌊n/2⌋ (LMS positions are never adjacent), and sort it recursively. The order of its suffixes is the order of the LMS suffixes.
3. Induce once more from the LMS suffixes in their true order.

The recurrence T(n) = T(n/2) + O(n + σ) gives O(n) on an integer alphabet. `suffix_array` implements it with the types in bytes, prefix doubling below 40 letters, and the naive sort below 10. It is checked against the naive sort on every binary word of length at most 12 and on random words of alphabets 1 to 256 up to 10⁵ letters. On repetitive words (Fibonacci, Thue–Morse), where the naive sort is quadratic, it is checked against prefix doubling (Manber and Myers 1993). On unary words it is checked against the closed form n − 1, …, 0.

### 2.3 Sets of Words

A set of words w₀, …, w_{m−1} is indexed as w₀ #₀ w₁ #₁ ⋯, with separators #₀ < #₁ < ⋯ below every letter. A comparison of two suffixes stops at or before the first separator, where the two suffixes either differ or meet distinct separators. The order is therefore the order of the suffixes up to their separators, with ties broken by word. Restricted to one word, it is that word's own suffix array (`concatenate`, `generalised_suffix_array`). The FMD-index of §5 is built this way.

## 3. The Burrows–Wheeler Transform as a Permutation

### 3.1 The Transform and the LF Map

**Definition 3.1.** Let $ be a letter below A. The *BWT* of w is the word L of length n + 1 whose i-th letter is the letter cyclically preceding the i-th smallest suffix of w$. Row 0 is the suffix $, so L[0] is the last letter of w.

**Proposition 3.2 (the LF property).** Let C[c] be the number of letters of L smaller than c. Then LF(i) = C[L[i]] + rank_{L[i]}(L, i) is the row of the rotation one step to the right of the rotation in row i.

*Proof.* The rows whose rotations end with c are, after moving that c to the front, exactly the rows beginning with c, which occupy [C[c], C[c + 1]). Moving the same letter to the front of two rotations preserves their relative order. So the k-th row ending with c becomes the k-th row beginning with c. ∎

*Corollary 3.3.* From row 0, LF walks the rotations of w$ to the right, reading the letters of L as w backwards. The transform is inverted in O(n + σ) (`inverse_bwt`).

### 3.2 The Standard Permutation and Necklaces

Forget the sentinel. For any word L of length N, let π_L be its *standard permutation*, the stable sort of the positions of L by letter, so that sorted(L)[i] = L[π_L(i)]. LF is its inverse. The cycles of LF are a partition of the positions. Reading L along a cycle backwards spells a word, defined up to rotation: a necklace.

The *ω-order* compares u^ω with v^ω, the infinite periodic words. By Theorem 1.6, u^ω and v^ω either agree everywhere, when u and v are powers of one word, or differ within their first |u| + |v| − gcd(|u|, |v|) letters. And u^ω < v^ω exactly when uv < vu. The conjugates of a multiset of words are sorted in ω-order by prefix doubling on the cyclic shifts (`sort_conjugates`). After ⌈log₂(2 max |w|)⌉ rounds, the ranks compare at least |u| + |v| letters and decide the order.

**Definition 3.4 (the extended BWT; Mantaci, Restivo, Rosone and Sciortino 2007).** For a multiset M of primitive words, sort every conjugate of every word in ω-order and take the last letter of each: eBWT(M).

Moving a last letter to the front preserves the ω-order among conjugates ending with that letter, because (ax)^ω = a(xa)^ω. So Proposition 3.2 holds for the eBWT, and LF sends each conjugate to its right rotation. Its cycles are the words of M.

**Theorem 3.5 (Gessel and Reutenauer 1993; Crochemore, Désarménien and Perrin 2005).** The map M ↦ eBWT(M) is a bijection from finite multisets of primitive necklaces onto words, preserving the multiset of letters. Its inverse sends L to the necklaces spelt by the cycles of its standard permutation.

*Proof.* Gessel and Reutenauer (1993); Crochemore, Désarménien and Perrin (2005), Theorem 1. That the cycles spell primitive necklaces is part of the theorem; that the eBWT of the necklaces of L is L again is the LF property read backwards. ∎

The checks verify both directions on every binary word of length at most 10 and every ternary word of length at most 6 (`necklaces`, `extended_bwt`). Each necklace is returned as its Lyndon word, the conjugate in the least row of its cycle.

### 3.3 The Bijective BWT

**Theorem 3.6 (Gil and Scott 2012; Kufleitner 2009).** Let w = ℓ₁ ⋯ ℓ_m be the Chen–Fox–Lyndon factorisation. The map w ↦ eBWT{ℓ₁, …, ℓ_m} is a bijection of Aⁿ, for every n.

*Proof.* It is the composite of three bijections:

1. w ↦ its multiset of Lyndon factors, a bijection onto multisets of Lyndon words of total length n, by Theorem 1.10. A multiset sorted non-increasingly is a valid factorisation, and its product is a word.
2. Lyndon words ↔ primitive necklaces, each necklace having exactly one Lyndon word.
3. The Gessel–Reutenauer bijection of Theorem 3.5. ∎

The inverse reads the cycles of LF as Lyndon words, sorts them non-increasingly and concatenates (`bijective_bwt`, `inverse_bijective_bwt`). Unlike the BWT of w$, it needs neither a sentinel nor an index. The checks show every binary word of length at most 12 and every ternary word of length at most 7 reached exactly once.

### 3.4 Sturmian Words, Runs and Move-to-Front

The BWT groups letters that precede equal contexts, so repetitive words have transforms with few *runs*, r. The extreme case is a theorem.

**Theorem 3.7 (Mantaci, Restivo and Sciortino 2003).** Over a < b, the rotation BWT of w is b^p a^q with gcd(p, q) = 1 exactly when w is a conjugate of a standard Sturmian word.

The standard words are s_{−1} = b, s₀ = a, s_{k+1} = s_k^{d_k} s_{k−1}, with d₀ ≥ 0 and d_k ≥ 1 afterwards (Lothaire 2002, §2.2.2). The checks generate every standard word of length at most 12 and verify both directions on every binary word of that length, and the form b^p a^q, r = 2, for the Fibonacci words up to 17 711 letters.

*Remark.* Move-to-front (Bentley, Sleator, Tarjan and Wei 1986) turns runs into runs of zeros. The pipeline BWT, then MTF, then an entropy coder is Burrows and Wheeler's compressor (`move_to_front`, `runs`).

## 4. Backward Search, Succinct Structures and Automata

### 4.1 Backward Search

**Proposition 4.1 (Ferragina and Manzini 2000).** If the rows of the BWT matrix of w$ beginning with Q form [b, e), those beginning with cQ form [C[c] + rank_c(L, b), C[c] + rank_c(L, e)).

*Proof.* A row begins with cQ exactly when it is the LF image of a row beginning with Q whose last letter is c. LF is increasing on the rows ending with c, by Proposition 3.2. ∎

Reading P from right to left gives its interval in |P| steps, each two ranks. The FM-index stores L in a wavelet matrix, so a rank costs O(log σ) and count costs O(|P| log σ). It also stores suffix-array values at the text positions divisible by s, marked in a bit vector. Locating a row follows LF until a marked row, at most s − 1 steps, and adds the steps back (`FmIndex`).

### 4.2 Rank, Select and the Wavelet Matrix

A bit vector answers rank in O(1) with o(n) extra bits (Jacobson 1989). The layout here is Vigna's rank9 (Vigna 2008). For each block of 512 bits it keeps a 64-bit absolute count and seven 9-bit relative counts in two adjacent words, so a rank is two loads and a popcount. Select is a binary search on the counts.

The *wavelet matrix* (Claude, Navarro and Ordóñez 2015) holds a sequence S of values below σ in ℓ = ⌈log₂ σ⌉ bit vectors. Level 0 holds the top bits. Each next level holds the next bits, after a stable partition of the values by the current bit, with z_l zeros first.

**Proposition 4.2.** At level l, the values whose first l bits agree occupy a contiguous range. The ranges are ordered by their prefixes read with the bits reversed.

*Proof.* By induction: a stable partition by bit l keeps each range's internal order and sends it to two ranges, one among the zeros and one among the ones. ∎

So rank(c, i) follows two positions down the levels: the start of c's range, and the image of i. Each moves by rank₀, or by z_l + rank₁, according to the bit of c. Their difference at the bottom is the answer. Select runs the same path upwards with select₀ and select₁ (`WaveletMatrix`).

### 4.3 The Aho–Corasick Automaton

For a finite set P of patterns, the trie of P is a deterministic automaton whose states are the prefixes of patterns. The *failure* f(u) is the longest proper suffix of u that is a state. The transition δ(u, a) is g(u, a) where the trie has the edge, and δ(f(u), a) otherwise.

**Proposition 4.3 (Aho and Corasick 1975).** δ*(root, x) is the longest suffix of x that is a prefix of a pattern. Hence the patterns ending at a position of the text are the suffixes of that state that are patterns, found along the failure chain. The automaton recognises A* P, a regular language, with one transition per letter.

*Proof.* Induction on x: the longest suffix of xa that is a state is ua, for u the longest suffix of x such that ua is a state; the failure chain lists the suffixes of x that are states, longest first. ∎

The code links each state to the nearest state on its failure chain that ends a pattern, so the occurrences at a position are listed in time proportional to their number (`AhoCorasick`).

## 5. The FMD-Index: Bidirectionality from Reverse Complementation

### 5.1 The Involution

Let ā = 3 − a be the complement on {A, C, G, T} = {0, 1, 2, 3}, and θ: A* → A* the reverse complement, θ(a₁ ⋯ a_n) = ā_n ⋯ ā₁.

θ is an *anti-automorphism* of the free monoid, θ(uv) = θ(v)θ(u), and an involution. A double strand is a word up to θ. The text of the FMD-index is

T = c₀ $₀ θ(c₀) $₁ c₁ $₂ θ(c₁) $₃ ⋯,

whose set of sequences is closed under θ. A pattern occurring on the forward strand of a contig at [p, p + m) occurs in θ(c_i) at |c_i| − p − m. So P and θ(P) have equal numbers of occurrences in T.

**Definition 5.1 (Li 2012).** The *bi-interval* of P is (k, l, s): [k, k + s) is the interval of P in the BWT of T, and [l, l + s) that of θ(P).

**Proposition 5.2 (backward extension; Li 2012, Algorithm 2).** Let a be a base. The interval of aP is given by Proposition 4.1. The interval of θ(aP) = θ(P)ā is the sub-interval of [l, l + s) of rows whose rotation continues with ā after θ(P). Ordering the rows of θ(P) by their next symbol, $ < A < C < G < T, which is the order T, G, C, A, $ of a, gives

l_$ = l,  l_T = l + s_$,  l_G = l_T + s_T,  l_C = l_G + s_G,  l_A = l_C + s_C,

where s_b is the size of the interval of bP and s_$ = s − Σ_b s_b.

*Proof.* The rows of θ(P) are sorted by what follows θ(P). They split into consecutive blocks by the next symbol, in the order of the alphabet. The block followed by ā holds the occurrences of θ(P)ā = θ(aP), whose number is that of aP, which is s_a. The rows followed by a sentinel are the occurrences of θ(P) at the end of a sequence, equivalently of P at its start, and they come first. ∎

*Corollary 5.3 (forward extension).* The bi-interval of Pa is the backward extension of θ(P) by ā, with k and l exchanged.

So the involution turns a unidirectional index into a bidirectional one at no extra cost. One index of T serves both directions, where a general text needs two indexes (Lam et al. 2009).

In the code, B is stored as BWA stores it. Each block of 128 symbols occupies a cache line: four 64-bit counts and four words of 2-bit symbols. A sentinel is stored as A, and the 2m rows of sentinels are kept apart and subtracted from Occ(A). Occ of the four bases is four popcounts a word (`occ4`).

### 5.2 Super-Maximal Exact Matches

**Definition 5.4.** For a query Q, an *exact match* is a factor Q[b, e) occurring in T. It is *maximal* (a MEM) when neither Q[b − 1, e) nor Q[b, e + 1) occurs. A MEM is *super-maximal* (an SMEM) when no other MEM contains it.

Let L(b) be the length of the longest match starting at b. The MEMs are exactly the intervals [b, b + L(b)) with L(b − 1) ≤ L(b). This is the definition the checks compute by brute force.

**Proposition 5.5 (Li 2012, Algorithm 5).** For a position x, extend Q[x] forwards and keep the bi-interval each time its size is about to drop. These are the right-maximal matches [x, e_t), shortest e_t first. Then extend all of them backwards in step, dropping one when it can go no further. A dropped match is an SMEM containing x exactly when no match with a later end survived at that step. Every SMEM containing x is found once.

*Proof.* Li (2012), §2.6. ∎

**Proposition 5.6.** Start with x₀ = 0, and let x_{j+1} be the end of the longest match starting at x_j. Then every SMEM of Q contains exactly one x_j, so the loop finds each SMEM once.

*Proof.* Let [b, e) be an SMEM, and x_j ≤ b < x_{j+1}. If x_j = b, it contains x_j. Otherwise x_j < b, and the longest match [x_j, x_{j+1}) from x_j satisfies x_{j+1} > b. If x_{j+1} ≥ e, then [x_j, x_{j+1}) strictly contains [b, e) and is itself contained in some MEM, contradicting super-maximality. So x_{j+1} < e, and the SMEM contains x_{j+1}.

Uniqueness: an SMEM containing x_j and x_{j+1} would start at or before x_j and end after x_{j+1}. Then the longest match from x_j would extend beyond x_{j+1}, a contradiction. ∎

The checks compare the SMEMs of mutated reads on both strands, with N and with min_occ 1 and 2, against the definition by brute force. On Escherichia coli, the longest SMEM of each of 4000 simulated reads lies on its true diagonal.

### 5.3 Bases That Stand for N

A pattern containing N occurs nowhere. A genome's runs of N are replaced, before θ is taken, by bases drawn from a stream keyed by the seed, the purpose FMD_N and the contig, as BWA's index replaces its holes. The runs are kept, and a hit can be tested against them. Taking θ after the replacement is what keeps T closed under θ, and so keeps Proposition 5.2 true.

## 6. Minimizers

A (w, k)-minimizer is a k-mer of least hash among w consecutive k-mers (Roberts et al. 2004; Schleimer, Wilkerson and Aiken 2003). The definitions, with ties and palindromes, are those of `genomics::sequence::minimizers`.

**Proposition 6.1 (density).** If the k-mers are distinct and their order is uniformly random, the expected fraction of positions selected is 2/(w + 1).

*Proof.* A new minimizer is selected as the window advances from [i, i + w) to [i + 1, i + w + 1) exactly when the least of these w + 1 k-mers is the first or the last one. By symmetry, this has probability 2/(w + 1). ∎

**Proposition 6.2 (chunks with margins).** Cut a contig into chunks [b, e) and compute the minimizers of each extended chunk [b − w + 1, e + w + k − 2). Keeping those that start in [b, e) gives exactly the contig's minimizers in [b, e).

*Proof.* A window containing a k-mer starting at q ∈ [b, e) consists of k-mers starting in [q − w + 1, q + w − 1), whose bases lie in the extended chunk. So every such window is a window of the extended chunk, with the same N-free runs.

It remains to consider the runs that the margins cut. A run piece at the left edge with fewer than w k-mers has all its k-mers starting before b. A piece at the right edge with fewer than w k-mers has all of them starting at or after e. So the margins add no minimizer inside [b, e). ∎

This is what makes the parallel construction exact and independent of the number of threads (`MinimizerIndex`). The chunks have a fixed length, and the merge is a stable sort into buckets by hash, at offsets counted in chunk order.

## 7. Sketches as Estimators

Let h be a hash of 64-bit words, modelled as a uniformly random injective function (`kmer_hash`, a bijection for each seed). Let A and B be the sets of canonical k-mers of two genomes.

### 7.1 Bottom-s MinHash

**Proposition 7.1 (Broder 1997).** P(argmin_{A∪B} h ∈ A ∩ B) = J(A, B) = |A ∩ B|/|A ∪ B|.

*Proof.* The minimum of h over A ∪ B is at a uniformly random element. ∎

Let S(X) be the s smallest hashes of X. Then the s smallest of S(A) ∪ S(B) are S(A ∪ B): an element among the s smallest of A ∪ B lies in A or in B and is among the s smallest there. So S is a morphism from (sets, ∪) to (sketches, merge-and-truncate).

S(A ∪ B) is a uniformly random s-subset of A ∪ B. The number of its elements in both sketches, which are exactly the elements in A ∩ B, is hypergeometric. Its mean is sJ, and its variance sJ(1 − J)(N − s)/(N − 1) ≤ sJ(1 − J), where N = |A ∪ B|. This is the binomial error √(J(1 − J)/s) of Ondov et al. (2016, eq. 3), independent of the genomes' sizes.

**The Mash distance (ibid., eq. 4).** Suppose each k-mer survives a mutation process with probability e^{−kD} (a Poisson process of rate D per base), and |A| = |B|. Then the shared fraction is w = |A ∩ B|/|A|, and J = w/(2 − w), so w = 2J/(1 + J) and

D = −(1/k) ln(2J/(1 + J)).

Under independent substitutions at rate p (Blanca et al. 2022; Hera et al. 2023), a k-mer survives with probability (1 − p)^k. The same formula then estimates −ln(1 − p) = p + p²/2 + ⋯. The checks show both:

- D tracks d under the Poisson model, and −ln(1 − p) under substitutions, for p ≤ 0.2;
- the sketches' D lies within Ondov et al.'s Table 1 bounds of the D of the exact sets.

**The p-value (ibid., eqs. 1, 5, 8).** A random k-mer, of K possible, lies in a random genome of n distinct k-mers with probability P = 1 − (1 − 1/K)ⁿ. Two random genomes then have the expected Jaccard index r = P_X P_Y/(P_X + P_Y − P_X P_Y). The p-value of x shared hashes among s is the binomial tail I_r(x, s − x + 1).

*Caveat.* Mash takes K = 4^k. But a sketch of canonical k-mers draws from the orbits of θ on A^k, of which there are (4^k + 4^{k/2})/2 for even k (the palindromes are fixed points) and 4^k/2 for odd k, by Burnside's lemma for the group {1, θ}. With 4^k the chance matches expected of two random megabases at k = 12 are 28, where 63 are observed (and 56 expected with the canonical K). Mash's p-value is therefore too small for canonical sketches. The code uses the canonical K (`canonical_kmer_space`), and reproduces Mash when given 4^k.

### 7.2 FracMinHash

FRAC_s(A) = {h(a) : h(a) ≤ sH} keeps every hash below a fixed threshold, so that |FRAC_s(A)| ~ Binom(|A|, s) (Irber et al. 2022; Hera et al. 2023).

*Remark (why containment).* FRAC_s preserves both lattice operations:

FRAC(A ∪ B) = FRAC(A) ∪ FRAC(B),  FRAC(A ∩ B) = FRAC(A) ∩ FRAC(B).

It is a lattice homomorphism from sets to sets of hashes, where bottom-s preserves only ∪. So FracMinHash estimates |A ∩ B| and |A| separately, and their ratio, the containment C(A, B) = |A ∩ B|/|A|. Bottom-s MinHash estimates only the Jaccard index, a ratio taken inside A ∪ B.

**Proposition 7.2 (Hera et al. 2023, Theorems 1, 3, 5).** Let Ĉ = |FRAC(A) ∩ FRAC(B)|/|FRAC(A)|.

1. E[Ĉ · 1_{FRAC(A) ≠ ∅}] = C(A, B)(1 − (1 − s)^{|A|}).
2. To first order, Var Ĉ ≈ mn(1 − s)/(s(m + n)³), with n = |A ∩ B| and m = |A ∖ B|.
3. Under independent substitutions at rate p, with no repeated k-mers, E[C] = (1 − p)^k.

So ANI = 1 − p is estimated by Ĉ^{1/k}. The checks verify each statement in standard errors (`containment`, `containment_debiased`, `ani`, `containment_variance`). The containment of a subset in its superset is 1 exactly.

### 7.3 HyperLogLog

A hash is split into a register index j (the top p bits) and the rest w. The register keeps the largest ρ(w), the position of the first one bit. A register that has seen ν elements holds about log₂ ν, so a harmonic mean of the 2^{M_j} estimates n/m.

**Theorem 7.3 (Flajolet, Fusy, Gandouet and Meunier 2007).** The estimate E = α_m m²/Σ 2^{−M_j} is asymptotically unbiased. Its relative standard error is β_m/√m, with β_m → √(3 ln 2 − 1) ≈ 1.03896.

Small cardinalities need linear counting, m ln(m/V) for V empty registers. With 64-bit hashes no large-range correction is needed (Heule, Nunkesser and Hall 2013). Ertl (2017) Poissonises the number of elements, so that the registers become independent, and corrects the contributions of empty and saturated registers by the series σ and τ. The resulting *improved raw estimator* (his eq. 32) is nearly unbiased over the whole range, with no switch between estimators and no empirical tables. It is the default (`estimate`), and the checks find its error over 300 sets of 10 to 3·10⁵ elements within 1.04/√m.

*Remark (the semilattice).* Registers under max form a join semilattice, and the sketch is a homomorphism from (sets, ∪) to it: HLL(A ∪ B) = max(HLL(A), HLL(B)), exactly. Intersections are not preserved. The Jaccard index is therefore estimated by inclusion and exclusion, |A ∩ B| = |A| + |B| − |A ∪ B|, from three estimates (Baker and Langmead 2019). Its error is about √3 (1 + J)/√m in absolute terms, larger than MinHash's when J is small.

### 7.4 Distance Matrices as Task Graphs

The n × n matrix is cut into blocks. Each block on or above the diagonal is a task of a `calx::scheduler::TaskGraph`, computing its pairs and writing both (i, j) and (j, i). No two tasks write one cell, and the distance is a pure function of two sketches. So the matrix is the same, bit for bit, under every schedule. The trace records where each block ran. This is the domain's rule of determinism (shards writing their own slots) in its simplest form, and the checks verify it for blocks of 1, 7 and 32 on one and eight threads.

## 8. Notation and Code

| Symbol | Meaning | Code |
|---|---|---|
| A*, φ, φ^ω | free monoid, morphism, fixed point | `Word`, `Morphism`, `fixed_point_prefix` |
| b[i], Z[i] | border array, Z-array | `border_array`, `z_array` |
| ℓ₁ ≥ ⋯ ≥ ℓ_m | Chen–Fox–Lyndon factors | `lyndon_factorisation`, `lyndon_factors` |
| L_k(n) | Lyndon words of length n | `lyndon_count` |
| SA, rank, LCP | suffix array, its inverse, longest common prefixes | `suffix_array`, `inverse_suffix_array`, `lcp_array` |
| p(ℓ) | factor complexity | `factor_complexity` |
| L, C, LF | BWT, letter starts, last-to-first map | `bwt`, `letter_starts`, `lf_map` |
| π_L | standard permutation | `standard_permutation` |
| eBWT, bBWT | extended and bijective BWT | `extended_bwt`, `bijective_bwt` |
| rank, select | on bit vectors and sequences | `BitVector`, `WaveletMatrix` |
| δ, f | Aho–Corasick transition and failure | `AhoCorasick::delta`, `failure` |
| θ | reverse complement | `sequence::reverse_complement` |
| (k, l, s) | bi-interval | `index::BiInterval` |
| S(A), FRAC_s(A) | bottom-s and FracMinHash sketches | `sketch::MinHash`, `sketch::FracMinHash` |
| J, D, C | Jaccard, Mash distance, containment | `compare`, `mash_distance`, `containment` |
| M_j, n̂ | HyperLogLog registers and estimate | `HyperLogLog::registers`, `estimate` |

## 9. References

- Aho, A. V. and Corasick, M. J. (1975), Efficient string matching: an aid to bibliographic search, *Comm. ACM* 18, 333–340.
- Allouche, J.-P. and Shallit, J. (2003), *Automatic Sequences*, Cambridge University Press.
- Baker, D. N. and Langmead, B. (2019), Dashing: fast and accurate genomic distances with HyperLogLog, *Genome Biology* 20, 265.
- Bentley, J. L., Sleator, D. D., Tarjan, R. E. and Wei, V. K. (1986), A locally adaptive data compression scheme, *Comm. ACM* 29, 320–330.
- Blanca, A., Harris, R. S., Koslicki, D. and Medvedev, P. (2022), The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches, *J. Comput. Biol.* 29, 155–168.
- Booth, K. S. (1980), Lexicographically least circular substrings, *Inf. Process. Lett.* 10, 240–242.
- Brlek, S. (1989), Enumeration of factors in the Thue–Morse word, *Discrete Appl. Math.* 24, 83–96.
- Broder, A. Z. (1997), On the resemblance and containment of documents, *Compression and Complexity of Sequences*, 21–29.
- Burrows, M. and Wheeler, D. J. (1994), A block-sorting lossless data compression algorithm, DEC SRC Research Report 124.
- Chen, K.-T., Fox, R. H. and Lyndon, R. C. (1958), Free differential calculus IV, *Ann. Math.* 68, 81–95.
- Claude, F., Navarro, G. and Ordóñez, A. (2015), The wavelet matrix: an efficient wavelet tree for large alphabets, *Information Systems* 47, 15–32.
- Crochemore, M., Désarménien, J. and Perrin, D. (2005), A note on the Burrows–Wheeler transformation, *Theoret. Comput. Sci.* 332, 567–572.
- de Luca, A. and Varricchio, S. (1989), Some combinatorial properties of the Thue–Morse sequence and a problem in semigroups, *Theoret. Comput. Sci.* 63, 333–348.
- Duval, J.-P. (1983), Factorizing words over an ordered alphabet, *J. Algorithms* 4, 363–381.
- Duval, J.-P. (1988), Génération d'une section des classes de conjugaison et arbre des mots de Lyndon de longueur bornée, *Theoret. Comput. Sci.* 60, 255–283.
- Ertl, O. (2017), New cardinality estimation algorithms for HyperLogLog sketches, arXiv:1702.01284.
- Ferragina, P. and Manzini, G. (2000), Opportunistic data structures with applications, *FOCS*, 390–398.
- Fine, N. J. and Wilf, H. S. (1965), Uniqueness theorems for periodic functions, *Proc. Amer. Math. Soc.* 16, 109–114.
- Flajolet, P., Fusy, É., Gandouet, O. and Meunier, F. (2007), HyperLogLog: the analysis of a near-optimal cardinality estimation algorithm, *AofA*, DMTCS Proceedings AH, 127–146.
- Gessel, I. M. and Reutenauer, C. (1993), Counting permutations with given cycle structure and descent set, *J. Combin. Theory A* 64, 189–215.
- Gil, J. Y. and Scott, D. A. (2012), A bijective string sorting transform, arXiv:1201.3077.
- Gusfield, D. (1997), *Algorithms on Strings, Trees and Sequences*, Cambridge University Press.
- Hera, M. R., Pierce-Ward, N. T. and Koslicki, D. (2023), Deriving confidence intervals for mutation rates across a wide range of evolutionary distances using FracMinHash, *Genome Research* 33, 1061–1068.
- Heule, S., Nunkesser, M. and Hall, A. (2013), HyperLogLog in practice, *EDBT*, 683–692.
- Irber, L. et al. (2022), Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers, bioRxiv 2022.01.11.475838.
- Jacobson, G. (1989), Space-efficient static trees and graphs, *FOCS*, 549–554.
- Kasai, T., Lee, G., Arimura, H., Arikawa, S. and Park, K. (2001), Linear-time longest-common-prefix computation in suffix arrays and its applications, *CPM*, LNCS 2089, 181–192.
- Knuth, D. E., Morris, J. H. and Pratt, V. R. (1977), Fast pattern matching in strings, *SIAM J. Comput.* 6, 323–350.
- Ko, P. and Aluru, S. (2003), Space efficient linear time construction of suffix arrays, *CPM*, LNCS 2676, 200–210.
- Kufleitner, M. (2009), On bijective variants of the Burrows–Wheeler transform, *Prague Stringology Conference*, 65–79.
- Lam, T. W., Li, R., Tam, A., Wong, S., Wu, E. and Yiu, S. M. (2009), High throughput short read alignment via bi-directional BWT, *BIBM*, 31–36.
- Li, H. (2012), Exploring single-sample SNP and INDEL calling with whole-genome de novo assembly, *Bioinformatics* 28, 1838–1844.
- Lothaire, M. (1983), *Combinatorics on Words*, Addison-Wesley; (2002), *Algebraic Combinatorics on Words*, Cambridge University Press.
- Lyndon, R. C. and Schützenberger, M.-P. (1962), The equation a^M = b^N c^P in a free group, *Michigan Math. J.* 9, 289–298.
- Manber, U. and Myers, G. (1993), Suffix arrays: a new method for on-line string searches, *SIAM J. Comput.* 22, 935–948.
- Mantaci, S., Restivo, A. and Sciortino, M. (2003), Burrows–Wheeler transform and Sturmian words, *Inf. Process. Lett.* 86, 241–246.
- Mantaci, S., Restivo, A., Rosone, G. and Sciortino, M. (2007), An extension of the Burrows–Wheeler transform, *Theoret. Comput. Sci.* 387, 298–312.
- Morse, M. and Hedlund, G. A. (1940), Symbolic dynamics II: Sturmian trajectories, *Amer. J. Math.* 62, 1–42.
- Nong, G., Zhang, S. and Chan, W. H. (2011), Two efficient algorithms for linear time suffix array construction, *IEEE Trans. Computers* 60, 1471–1484.
- Ondov, B. D. et al. (2016), Mash: fast genome and metagenome distance estimation using MinHash, *Genome Biology* 17, 132.
- Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M. and Yorke, J. A. (2004), Reducing storage requirements for biological sequence comparison, *Bioinformatics* 20, 3363–3369.
- Schleimer, S., Wilkerson, D. S. and Aiken, A. (2003), Winnowing: local algorithms for document fingerprinting, *SIGMOD*, 76–85.
- Thue, A. (1912), Über die gegenseitige Lage gleicher Teile gewisser Zeichenreihen, *Norske Vid. Selsk. Skr. I Mat.-Nat. Kl.* 1, 1–67.
- Vigna, S. (2008), Broadword implementation of rank/select queries, *WEA*, LNCS 5038, 154–168.
