# index: The Indexes of a Genome

*Status: draft (2026-10-24).* The FMD-index of a genome and its reverse complement, with bidirectional extension and super-maximal exact matches; the minimizer index of a genome; patterns of IUPAC letters found on both strands by one Aho–Corasick automaton. k-mer counting, Bloom filters and count-min sketches are planned (the design note's §4.2).

## 1. Purpose

A genome is searched by its indexes. Two are here, and they serve the two families of mappers:

- **The FMD-index** finds every exact occurrence of a pattern on both strands in time proportional to the pattern, and the super-maximal exact matches of a read: the seeds of BWA-MEM. It is an FM-index of the genome's contigs and their reverse complements together. That set of words is closed under reverse complementation, so the index extends a match to the right as easily as to the left (Li 2012).
- **The minimizer index** stores the (w, k)-minimizers of the genome with where they occur, and turns a read's minimizers into anchors: the seeds of a minimap2-like mapper (Li 2018).

The combinatorics of words they rest on is morphe's: SA-IS, the Burrows–Wheeler transform, bit vectors with rank, and the FM-index over any alphabet live in `morphe::formal`. Only what knows the double strand is here. The theory is [theory/indexes.md](../theory/indexes.md).

## 2. Content

| Header | Content |
|---|---|
| `fmd.hpp` | `FmdIndex(contigs, FmdOptions{sample, seed})`. The text c₀ $₀ c̄₀ $₁ … is suffix-sorted by morphe's SA-IS, and each N is replaced by a base drawn from `Stream(seed, FMD_N, contig, 0)` and recorded (`ambiguous()`, `overlaps_ambiguous`). The BWT over {$, A, C, G, T} is kept in BWA's interleaved blocks of 128 symbols (counts and 2-bit symbols in one cache line): `occ(c, i)`, `occ4(i)`, `C(c)`, `bwt(i)`, `lf(i)`. `BiInterval{k, l, s}` with `init`, `backward`, `forward` (Li 2012, Algorithms 2–3); `search`, `count`, `locate` → `Hit{contig, position, reverse}`, `suffix(row)` by the suffix array sampled at text positions and sequence starts; `smems(query, min_length, min_occ)` → `Smem{begin, end, interval}` (Algorithm 5, BWA-MEM's `bwt_smem1` without re-seeding), `smems_at`, and `smems(pool, queries, …)` |
| `minimizer_index.hpp` | `MinimizerIndex(pool, contigs, MinimizerOptions{k, w, ignore, chunk})`. It is built from fixed chunks with margins of w − 1 and w + k − 2 bases, which give exactly the minimizers of the whole contigs, then bucketed by hash and sorted on the pool, deterministically. It holds a flat array sorted by (hash, contig, position, strand) and an open-addressing table of hash ranges: `lookup(hash)`, `decode` → `Occurrence{contig, position, reverse}`, `max_occurrences()` (minimap2's `-f` cut-off), `anchors(query)` → `Anchor{contig, target, query, reverse}`, sorted for chaining |
| `patterns.hpp` | `iupac_words(pattern, max_words)`: the words a pattern of IUPAC letters stands for, the product of its letters' sets, in order. `PatternSet(patterns, both_strands, max_words)`: every pattern's words, and the reverse complements of them that are not its own words, in one `morphe::formal::AhoCorasick` automaton; `find(text)` → `PatternHit{begin, end, pattern, reverse}` in one pass over the forward strand, a palindromic site once, forward; a window with N matches nothing |

## 3. Dependencies

`genomics::sequence` (codes, k-mers, minimizers); `morphe::formal` (SA-IS, bit vectors with rank); `calx::scheduler` (the pool, parallel loops); `calx::random` (the bases for N). The checks also use `genomics::simulate`.

## 4. Checks

| Check | What it verifies |
|---|---|
| `fmd` | The BWT and Occ, at every row, against a BWT from an independent suffix array (prefix doubling). The text closed under reverse complement, and the runs of N recorded. Counts and positions of 3000 patterns on both strands against brute force across five contigs. Bidirectional extension (forwards, backwards, and mixed from any base) against backward search of P and P̄. SMEMs against their definition by brute force on mutated reads with N, for min_occ 1 and 2, with the length filter. SMEMs the same on one and four threads. Escherichia coli K-12 (NC_000913.3): build time and memory, 1000 sampled 32-mers located at their true places, and the longest SMEM of 4000 simulated 150-base reads on its true diagonal |
| `minimizer_index` | Chunks of 997 bases with margins give exactly the minimizers of each whole contig, for four (k, w), across runs of N and contigs shorter than a window. The index is the same on one and eight threads and for any chunk length. Lookups of absent hashes find nothing. The cut-off is minimap2's quantile plus one. Anchors equal their brute-force definition. Exact reads of either strand anchor on their true diagonal. Escherichia coli: size, time on one and eight threads, density against 2/(w + 1) |
| `patterns` | The words of GAATTC, GGWCC, CCNGG and RYBN, in order, each reading its pattern; empty patterns, letters outside IUPAC's and too many words refused. Nine patterns on both strands of 300 kb of random sequence with runs of N against their definition by brute force (forward where the window reads the pattern, reverse where its reverse complement does and it does not); palindromic sites once; one strand alone |

Measured on the development machine (i7-11370H, WSL, other builds running, so the ranges are of several runs), 2026-10-24:

- the FMD-index of Escherichia coli builds in 2.2 to 3.0 s and takes 7.3 MB, with a peak of 207 MB while building;
- a 32-mer is located in 6 µs;
- the SMEMs of a 150-base read take 7 to 16 µs on eight threads;
- the (10, 15)-minimizer index (867 206 occurrences of 838 542 minimizers, 26 MB) builds in 0.23 to 0.31 s on one thread and 0.11 to 0.14 s on eight.

## 5. References

- Aho, A. V. and Corasick, M. J. (1975), Efficient string matching: an aid to bibliographic search, *Communications of the ACM* 18, 333–340.
- Cornish-Bowden, A. (1985), Nomenclature for incompletely specified bases in nucleic acid sequences, *Nucleic Acids Research* 13, 3021–3030.
- Li, H. (2012), Exploring single-sample SNP and INDEL calling with whole-genome de novo assembly, *Bioinformatics* 28, 1838–1844.
- Li, H. (2013), Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM, arXiv:1303.3997.
- Li, H. and Durbin, R. (2009), Fast and accurate short read alignment with Burrows–Wheeler transform, *Bioinformatics* 25, 1754–1760.
- Li, H. (2018), Minimap2: pairwise alignment for nucleotide sequences, *Bioinformatics* 34, 3094–3100.
- Ferragina, P. and Manzini, G. (2000), Opportunistic data structures with applications, *FOCS*, 390–398.
- Nong, G., Zhang, S. and Chan, W. H. (2011), Two efficient algorithms for linear time suffix array construction, *IEEE Trans. Computers* 60, 1471–1484.
- Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M. and Yorke, J. A. (2004), Reducing storage requirements for biological sequence comparison, *Bioinformatics* 20, 3363–3369.
- Lemire, D. (2019), Fast random integer generation in an interval, *ACM TOMACS* 29, 3 (the multiply-shift reduction).

## 6. Open Questions

1. **Genomes beyond 2^32 − 1 letters of text** (2.1 Gb on both strands) need 64-bit suffix arrays and an external-memory or blockwise construction (BWA's BCR or ropebwt); SA-IS here holds the whole suffix array.
2. **Re-seeding** of long SMEMs and the third round of BWA-MEM's seeding belong with the mapper.
3. **The minimizer index of a large genome** would split its table by hash buckets, as minimap2 does, and could store singletons in the table itself.
