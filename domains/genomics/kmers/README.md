# kmers

*Status: draft.* The k-mers of a sample of reads, counted in parallel and deterministically; the k-mer spectrum and what it says of the genome sequenced (the error threshold, the coverage, the genome size and the heterozygosity); approximate membership and counts in fixed memory (Bloom filters and count-min sketches); and strobemers, seeds that span more than they read and survive an indel between their strobes. It serves use, as the first step of assembly and of reference-free profiling, and learning, through the spectrum's model. The theory is [theory/assembly.md](../theory/assembly.md), §1–§2.

## 1. Concepts

- **A canonical k-mer** names a piece of a double strand: a k-mer and its reverse complement are counted under the smaller of the two (sequence/kmer.hpp), k ≤ 32.
- **Counting by partitions.** The k-mers are partitioned by the high bits of a 64-bit mix (SplitMix64's finaliser), and each partition is counted by one task in an open-addressing table of its own, so no locks are needed (the scheme of Jellyfish and KMC in spirit). Reads are cut into chunks that depend only on the input and the options; each chunk scatters its k-mers by a counting sort, and each partition takes them from the chunks in chunk order. Every table therefore receives the same keys in the same order whatever the number of threads, and ends in the same state.
- **The spectrum** h(c), the number of distinct k-mers seen c times, is a descent of errors, a first minimum, a peak at the k-mer coverage λ = c (L − k + 1)/L, and peaks at λ/2 (heterozygous k-mers) and at jλ (repeats). It is fitted as a mixture of Poisson laws at λ/2, λ, 2λ, 3λ and 4λ by EM, a stated simplification of GenomeScope; the genome size is the k-mers above the threshold, counted with multiplicity, over λ, and the heterozygosity follows from the weight of the λ/2 component.
- **Strobemers** (Sahlin 2021) join n strobes of ℓ bases, the first at a position and each other chosen from a window downstream: the ℓ-mer of least hash (minstrobes) or the one that minimises a hash conditioned on the strobes before (randstrobes). A one-base indel between the first strobe and its window leaves the seed whole with probability 1 − 2/(w + 1) for a window of w positions, where a k-mer of the same span is always lost. Their sibling seeds, closed and open syncmers, are `sequence`'s (`sequence/kmer.hpp`); their density is checked here against its exact value. The theory is [theory/signals-classification-and-comparison.md](../theory/signals-classification-and-comparison.md), §3.
- **Filters.** A Bloom filter answers membership with no false negatives and a false-positive rate (1 − e^{−hn/m})^h; a count-min sketch never under-counts and over-counts by more than ε‖a‖₁ with probability at most δ. Both are functions of a seed, so shards merge.

## 2. Content

| Header | Content |
|---|---|
| `counter.hpp` | `mix64`; `KmerCount{kmer, count}`; `KmerMap` (open addressing, linear probing: `at`, `find`, `prefetch`, `for_each`); `CounterOptions{k, partition_bits, chunk_bases, batch_bases}`; `KmerCounter` (`add(n, get, pool)`, `add(sequences, pool)`, `count`, `distinct`, `total`, `for_each`, `sorted(min_count, pool)`, `table`, `memory_bytes`, `peak_bytes`, `release_buffers`); `count_kmers(sequences, pool, options)` |
| `spectrum.hpp` | `Spectrum{k, histogram, distinct, total}`; `spectrum(counts, pool, max_count)`; `error_threshold`; `coverage_peak`; `kmer_coverage`, `base_coverage`; `FitOptions{threshold, coverage, max_iterations}`; `fit_spectrum(spectrum, options)` → `SpectrumFit{threshold, peak, coverage, kmers[5], log_likelihood, genome_size, genome_size_at_peak, heterozygous_windows, heterozygosity, error_kmers, solid_kmers, iterations, converged}` |
| `filters.hpp` | `bloom_false_positive_rate(n, m, h)`, `bloom_optimal_hashes(n, m)`, `bloom_bits(n, p)`; `BloomFilter(bits, hashes, seed)` with `for_capacity`, `insert`, `contains`, `fill`, `merge`; `CountMinSketch(width, depth, seed)` with `for_error(ε, δ)`, `add`, `estimate`, `epsilon`, `delta`, `merge` |
| `strobemers.hpp` | `StrobemerShape{order, length, w_min, w_max}`; `Strobemer{hash, position, strobes, order}` with `end(length)`; `strobemer_positions(L, shape)`; `minstrobes(x, shape)` (one sliding minimum, O(L)), `randstrobes(x, shape)` (O(L n w)), of orders 2 to 4, strobes of 1 to 31 bases |

## 3. Performance

Counting the 150 million k-mers of 30× reads of a 5 Mb genome (24 million distinct) takes 1.8 to 2.6 s on 8 threads and 4.9 to 11 s on one, on the development machine (i7-11370H, shared with other builds at a load of 2.5 to 6.5), in 871 MB: 12 bytes a slot at a load between 5/16 and 5/8, 19 to 38 bytes for each distinct k-mer, and 8 bytes for each k-mer of a batch. The counts, their order and the memory are the same on any number of threads. The solid k-mers in canonical order (`sorted`) take 0.25 s on 8 threads, 0.9 s on one.

## 4. Dependencies

`genomics::sequence` (codes, canonical k-mers, `for_each_kmer`, hash64, syncmers); calx `scheduler` (the pool and `parallel_for`) and `random` (the seeds of the filters). The checks use `genomics::simulate`.

## 5. Checks

| Check | What it verifies |
|---|---|
| counter | Counts of 4000 simulated reads with errors and Ns against a `std::unordered_map`, for k = 1, 15, 21, 31, 32; the same tables, slot order and memory on 1, 2 and 8 threads; a long sequence cut into chunks of 1 to 10⁵ bases overlapping by k − 1; 1 to 4096 partitions and two calls of `add` against one; the canonical order of `sorted`; absent k-mers; the memory against its statement |
| spectrum | Error-free reads at 30×: the peak and λ at c (L − k + 1)/L, the genome size within 1 %; reads with errors: λ thinned by error-free k-mers within 2 %, the threshold separating errors from the genome to 10⁻³, the genome size with repeats within 2 %, the 4-copy repeats in the 4λ component; a diploid genome: the heterozygous windows and the heterozygosity against the truth of its mutations, the haploid size; at 2 % heterozygosity the mode at λ/2, resolved by the coverage given |
| filters | Bloom: no false negatives; the false-positive rate within 4σ of (1 − e^{−hn/m})^h for h = 1 to 12 and of ρ^h given the fill; the fill against 1 − e^{−hn/m}; filters sized for a rate; shards merged. Count-min: no under-count on a Zipf stream; over-counts beyond ε‖a‖₁ below δ; one row's mean excess against Σ_{y≠x} a(y)/w; shards merged |
| strobemers (8) | Minstrobes and randstrobes of orders 2 to 4 against their definitions by brute force on 60 random sequences with bases N (101 326 seeds); constructed: a minstrobe surviving a deletion between its strobes with the same strobes and hash, where the 85-base k-mer of its span occurs nowhere after it; over 3000 random one-base indels before the window (ℓ = 15, w = 31), minstrobes kept 0.942 and randstrobes 0.937 of the time (at least 0.9 asserted; 1 − 2/(w + 1) = 0.938), k-mers of the same span never; the density of closed syncmers at (k, s) = (9, 4), (11, 7), (12, 8) on a random 2 Mb against its exact value by enumeration of every k-mer, within 5σ (the exact values 0.33408, 0.40376, 0.39952 beside Edgar's 1/3, 2/5, 2/5) |

## 6. References

- Marçais, G. and Kingsford, C. (2011), A fast, lock-free approach for efficient parallel counting of occurrences of k-mers, *Bioinformatics* 27, 764–770.
- Kokot, M., Długosz, M. and Deorowicz, S. (2017), KMC 3: counting and manipulating k-mer statistics, *Bioinformatics* 33, 2759–2761.
- Steele, G. L., Lea, D. and Flood, C. H. (2014), Fast splittable pseudorandom number generators, *OOPSLA*.
- Knuth, D. E. (1998), *The Art of Computer Programming* 3, 2nd ed., §6.4.
- Lander, E. S. and Waterman, M. S. (1988), Genomic mapping by fingerprinting random clones: a mathematical analysis, *Genomics* 2, 231–239.
- Li, X. and Waterman, M. S. (2003), Estimating the repeat structure and length of DNA sequences using ℓ-tuples, *Genome Research* 13, 1916–1922.
- Vurture, G. W. et al. (2017), GenomeScope: fast reference-free genome profiling from short reads, *Bioinformatics* 33, 2202–2204.
- Dempster, A. P., Laird, N. M. and Rubin, D. B. (1977), Maximum likelihood from incomplete data via the EM algorithm, *J. R. Stat. Soc. B* 39, 1–38.
- Bloom, B. H. (1970), Space/time trade-offs in hash coding with allowable errors, *Comm. ACM* 13, 422–426.
- Kirsch, A. and Mitzenmacher, M. (2006), Less hashing, same performance: building a better Bloom filter, *ESA*, LNCS 4168, 456–467.
- Broder, A. and Mitzenmacher, M. (2004), Network applications of Bloom filters: a survey, *Internet Mathematics* 1, 485–509.
- Cormode, G. and Muthukrishnan, S. (2005), An improved data stream summary: the count-min sketch and its applications, *J. Algorithms* 55, 58–75.
- Sahlin, K. (2021), Effective sequence similarity detection with strobemers, *Genome Research* 31, 2080–2094.
- Edgar, R. (2021), Syncmers are more sensitive than minimizers for selecting conserved k-mers in biological sequences, *PeerJ* 9, e10805.
- Carter, J. L. and Wegman, M. N. (1979), Universal classes of hash functions, *J. Comput. Syst. Sci.* 18, 143–154.

## 7. Open Questions

1. **The spectrum's model.** Negative binomial components and GenomeScope's parametrisation of the weights by heterozygosity and duplication, once real reads are fitted; the general mixture EM belongs to the inference foundation when it exists (CLAUDE.md §8, open decision 3).
2. **Counting beyond memory.** KMC's disk bins, or a Bloom filter in front of the tables to keep singletons out (BFCounter), for genomes whose distinct k-mers do not fit.
3. **k above 32**, in two words, for long-read and pan-genome work.
4. **`sequence::for_each_kmer` on a span.** The counter copies the pieces of sequences longer than a chunk to roll them; a span overload in `sequence` would avoid the copy.
