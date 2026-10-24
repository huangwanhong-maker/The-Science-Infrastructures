# assembly

*Status: draft.* De novo assembly of short reads on the de Bruijn graph, and assembly as an Eulerian path. The compacted de Bruijn graph of canonical k-mers as a bidirected graph, its unitigs, its cleaning of the tips and bubbles that errors and heterozygous sites make, recorded step by step, the contigs and their statistics; and, for learning, the directed de Bruijn graph of a k-mer composition, the strings it spells by Hierholzer's algorithm, their number by the BEST theorem, and the de Bruijn sequences. It serves use (contigs from reads) and learning (why repeats make assembly ambiguous, and by how much). The theory is [theory/assembly.md](../theory/assembly.md).

## 1. Concepts

- **The graph of a double strand.** A node is a canonical k-mer, read in two orientations; an arc joins two oriented nodes that overlap by k − 1 bases, and comes paired with the arc that reads the same overlap on the other strand. k is odd, so no k-mer is its own reverse complement. The structure is a skew-symmetric graph: a digraph on oriented nodes with the involution o ↦ o⁻ reversing its arcs.
- **Unitigs** are the maximal non-branching paths. They partition the nodes, and the arcs between them join only their ends, so the compacted graph spells the same paths. It is put in a canonical form (each unitig spelled from its smaller end, unitigs numbered by their least k-mer), so that it is a function of the k-mer set alone: compacting after cleaning gives the graph that building from the remaining k-mers gives.
- **Cleaning** removes tips (dead ends shorter than 2k, less covered than the branch they leave) and short islands, pops bubbles (two or more single-unitig arms between the same two sides, and the crossed bubble two heterozygous sites exactly k apart make), and removes low-coverage unitigs between better covered neighbours. Each step is recorded with what was removed, why, and which unitigs were merged after it, so an application can play the cleaning.
- **Assembly as an Eulerian path.** The strings with a given k-mer composition are the Eulerian trails of its directed de Bruijn graph, read by their labels; their number is the BEST count of trails divided by ∏ m(x)! for the multiplicities of repeated k-mers. A repeat of k − 1 bases or more is a vertex the path passes twice, and the number of strings is the measure of the ambiguity.

## 2. Content

| Header | Content |
|---|---|
| `debruijn.hpp` | `oriented`, `element`, `is_reverse`, `flip`; `Adjacency`; `Unitig{sequence, count_sum}` with `coverage`; `Link{from, to, from_reverse, to_reverse}`; `Graph{k, unitigs, next}` with `successors`, `predecessors`, degrees, `sequence(o)`, `links`, `gfa`; `kmer_successors`; `build_graph(solid, k, pool, report)` → `Graph`; `compact(graph, removed, merged_from)` → `Graph` |
| `clean.hpp` | `CleaningOptions{tip_length, tip_ratio, bubble_length, bubble_indel, low_length, low_ratio, rounds, keep_graphs}`; `StepKind`; `Removal{unitig, length, coverage, reference, kept, from, to}`; `CleaningStep{kind, round, removed, unitigs_before, unitigs_after, merged_from, after}`; `find_tips`, `find_bubbles`, `find_low_coverage`; `clean(graph, options)` → `Cleaning{graph, steps}` |
| `assemble.hpp` | `AssemblyStats{contigs, total, largest, n50, l50, ng50, lg50}`; `statistics(lengths or contigs, genome_size)`; `lander_waterman_contigs`, `lander_waterman_contig_length`; `AssemblyOptions{k, min_count, counter, cleaning, clean, min_contig}`; `assemble(counts, pool, options, genome_size)`, `assemble(n, get, pool, …)`, `assemble(reads, pool, …)` → `Assembly{spectrum, fit, threshold, solid_kmers, raw, cleaning, contigs, stats, seconds, count_bytes, graph_bytes}` |
| `eulerian.hpp` | `Composition{k, kmers}`, `composition(s, k)`; `KmerGraph{k, nodes, edge_kmers, graph}`, `kmer_graph(composition)`; `spell(graph, walk)`; `reconstruct(graph)` (Hierholzer); `count_strings(graph)`; `complete_de_bruijn(σ, n)`; `de_bruijn_sequences(σ, n)` (BEST) and `de_bruijn_sequences_formula(σ, n)`; `de_bruijn_sequence(σ, n)` |

The graph theory (Hierholzer, the matrix-tree and BEST theorems) and the exact determinant are morphe's: `morphe/combinatorics/eulerian.hpp` and `morphe/algebra/determinant.hpp`.

## 3. Performance

A 5 Mb genome with repeats at 30× (one million reads of 150 bases with errors) is counted and assembled in 4.8 s on 8 threads (counting 1.8 s, the spectrum 0.1 s, the graph of 4.9 million solid k-mers 2.8 s, the cleaning 0.1 s) and in 14.8 s on one (6.3, 0.3, 8.0 and 0.1 s), on the development machine (i7-11370H, shared with other builds at a load near 3; runs at higher loads took up to 5.4 and 22 s). The counter holds 871 MB and the graph build 357 MB; the process peaks at 2.2 GB with the reads (166 MB) and both runs. The contigs are the same on 1 and 8 threads: 2146 of them, N50 92 566.

## 4. Dependencies

`genomics::kmers` and `genomics::sequence`; morphe `combinatorics` (Eulerian circuits, BEST), `algebra` (the determinant) and `arithmetic` (exact integers); calx `scheduler` (through kmers). The checks use `genomics::simulate`.

## 5. Checks

| Check | What it verifies |
|---|---|
| debruijn | A genome without repeats is one unitig; a R b R c makes the four unitigs and four links predicted, R at coverage 2; on error-free reads of 200 kb with repeats longer and shorter than k, every unitig a substring of the genome or its reverse complement, every solid k-mer in exactly one unitig, every unitig maximal, every link an overlap of k − 1, coverage the mean count; the same graph on 1 and 8 threads; the canonical form; compaction after removals equal to the graph of the remaining k-mers; GFA |
| clean | A planted tip and bubble found and removed, and nothing else; on reads with errors, tips and bubbles before and none after, the removed unitigs erroneous, the erroneous k-mers left below 10⁻³ and the genome's k-mers kept to 10⁻³, the cleaned graph equal to the graph of its k-mers, every step's merges accounted for; a diploid genome's heterozygous SNVs make one arm to pop for each cluster closer than k (crossed bubbles included), and the cleaning leaves one contig, a mosaic of the haplotypes |
| assemble | N50, L50, NG50 and LG50 on a worked example; contigs of error-free reads at 2×, 3× and 4× against Lander and Waterman's N e^{−cσ} and the exact expectation for uniform reads; 30× reads with errors of 200 kb with repeats: every contig exactly on the genome, at least 99.5 % of the genome's k-mers in contigs (99.9 % measured), N50 above the spacing of the repeats longer than k, the 20-base repeats inside contigs, the genome size from the spectrum, a cleaned graph without tips or bubbles, the same contigs on 1 and 8 threads |
| eulerian | Hierholzer reconstructs random strings whose composition has one string; the count of strings against all 4⁹ strings of length 9 (k = 2 to 5) and against the distinct spellings of every trail of TAATGCCATGGGATGTT; a balanced composition against all 4¹⁰ strings; the binary de Bruijn sequences 2^(2^(n−1) − n) and (σ!)^(σ^(n−1))/σ^n for σ ≤ 4; a de Bruijn sequence spelled by Hierholzer holds every word once |
| speed | A 5 Mb genome at 30×, counted and assembled on 1 and 8 threads (`GENOMICS_SPEED_MB` changes the size): times and memory reported; the same contigs; the genome size within 3 %; every contig k-mer the genome's |

## 6. References

- Pevzner, P. A. (1989), l-tuple DNA sequencing: computer analysis, *J. Biomol. Struct. Dyn.* 7, 63–73.
- Idury, R. M. and Waterman, M. S. (1995), A new algorithm for DNA sequence assembly, *J. Comput. Biol.* 2, 291–306.
- Pevzner, P. A., Tang, H. and Waterman, M. S. (2001), An Eulerian path approach to DNA fragment assembly, *PNAS* 98, 9748–9753.
- Medvedev, P., Georgiou, K., Myers, G. and Brudno, M. (2007), Computability of models for sequence assembly, *WABI*, LNBI 4645, 289–301.
- Zerbino, D. R. and Birney, E. (2008), Velvet: algorithms for de novo short read assembly using de Bruijn graphs, *Genome Research* 18, 821–829.
- Compeau, P. E. C., Pevzner, P. A. and Tesler, G. (2011), How to apply de Bruijn graphs to genome assembly, *Nature Biotechnology* 29, 987–991.
- Compeau, P. and Pevzner, P. (2015), *Bioinformatics Algorithms*, Active Learning Publishers, ch. 3.
- Chikhi, R., Limasset, A. and Medvedev, P. (2016), Compacting de Bruijn graphs from sequencing data quickly and in low memory, *Bioinformatics* 32, i201–i208.
- Goldberg, A. V. and Karzanov, A. V. (1996), Path problems in skew-symmetric graphs, *Combinatorica* 16, 353–382.
- Lander, E. S. and Waterman, M. S. (1988), Genomic mapping by fingerprinting random clones: a mathematical analysis, *Genomics* 2, 231–239.
- Miller, J. R., Koren, S. and Sutton, G. (2010), Assembly algorithms for next-generation sequencing data, *Genomics* 95, 315–327; Gurevich, A. et al. (2013), QUAST, *Bioinformatics* 29, 1072–1075.
- Hierholzer, C. (1873), Über die Möglichkeit, einen Linienzug ohne Wiederholung und ohne Unterbrechung zu umfahren, *Math. Ann.* 6, 30–32.
- de Bruijn, N. G. (1946), A combinatorial problem, *Proc. Ned. Akad. Wet.* 49, 758–764.
- van Aardenne-Ehrenfest, T. and de Bruijn, N. G. (1951), Circuits and trees in oriented linear graphs, *Simon Stevin* 28, 203–217.
- GFA Format Specification Working Group (2016), *Graphical Fragment Assembly (GFA) Format Specification*, version 1.

## 7. Open Questions

1. **Tour Bus in full.** Bubbles whose arms are paths of several unitigs, compared by their sequences, and the merging of an arm into the one kept (Velvet), rather than single-unitig arms and dropping.
2. **Read threading and pairs.** Contigs end at every repeat longer than k − 1; resolving repeats needs reads threaded through the graph, read pairs, or long reads (the overlap–layout–consensus path the design note plans).
3. **Iterative k.** Assemblies at several k (IDBA, SPAdes) combine the contiguity of large k with the sensitivity of small.
4. **Diploid assembly.** The cleaning gives a mosaic of the haplotypes; keeping the bubbles as variants, with their phase, is the variant module's and a haplotype-resolved assembler's.
