# genomics

*Status: draft. Field, height 0.* Genomes and the sequences they are read as: their alphabets and formats, their composition, the alignment of one sequence to another, the indexes and sketches that make a genome searchable, the mapping and assembly of reads, the variants of individuals, hidden Markov models and multiple alignment, expression and annotation ([ADR 0043](../../development/decisions/0043-the-genome-infrastructure.md)). Drawn by the [Genome Atlas](../../applications/genome-atlas/). Serves use, learning and research.

## 1. Concepts

- **A sequence is a word; a genome is a word read on two strands.** Everything that knows only words (suffix arrays, the Burrows–Wheeler transform, Lyndon words, FM-indexes, automata, the alignment of words over semirings, the polytope algebra) is morphe's `formal` module, and this domain uses it. What knows the double strand is here: complements, canonical k-mers, the FMD-index of a genome and its reverse complement, the genetic codes.
- **An alignment is a path; a scoring scheme is a weight in a semiring; an algorithm is a schedule.** Every pairwise alignment is a sum over the paths of an edit graph. The tropical semiring gives the best score, the counting semiring the number of alignments, the probability semiring a pair HMM, the polytope algebra every optimal alignment of every scoring at once (Pachter and Sturmfels 2004). morphe's generic dynamic programme defines the answers; the fast algorithms here (Gotoh, Hirschberg, banded, Myers' bit vectors, Farrar's stripes on SIMD lanes, the wavefront algorithm, X-drop, partial orders) are schedules of the same sum, each checked against it ([theory](theory/alignment.md)).
- **Speed is a property of a schedule, never of a definition.** Vector lanes come from calx's CPU backend (`calx::cpu::Lanes`, compiled to AVX-512, AVX2 or SSE as the build allows), threads from calx's scheduler (`ThreadPool`, `parallel_for`, `TaskGraph`); a parallel computation writes its shards to their own slots and merges them in order, so its answer never depends on the number of threads, processes or machines.
- **Truth comes from simulation.** A real genome does not say where its reads came from. Seeded simulation does: genomes with repeats, individuals with their variants, reads with their origins and errors, families with their true alignment. A shard of a simulated sample is two numbers, which is how work travels to other machines.
- **Real genomes are held to their own records.** The GenBank reader translates every CDS of fourteen RefSeq genomes exactly as RefSeq's `/translation` (4,522 proteins, selenoproteins included); composition finds *Escherichia coli*'s origin of replication from the cumulative GC skew within two hundred bases of oriC.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `sequence` | Alphabets and codes; packed DNA; the 27 NCBI genetic codes, translation, ORFs; FASTA, FASTQ and GenBank with its location grammar; k-mers, hash64, ntHash, minimizers, syncmers; composition, GC skew, CpG, codon usage ([`sequence/`](sequence/include/genomics/sequence/)) | draft |
| `simulate` | Genomes with repeats; diploid mutation with its truth and liftovers; short and long reads with their origins and errors; families of sequences with their true alignment ([`simulate/`](simulate/include/genomics/simulate/)) | draft |
| `kmers` | Canonical k-mers counted in parallel and deterministically by hash partitions; the k-mer spectrum, its error threshold, coverage, genome size and heterozygosity (a simplified GenomeScope); Bloom filters and count-min sketches ([`kmers/`](kmers/include/genomics/kmers/)) | draft |
| `align` | Pairwise alignment on morphe's generic sum over the edit graph: NCBI matrices and DNA presets; Gotoh with traceback and the field for drawing; Myers–Miller linear space; banded with a certificate; X-drop and Z-drop; Karlin–Altschul λ, K, H, E-values and BLAST's gapped tables; parametric alignment by polytope propagation; the pair HMM (Viterbi, forward and backward, posteriors, MEA); Myers' bit vectors, Farrar's striped and inter-sequence Smith–Waterman on lanes, the wavefront algorithm and BiWFA, partial-order alignment ([`align/`](align/include/genomics/align/)) | draft |
| `index` | The FMD-index of a genome and its reverse complement on morphe's SA-IS: BWA's interleaved occurrence blocks, a sampled suffix array, bi-intervals with forward and backward extension, super-maximal exact matches (Li 2012); the minimizer index, built in parallel and deterministically, with minimap2's frequency cut-off and anchors ([`index/`](index/include/genomics/index/)) | draft |
| `sketch` | Bottom-s MinHash and FracMinHash sketches of canonical k-mer sets: Jaccard, the Mash distance and p-value (with a null model corrected for canonical k-mers), containment and ANI; HyperLogLog with Ertl's estimator, unions and Jaccard; all-against-all distance matrices as task graphs ([`sketch/`](sketch/include/genomics/sketch/)) | draft |
| `mapping` | A minimap2-like mapper: minimizer or SMEM seeds, chaining by Li's DP and gap cost, primary and secondary chains and their mapping quality, base-level alignment (certified ungapped gaps, the wavefront algorithm, certified bands, X-drop and Z-drop ends, a clipping penalty), SAM, accuracy against the simulator's truth ([`mapping/`](mapping/include/genomics/mapping/)) | draft |
| `assembly` | The bidirected de Bruijn graph of canonical k-mers and its unitigs in canonical form; cleaning of tips, islands, bubbles (crossed ones included) and low-coverage unitigs, recorded step by step; contigs, N50/NG50/L50 and Lander–Waterman; the directed graph of a k-mer composition, Hierholzer, the number of strings and of de Bruijn sequences by the BEST theorem ([`assembly/`](assembly/include/genomics/assembly/)) | draft |
| `variant` | Pileups with left-aligned indels, Li's genotype likelihoods and their posterior, SNV and indel calling with depth and strand-bias filters, VCF 4.2, Tan et al.'s normalisation of indels, evaluation against a truth set ([`variant/`](variant/include/genomics/variant/)) | draft |
| `hmm` | Hidden Markov models: Viterbi, forward, backward, posteriors over the tropical and log semirings, Baum–Welch with pseudocounts in parallel, sampling; the CpG-island model of Durbin et al.; profile HMMs from alignments, scored and aligned in global, glocal and local modes ([`hmm/`](hmm/include/genomics/hmm/)) | draft |
| `msa` | The multiple alignment and its formats (aligned FASTA, Clustal, ClustalW's conservation marks), consensus and column entropy; guide trees from MUSCLE's k-mer distances and corrected identities by UPGMA (nearest-neighbour chain), CLUSTAL W weights; progressive alignment of profiles with sum-of-pairs scores and occupancy-scaled gaps, MUSCLE's re-estimated tree and tree-dependent refinement; the centre star; alignment by partial orders; the SP score, and BAliBASE's Q and TC against a reference ([`msa/`](msa/include/genomics/msa/)) | draft |
| `expression` | Transcriptomes and simulated RNA-seq with their truth; pseudoalignment by k-mer equivalence classes (kallisto's, k = 31); abundances by EM with effective lengths, TPM and the bootstrap; size factors by the median of ratios, negative binomial GLMs by IRLS, Cox–Reid dispersions shrunk toward a trend, Wald tests and Benjamini–Hochberg ([`expression/`](expression/include/genomics/expression/)) | draft |
| `annotation` | Features of GenBank and GFF3 as tracks (across the origin, translations checked); PWMs, information, exact score distributions, scanning, Gibbs sampling and MEME's EM; a self-trained prokaryotic gene finder; maximal and supermaximal repeats on both strands by morphe's suffix array ([`annotation/`](annotation/include/genomics/annotation/)) | draft |

## 3. Conventions

- **Codes.** DNA is A = 0, C = 1, G = 2, T = 3, N = 4; protein is the order of the NCBI matrices `ARNDCQEGHILKMFPSTWYVBZX*`; `GAP` = 0xFF in aligned rows (`sequence/alphabet.hpp`).
- **Coordinates** are 0-based and half-open inside, 1-based only in the formats that say so (SAM, VCF, GenBank, GFF3). A strand is `bool reverse`.
- **Scores and penalties.** A scoring is maximised: a substitution matrix and affine gaps, a gap of length L costing open + L · extend. Penalties (the wavefront algorithm's) are minimised. CIGARs use the SAM letters, I for a base of the query against a gap and D for a base of the target against one.
- **Readers never throw on content**; they return problems with their lines. **Randomness** is a `calx::random::Stream` keyed by purpose and object. **Parallel** computation is deterministic.

## 4. Dependencies

- **Domains:** none, so genomics has height 0. Evolution stands on it.
- **Foundations:** morphe (`formal`: words, suffix arrays, BWT, FM-indexes, alignment over semirings, the polytope algebra; `arithmetic` and `algebra` for exact counts), calx (`cpu` lanes, `scheduler`, `random`, `linalg`, `special_functions`, and the report of the checks).
- **The data library** fetches NCBI's genomes into the cache (`data/registry/ncbi-genomes.json`); the domain reads only the files whose paths it is given.

## 5. Code and Checks

The C++ is header-only, C++20, built with CMake (`CMakeLists.txt`), each module a directory with its own `CMakeLists.txt` and the target `genomics::<module>`. The checks are built for the machine's widest lanes (`-march=native`) and run by

```
node domains/genomics/tests/run-tests.mjs
```

Each module's README lists its checks and what they verify.

The theory, with statements and proofs, is in `theory/`: [alignment](theory/alignment.md) (the edit graph over semirings, the parametric polygon, the pair HMM), [fast alignment](theory/fast-alignment.md) (bit vectors, stripes, wavefronts and BiWFA, partial orders, as schedules of one sum), [indexes](theory/indexes.md) (the FMD-index, SMEMs, minimizers, sketches), [assembly](theory/assembly.md) (the spectrum, the bidirected de Bruijn graph, Hierholzer, the matrix-tree and BEST theorems), [mapping and variants](theory/mapping-and-variants.md) (chaining, mapping quality, genotype likelihoods, the canonical form of indels), [HMMs](theory/hmm.md), [annotation](theory/annotation.md) and [multiple alignment and expression](theory/msa-and-expression.md) (sum-of-pairs, UPGMA, Gusfield's bound; pseudoalignment, the EM, the negative binomial).

## 6. Open Questions

1. **Where general sketches belong.** MinHash and HyperLogLog sketch any set; Bloom filters filter any keys. They are here because genomes need them first; when a second domain does, they move to a foundation.
2. **HMMs and the negative binomial model** are statistics, and move to the inference foundation when it exists (CLAUDE.md §8, open decision 3).
3. **Trees, and the methods that build them from distances.** `msa` keeps a guide tree (UPGMA, distances) of its own because genomics cannot depend on evolution, whose `phylogeny` holds the full ones; both belong in a foundation (ADR 0043, open).
4. **Calling from long reads.** The mapper places long reads well, but the caller's genotype likelihoods (Li 2011) treat each base's error as a substitution of known probability; long reads' errors, mostly insertions and deletions at several percent, break that model, and most calls from them are false. A model of their errors (a pair HMM over the read against each haplotype, or a consensus of the reads first, as the partial order already gives) is needed.
5. **Large genomes.** Everything is in memory: a bacterial genome is instant, a human chromosome comfortable, a whole human genome heavy. Memory-mapped indexes would carry it further.

## 7. References

Each module's README gives the sources of its algorithms. The standard references of the domain:

- Durbin, R., Eddy, S. R., Krogh, A. and Mitchison, G. (1998), *Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids*, Cambridge University Press.
- Gusfield, D. (1997), *Algorithms on Strings, Trees, and Sequences*, Cambridge University Press.
- Pachter, L. and Sturmfels, B. (eds.) (2005), *Algebraic Statistics for Computational Biology*, Cambridge University Press.
- Mäkinen, V., Belazzougui, D., Cunial, F. and Tomescu, A. I. (2015), *Genome-Scale Algorithm Design*, Cambridge University Press.
- Compeau, P. and Pevzner, P. (2015), *Bioinformatics Algorithms: An Active Learning Approach*, Active Learning Publishers.
