---
status: seed
---

# The Genome Infrastructure: Domains, Foundations and the Genome Atlas

*Status: seed, 2026-10-24; promoted into [ADR 0043](../decisions/0043-the-genome-infrastructure.md), the READMEs of the genomics and evolution domains and their theory notes, where its substance now lives and is kept current.* The design of the infrastructure's genome analysis, asked for by the user on 2026-10-24: "a genome analysis infrastructure ... with visualization and interactive well-designed application ... various algorithms implemented at low level, classical and advanced, high speed and high quality ... distributed parallel structures on heterogeneous devices". It is written before the code, as the contract the modules are built to, and is promoted into ADR 0043, the domains' READMEs and their theory notes as the work lands.

## 1. Principles

1. **An alignment is a path; a scoring scheme is a weight in a semiring; an algorithm is a schedule.** Every pairwise alignment problem here is a sum over the paths of an edit graph (a DAG) of products of edge weights in a semiring. The tropical semiring (max, +) gives the optimal score; (min, +) the edit distance; (+, ×) the probability of a pair HMM; the tropical semiring with multiplicities the number of co-optimal alignments; the polytope algebra the whole parametric family (Pachter and Sturmfels 2004). One generic, slow, obviously correct dynamic programme over the edit graph defines the answers; it names no object in nature, so it is morphe's (`formal/alignment.hpp`, with the polytope algebra in `formal/polytope.hpp`), and the user's rule that domains reuse the foundations and never reimplement them holds it there. The fast algorithms (Gotoh's recurrences, Hirschberg's linear space, Myers' bit vectors, Farrar's striped lanes, the wavefront algorithm, banded and X-drop extension) are *schedules* of the same sum for one semiring, each checked against the generic one. This is calx's separation of what from how (calx architecture §1.2), in the domain's own terms.
2. **Words are mathematics; genomes are nature.** The combinatorics of words and their indexes (suffix arrays, the Burrows–Wheeler transform, Lyndon words, wavelet matrices, FM-indexes over any alphabet, Aho–Corasick automata) name no object in nature, and linguistics, the next domain, needs them too: they are the first slice of morphe's `formal` module (docs/mathematics.md §3). What knows the double strand (the FMD-index of a genome and its reverse complement, canonical k-mers, the genetic code) is the genomics domain's.
3. **Pure, deterministic, parallel.** The domains read only what they are given. Every random draw is a `calx::random::Stream` keyed by what it is for, so a simulation is the same on any number of threads, processes or hosts. Parallel work is expressed as shards with an associative merge in a fixed order, so the answer never depends on the schedule.
4. **Heterogeneous by target description, not by device code in the domains.** A kernel that uses vector lanes is written once against `calx::cpu::Lanes<T, W>` (GCC's vector extensions), which the compiler realises as SSE2, AVX2 or AVX-512 for the target it is built for; `calx::cpu::host()` describes the processor. Across processes and hosts the edge (the application's cluster) dispatches shards to workers whose targets differ, weighted by measured throughput. A GPU enters through the browser (WebGPU) for kernels the page implements, checked against the C++.
5. **Checked against references.** Textbook examples with published answers (Durbin et al. 1998; Gusfield 1997; Compeau and Pevzner), closed forms (expected minimizer density 2/(w + 1), the coalescent's E[T_MRCA] = 2(1 − 1/n), Witt's formula for Lyndon words, the BEST theorem), published constants (BLOSUM62's λ = 0.3176, K = 0.134), brute force on small inputs, and truth sets from our own seeded simulations.

## 2. The Order

```
height 1   evolution      phylogeny · population genetics · selection          (stands on genomics)
height 0   genomics       sequences · alignment · indexes · sketches · mapping · assembly · variants ·
                          hidden Markov models · multiple alignment · expression · annotation
foundations   morphe::formal (words, suffix arrays, BWT, succinct structures, FM-index, Aho–Corasick)
              calx::cpu (target, lanes) · calx::scheduler (task graphs, dynamic loops) ·
              calx::linalg (symmetric eigenproblems) · calx::special_functions (gamma, beta)
```

Both domains are fields. Evolution is the README's 3.2, brought forward because the genome analysis the user asked for reaches it (phylogenies of genomes, populations of variants); its scope is the README's: phylogeny, population genetics, selection.

## 3. The Foundations

### 3.1 morphe::formal: Words and Their Indexes

Words over an integer alphabet [0, σ) as `std::vector<std::uint32_t>` (`morphe::formal::Word`); the algorithms that index words take any contiguous range of unsigned integers, so that the genomics domain passes its byte codes without copying.

| Header | Content | Sources |
|---|---|---|
| `words.hpp` | The free monoid and its morphisms (substitutions, iterated: Thue–Morse, Fibonacci); borders and periods (the failure function of Knuth, Morris and Pratt), the Z-array, KMP search; primitive roots; Lyndon words and the Chen–Fox–Lyndon factorisation by Duval's algorithm; least rotations | Lothaire 1983/2002; Knuth, Morris and Pratt 1977; Gusfield 1997 §1.3; Duval 1983; Booth 1980 |
| `suffix_array.hpp` | SA-IS in linear time; prefix doubling (Manber–Myers) and a naive sort as references; the inverse; the LCP array by Kasai's algorithm; counts of distinct substrings | Nong, Zhang and Chan 2011; Manber and Myers 1993; Kasai et al. 2001 |
| `bwt.hpp` | The Burrows–Wheeler transform from the suffix array and its inverse by the LF map; the bijective BWT through the Lyndon factorisation (Gil–Scott); the extended BWT of a multiset of words; runs and move-to-front | Burrows and Wheeler 1994; Gil and Scott 2012; Kufleitner 2009; Mantaci et al. 2007; Gessel and Reutenauer 1993 |
| `succinct.hpp` | Bit vectors with rank in O(1) and select; the wavelet matrix: access, rank and select on integer sequences | Jacobson 1989; Claude, Navarro and Ordóñez 2015 |
| `fm_index.hpp` | The FM-index of a word over any alphabet: count, locate (sampled suffix array), extract | Ferragina and Manzini 2000 |
| `aho_corasick.hpp` | The Aho–Corasick automaton of a dictionary: goto, failure and output functions; all occurrences in one pass | Aho and Corasick 1975 |
| `semiring.hpp` | Semirings as a concept, and the exact ones: Boolean, tropical (max, +) and (min, +) over the integers, counting over ℕ (morphe's integers), the tropical semiring with multiplicities (score and number of optimal paths) | Golan 1999; Mohri 2002; Pachter and Sturmfels 2005, ch. 2 |
| `alignment.hpp` | The edit graph of two words and the sum over its paths in any semiring: global, local, semi-global and overlap; linear and affine gaps (three states); weights from a function of the two letters; traceback in idempotent semirings; the reference every fast alignment is checked against | Needleman and Wunsch 1970; Smith and Waterman 1981; Gotoh 1982; Mohri 2002 |
| `polytope.hpp` | The polytope algebra of lattice polygons (⊕ the convex hull of the union, ⊗ the Minkowski sum), exact on integer points; polytope propagation over the edit graph, the alignment polygon, its normal fan and the optimal alignment of each cone (parametric alignment) | Pachter and Sturmfels 2004 (PNAS 101, 16138); Gusfield, Balasubramanian and Naor 1994; Waterman, Eggert and Lander 1992 |

Checks: the semiring axioms on samples; the generic DP against enumeration of every alignment of short words (scores, counts, polygons); the polygon's vertices each attained by the alignment its cone selects; SA-IS against the naive sort on random words, unary words, Fibonacci and Thue–Morse words; the BWT inverted; the bijective BWT a bijection on all binary words up to length 12; Witt's count of Lyndon words; the BWT of a standard Sturmian word of the form b^p a^q (Mantaci, Restivo and Sciortino 2003); distinct substrings n(n + 1)/2 − Σ LCP; FM counts, wavelet ranks and Aho–Corasick matches against brute force.

### 3.2 calx

| Where | Content |
|---|---|
| `compute/backends/cpu` → `calx::cpu` | `target.hpp`: the target description of the host (CPUID: SSE4.1, AVX2, AVX-512, POPCNT, BMI2; threads) and of the build (`LANE_BYTES`, the widest vector the compiler was allowed). `lanes.hpp`: `Lanes<T, W>`, W lanes of T on GCC's vector extensions, with splat, load, store, arithmetic, max, min, compare, select, saturating add and subtract, shift by a lane, horizontal maximum. The first content of the CPU backend's native path. |
| `compute/scheduler` | `task_graph.hpp`: tasks with dependencies and costs, run on the `ThreadPool` by priority of upward rank (the critical path by cost, HEFT), with a trace (task, thread, start, end). `parallel_for.hpp`: a loop of irregular iterations in chunks taken from an atomic counter. |
| `numerics/linalg` | `eigen.hpp`: the symmetric eigenproblem by Householder tridiagonalisation and the implicit QL algorithm (PCA of genotypes, the spectral decomposition of reversible rate matrices). |
| `numerics/special_functions` | `gamma.hpp`: log Γ (Lanczos), the regularised incomplete gamma functions and their inverse (χ² tails and quantiles, discrete-Γ rates); the regularised incomplete beta function. |

## 4. The Genomics Domain

C++20, header-only, each module a directory `<module>/` with `README.md`, `CMakeLists.txt` (the target `genomics::<module>` and its checks), `include/genomics/<module>/`, `tests/`. Namespace `genomics::<module>`.

### 4.1 Conventions Shared by Every Module

- **Codes.** `genomics::sequence::Code` is `std::uint8_t`; DNA is A = 0, C = 1, G = 2, T = 3, N = 4 (any other letter); the complement of x < 4 is 3 − x. Protein is the order of the NCBI matrices, `ARNDCQEGHILKMFPSTWYVBZX*` (0 to 23). In aligned rows, `GAP` = 0xFF.
- **Coordinates.** 0-based, half-open [begin, end) everywhere inside; 1-based only at the boundaries of formats that say so (SAM, VCF, GenBank, GFF3). A strand is `bool reverse`.
- **Scores and penalties.** A scoring is a substitution matrix and affine gaps: a gap of length L costs `open + L · extend` (the convention of BLAST and minimap2), so `open = 0` is a linear gap. Scores are maximised; penalties (the wavefront algorithm's) are minimised, and the domain converts between them.
- **Reading formats never throws on content.** A reader returns what it could read and a list of problems, each with its line and why. Invalid arguments throw `std::invalid_argument`.
- **Randomness** is a `calx::random::Stream(seed, a, b, c)`: a names the purpose (a constant per use), b and c the object (read index, replicate, chromosome). Never a global generator.
- **Parallelism** takes a `calx::scheduler::ThreadPool&` and is deterministic: shards write to their own slots and are merged in index order.
- **Comments and documentation.** Each header opens with what it computes, the algorithm and its source (paper, equation or algorithm number), complexity, and the limits of its accuracy. Each module's README gives its purpose, status, content by header, dependencies, checks and references. Headings are declarative, numbered, academic.
- **Checks** are programs on `calx::verification::Suite` (`require` for conditions, `check` for an error against a limit), named by the sentence they verify, registered by `calx_add_check` with the name `check_genomics_<module>[_<topic>]`.

### 4.2 The Modules, in Their Order

| Module | Content | Depends on |
|---|---|---|
| `sequence` | Alphabets (DNA, IUPAC, protein); packed 2-bit DNA; reverse complement; the NCBI genetic codes, translation, open reading frames; FASTA, FASTQ and GenBank (locations as a small grammar: join, order, complement, partial ends) read from text; k-mers (canonical, rolling), ntHash, minimizers (w, k) and syncmers; composition, GC and its skew (the origin of replication), CpG observed/expected, entropy, codon usage | calx::verification (checks) |
| `simulate` | Genomes (i.i.d. or with repeats), diploid mutation with a truth set of SNVs and indels, short reads with position-dependent substitution errors and Phred qualities, long reads with indel-dominated errors and gamma lengths; evolution of a root sequence with indels (truth alignments) | sequence, calx::random |
| `align` | Biological scoring on morphe's generic alignment (the reference): Gotoh with traceback; Hirschberg; banded; Myers' bit-vector edit distance; Farrar's striped Smith–Waterman on lanes; the wavefront algorithm (gap-affine, with traceback); X-drop extension; partial-order alignment (a sequence to a DAG, and the graph's consensus); substitution matrices (BLOSUM, PAM); Karlin–Altschul statistics (λ, K, H, E-values); parametric alignment of DNA and protein on morphe's polytope algebra; the pair HMM in a log-probability semiring (forward, backward, posterior match probabilities, maximum expected accuracy) | sequence, morphe::formal, calx::cpu |
| `index` | The FMD-index of a genome and its reverse complement on morphe's SA-IS: interleaved occurrence blocks, sampled suffix array, bidirectional extension, super-maximal exact matches; the minimizer index; k-mer counting by partitioned hash tables, the k-mer spectrum and genome size; Bloom filters and count-min sketches | sequence, morphe::formal, calx::scheduler |
| `sketch` | MinHash bottom-k and FracMinHash sketches, Jaccard, containment, the Mash distance and ANI; HyperLogLog; all-against-all distance matrices in parallel | sequence, calx::scheduler |
| `mapping` | A minimap2-like mapper: minimizer anchors, chaining by dynamic programming with Li's gap cost, primary and secondary chains, mapping quality, base-level alignment of the gaps between anchors (wavefront or banded Gotoh) and of the ends (X-drop / Z-drop); SAM; presets for short and long reads; evaluation against the simulator's truth | sequence, index, align, calx::scheduler |
| `assembly` | The de Bruijn graph of canonical k-mers (bidirected), compaction into unitigs, cleaning (tips, bubbles, low coverage), Eulerian paths (Hierholzer) and the number of Eulerian circuits (the BEST theorem, exact by morphe's integers), contigs and their statistics (N50, NG50, L50); an overlap–layout–consensus path for long reads (minimizer overlaps, the string graph and its transitive reduction, POA consensus) | sequence, index, align, morphe::arithmetic |
| `variant` | Pileups from alignments; genotype likelihoods and the Bayesian caller (Li 2011); VCF; normalisation of indels (left-align and trim); read likelihoods by the pair HMM; evaluation against a truth set | sequence, align |
| `hmm` | Hidden Markov models over semirings: Viterbi, forward, backward, posterior decoding, Baum–Welch; the CpG-island model; profile HMMs from alignments, with Viterbi and forward scoring against a null model | sequence |
| `msa` | Progressive multiple alignment: k-mer distances, a guide tree by UPGMA, profile–profile alignment with affine gaps, refinement; the centre-star method; POA alignment; sum-of-pairs and column scores against a truth | sequence, align |
| `expression` | Pseudoalignment to a transcriptome by k-mer equivalence classes; the EM of abundances, effective lengths, TPM; size factors by the median of ratios; the negative binomial GLM, Wald tests, Benjamini–Hochberg | sequence, index, calx::special_functions |
| `annotation` | Features from GenBank and GFF3, extracted and translated (checked against the record's own translations); position weight matrices, information content, scanning; motif discovery by Gibbs sampling and by MEME's EM; maximal repeats by the suffix array; a self-trained gene finder (periodic Markov chains) | sequence, hmm, morphe::formal |

## 5. The Evolution Domain

| Module | Content | Depends on |
|---|---|---|
| `phylogeny` | Trees (rooted, unrooted), Newick, splits, Robinson–Foulds, majority-rule consensus; distances (p, JC69, K80, F81, TN93, with Γ); UPGMA and neighbour joining; substitution models JC69 to GTR with +Γ (Yang's discrete rates) and +I, P(t) by spectral decomposition; Felsenstein's pruning, branch lengths by Newton's method, NNI search; Fitch and Sankoff parsimony; the bootstrap in parallel; sequences simulated along a tree; root-to-tip regression of dated tips | genomics::sequence, genomics::msa, calx::linalg, calx::special_functions, calx::random, calx::scheduler |
| `population` | Genotype matrices from VCF; allele frequencies; Hardy–Weinberg (χ² and Wigginton's exact test); linkage disequilibrium (D, D′, r²); diversity (π, θ_W, Tajima's D, Fay and Wu's H), site frequency spectra, F_ST (Hudson, Weir–Cockerham); the Kingman coalescent with mutation, growth and recombination (Hudson's algorithm), with its genealogies; Wright–Fisher with selection, mutation and migration (Kimura's fixation probability); PCA (Patterson, Price and Reich) and admixture by EM | genomics::variant, calx::linalg, calx::random |
| `selection` | dN/dS by Nei–Gojobori; the McDonald–Kreitman test | genomics::sequence |

## 6. The Genome Atlas

`applications/genome-atlas/`, port 8890. Main colour: mulberry (#A8577E) on a blush paper, provisionally; the four bases in softened conventional colours (A sage, C dusty blue, G ochre, T rose).

### 6.1 Views

| View | What it shows |
|---|---|
| Genome | A genome as a circle and as linear tracks: genes and features from GenBank, ORFs in six frames, GC content and its cumulative skew (the origin and terminus of replication), CpG islands by the HMM, repeats, motif hits; sequence search by the FMD-index |
| Align | The DP of two sequences as a field: the matrices of Gotoh, the traceback, the band, the wavefronts of WFA rising by score, Myers' bit vectors, Farrar's stripes; the semiring changed (score, count, posterior of the pair HMM); the parametric polygon and its normal fan, a click on a region giving its optimal alignment |
| Index | Suffix array, LCP, the BWT matrix of rotations, the LF map, backward search narrowing its interval, SMEMs of a read; Lyndon factorisation and the bijective BWT |
| Reads | Reads simulated from a genome (or given), mapped: anchors and chains as a dot plot, coverage, pileup, called variants against the truth, mapping accuracy by MAPQ |
| Assembly | The k-mer spectrum, the de Bruijn graph drawn by force, unitigs, cleaning step by step, Eulerian paths and their number, contigs against the reference, N50 |
| Trees | Sequences aligned (MSA view), distances, trees by UPGMA, NJ and maximum likelihood with bootstrap support; trees of genomes from sketches; a tanglegram of two trees |
| Populations | The coalescent genealogy with its mutations, the site frequency spectrum and Tajima's D; Wright–Fisher drift; PCA and admixture of a structured sample |
| Run | The plan of a job as a task graph; shards on workers and threads on a timeline; targets and lanes; throughput by thread count; the striped kernel scalar against SIMD; WebGPU against the CPU |

### 6.2 The Engine and the Cluster

- **The engine** (`engine/genome.cpp`, built by the recipe into `data/cache/engines/genome/`) speaks JSON lines: `--serve` reads one request per line and answers one line per request, so the server keeps workers alive between requests.
- **Shardable operations** have three parts: `plan` (the shards of a job), `shard` (one shard's partial result), `merge` (partial results combined in shard order). Simulated inputs are regenerated from their seed by each worker, so a shard is a few numbers, not data.
- **The cluster** (`cluster.py`, an edge) holds workers: local engine processes, and remote ones by a command (`ssh host engine --serve`) from a configuration the user writes. Workers report their targets (`op: target`) and a calibration (`op: benchmark`); shards are pulled by free workers, sized by measured throughput, and every shard's start, end, worker and inner trace is kept for the Run view.

## 7. Data

Fetched into `data/cache/` by manifests (NCBI's E-utilities, sequences in the public domain as NCBI distributes them): phiX174 (NC_001422.1), phage lambda (NC_001416.1), the human mitochondrion (NC_012920.1), SARS-CoV-2 Wuhan-Hu-1 (NC_045512.2), *Escherichia coli* K-12 MG1655 (NC_000913.3), and the mitochondrial genomes of the great apes for the trees. Everything else the atlas shows can be simulated from a seed, so that it runs without the network.

## 8. Open Questions

1. The application's name and colour are provisional (the user proposed "genom-lab"; the infrastructure retires the word laboratory).
2. Whether the general HMM machinery and the negative binomial GLM move to the inference foundation when it exists (CLAUDE.md §8, open decision 3).
3. Real population data (1000 Genomes) by range requests on BGZF and tabix, and real reads from ENA.
4. A CUDA backend for calx when a toolkit is installed; the WebGPU kernels of the page as the first accelerator path.
