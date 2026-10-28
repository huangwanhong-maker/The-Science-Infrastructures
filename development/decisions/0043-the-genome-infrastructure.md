# 0043. The Genome Infrastructure: the Genomics and Evolution Domains and the Genome Atlas

- **Status:** Accepted; amended by [0045](0045-the-genome-workbench.md) (the Genome Atlas renamed the Genome Workbench, with a Toolbox and the user's files); amended by [0047](0047-group-grammars-and-advanced-genome-analysis.md) (group grammars, and the survey's new modules)
- **Date:** 2026-10-24
- **Decided by:** the user asked for "a tremendous new infrastructure and application" for genome analysis, proposing the name "genom-lab": "various algorithms implemented at low level ... classical algorithms ... also advanced high speed and high quality algorithms ... distributed parallel structures ... as some high speed algorithms support distributed parallel computing in heterogeneous devices", and left Claude "fully free to design and implement", asking for high quality, creativity and effort. While the work went on they told another session that domains must reuse the foundations and never reimplement them. Claude designed the domains, the foundations' parts, the engine and the application, wrote the core and the application, had seven builders write modules to one contract in parallel, and wrote this record.

## 1. Context

Genomics is the README's 1.5, a field of height 0: "sequence alignment, k-mers, genome assembly on de Bruijn graphs, variant analysis". Evolution is its 3.2, standing on genomics: "phylogeny, population genetics, selection". The analysis of genomes the user asked for reaches both: alignment and assembly are genomics, but the trees of genomes and the genealogies of populations are evolution.

Much of what genome analysis needs names no object in nature: the combinatorics of words and their indexes (suffix arrays, the Burrows–Wheeler transform, FM-indexes), the alignment of words over semirings, the polytope algebra of parametric alignment, vector lanes, task graphs, symmetric eigenproblems, the gamma family of functions and its variates. The user's rule, and CLAUDE.md §3.7, put these in the foundations.

The user's word "lab" is retired with "laboratory" (CLAUDE.md §1; ADR 0026).

## 2. Decision

1. **Two domains.** `domains/genomics/` (a field of height 0, on calx and morphe) and `domains/evolution/` (a field of height 1, on genomics), C++20 and header-only from the start, each module a directory with its CMake target `genomics::<module>` or `evolution::<module>`:
   - genomics: `sequence`, `simulate`, `kmers`, `align`, `index`, `sketch`, `mapping`, `assembly`, `variant`, `hmm`, `msa`, `expression`, `annotation`;
   - evolution: `phylogeny`, `population`, `selection`.
2. **An alignment is a path; a scoring scheme is a weight in a semiring; an algorithm is a schedule.** morphe's generic dynamic programme over the edit graph defines every alignment, in any semiring (tropical for scores, counting for numbers of alignments, with multiplicities for numbers of optimal ones, the polytope algebra for the parametric family). The domain's fast algorithms (Gotoh, Hirschberg, banded, X-drop, Myers' bit vectors, Farrar's striped lanes, the wavefront algorithm and BiWFA, partial orders) are schedules of the same sum, checked against it. This is calx's separation of what from how, in sequence analysis.
3. **The foundations' parts** (each recorded in its own log):
   - **morphe `formal`** (its decision 0005), its first slice: words and their morphisms, borders, Lyndon factorisation; suffix arrays by SA-IS and LCP; the BWT, bijective and extended; rank and select, wavelet matrices; FM-indexes; Aho–Corasick automata; semirings; the alignment of words over semirings; the polytope algebra. And in `algebra`, exact determinants (Bareiss); in `combinatorics`, Hierholzer's circuits, the matrix-tree and BEST theorems, which joins `combinatorics` to `algebra` in morphe's order.
   - **calx** (its decision 0006): the CPU backend's target description and vector lanes (`calx::cpu`), task graphs by rank and dynamic loops in the scheduler, the symmetric eigenproblem, LU and determinants, the gamma family, gamma, beta and binomial variates, and the first content of `optimisation`, Brent's minimiser.
4. **Truth by simulation.** The domain's `simulate` module makes genomes with repeats, diploid individuals with their variants and liftovers, reads with their origins and errors, and families with their true alignment, each draw keyed by purpose and object, so that read i is the same whichever process makes it.
5. **Real genomes** come through the data library: fourteen RefSeq genomes by NCBI's E-utilities (`data/registry/ncbi-genomes.json`), fetched into the cache. The GenBank reader is held to RefSeq's own `/translation` on every one of their 4,522 CDSs.
6. **The Genome Atlas**, `applications/genome-atlas/`, on port 8890, in the pattern of the other applications, with three new parts:
   - **an engine that stays alive.** `engine/` on both domains, built by the recipe with -O3 -march=native into `data/cache/engines/genome/`, a unit for each group of operations, compiled one at a time and only when its includes change (one unit holding both domains needed more memory than WSL gave without failing), speaks JSON lines; with --serve it is a worker answering requests until its input ends, keeping its genomes in memory. Its staleness is judged by the files the compiler reports it read.
   - **a cluster** (`cluster.py`, an edge): local workers, and remote ones started by a command the user writes (ssh), requests to the least busy (or to the one holding the genome), and **jobs** in three parts, `plan`, `shard` and `merge`: shards are a few numbers, pulled by free workers, merged in their order, so that the answer is the same however the job is cut and wherever it runs; every shard's worker, times and inner task-graph trace are kept.
   - **views:**
     - **Genome:** the circular map and a zoomable strip, features, ORFs in six frames, composition, the cumulative skew and the origin of replication, a card per feature with its protein checked against RefSeq's; genes predicted, CpG islands and repeats found from the sequence alone, the repeats drawn as chords.
     - **Align:** the edit graph with the footprint of each schedule (the matrix, the band, the WFA's wavefronts), every schedule timed and held to the reference, the other semirings, the parametric polygon and its fan, and many sequences by partial order.
     - **Index:** suffix arrays, the BWT and the FMD-index, SMEMs.
     - **Reads:** a diploid read, mapped and called, the calls against the truth, a pileup browser, mapping quality measured on a job across the workers, SAM and VCF.
     - **Assembly:** the k-mer spectrum, the de Bruijn graph cleaned step by step, the contigs on the genome, Eulerian trails and the BEST theorem.
     - **Expression:** pseudoalignment and the EM against the truth, and a study of two conditions (MA plot, volcano, dispersions, p-values).
     - **Trees:** a gene of chosen genomes aligned and its tree by distances or likelihood with the bootstrap, dN/dS between every pair, and whole genomes by their sketches.
     - **Populations:** the coalescent, Wright–Fisher, PCA and admixture.
     - **Run:** the machine, the workers, jobs on a timeline, threads, and the GPU of the browser (a WebGPU kernel of Smith–Waterman checked against the CPU's three schedules).

   Its main colour is mulberry (#A8577E) on blush paper, provisionally, until the user chooses one; the four bases keep softened conventional colours.
7. **The name.** The application is the Genome Atlas, in the family of the infrastructure's atlases, and not "genom-lab", the word laboratory being retired; the user may rename it.

## 3. Alternatives Considered

- **Phylogeny and population genetics inside genomics.** Rejected: the README's map gives them to evolution, which stands on genomics; a domain of their own keeps the order and lets evolution grow toward selection and systematics.
- **The general algorithms in the domain** (suffix arrays, alignment over semirings, SIMD lanes, task graphs). Rejected by the user's rule and CLAUDE.md §3.7: linguistics, the next domain, needs the words and their indexes too.
- **Intrinsics in the domain.** Rejected: the domains never address a device; calx's lanes are compiled for whichever device the build targets.
- **One engine process per request** (the History Complex's pattern). Rejected for this application: a genome parsed once serves the many requests of a view, and a worker kept alive is also what a cluster needs; the protocol is the same pipes, local or over ssh.
- **Distributed execution inside calx.** Rejected: calx is pure and has no network. The core gives deterministic shards and traces; the edge carries them.
- **A borrowed aligner, index or assembler** (minimap2, BWA, WFA2-lib, SPAdes). Licence-compatible ones exist (MIT), and they would be faster to adopt; but writing the methods is how they are understood (the map of the mathematics, §10.3), and the infrastructure's own interfaces, checks against references and visual explanations are the point of the work. They remain references the checks can be held to.

## 4. Consequences

- **The order of the domains:** genomics at height 0, evolution at height 1. The root `CMakeLists.txt` builds them after history.
- **The applications are seven.** CLAUDE.md §1, §3, §5, §7 and §9, the README, `docs/architecture.md`, `docs/mathematics.md`, the domains' and the applications' indexes name the new domains and application.
- **The design note** `development/ideas/2026-10-24-the-genome-infrastructure.md` is promoted into this record, the domains' READMEs and their theory notes.
- **What the builders found.** Each module was held to published examples and to brute force, and its README states where it departs from the method it follows; the departures that change the method are recorded there, not here (for example the integer gap cost of chaining, the peeled gap breakpoints of BiWFA, the Poisson-mixture fit of the k-mer spectrum).
- **Open:**
  - the application's name and colour, the user's to choose;
  - where general sketches (MinHash, HyperLogLog, Bloom filters) and the HMMs belong when a second domain needs them (the inference foundation, CLAUDE.md §8, open decision 3);
  - **trees and the methods that build them from distances.** The `Tree` structure, splits, UPGMA and neighbour joining name no object in nature, and two modules need them: evolution's `phylogeny`, and genomics' `msa`, which cannot depend on evolution and so keeps a guide tree of its own. They belong in a foundation (morphe's combinatorics for trees and splits; the inference foundation for clustering), to be moved there with the user's agreement;
  - the statistics written in the domain because no foundation holds them yet: the EM of mixtures, the negative binomial GLM, the shrinkage of dispersions, Benjamini–Hochberg (expression), Fisher's exact test (variant, selection), the hypergeometric law; they move to the inference foundation when it exists;
  - waking calx's pool costs about a millisecond under WSL, so loops of many short rounds run serially below a threshold (the EM of expression); a pool that spins longer, or adapts its spin to the gaps it sees, would serve every caller;
  - real population data (1000 Genomes) by range requests on BGZF and tabix, and real reads from ENA;
  - a CUDA backend in calx when a toolkit is installed; the browser's WebGPU as the first accelerator of the atlas.
