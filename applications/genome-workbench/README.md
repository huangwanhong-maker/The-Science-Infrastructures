# Genome Workbench

*Status: draft (2026-10-24). Serves learning, use and research.* Genomes drawn and computed: a genome's map and composition; the alignment of two sequences as a path in a semiring, with the footprint of every algorithm that finds it; the indexes that make a genome searchable; reads simulated, mapped, called and assembled; transcripts quantified and tested; the trees of genomes, the selection on their genes and the genealogies of populations; and the engine's own work across processes and machines ([ADR 0043](../../development/decisions/0043-the-genome-infrastructure.md)); and a toolbox of instruments for the user's own sequences and files ([ADR 0045](../../development/decisions/0045-the-genome-workbench.md)), the Genome Atlas renamed. It composes the [genomics](../../domains/genomics/) and [evolution](../../domains/evolution/) domains, on calx and morphe.

## 1. Concepts

- **A shell built from registries.** Views, tools and commands are declared (`src/core/registry.js`), and the shell (`src/shell.js`) builds itself from them: the rail of views, the command palette (Ctrl K, or /), which finds any view, tool, genome, file or command and reads `1200-4800` as a position, a word of DNA as a motif and a name as a gene; the dataset switcher; the activity tray of requests and jobs in flight; the status bar; day and night themes. A new view, tool or command is one declaration.
- **A kit, not raw elements.** Every control and surface comes from `src/ui/kit.js` (buttons, segmented choices, selects, number fields, sliders, switches, the sequence editor, drop zones, cards, tables, metrics, popovers, menus, dialogs, toasts), drawn from the tokens of `styles.css`, so the views read as one instrument.
- **The Toolbox: a tool is a declaration.** `{id, title, group, icon, summary, refs, inputs, example, run, render}` in `src/tools/`: its `inputs` are a typed schema from which `src/ui/form.js` builds the form, its `run` calls the engine, its `render` uses the shared renderers of `src/ui/results.js` (alignments, tables with CSV, trees with Newick, figures, downloads). Its examples are real sequences, read on request from the RefSeq genomes the recipe fetched.
- **Your files.** Files dropped anywhere on the page (FASTA, FASTQ, GenBank, VCF, Newick, SAM, plain or gzipped) are kept in `data/cache/genome-workbench/files/` under their content's hash, recognised and summarised (`files.py`); FASTA and GenBank join the catalogue as genomes, and every file can be handed to the tools that take its kind. A request names a file as `file:<id>`, and only the server turns it into a path for the engine.
- **One engine, kept alive.** Everything the page shows is computed by the engine: the two domains in C++, compiled for this processor (-O3 -march=native, so calx's vector lanes are as wide as it allows). The server keeps engines alive as **workers** between requests, so a genome is parsed once and then serves every request of a view.
- **Jobs travel as a few numbers.** A large computation is a **job** in three parts: a *plan* cuts it into shards; any worker computes a *shard* (simulated reads are regenerated from their seed, so a shard is "reads 40,000 to 49,999"); the *merge* combines the shards in their order. The answer never depends on how the job was cut or which worker took which shard; the Engine view draws who did what, when, and on which threads.
- **Workers can be elsewhere.** A worker is a process speaking JSON lines on its pipes; started through ssh on another machine, it is a remote worker, with its own processor, threads and lanes. Free workers pull shards, so a faster machine takes more.
- **The browser's GPU is a device too.** The Engine view's Smith–Waterman kernel in WGSL scores thousands of alignments at once on the GPU through WebGPU, and the page checks its scores against the engine's, which must agree exactly.
- **Truth is shown beside the answer.** Simulated genomes carry their repeats, individuals their variants, reads their origins; real genomes are held to their records (every protein against RefSeq's own translation).

## 2. Views

| View | What it shows |
|---|---|
| Genome | The circular map (features by strand, GC against its mean, GC skew, the cumulative skew whose extremes mark a bacterial chromosome's origin and terminus) and a zoomable strip (features in rows, ORFs in six frames, composition, the letters at base resolution); search by sequence or gene; a card for the genome and for each feature, its protein checked against RefSeq's and its codons; from the sequence alone, genes predicted by a self-trained finder (with their agreement with the record), CpG islands by the HMM and repeats by the suffix array, drawn on the strip and as chords on the circle |
| Align | Two sequences aligned by any schedule (Gotoh, Hirschberg, banded, WFA, BiWFA, Myers, striped) in four modes; the edit graph as a field with its path and the schedule's footprint (the band, the WFA's wavefronts); every schedule's value, work and speed on the pair, held to Gotoh's and to morphe's generic sum; other semirings (the number of alignments, the pair HMM's posterior); the parametric polygon and its normal fan; many noisy reads fused into a partial order, its graph, consensus and multiple alignment |
| Index | Suffix arrays and the BWT of a text, the LF map and backward search; the FMD-index of a genome and the super-maximal exact matches of a read; Lyndon words and the bijective transform |
| Reads | A diploid made from a genome, read, mapped and called: the calls against the truth by type (precision, recall, genotypes), a pileup browser over any call or missed variant, Li's curve of mapping accuracy against mapping quality from a job across the workers, the SAM and VCF as written |
| Assembly | The k-mer spectrum with its error cut and genome size, the de Bruijn graph cleaned step by step (tips, bubbles, low coverage), the contigs placed on the genome beside its repeats; the de Bruijn graph of a text, an Eulerian trail and the number of texts by the BEST theorem |
| Expression | A simulated transcriptome's RNA-seq sample pseudoaligned and quantified by the EM, against the truth and an equal split, with equivalence classes and the bootstrap; a study of two conditions: MA plot, volcano, dispersions shrunk toward their trend, p-values, power and false discoveries |
| Trees | A gene of chosen genomes aligned (progressive), distances, its tree by UPGMA, neighbour joining or maximum likelihood (JC69 to GTR, +Γ) with the bootstrap; dN/dS between every pair on codon alignments; whole genomes by their sketches |
| Populations | The coalescent's genealogies, the site frequency spectrum, Tajima's D; Wright–Fisher; principal components and admixture |
| Toolbox | Every instrument on the user's own sequences and files, by group (sequences, alignment, annotation, indexes, reads, assembly, expression, phylogeny, populations, statistics): a form built from the tool's inputs, its results, and its earlier runs |
| Your data | The user's files, recognised and summarised: opened as genomes, handed to a tool, previewed, removed; the catalogue's genomes |
| Engine | The machine and what the engine was compiled for; the workers; a job on a timeline with each shard's threads; speed-up by threads; the GPU against the CPU's three schedules of Smith–Waterman |

## 3. The Protocol

The server (`serve.py`, port 8890) answers:

| Request | Answer |
|---|---|
| `GET /api/engine` | The engine's build (current or stale, the target it was compiled for) and the cluster's workers with their targets and work |
| `GET /api/genomes` | NCBI's genomes of `data/registry/ncbi-genomes.json` (fetched or not) and the simulated presets |
| `POST /api/run` `{op, genome, ...}` | One operation on one worker; a genome is `ncbi:<accession>`, a simulated preset `sim:<name>`, `{"simulated": spec}` or `{"text", "name"}` |
| `POST /api/job` `{op, genome, ...}` | A job: `<op>.plan`, `<op>.shard` on the workers, `<op>.merge`; the answer carries `_job`, its trace |
| `POST /api/restart` | Stop the workers; they start again from the present build |
| `GET /api/files`, `POST /api/files?name=…`, `DELETE /api/files/<id>` | The user's files: their list with kinds and summaries; one added (the body is the file, gzipped or not); one removed |
| `GET /api/files/<id>/head`, `GET /api/files/<id>/text` | A file's first lines, or its text (to 16 MB) |
| `GET /setup`, `POST /setup/run` | The recipe: the engine and the genomes |

The engine (`engine/`: the loop in `genome.cpp`, each group of operations in a unit of its own, `ops_*.cpp`, compiled apart and linked; the helpers they share in `shared.hpp`, the cores that a view's operation and a tool's share in `shared_*.hpp`, and the reader of a tool's inputs in `inputs.hpp`) reads one JSON request a line and writes one answer a line: `{"id", "ok": true, "seconds", "result"}` or `{"id", "ok": false, "error"}`. Its operations:

| Operation | Computes |
|---|---|
| `target` | The host's processor, threads and instruction sets; what the build was compiled for; the operations |
| `genome.summary`, `genome.tracks`, `genome.feature`, `genome.gene`, `genome.search` | A genome's records; its tracks over a window in bins; one feature in full; a gene by name; a pattern on both strands |
| `genome.annotate` | Genes predicted, CpG islands or repeats of a contig, from the sequence alone (kept by the worker) |
| `align.random`, `align.pair`, `align.compare`, `align.semiring`, `align.parametric` | Two relatives; one schedule's alignment with its field and footprint; every schedule on the pair; the other semirings; the parametric polygon and its fan |
| `align.batch_input`, `align.batch` | Reads to score against a query; their local scores by scalar rows, stripes or one target to a lane, on the worker's threads (the CPU's side of the GPU's comparison) |
| `align.poa` | Noisy long reads of a sequence fused into a partial-order graph, its consensus and multiple alignment |
| `index.word`, `index.count`, `index.sample_read`, `index.smems` | A text's suffix array and BWT; counts and places in a genome's FMD-index; a read and its SMEMs |
| `reads.sample.plan`, `.shard`, `.merge` | A job of simulated reads and their measured errors and qualities |
| `reads.map.plan`, `.shard`, `.merge` | A job of reads simulated and mapped, counted by mapping quality against their truth |
| `reads.call`, `reads.window` | A diploid read, mapped and called, against its truth (kept by the worker); the reads, calls and truth of a stretch |
| `assembly.run`, `assembly.eulerian` | Reads assembled on the de Bruijn graph with every step of cleaning; a text's k-mer graph, a trail and the BEST counts |
| `expression.quantify`, `expression.differential` | An RNA-seq sample quantified against its truth; a study of two conditions tested gene by gene |
| `trees.gene`, `trees.dnds`, `trees.sketches` | A gene's alignment and tree with the bootstrap; dN/dS between pairs; whole genomes by Mash distances and neighbour joining |
| `pop.coalescent`, `pop.wright_fisher`, `pop.structure` | Genealogies and their statistics; allele trajectories; principal components and admixture |
| `run.scaling` | One computation on 1, 2, 4 … threads, with each run's task-graph trace |
| `tool.seq.stats`, `.revcomp`, `.translate`, `.orfs`, `.find`, `.kmers`, `.codons` | The toolbox on the user's sequences (`text`, or `fasta_file`: FASTA, FASTQ or GenBank, told by their first letters): composition; the other strand; translations in 1, 3 or 6 frames; ORFs; IUPAC motifs on both strands by one Aho–Corasick automaton; the k-mer spectrum; codon usage, RSCU and N_c |
| `tool.annot.genes`, `.cpg`, `.repeats` | Genes, CpG islands or maximal repeats of the first sequence, from it alone |
| `tool.align.msa`, `tool.align.search`, `tool.align.dnds`, `tool.stats.evalue` | Many sequences aligned (progressive, centre star, partial order); a query against many, one target to a lane, with E-values; dN/dS of two coding sequences; λ, K, H, bits and E-value of a score |
| `tool.reads.qc`, `.map`, `.call`, `tool.assembly.reads`, `tool.expr.quant` | The user's reads (`reads_file`, read a record at a time): checked; mapped, with coverage and SAM; called, with VCF; assembled; quantified against `transcripts_file` |
| `tool.tree.build`, `tool.tree.read`, `tool.pop.vcf` | A tree of the user's sequences with the bootstrap; a Newick tree read and rooted; a VCF's samples on their principal components, with θ_W, π and Tajima's D |

Every tool reports in `reading` what it read of its input: from where, in what format, how many records and bytes, the reader's problems by line, and whether a bound stopped it.

## 4. Running It

```
python data/acquisition/recipe.py applications/genome-workbench/recipe.json   # the engine and the genomes, once
python applications/genome-workbench/serve.py [--workers 2]
```

then open http://127.0.0.1:8890/. The recipe can also be run from the page (the Engine view, or the Setup shown when the engine is missing). The engine needs g++ 11 or later, natively or in WSL; the genomes need the network once.

**Remote workers.** Build the engine on the other machine from the same sources, then write `data/cache/engines/genome/cluster.json`:

```json
{"workers": [{"name": "eu_region_server", "command": ["ssh", "eu_region_server", "~/The-Science-Infrastructures/data/cache/engines/genome/genome --serve"]}]}
```

and restart the workers (the Engine view). Remote workers take the shards of jobs on simulated or inline genomes; genomes read from this computer's files stay with the local workers.

## 5. Checks

```
node applications/genome-workbench/tests/run-tests.mjs
```

What is the workbench's own: how it writes lengths, positions, ticks and its address; that every tool is declared whole (a known group, inputs of known types, run and render), takes only kinds of file the library recognises, and is found by the files it takes; how the palette matches; and, when the engine is built, its answers: phiX174's summary and tracks, a search on both strands against brute force, a protein as RefSeq translates it, *Escherichia coli*'s origin of replication found from the cumulative skew within 500 bases of oriC and its genes found by the self-trained finder (over 85% of RefSeq's CDSs, over 95% of the predictions among them), every schedule of an alignment at Gotoh's optimum, the three batch kernels giving the same integers, a partial-order consensus within 2 edits of its truth, and a job whose answer is the same in 5 shards and in 13; and the toolbox's operations on phiX174's file and on reads drawn from it: its length and GC as counted on its letters, the reverse complement an involution, a CDS translated as its /translation, a motif counted on both strands as brute force counts it, every ORF whole codons from ATG to its stop, L − k + 1 k-mers, BLAST's gapped λ and K with an E-value that falls as the score rises, a tree's length kept by its midpoint root, and the reads checked and mapped back at their depth. The science is checked in the domains.

## 6. Design

Blush paper and plum ink by day, plum by night, mulberry (#A8577E) for choices, provisionally until the user chooses; the four bases in softened conventional colours (A sage, C dusty blue, G ochre, T rose); Source Serif 4; hairlines rather than boxes.
