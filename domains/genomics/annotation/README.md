# annotation: Features, Motifs, Genes and Repeats

*Status: draft (2026-10-24).* What a genome's sequence says about its parts: the features of its annotation as tracks to draw (from GenBank and GFF3), the motifs of DNA (position weight matrices, their information and the exact statistics of their scores, and their discovery by Gibbs sampling and by EM), the protein-coding genes of a prokaryote found from the statistics of codons the genome itself teaches, and its exact repeats on both strands. Serves use and learning; the theory is [theory/annotation.md](../theory/annotation.md).

## 1. Concepts

- **A feature is a location and a name.** A GenBank location is a word of a small grammar (join, complement, order, partial ends) that sequence::genbank reads; here its parts are put in the forward strand's order, a feature of a circular genome that reads past the end continues at the start, and its extent then wraps. GFF3 says the same with lines: lines that share an ID are one feature, Parent is the part-of relation, and the hierarchy is a directed acyclic graph (an exon may be part of several transcripts) whose cycles are faults.
- **A motif is a distribution, its score a log-odds.** A PWM's weight is log₂ of the motif's probability over the background's; the information of a column is its relative entropy, the mean score of a site. The distribution of the score of a random word is computed exactly on a lattice, so p-values are exact for the rounded matrix and bracketed for the real one.
- **Discovery is inference of a hidden alignment.** The Gibbs sampler draws it, site by site, from its conditional distribution; MEME's EM computes the posterior of every window and re-estimates the matrix. Both are deterministic here: the sampler's streams are keyed by its seed and restart, and EM has no randomness.
- **Coding has a period of three.** A three-periodic Markov chain against a homogeneous one tells a coding frame; the finder learns both from the genome's long ORFs, annotates, and re-learns from its own annotation, as GeneMarkS does in its far richer way.
- **Repeats are nodes of the suffix tree.** Maximal repeats are the left-diverse internal nodes, read from morphe's suffix array and LCP array as lcp-intervals; on both strands the text is the genome and its reverse complement.

## 2. Content

| Header | Content |
|---|---|
| [`features.hpp`](include/genomics/annotation/features.hpp) | `feature_tracks` of a `GenBank` record: `TrackFeature` (kind, key, name, strand, segments, extent wrapping the origin, translation checked, lane); `read_gff3` (`Gff`: features, entities with parents and children, roots, sequence regions, circular landmarks, FASTA), `percent_decode`, `segments` split at the origin |
| [`motifs.hpp`](include/genomics/annotation/motifs.hpp) | `Pwm` from counts or sites, `information_content`, `small_sample_correction`, `logo`; `score_distribution` (exact on a lattice: tails, bounds, thresholds, p-values of words); `scan` of both strands; `gibbs_motif` (OOPS, phase shifts, restarts on a pool); `meme_motif` (OOPS, ZOOPS, MEME's starting points), `motif_log_likelihood` |
| [`genes.hpp`](include/genomics/annotation/genes.hpp) | `PeriodicChain` (order up to 8, phases, Dirichlet back-off); `find_genes`: candidates on both strands with ATG, GTG and TTG starts, long-ORF training, rounds of re-training, starts by likelihood ratio and prior, overlaps resolved by weighted interval scheduling; `Gene`, `GeneReport` |
| [`repeats.hpp`](include/genomics/annotation/repeats.hpp) | `maximal_repeats` (with supermaximality), `supermaximal_repeats`, `longest_repeat`, on one strand or both, repeats with N left out |

## 3. Dependencies

`genomics::sequence` (codes, GenBank, FASTA, genetic codes, ORFs), `genomics::hmm` (declared for the HMM-based annotators to come; nothing here uses it yet), `morphe::formal` (the suffix array and LCP array), calx (`random`, `scheduler`, `special_functions`). The checks also use `genomics::simulate`.

## 4. Checks

`check_genomics_annotation_features`, `_motifs`, `_genes`, `_repeats`:

| Check | What it verifies |
|---|---|
| features | A GenBank record written in the check: a gene and CDS across the origin, a reverse-strand feature across it, a complement CDS, a spliced CDS with a wrong /translation (found at its line), an RNA, repeat and mobile element of each kind, named, segmented and laned without overlap; the GFF3 specification's canonical gene (fixture) with its hierarchy, discontinuous CDSs and phases; percent-encoding; the specification's circular example split at the origin; ##FASTA; eight kinds of fault at their lines; every NCBI genome in the cache drawn, all checked CDSs agreeing (4294 of 4294 in E. coli) |
| motifs | Information of known columns (2, 1, 0 bits); e(1) = 2, e(2) = 5/4, e(150) against Miller's formula; the lattice distribution against enumeration of every word, point by point, and real tails within its bounds; a reversed site found by scanning with p-value 4⁻¹²; Gibbs and MEME (OOPS) placing a planted motif exactly in all of 25 sequences, on one strand and both; the sampler equal on one thread and four; EM never lowering its log posterior; ZOOPS finding every planted site at a likelihood above the truth's |
| genes | A periodic chain recovered from its samples; a simulated genome of 554 kb with 500 genes from a known codon usage: stops with Sn 1.000, Sp 1.000, starts 0.95; E. coli K-12 against its 4300 CDSs: Sn 0.877 (0.954 of the 3952 of 80 codons or more), Sp 0.975, starts 0.74, in about a second |
| repeats | Maximal, supermaximal and longest repeats against brute force on 400 words, one strand and both, with N; planted exact copies on random strands found as one repeat each; E. coli's 237 maximal repeats of 100 bases or more on both strands (the longest 3027 bases), in 3.9 s |

## 5. References

- The INSDC Feature Table Definition, version 11; The Sequence Ontology, *Generic Feature Format Version 3*, version 1.26 (L. Stein, 2020).
- Stormo, G. D. (2000), DNA binding sites: representation and discovery, *Bioinformatics* 16, 16–23.
- Schneider, T. D., Stormo, G. D., Gold, L. and Ehrenfeucht, A. (1986), Information content of binding sites on nucleotide sequences, *Journal of Molecular Biology* 188, 415–431; Schneider, T. D. and Stephens, R. M. (1990), Sequence logos: a new way to display consensus sequences, *Nucleic Acids Research* 18, 6097–6100.
- Staden, R. (1989), Methods for calculating the probabilities of finding patterns in sequences, *Computer Applications in the Biosciences* 5, 89–96; Touzet, H. and Varré, J.-S. (2007), Efficient and accurate P-value computation for position weight matrices, *Algorithms for Molecular Biology* 2, 15.
- Lawrence, C. E., Altschul, S. F., Boguski, M. S., Liu, J. S., Neuwald, A. F. and Wootton, J. C. (1993), Detecting subtle sequence signals: a Gibbs sampling strategy for multiple alignment, *Science* 262, 208–214.
- Bailey, T. L. and Elkan, C. (1994), Fitting a mixture model by expectation maximization to discover motifs in biopolymers, *Proceedings of ISMB* 2, 28–36; (1995), Unsupervised learning of multiple motifs in biopolymers using expectation maximization, *Machine Learning* 21, 51–80.
- Borodovsky, M. and McIninch, J. (1993), GENMARK: parallel gene recognition for both DNA strands, *Computers & Chemistry* 17, 123–133; Besemer, J., Lomsadze, A. and Borodovsky, M. (2001), GeneMarkS: a self-training method for prediction of gene starts in microbial genomes, *Nucleic Acids Research* 29, 2607–2618; Delcher, A. L., Harmon, D., Kasif, S., White, O. and Salzberg, S. L. (1999), Improved microbial gene identification with GLIMMER, *Nucleic Acids Research* 27, 4636–4641.
- Gusfield, D. (1997), *Algorithms on Strings, Trees, and Sequences*, Cambridge University Press, §7.12; Abouelhoda, M. I., Kurtz, S. and Ohlebusch, E. (2004), Replacing suffix trees with enhanced suffix arrays, *Journal of Discrete Algorithms* 2, 53–86.

## 6. Open Questions

1. **Starts.** The finder chooses starts by likelihood and a start-codon prior (74% right on E. coli); a ribosome-binding-site motif learnt by `gibbs_motif` upstream of the predicted starts, as GeneMarkS does, is the natural next step.
2. **A gene HMM.** A hidden semi-Markov model of genes and intergenic regions with length distributions (GeneMark.hmm) on `genomics::hmm`, which is why the module is declared.
3. **The number of sites.** MEME's ZOOPS by maximum likelihood takes weak sites where none was planted (precision 0.70 in the check); MEME chooses the number of sites by the E-value of the motif, and erases found sites to look for the next motif; neither is here.
4. **Genes across the origin** of circular genomes, and eukaryotic genes (splicing) are out of reach of this finder.
5. **Approximate repeats** (diverged copies, as `simulate` plants them) need seeds and extension (the align and index modules), beyond exact repeats.
