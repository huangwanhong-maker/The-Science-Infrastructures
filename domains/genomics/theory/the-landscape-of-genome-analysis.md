# The Landscape of Genome Analysis: a Survey of Its Algorithms, and What the Domains Hold

*Status: draft, 2026-10-24; revised 2026-11-01 for [ADR 0064](../../../development/decisions/0064-abnormal-genomes-their-description-and-expression.md).*

This survey serves learning and planning. It maps the problems of genome analysis to the algorithms that solve them, classical and advanced, with their primary sources. It says what the genomics and evolution domains and the foundations under them already hold, what was added in response to the user's request of 24 October 2026 ([ADR 0047](../../../development/decisions/0047-group-grammars-and-advanced-genome-analysis.md)) and in response to their request of 1 November 2026 for abnormal genomes, their international description and the prediction of expression (ADR 0064, sections 10–12), and what is planned.

**Status words.** Every entry has one of these:
- **held**: implemented and checked;
- **added**: implemented in the round of ADR 0047;
- **added (0064)**: implemented in the round of ADR 0064;
- **planned**: named here, not yet written.

**Placement.** Follows the repository's reuse rule:
- a method that names no object in nature lives in a foundation, morphe (exact) or calx (numerical);
- the domains use it.

**Citations.** These are the sources we are sure of. Where a detail is not certain, the author and year are given alone.

## 1. Sequences and Their Statistics

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Formats and codes | FASTA, FASTQ (a line at a time), GenBank with its location grammar; NCBI's 25 genetic codes | NCBI | genomics `sequence` | held; line-at-a-time readers added |
| Composition | GC, skew, CpG observed over expected, entropy; the cumulative skew and the origin of replication | Lobry 1996; Grigoriev 1998 | genomics `sequence`, the Workbench's engine | held |
| k-mers | counting, spectra, canonical k-mers; their number up to the strand group by Burnside's lemma | Burnside | genomics `kmers`, `grammar` | held; Burnside counts added |
| Sampled seeds | minimizers | Roberts et al. 2004 | genomics `index` | held |
| | syncmers | Edgar 2021 | genomics `sequence` | held; their density checked exactly here |
| | strobemers | Sahlin 2021 | genomics `kmers` | added |
| Spectral analysis | the Fourier transform of the indicator sequences (Voss), period-3 of coding DNA | Voss 1992; Tiwari et al. 1997; Anastassiou 2001 | genomics `signal` on calx's FFT | added |
| Symmetry of the statistics | which elements of the bases' and strands' group a genome's k-mer statistics keep (Chargaff's second rule) | Rudner, Karkas and Chargaff 1968 | genomics `grammar` | added |

## 2. Pairwise Alignment

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Global and local alignment | dynamic programming over the edit graph, defined once in any semiring | Needleman and Wunsch 1970; Smith and Waterman 1981; Gotoh 1982 | morphe `formal/alignment`; genomics `align` | held |
| Linear space | divide and conquer | Hirschberg 1975; Myers and Miller 1988 | genomics `align` | held |
| Restricted search | bands, X-drop | Chao, Pearson and Miller 1992; Zhang et al. 2000 | genomics `align` | held |
| Bit parallelism | Myers' bit vectors for edit distance | Myers 1999 | genomics `align` | held |
| Vector lanes | striped Smith–Waterman; one target per lane | Farrar 2007; Rognes 2011 | genomics `align` on calx's lanes | held |
| Wavefronts | WFA; BiWFA in memory proportional to the score | Marco-Sola et al. 2021, 2023 | genomics `align` | held |
| Parametric alignment | the polytope of all optimal alignments, its normal fan | Gusfield, Balasubramanian and Naor 1994; Pachter and Sturmfels 2004 | morphe `formal/polytope` | held |
| Statistics of scores | Karlin–Altschul λ, K, E-values | Karlin and Altschul 1990 | genomics `align` | held |
| Pair HMMs | forward, posterior | Durbin et al. 1998 | genomics `align` | held |

## 3. Indexes and Search

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Suffix structures | SA-IS, LCP | Nong, Zhang and Chan 2011; Kasai et al. 2001 | morphe `formal` | held |
| Compressed search | the BWT and FM-index | Burrows and Wheeler 1994; Ferragina and Manzini 2000 | morphe `formal` | held |
| | both strands at once, the FMD-index and super-maximal exact matches | Li 2012 | genomics `index` | held |
| Repetitive collections | the run-length FM-index and the r-index, space proportional to the BWT's runs | Mäkinen and Navarro 2005; Gagie, Navarro and Prezza 2020 | morphe `formal` | added |
| Haplotype panels | the positional BWT, long and set-maximal matches | Durbin 2014 | morphe `formal` | added |
| Dictionaries | Aho–Corasick; IUPAC patterns on both strands in one pass | Aho and Corasick 1975 | morphe `formal`; genomics `index` | held; IUPAC patterns added |
| Sketches | MinHash, Mash distances; HyperLogLog | Broder 1997; Ondov et al. 2016; Flajolet et al. 2007 | genomics `sketch` | held |
| Pangenome graphs | the GBWT, variation-graph indexes | Sirén et al. 2020 | — | planned |

## 4. Mapping, Variants and Haplotypes

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Read mapping | minimizer seeds, chaining, base-level extension, mapping qualities | Li 2018 | genomics `mapping` | held |
| | mapping to graphs (sequence-to-graph alignment) | Jain et al. 2020 | — | planned |
| Small variants | pileups, genotype likelihoods of a diploid, normalisation of indels | Li 2011; Tan, Abecasis and Kang 2015 | genomics `variant` | held |
| Structural variants, copy number | discordant pairs and split reads clustered into breakpoints with LUMPY's probability intervals; SVTyper's genotypes; read depth corrected for GC and mappability, segmented by PELT; runs of homozygosity by bcftools roh's HMM | Rausch et al. 2012; Layer et al. 2014; Chiang et al. 2015; Benjamini and Speed 2012; Killick et al. 2012; Narasimhan et al. 2016 | genomics `structural`; tyche `time_series` for the segmentation | added (0064) |
| Haplotypes | Li and Stephens's copying model: forward–backward in O(NK), imputation, phasing | Li and Stephens 2003 | evolution `haplotypes` | added |
| Selection from haplotypes | EHH and iHS | Sabeti et al. 2002; Voight et al. 2006 | evolution `haplotypes` | added |

## 5. Assembly and Comparison

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Short reads | de Bruijn graphs cleaned of tips, bubbles and weak links | Pevzner, Tang and Waterman 2001; Zerbino and Birney 2008 | genomics `assembly` | held |
| | Eulerian circuits counted by the BEST theorem | van Aardenne-Ehrenfest and de Bruijn 1951 | morphe `combinatorics` | held |
| Long reads | string graphs, overlap–layout–consensus; minimizer-space de Bruijn graphs | Myers 2005; Ekim, Berger and Chikhi 2021 | — | planned |
| Whole-genome comparison | maximal unique matches and their longest increasing chain (MUMmer) | Delcher et al. 1999 | genomics `comparative`; morphe `combinatorics` for the chain | added |

## 6. Structure: Grammars, RNA, Symmetry

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Context-free grammars | Earley's parser; inside–outside over semirings; Viterbi parses | Earley 1970; Aycock and Horspool 2002; Goodman 1999; Lari and Young 1990 | morphe `formal` | added |
| Grammars with symmetry | the user's GE-CSG, made exact: rule orbits, the equivariant Earley parser, spans shared across orbits, orbit-canonical search of context-sensitive grammars, breakings, coarse-to-fine parsing | Huang 2026; [group-grammars.md](group-grammars.md) | morphe `formal`; genomics `grammar` | added |
| Representations of finite groups | exact character tables (Dixon–Schneider), group algebras and their Wedderburn data, Young's seminormal form, Fourier transforms on S_n and abelian groups, Walsh–Hadamard | Dixon 1967; Schneider 1990; Serre 1977; Diaconis 1988 | morphe `algebra`, `arithmetic` | added |
| The nucleotides' group | the Klein group of complementarity and transition: IUPAC's codes as cosets, the classifications as characters | Cornish-Bowden 1985 | genomics `grammar` | added |
| The genetic code's symmetry | the code as a function on (C₂)⁶, its stabiliser and invariance profile | NCBI's tables | genomics `grammar` | added |
| RNA secondary structure | Nussinov's maximal pairing; stochastic grammars trained by inside–outside, with orbit-tied (equivariant) EM; pair probabilities and maximum expected accuracy structures | Nussinov et al. 1978; Knudsen and Hein 1999; Do, Woods and Batzoglou 2006 | genomics `rna` | added |
| | free-energy minimisation with nearest-neighbour parameters | Zuker and Stiegler 1981; McCaskill 1990 | — | planned (needs the parameter tables as data) |
| | covariance models for families | Eddy and Durbin 1994 | — | planned |

## 7. Annotation, Expression, Classification

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Gene finding | periodic Markov chains trained on the genome itself | Besemer, Lomsadze and Borodovsky 2001 | genomics `annotation` | held |
| Islands, repeats, motifs | the CpG HMM; maximal repeats; PWMs, motif discovery by Gibbs sampling and EM | Durbin et al. 1998; Gusfield 1997; Lawrence et al. 1993; Bailey and Elkan 1994; Stormo 2000 | genomics `hmm`, `annotation` | held |
| Profile HMMs | — | Krogh et al. 1994; Eddy 1998 | genomics `hmm` | held |
| Expression | pseudoalignment by equivalence classes and EM; differential expression | Bray et al. 2016; Li and Dewey 2011 | genomics `expression` | held |
| Expression from codons | CAI, tAI and stAI, ENC (Wright; Sun et al.), Fop, MILC, Karlin's E(g) and PHX genes, dos Reis's S value; held against measured protein abundance (PaxDb) | Sharp and Li 1987; dos Reis et al. 2004; Sabi and Tuller 2014; Wright 1990; Sun, Yang and Xia 2013; Karlin and Mrázek 2000 | genomics `codon` | added (0064) |
| Expression from regulatory sequence | σ70's two boxes and spacer learnt by EM from measured start sites; thermodynamic occupancy (Bintu et al.); k-mer and gapped k-mer elastic nets; a convolutional network whose filters are read back as motifs, on the DREAM Challenge's promoters | Shea and Ackers 1985; Bintu et al. 2005; Ghandi et al. 2014; Alipanahi et al. 2015; Rafi et al. 2025 | genomics `regulation`; tyche `regression` and `learning` | added (0064) |
| | sequence-to-function models of whole loci (Enformer, Borzoi, AlphaGenome) | Avsec et al. 2021, 2025; Linder et al. 2025 | — | planned (their scale is beyond the machine) |
| Taxonomic classification | k-mers to the lowest common ancestor (Kraken) | Wood and Salzberg 2014 | genomics `classify` | added |

## 8. Evolution

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Trees | distances, UPGMA, neighbour joining; maximum likelihood under JC69 to GTR + Γ; the bootstrap; parsimony | Saitou and Nei 1987; Felsenstein 1981, 1985; Fitch 1971 | evolution `phylogeny` | held |
| | Bayesian trees by MCMC | Yang and Rannala 1997 | — | planned |
| Populations | the coalescent with recombination; Wright–Fisher; PCA and admixture | Kingman 1982; Hudson 1983; Patterson, Price and Reich 2006 | evolution `population` | held |
| | inference of ancestral recombination graphs | — | — | planned |
| Selection | dN/dS; McDonald–Kreitman | Nei and Gojobori 1986; McDonald and Kreitman 1991 | evolution `selection` | held |
| Association | scans with principal components and genomic control; clumping; cis-eQTLs with FastQTL's beta approximation and q-values; fine-mapping by ABF and SuSiE; colocalisation; expression predicted from genotypes and its association (PrediXcan, S-PrediXcan); heritability by Haseman–Elston and REML; polygenic scores | Devlin and Roeder 1999; Ongen et al. 2016; Storey and Tibshirani 2003; Wakefield 2009; Wang et al. 2020; Giambartolomei et al. 2014; Gamazon et al. 2015; Barbeira et al. 2018; Yang et al. 2011; Vilhjálmsson et al. 2015 | evolution `association`; tyche `regression`, `testing` | added (0064) |
| | haplotype scans, iHS | Voight et al. 2006 | evolution `haplotypes` | added |

## 9. Learning from Sequences

The draft's §8 proposes neural encoders that respect the group, built from the group algebra's blocks.

**What exists.**
- **The exact parts.** The isotypic decomposition of a sequence's one-hot encoding, and equivariant linear maps, block-diagonal by Schur's lemma. These are morphe's `algebra`, through its characters and group algebras.
- **The models.** These would run in the repository's Python environment (CLAUDE.md §7).

**Status.** Planned.

## 10. The Description of Variants and Chromosomes

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| A notation as a language | grammars written in BNF over bytes, parsed by Earley's parser, with exact diagnostics from its chart | Earley 1970; Aycock and Horspool 2002 | morphe `formal` (decision 0009) | added (0064) |
| Sequence variants | HGVS 21.1 parsed by its grammar and written back; the three normal forms of an indel (VCF's left, HGVS's 3′, SPDI's hull, VOCA); projection through transcripts; protein consequences; Sequence Ontology terms in VEP's order | den Dunnen et al. 2016; Hart et al. 2024; Tan et al. 2015; Holmes et al. 2020; McLaren et al. 2016 | genomics `nomenclature` | added (0064) |
| Structural records | VCF 4.3's symbolic alleles and breakend notation as adjacencies; HGVS of rearrangements | the VCF specification §5 | genomics `nomenclature` | added (0064) |
| Chromosomes | ISCN 2020/2024 parsed; karyotypes interpreted as derivative chromosomes of oriented band segments, copy number per band, the loss–gain–fusion model; the detailed system; `arr` and `seq` written from results | ISCN 2020, 2024; Ordulu et al. 2014; Abrams et al. 2015, 2019 | genomics `karyotype` | added (0064) |
| Decomposition of a delins into variants | the description extractor | Vis et al. 2015 | — | planned |

## 11. The Somatic Genome

| Problem | Algorithms | Sources | Where | Status |
|---|---|---|---|---|
| Somatic mutations | MuTect's tumour and normal log odds, its filters and power | Cibulskis et al. 2013 | genomics `somatic` | added (0064) |
| Purity, ploidy, allele-specific copy number | ASPCF on log R and B-allele fraction; ASCAT's fit over ρ and ψ | Van Loo et al. 2010 | genomics `somatic`; tyche `time_series` | added (0064) |
| Clonality | cancer cell fractions with multiplicity; binomial mixtures by EM and BIC; clone trees by the sum rule | McGranahan et al. 2015; Dentro et al. 2017 | genomics `somatic`; tyche `multivariate` | added (0064) |
| Mutational signatures | SBS-96 and DBS-78 catalogues; NMF with SigProfiler's stability; refitting by NNLS; bootstrap | Alexandrov et al. 2013, 2020; Islam et al. 2022 | genomics `somatic`; tyche `multivariate`, `regression` | added (0064) |
| Microsatellite instability | MSIsensor's χ² of repeat lengths, tumour against normal | Niu et al. 2014 | genomics `somatic` | added (0064) |
| | haplotype realignment (Mutect2, Strelka2), subclonal copy number (Battenberg), PyClone-VI | Benjamin et al. 2019; Kim et al. 2018; Gillis and Roth 2020 | — | planned |

## 12. Foundations Added for These

- **morphe** `formal`: BNF grammars and Earley's diagnostics (decision 0009).
- **tyche**:
  - change points: optimal partitioning, PELT, binary, wild binary and circular binary segmentation;
  - penalised regression: ridge, lasso, elastic net, cross-validation;
  - NNLS; NMF and its stability; binomial mixtures (decision 0006);
  - SuSiE, variance components and multiple testing (decision 0007);
  - networks and their training, and ranks (decision 0008).
- **calx**: the Hungarian method of assignment, and multinomial variates (decision 0014).
