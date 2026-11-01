# 0064. Abnormal Genomes: Their Description in the International Notations, Their Detection, and Expression Read from Sequence

- **Status:** Accepted, provisionally; amends 0043 and 0047
- **Date:** 2026-11-01
- **Decided by:**
  - **The request.** The user asked Claude to "continue implement traditional and advanced genome analysis algorithms, especially abnormal detection, expression prediction". They added that abnormalities have "an international representation description system … from where to where is lack, which part swap to which part", and asked for it too. Then: "you can implement various algorithms, also can see STOA papers. Remain UI elegance. Be creative."
  - **What Claude did.**
    - It placed the work in the order and wrote the contracts between the parts.
    - It wrote morphe's BNF grammars and Earley's diagnostics (morphe decision 0009), the Workbench's five views, the server's naming of data, and this record.
  - **What the agents did.** Seven agents wrote the modules to those contracts, each checked against references and simulated truth:
    - tyche's change points, penalised regression, NNLS and NMF;
    - HGVS;
    - ISCN;
    - structural variants;
    - the somatic genome;
    - codons and regulation;
    - association.

## 1. Context

Two families were missing from the genomics and evolution domains (ADR 0043, ADR 0047). Both were listed as planned in the survey `domains/genomics/theory/the-landscape-of-genome-analysis.md`.

- **Abnormal genomes.** The domain called small variants, but nothing larger:
  - structural variants and copy number;
  - losses of heterozygosity and aneuploidy;
  - the somatic genome of a tumour.
- **Expression predicted from sequence or genotype.** The domain quantified expression from reads, but did not predict it.

What the user calls the international description system is two standards:

- **HGVS** (the HGVS Nomenclature, now 21.1) describes a change of sequence against a named reference.
- **ISCN** (the International System for Human Cytogenomic Nomenclature, 2020 and 2024) describes chromosomes:
  - a karyotype seen under the microscope;
  - an array's copy number (`arr`);
  - the junctions found by sequencing (`seq`).

VCF's structural records and NCBI's SPDI are their machine forms. Both standards are languages with grammars, and the user's field is formal languages, so the notations are treated as languages: parsed by their grammars, and given their meaning by computation.

## 2. Decision

### 2.1 The Notations as Languages

- **morphe `formal`** (decision 0009):
  - grammars written in BNF over bytes, compiled over the atoms of the byte sets they name;
  - parsed by Earley's parser into named trees;
  - exact diagnostics from the chart: the furthest viable prefix, and what may follow it.
- **genomics `nomenclature`**:
  - **The grammar.** HGVS 21.1 as a grammar of 206 rules. Each of the standard's 702 examples parses with one derivation and writes back to itself.
  - **The three normal forms.** An indel's minimal representations form an interval of placements (the theory note's Theorem 1, proved and checked by brute force). VCF's left form (Tan et al.), HGVS's 3′ form and SPDI's hull (VOCA) are its two ends and its hull.
  - **Transcripts and proteins.** The projection through RefSeq transcripts, numbered as VEP numbers them. The protein consequence, and the Sequence Ontology's terms in VEP's order.
  - **Structural records.** VCF 4.3's symbolic and breakend records as adjacencies, and their HGVS.
- **genomics `karyotype`**:
  - **The grammar and the bands.** ISCN's grammar, and the human ideogram at 400, 550 and 850 bands as a tree in which a band is a word and the hierarchy is its prefix order.
  - **The interpretation.** Every chromosome of a karyotype is rebuilt as a word of oriented segments of the normal chromosomes. Copy number band by band, losses, gains and fusions from where to where, CytoGPS's loss–gain–fusion vectors, and ISCN's detailed system all follow from these words.
  - **Description.** `arr` and `seq` are written from detected results.
- **The rules' sources.** ISCN 2020 and 2024 are not open, so the rules are taken from open summaries and examples, each cited in the theory note. Where only secondary sources exist, the note says so.

### 2.2 Detection

- **genomics `structural`**, with `simulate`'s structural samples and paired reads:
  - discordant pairs and split reads, combined into breakpoints by LUMPY's probability intervals;
  - events read from their junctions' signatures, and genotypes after SVTyper;
  - read depth corrected for GC and mappability, segmented by tyche's PELT into integer copy numbers;
  - B-allele fractions and runs of homozygosity by bcftools roh's HMM. With depth, they tell copy-neutral LOH from a deletion;
  - aneuploidy;
  - evaluation as Truvari evaluates.
  - Every call is written as VCF, HGVS and ISCN.
- **genomics `somatic`**, with `simulate`'s tumours:
  - MuTect's tumour and normal tests;
  - ASCAT's purity, ploidy and allele-specific copy number;
  - cancer cell fractions, clones by binomial mixtures, and clone trees by the sum rule;
  - SBS-96 and DBS-78 catalogues, and signatures by NMF with SigProfiler's stability, refitted by NNLS;
  - MSIsensor's microsatellite instability.
- **Smaller additions to existing modules:**
  - `hmm` gained chains with per-step emissions;
  - `sequence` gained the trinucleotide and dinucleotide contexts;
  - `mapping`'s choice of supplementary parts now compares query overlap in the read's own coordinates across strands. Before, the two parts of a read across an inversion were taken as rivals.

### 2.3 Expression

- **genomics `codon`:**
  - CAI, tAI and stAI from a genome's tRNA genes;
  - ENC by Wright's and by Sun's estimators, Fop, MILC, Karlin's E(g), and dos Reis's S value;
  - all held against measured protein abundance.
- **genomics `regulation`:**
  - σ70's promoter learnt by EM from measured start sites;
  - thermodynamic occupancy;
  - expression learnt from measured promoters, by k-mer elastic nets and by a convolutional network whose filters are read back as motifs.
- **evolution `association`:** expression and traits from genotypes:
  - scans with principal components and genomic control;
  - cis-eQTLs, fine-mapping and colocalisation;
  - PrediXcan and S-PrediXcan;
  - heritability and polygenic scores.

### 2.4 Foundations

Every method that names no object in nature went to a foundation:

- **morphe** decision 0009: BNF grammars and Earley's diagnostics.
- **tyche** decision 0006:
  - change points: optimal partitioning, PELT, binary, wild binary and circular binary segmentation;
  - penalised regression: ridge, lasso, elastic net, cross-validation;
  - NNLS; NMF and its stability; binomial mixtures.
- **tyche** decision 0007:
  - the sum of single effects (SuSiE) on individual data, sufficient statistics and summary statistics, as susieR fits it;
  - variance components by REML (exact over one kernel, AI-REML over several) and Haseman–Elston;
  - Reinsch's smoothing spline;
  - a module `testing` (Bonferroni, Holm, Benjamini–Hochberg, π₀ and q-values, permutation p-values, and FastQTL's beta approximation and effective degrees of freedom);
  - a correction of the Student t law's upper tail, which lost every digit below about 10⁻¹⁶ for large degrees of freedom, as GWAS p-values need.
- **tyche** decision 0008: networks and their training, and ranks.
- **calx** decision 0014: the Hungarian method, and multinomial variates.

genomics and evolution now declare physis and tyche among their foundations.

### 2.5 The Genome Workbench

Five views:

- **Describe.** HGVS and ISCN, read by their grammars, their derivations drawn over the text.
  - HGVS: a variant on its reference in its three forms, over the letters of the interval of its placements, carried through its transcripts to the protein.
  - ISCN: a karyotype's chromosomes rebuilt and painted by the chromosome each piece came from. Then gains and losses on the ideogram, and the changes from band to band in GRCh38's bases.
- **Rearrangements.** Copies per bin, B-allele fractions, allelic states and calls along the genome, joined by arcs to their other ends.
  - Each call's pairs and splits are drawn on a broken axis, with its three descriptions.
  - The individual's planted chromosomes are painted as pieces of the reference.
- **Tumour.**
  - the allele-specific copy-number profile, and ASCAT's sunrise with its fit and the truth;
  - clones, their fractions and their tree;
  - the catalogue and its refit, microsatellites, and signatures found from a cohort's catalogues alone.
- **Regulation.**
  - codon indices round the genome and against abundance, Wright's Nc curve, the code table of weights;
  - σ70's boxes as logos with their spacer;
  - promoters along the genome;
  - the learnt models against held-out measurement, with their filters as logos matched to JASPAR.
- **Association.** The study's scan, locus, fine-mapping, eQTLs, predicted expression, heritability and scores, each beside its truth.

The server names a dataset of the data library, or one of its files, as `data:<id>[/<file>]`, and gives an operation the data it needs (the ideogram for the karyotype's operations). Simulated samples are labelled so on every page; their references are real.

### 2.6 Data

New manifests, each with its licence:

- `ncbi-ideogram`: NCBI's Genome Decoration Page and UCSC's cytoBand;
- `ncbi-human-genes`: RefSeq regions and transcripts of HBB, CFTR and BRCA1;
- `ncbi-grch38-chromosomes`: GRCh38's chromosomes 21 and 22;
- `paxdb-ecoli`: PaxDb 6.0, CC BY 4.0;
- `ecoli-tss-cappable`: Ettwiller et al. 2016, CC BY 4.0;
- `dream-promoters`: CC BY 4.0;
- `jaspar-fungi`: JASPAR 2024, CC BY 4.0.

COSMIC's signatures are not registered. Their terms forbid redistribution and display on free platforms, so the signatures of the checks and the Tumour view are built from their mechanisms and labelled synthetic. A licensed user can read COSMIC's matrix with `somatic::read_profiles`.

## 3. Alternatives Considered

- **Hand-written parsers for HGVS and ISCN.** Rejected. A grammar in a foundation is the notation's definition, checked on the standard's examples, unambiguous by count, and diagnosable from its chart. A parser in the domain would hide it.
- **The reference tools run as services** (Mutalyzer, VEP, CytoGPS). Used to make fixtures, and not as dependencies. Their answers are the references the domain is held to, each with its version and date.
- **COSMIC's signatures in the cache.** Rejected for the terms above.
- **Sequence-to-function models of whole loci** (Enformer, Borzoi, AlphaGenome). Named in the theory note and not attempted: their training needs data and devices beyond the machine. The regulatory models here are trained on 8,500 promoters, where k-mers beat the network, as measured.
- **ClinGen's dosage regions.** Deferred until their licence is checked.

## 4. Consequences

- **The checks.**
  - **nomenclature (90):** HGVS 21.1's 702 examples; Mutalyzer's extractor (92 of 92); 53 real variants against Mutalyzer, NCBI's Variation Services and VEP; 352 positions against VEP.
  - **karyotype (173):** ISCN 2009's 181 examples, with their detailed forms; CytoGPS's computed clones; Ordulu et al.'s and Conlin et al.'s published strings.
  - **structural:** precision 1.000 and recall 0.97–1.00 on simulated individuals on E. coli and GRCh38; SVTyper's genotype likelihoods reproduced to 10⁻⁹.
  - **somatic (104):** purity exact and ploidy within 0.003 on a simulated tumour; signatures within 4%.
  - **codon and regulation (110):** CAI equal to Biopython's to 5 × 10⁻¹⁶; ρ 0.58 between CAI and measured abundance; held-out r 0.66 on the DREAM promoters.
  - **association and tyche 0007 (182):** SuSiE reproduces susieR 0.14.2, its credible sets covering 99.6% at 95% over 300 simulated regions; REML matches R's regress; coloc's posteriors match R's coloc to 10⁻¹⁵; the beta approximation matches 20,000 direct permutations; S-PrediXcan matches the individual-level TWAS within 1.9%.
  - **tyche (650 of 653 for decision 0006):** the three failures are Hodrick–Prescott round-off checks that predate this round.
  - **The suites:**
    - genomics 1696 of 1696 and evolution 514 of 514, each also alone with its down-set;
    - morphe 713 of 713, calx 902 of 902;
    - the whole run 9515 of 9520 in 38 suites. The five failures are not of this round: two are in the nomos foundation, then in progress, and three are tyche's Hodrick–Prescott round-off checks.
- **Planned:**
  - the decomposition of a delins into variants (Mutalyzer's extractor), splicing and NMD;
  - ish and rsa in ISCN;
  - local assembly of breakends; long-read SVs;
  - haplotype realignment and subclonal copy number for tumours;
  - an E. coli reporter assay for the thermodynamic model;
  - a faster PrediXcan: training an elastic net per gene with nested cross-validation takes minutes, most of the Association view's wait;
  - LD score regression, mixed-model scans and real genotype panels;
  - the reference tools' licences and data where they allow.
