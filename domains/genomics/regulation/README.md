# regulation: Expression from Regulatory Sequence

*Status: draft (2026-11-01).* What a promoter's sequence says of its activity: the σ⁷⁰ promoter of bacteria as two boxes and a spacer, learnt by expectation maximisation from *Escherichia coli*'s transcription start sites and scanned along genomes; the thermodynamic model of polymerase occupancy with activators and repressors (Shea and Ackers 1985; Bintu et al. 2005); and models learnt from measured promoter activity, a penalised linear model of k-mer and gapped k-mer counts and a convolutional network whose filters are read as motifs and matched to JASPAR. Serves use (the Genome Workbench draws promoters, logos and predictions), learning and research (the theory note [regulatory-sequence.md](../theory/regulatory-sequence.md)).

## 1. Concepts

- **A promoter is a mixture model.** A window about a start site holds a promoter or not; if it does, its −10 and −35 boxes lie at a variable offset and spacer. EM learns the boxes, the spacer's law and the fraction of windows with a promoter from start sites alone, started from MEME's sites (the annotation module) moved by a few bases either way, the most probable model kept.
- **Expression is occupancy.** In equilibrium the probability that polymerase occupies its site follows from the weights of all states of a promoter's sites; regulators change it by a regulation factor, which is the fold change for a weak promoter. Sequence enters through binding energies, which weight-matrix scores approximate.
- **Learnt models are judged on families they never saw.** The DREAM Challenge's measured promoters include variants of one another; held-out sets keep each family on one side.
- **A filter is a weight matrix.** On a one-hot sequence a convolution's filter scores each window as a PWM does; the windows where it fires most make its motif, which is matched to known profiles.

## 2. Content

| Header | Content |
|---|---|
| `encoding.hpp` | `one_hot<T>(codes)`, `one_hot_into`; `KmerFeatures(k_min, k_max, canonical)`, `GappedKmerFeatures(l, k, canonical)` with `columns()`, `names()`, `count(x, row)` |
| `sigma70.hpp` | `StartSite`, `read_start_sites(gtf)` (Ettwiller et al.'s Cappable-seq sites), `site_window(genome, site, upstream, downstream)`; `Sigma70{minus10, minus35, pi_r, pi_s, gamma, background, trace}`, `Sigma70Options`, `learn_sigma70(windows, options, pool)`, `window_posterior`, `promoter_score(model, window)` (log₂ likelihood ratio); `PromoterHit`, `scan_promoters(model, x, threshold)` (both strands, forward coordinates, the predicted +1), `promoter_threshold(model, rate)`; `consensus(pwm)` |
| `thermodynamic.hpp` | `Species`, `BindingSite`, `SiteInteraction`, `Architecture`; `occupancy(architecture)` → `Occupancy{z, p_bound, site, states, weights}` (exact enumeration); `fold_change`; Bintu et al.'s closed forms `p_bound_basal`, `p_bound_regulated`, `regulation_activation`, `regulation_repression`; `site_energy(bits, λ, μ)` |
| `learn.hpp` | `Measured`, `read_measured(text)`, `subset`, `grouped_split(data, keys, validation, test, seed)`, `agreement(predicted, measured)` (Pearson, Spearman, R²); `KmerBaselineOptions`, `fit_kmer_baseline(train, options)` (elastic net, λ by cross-validation), `predict`; `CnnOptions`, `SequenceCnn`, `train_cnn(train, validation, options, pool)`, `predict`; `FilterMotif`, `filter_motifs(cnn, sequences)`; `KnownMotif`, `read_meme_motifs(text)`, `match_motif(query, width, known, min_overlap, shuffles, seed)` → `MotifMatch`, `informative_core` |

## 3. Dependencies

`genomics::annotation` (PWMs, logos, MEME's EM), `genomics::sequence`; tyche `learning` (networks, Adam, training; decision 0008), `regression` (the elastic net; decision 0006), `descriptive` (Pearson, Spearman); calx `random`, `scheduler`. Data, by the data library: `ecoli-tss-cappable` (Ettwiller et al. 2016, CC BY 4.0, with U00096.2), `dream-promoters` (Rafi et al. 2025, CC BY 4.0), `jaspar-fungi` (JASPAR 2024, CC BY 4.0).

## 4. Checks

`node domains/genomics/tests/run-tests.mjs`, or one program at a time (`check_genomics_regulation_<topic>`).

| Check | What it verifies |
|---|---|
| `encoding` | One-hot by hand; k-mer counts (k = 2 … 4, plain and canonical) and gapped 4-mers in windows of 6 (plain and canonical) against brute force, with the number of canonical classes in closed form; occupancy against Bintu et al.'s eqs. 4, 5, 8, 9 and 10 and the dual repression of their Table 1 to 10⁻¹⁵, the fold change of a weak promoter equal to F_reg |
| `sigma70` | A planted two-box model recovered from 1500 simulated windows (boxes GATAAT and CTGACA exactly, the spacer law within 0.02 in total variation, γ within 0.02), EM's log posterior never falling; a planted promoter found on both strands at its coordinates and +1; Ettwiller et al.'s 16,359 sites and U00096.2 read; learnt from *E. coli*'s 2000 strongest sites: −10 box TATAAT (Harley and Reynolds's consensus), −35 box CTTGCC holding the conserved TTG one base upstream of the hexamer, the commonest spacer 17 from the canonical hexamer (18 from the learnt box), the −10 box most often at −12; held-out sites against random genomic windows, area under the ROC curve 0.93; held-out promoter scores against measured strength (RRS), Spearman ρ = 0.22 |
| `learn` | MEME's format; the matcher (a word against itself, similarity 1 at its offset; reverse strands); a planted code (two words and GC content, 6000 sequences, families held out): the k-mer baseline held-out R² 0.85, the CNN 0.84, both words among the CNN's filter motifs; DREAM promoters (8500 to train, 3556 held-out sequences of unseen families): the k-mer baseline Pearson r 0.66 (ρ 0.67), the CNN r 0.56 (ρ 0.56); JASPAR's 178 fungal profiles read, at least one filter matched (SPT23, column correlation 0.92, empirical p ≤ 0.01) |

## 5. Performance

Measured in WSL on the development machine (i7-11370H), Release, four threads: the σ⁷⁰ model learnt from 2000 windows (25 EM runs from shifted starts) in about 6 s; the CNN (32 filters of 15, 8 epochs) on 8500 sequences of 110 bases in 5 s; the k-mer baseline (340 columns, 8500 rows) by 3-fold cross-validation over 12 penalties in 31 s (360 s with 5 folds, 30 penalties and a path to 10⁻⁴ λ_max: coordinate descent on dense counts is the slow part).

## 6. References

- Hawley, D. K. and McClure, W. R. (1983), Compilation and analysis of *Escherichia coli* promoter DNA sequences, *Nucleic Acids Res.* 11, 2237–2255.
- Harley, C. B. and Reynolds, R. P. (1987), Analysis of *E. coli* promoter sequences, *Nucleic Acids Res.* 15, 2343–2361.
- Ross, W., Gosink, K. K., Salomon, J., Igarashi, K., Zou, C., Ishihama, A., Severinov, K. and Gourse, R. L. (1993), A third recognition element in bacterial promoters: DNA binding by the α subunit of RNA polymerase, *Science* 262, 1407–1413.
- Shea, M. A. and Ackers, G. K. (1985), The O_R control system of bacteriophage lambda: a physical-chemical model for gene regulation, *J. Mol. Biol.* 181, 211–230.
- Berg, O. G. and von Hippel, P. H. (1987), Selection of DNA binding sites by regulatory proteins, *J. Mol. Biol.* 193, 723–750.
- Bintu, L., Buchler, N. E., Garcia, H. G., Gerland, U., Hwa, T., Kondev, J. and Phillips, R. (2005), Transcriptional regulation by the numbers: models; … applications, *Curr. Opin. Genet. Dev.* 15, 116–124 and 125–135.
- Kinney, J. B., Murugan, A., Callan, C. G. and Cox, E. C. (2010), Using deep sequencing to characterize the biophysical mechanism of a transcriptional regulatory sequence, *PNAS* 107, 9158–9163.
- Ettwiller, L., Buswell, J., Yigit, E. and Schildkraut, I. (2016), A novel enrichment strategy reveals unprecedented number of novel transcription start sites at single base resolution in a model prokaryote and the gut microbiome, *BMC Genomics* 17, 199.
- Bailey, T. L. and Elkan, C. (1994), Fitting a mixture model by expectation maximization to discover motifs in biopolymers, *Proc. ISMB* 2, 28–36; Liu, X., Brutlag, D. L. and Liu, J. S. (2001), BioProspector, *Pac. Symp. Biocomput.* 6, 127–138.
- Ghandi, M., Lee, D., Mohammad-Noori, M. and Beer, M. A. (2014), Enhanced regulatory sequence prediction using gapped k-mer features, *PLoS Comput. Biol.* 10, e1003711.
- Alipanahi, B., Delong, A., Weirauch, M. T. and Frey, B. J. (2015), Predicting the sequence specificities of DNA- and RNA-binding proteins by deep learning, *Nat. Biotechnol.* 33, 831–838.
- Kelley, D. R., Snoek, J. and Rinn, J. L. (2016), Basset, *Genome Res.* 26, 990–999.
- de Boer, C. G. et al. (2020), Deciphering eukaryotic gene-regulatory logic with 100 million random promoters, *Nat. Biotechnol.* 38, 56–65; Vaishnav, E. D. et al. (2022), The evolution, evolvability and engineering of gene regulatory DNA, *Nature* 603, 455–463.
- Agarwal, V. and Shendure, J. (2020), Predicting mRNA abundance directly from genomic sequence using deep convolutional neural networks, *Cell Rep.* 31, 107663.
- Avsec, Ž. et al. (2021), Effective gene expression prediction from sequence by integrating long-range interactions (Enformer), *Nat. Methods* 18, 1196–1203; Linder, J. et al. (2025), Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation (Borzoi), *Nat. Genet.* 57; Avsec, Ž. et al. (2025), AlphaGenome.
- Rafi, A. M. et al. (2025), A community effort to optimize sequence-based deep learning models of gene regulation, *Nat. Biotechnol.* 43, 1373–1383.
- Gupta, S., Stamatoyannopoulos, J. A., Bailey, T. L. and Noble, W. S. (2007), Quantifying similarity between motifs, *Genome Biol.* 8, R24; Rauluseviciute, I. et al. (2024), JASPAR 2024, *Nucleic Acids Res.* 52, D174–D182.

## 7. Open Questions

1. **An *E. coli* reporter assay.** Urtecho et al.'s σ⁷⁰ library (GEO GSE108535) would let the thermodynamic model's energies be fitted to measured activity; NCBI did not serve its files on 2026-11-01. Registering it (or LaFleur et al. 2022's) is planned.
2. **The extended −10, the UP element and the discriminator** are not in the σ⁷⁰ model; nor are other σ factors.
3. **Larger networks.** Dilated convolutions, attention and training on the challenge's 6.7 million sequences need calx's faster kernels; the module's networks are small by design.
4. **Variant effects.** The difference of predictions between alleles (in silico mutagenesis) is a few lines on top of `predict` and is planned for the Workbench.
