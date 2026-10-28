# evolution

*Status: draft. Field, height 1, on genomics.* The history of genomes: the trees of their descent, inferred from their sequences by distances, likelihood and parsimony; the genealogies and frequencies of variants in populations, by the coalescent and by Wright–Fisher; haplotypes read against a reference panel, as mosaics of it, to phase and impute them; and the signatures of selection, in coding sequences and in the length of haplotypes ([ADR 0043](../../development/decisions/0043-the-genome-infrastructure.md)). Drawn by the [Genome Workbench](../../applications/genome-workbench/)'s Trees and Populations views. Serves use, learning and research.

## 1. Concepts

- **A phylogeny is a tree metric, and a probability model on it.** A tree with lengths on its branches is a metric on its leaves satisfying the four-point condition; distance methods (UPGMA, neighbour joining) recover it from estimated distances. A continuous-time Markov chain on its branches makes it a probability model of the sequences at its leaves, whose likelihood Felsenstein's pruning computes as a sum–product, and parsimony as its tropical shadow.
- **A population is a process whose genealogy is a random tree.** The Wright–Fisher model moves allele frequencies forward by binomial sampling; Kingman's coalescent runs its genealogy backward, and mutations on the genealogy's branches give the variation a sample shows. Estimators of θ, Tajima's D and the site frequency spectrum are statistics of that tree.
- **A haplotype is a mosaic of others.** Li and Stephens's copying model reads a haplotype as copied, piece by piece, from a panel, switching at the rate of recombination: a hidden Markov model over the panel whose uniform jumps make each step linear in the panel's size. Its posterior imputes the sites a target was not typed at, and a pair of such processes phases a diploid's genotypes.
- **Selection leaves ratios, and long haplotypes.** Synonymous and non-synonymous change (dN/dS), and fixed against polymorphic differences (the McDonald–Kreitman test), compare what selection sees with what it does not. An allele that rose fast carries long haplotypes for its frequency, which EHH and iHS measure.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `phylogeny` | Trees and Newick, splits, Robinson–Foulds, majority-rule consensus, layouts; distances (p, JC69, K80, F81, TN93 with Γ, LogDet); UPGMA, neighbour joining, BIONJ; JC69 to GTR with discrete Γ and +I by spectral decomposition; Felsenstein's pruning on lanes, branch lengths by Newton's method, NNI search, κ, α and invariant sites; Fitch, Hartigan and Sankoff parsimony, stepwise addition and branch and bound; the bootstrap on task graphs; simulation along trees, Yule, coalescent and uniform trees; root-to-tip dating ([`phylogeny/`](phylogeny/include/evolution/phylogeny/)) | draft |
| `population` | Genotype matrices and allele counts; Hardy–Weinberg (χ², Wigginton's exact test); linkage disequilibrium (D, D′, r², EM for unphased, decay); diversity (S, θ_W, π, Tajima's D, Fay and Wu's H), spectra and F_ST (Hudson, Weir–Cockerham) in windows; the coalescent with growth, bottlenecks and recombination (Hudson's algorithm, marginal trees); Wright–Fisher with selection, mutation and migration; a selective sweep forward in time from the coalescent's standing variation; PCA (Patterson) and admixture by EM ([`population/`](population/include/evolution/population/)) | draft |
| `haplotypes` | Li and Stephens's copying model on a reference panel: forward, backward, posteriors and Viterbi in O(N K) by the uniform jump; imputation of untyped sites from the posterior; the pair model of a diploid in O(N K²), its posterior, imputation of unphased genotypes and phasing as the joint maximum, switch errors; EHH by allele, iHH by the trapezoid rule over the map, iHS standardised within bins of frequency ([`haplotypes/`](haplotypes/include/evolution/haplotypes/)) | draft |
| `selection` | Synonymous and non-synonymous sites and differences by Nei and Gojobori, dN/dS under any NCBI genetic code; the McDonald–Kreitman table, neutrality index, α, direction of selection, Fisher's exact and Williams-corrected G-tests ([`selection/`](selection/include/evolution/selection/)) | draft |

## 3. Dependencies

- **Domains:** genomics (sequences, genetic codes, alignments, variants), so evolution has height 1.
- **Foundations:** calx (`linalg` for spectral decompositions and principal components, `special_functions` for Γ rates and tests, `random` for simulation, `scheduler` for bootstrap replicates and many simulations) and morphe.

## 4. Code and Checks

The C++ is header-only, C++20, built with CMake from this directory (which adds the genomics domain beneath it), each module a target `evolution::<module>`. The checks are run by

```
node domains/evolution/tests/run-tests.mjs
```

The theory, with statements and proofs, is in `theory/`: [phylogeny](theory/phylogeny.md) (splits and Buneman's theorem, the four-point condition and neighbour joining, reversible Markov chains, pruning as a sum–product with parsimony its tropical case, Γ rates, the bootstrap, root-to-tip dating), [population](theory/population.md) (the coalescent, Wright–Fisher and its diffusion, diversity, linkage, structure), [haplotypes](theory/haplotypes.md) (the copying chain as a semigroup and its linear steps, imputation as the posterior predictive, the pair model and phasing, EHH and iHS, the sweep simulated forward) and [selection](theory/selection.md) (the genetic code on the Hamming graph, synonymous sites, ω, the McDonald–Kreitman test).

## 5. References

- Felsenstein, J. (2004), *Inferring Phylogenies*, Sinauer.
- Yang, Z. (2014), *Molecular Evolution: A Statistical Approach*, Oxford University Press.
- Semple, C. and Steel, M. (2003), *Phylogenetics*, Oxford University Press.
- Hein, J., Schierup, M. H. and Wiuf, C. (2005), *Gene Genealogies, Variation and Evolution: A Primer in Coalescent Theory*, Oxford University Press.
- Wakeley, J. (2009), *Coalescent Theory: An Introduction*, Roberts and Company.
- Nei, M. and Kumar, S. (2000), *Molecular Evolution and Phylogenetics*, Oxford University Press.
- Li, N. and Stephens, M. (2003), Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data, *Genetics* 165, 2213–2233; Voight, B. F., Kudaravalli, S., Wen, X. and Pritchard, J. K. (2006), A map of recent positive selection in the human genome, *PLoS Biology* 4, e72 (the haplotypes module's README has the rest).
