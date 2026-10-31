# microbiology

*Status: planned. Field, height 2, on genomics, evolution and cell biology.* Bacteria, archaea, viruses and microbial eukaryotes: how they grow, what their genomes share and exchange, how they are classified, how they live together, how phages and their hosts contend, and how microbes infect. Bacteriology, virology, the microbiome and the systematics of microbes are its modules. Placed in the order with the other domains of life at the user's request (2026-10-31, [ADR 0061](../../development/decisions/0061-the-domains-of-life.md)); its work is in tier 4 of the [README's order](../../README.md#46-the-order-of-implementation).

## 1. Purpose

- **Growth is a law of the cell's economy.** A population in a chemostat grows at Monod's rate μ(s) = μ_max s / (K_s + s) in its limiting substrate s; at steady state μ equals the dilution rate, unless the population washes out. Within the cell, the growth rate fixes how the proteome is shared between ribosomes and the rest, linearly (Scott et al. 2010).
- **Mutants arise before selection.** The number of resistant mutants in parallel cultures follows the Luria–Delbrück distribution, not a Poisson distribution, because mutations arise during growth and are inherited by their descendants (Luria and Delbrück 1943).
- **A species is its pangenome.** The genes of a species' genomes divide into a core that all carry and an accessory set that some carry; the pangenome is open when every new genome keeps adding genes (Heaps' law; Tettelin et al. 2008). Genes pass between lineages by horizontal transfer, so a microbe's history is a network with reticulations, which evolution's trees approximate.
- **Microbes are classified by their genomes.** An average nucleotide identity near 95% bounds a species (Jain et al. 2018), and the tree of bacteria and archaea is inferred from concatenated marker proteins (Parks et al. 2018). The systematics of microbes needs evolution's phylogeny; the systematics of plants, animals and fungi is evolution's own (README §4.5).
- **Phages and their hosts contend.** Phage and bacteria form a predator–prey system (Levin, Stewart and Chao 1977); phage λ chooses between lysis and lysogeny by a stochastic switch (Arkin, Ross and McAdams 1998); CRISPR is an adaptive immunity written into the host's genome.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `growth` | Growth curves and their phases, Monod kinetics, the chemostat and washout, the growth laws of bacteria, persistence, the Luria–Delbrück distribution | planned |
| `genomes` | Pangenomes (core and accessory genes, Heaps' law), horizontal transfer, plasmids and prophages, 16S rRNA; on genomics's annotation and sketches | planned |
| `systematics` | The classification of bacteria, archaea and viruses: species by ANI, trees of marker genes on evolution's phylogeny, the taxonomies of GTDB and ICTV | planned |
| `communities` | Microbiomes: composition from reads (genomics's `classify`), compositional data in Aitchison's geometry, diversity, consumer–resource models and cross-feeding | planned |
| `phages` | The dynamics of phages and their hosts, lysis and lysogeny in λ, CRISPR immunity | planned |
| `infection` | Infection within a host: target-cell-limited models of viral infection (Nowak and May 2000), bacterial growth under antibiotics, the rise of resistance | planned |

## 3. Dependencies

| Source | What microbiology uses from it |
|---|---|
| genomics | Genomes and their annotation, sketches and ANI, the taxonomic classification of reads |
| evolution | Phylogeny; population genetics (mutation, selection, the rise of resistance) |
| cell-biology | Regulation (the λ switch), metabolism (growth, cross-feeding) |
| physis | `chemistry` (mass action and its exact Markov chains), `statistical_mechanics` (master equations) |
| tyche | Regression; compositional data by Aitchison's log-ratios (to be written there) |
| calx | Integrators, random numbers, optimisation |
| morphe | Graphs (networks of descent, of cross-feeding) |

## 4. Used By

No work uses it yet. The Genome Workbench already holds *Escherichia coli* and phage λ among its genomes, and classifies reads by Kraken's method (genomics's `classify`).

## 5. Open Questions

- Infection between hosts is epidemiology's; within a host, the pathogen's side is here and the host's response is physiology's `immunity`. Which bridge joins them.
- Whether a network of descent with reticulations belongs to evolution, beside its trees.
- Microbial ecology: here (`communities`) for the microbiome, or ecology's once ecology declares this domain, which it can do without a cycle.

## 6. References

- Madigan, M. T., Bender, K. S., Buckley, D. H., Sattley, W. M. and Stahl, D. A. (2018), *Brock Biology of Microorganisms*, 15th ed., Pearson.
- Monod, J. (1949), The growth of bacterial cultures, *Annual Review of Microbiology* 3, 371–394.
- Luria, S. E. and Delbrück, M. (1943), Mutations of bacteria from virus sensitivity to virus resistance, *Genetics* 28, 491–511.
- Scott, M., Gunderson, C. W., Mateescu, E. M., Zhang, Z. and Hwa, T. (2010), Interdependence of cell growth and gene expression: origins and consequences, *Science* 330, 1099–1102.
- Tettelin, H. et al. (2005), Genome analysis of multiple pathogenic isolates of *Streptococcus agalactiae*: implications for the microbial "pan-genome", *PNAS* 102, 13950–13955.
- Tettelin, H., Riley, D., Cattuto, C. and Medini, D. (2008), Comparative genomics: the bacterial pan-genome, *Current Opinion in Microbiology* 11, 472–477.
- Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T. and Aluru, S. (2018), High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries, *Nature Communications* 9, 5114.
- Parks, D. H. et al. (2018), A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life, *Nature Biotechnology* 36, 996–1004.
- Levin, B. R., Stewart, F. M. and Chao, L. (1977), Resource-limited growth, competition, and predation: a model and experimental studies with bacteria and bacteriophage, *American Naturalist* 111, 3–24.
- Arkin, A., Ross, J. and McAdams, H. H. (1998), Stochastic kinetic analysis of developmental pathway bifurcation in phage λ-infected *Escherichia coli* cells, *Genetics* 149, 1633–1648.
- Aitchison, J. (1986), *The Statistical Analysis of Compositional Data*, Chapman and Hall.
- Nowak, M. A. and May, R. M. (2000), *Virus Dynamics: Mathematical Principles of Immunology and Virology*, Oxford University Press.
