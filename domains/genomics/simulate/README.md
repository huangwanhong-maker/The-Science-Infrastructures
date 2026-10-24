# simulate: Genomes, Individuals and Reads with Their Truth

*Status: draft (2026-10-24).* Seeded simulation of what the genome's algorithms are checked on: genomes with repeats, diploid individuals mutated from a reference, reads from them with a sequencer's errors, and families of sequences with their true alignment. A simulation always comes with its truth, which is why it is here: a real genome tells no one where its reads came from.

## 1. Concepts

- **A simulation is a function of its seed.** Every draw comes from a `calx::random::Stream` keyed by the seed, a purpose and an object; read i comes from its own stream. Reads a to b are therefore the same whether drawn alone or within reads 0 to n, so that a shard of a simulated sample is two numbers, and a worker on another machine regenerates it exactly.
- **Truth travels with the data.** A diploid individual carries its variants (VCF-style, indels anchored on the base before them) and, for each haplotype, a liftover to the reference; a read carries its contig, its reference interval, its strand, its haplotype and its number of errors; a family carries its true alignment and its tree.

## 2. Content

| Header | Content |
|---|---|
| [`genome.hpp`](include/genomics/simulate/genome.hpp) | Random sequence of a GC content; genomes with interspersed repeats (diverged copies on both strands) and tandem arrays, with where each copy lies |
| [`mutate.hpp`](include/genomics/simulate/mutate.hpp) | A diploid individual: SNVs and indels at given rates, heterozygous or homozygous, never overlapping; its two haplotypes; liftovers to the reference |
| [`reads.hpp`](include/genomics/simulate/reads.hpp) | Short reads with substitution errors rising along the read, a few one-base indels, informative Phred qualities; long reads with gamma lengths and indel-dominated errors; from one or two copies of several contigs; as FASTQ |
| [`family.hpp`](include/genomics/simulate/family.hpp) | Sequences evolved from a root along a balanced tree with substitutions and indels, with their true alignment (every base keeps the column it was born in) and the tree in Newick |

## 3. Dependencies

`genomics::sequence`; calx (`random`).

## 4. Checks

`check_genomics_simulate`:

| Check | What it verifies |
|---|---|
| genomes | A function of the seed; the GC asked for; the copies placed without overlap; copies of an element differing at 2p − 4p²/3 |
| individuals | The haplotypes are the reference with exactly their variants applied; SNV and indel rates and zygosity as asked; the liftover sends unmutated bases home |
| short reads | A shard equals the same reads of the whole sample; lengths and qualities; error rates at the first and last bases as the model says; qualities lower on erroneous bases |
| long reads | The mean length and the error rate per template base asked |
| families | The true alignment without its gaps gives back each member; no empty column; wider than any member; the tree in Newick |

## 5. References

- Dohm, J. C., Lottaz, C., Borodina, T. and Himmelbauer, H. (2008), Substantial biases in ultra-short read data sets from high-throughput DNA sequencing, *Nucleic Acids Research* 36, e105.
- Danecek, P. et al. (2011), The variant call format and VCFtools, *Bioinformatics* 27, 2156–2158.
- Thompson, J. D., Plewniak, F. and Poch, O. (1999), A comprehensive comparison of multiple sequence alignment programs, *Nucleic Acids Research* 27, 2682–2690.
- Weirather, J. L. et al. (2017), Comprehensive comparison of Pacific Biosciences and Oxford Nanopore Technologies and their applications to transcriptome analysis, *F1000Research* 6, 100.

## 6. Open Questions

1. **Error profiles from data.** The short-read model is linear in position; profiles learnt from real runs (ART's, or a sequencer's own quality calibration) would be more faithful.
2. **Paired reads** and their fragment lengths, for the mapper's pairing, are not yet simulated.
