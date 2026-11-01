# simulate: Genomes, Individuals and Reads with Their Truth

*Status: draft (2026-10-24).* Seeded simulation of what the genome's algorithms are checked on: genomes with repeats, diploid individuals mutated from a reference, reads from them with a sequencer's errors, and families of sequences with their true alignment. A simulation always comes with its truth, which is why it is here: a real genome tells no one where its reads came from.

## 1. Concepts

- **A simulation is a function of its seed.** Every draw comes from a `calx::random::Stream` keyed by the seed, a purpose and an object; read i comes from its own stream. Reads a to b are therefore the same whether drawn alone or within reads 0 to n, so that a shard of a simulated sample is two numbers, and a worker on another machine regenerates it exactly.
- **Truth travels with the data.** A diploid individual carries its variants (VCF-style, indels anchored on the base before them) and, for each haplotype, a liftover to the reference; a read carries its contig, its reference interval, its strand, its haplotype and its number of errors; a family carries its true alignment and its tree.

## 2. Content

| Header | Content |
|---|---|
| [`genome.hpp`](include/genomics/simulate/genome.hpp) | Random sequence of a GC content; genomes with interspersed repeats (diverged copies on both strands) and tandem arrays, with where each copy lies |
| [`mutate.hpp`](include/genomics/simulate/mutate.hpp) | A diploid individual: SNVs and indels at given rates, heterozygous or homozygous, never overlapping; its two haplotypes; liftovers to the reference (`Liftover` of blocks, each with its contig and strand, and `locate(h)` → `Locus{contig, position, reverse, inserted}`) |
| [`reads.hpp`](include/genomics/simulate/reads.hpp) | Short reads with substitution errors rising along the read, a few one-base indels, informative Phred qualities; long reads with gamma lengths and indel-dominated errors; from one or two copies of several contigs; as FASTQ |
| [`structural.hpp`](include/genomics/simulate/structural.hpp) | A diploid with structural variants (ADR 0064): `SvType` (deletion, tandem and dispersed duplication, inversion, novel insertion, balanced reciprocal translocation, aneuploidy, copy-neutral LOH), `StructuralVariant` heterozygous or homozygous; `make_structural_sample(contigs, variants, rates, seed, names)` → `StructuralSample{diploids, truth, chromosomes, junctions, copy_number, novel}`: chromosomes as words of oriented `Piece`s of the haplotypes (with mutate.hpp's small variants), their sequences and liftovers of oriented blocks of any contig, their `Junction`s (breakends and novel bases), the copy number of every reference base by haplotype (`copies_at`); `validate_structural`; `random_structural(contigs, StructuralSpec, seed)`, kept apart from one another, from the ends and from runs of N |
| [`pairs.hpp`](include/genomics/simulate/pairs.hpp) | Paired-end reads from rearranged chromosomes: fragments of a normal length law, FR orientation, an optional GC bias of whole fragments (Benjamini and Speed 2012), none touching N; each mate a short read of reads.hpp's model with its chromosome interval, strand, reference loci of its ends and whether it crosses a junction; `pairs_for_coverage` |
| [`family.hpp`](include/genomics/simulate/family.hpp) | Sequences evolved from a root along a balanced tree with substitutions and indels, with their true alignment (every base keeps the column it was born in) and the tree in Newick |
| [`tumour.hpp`](include/genomics/simulate/tumour.hpp) | A tumour derived from a normal diploid (ADR 0064): a trunk and subclones with their cancer cell fractions (the sum rule enforced); substitutions drawn from SBS-96 signature profiles (`synthetic_signatures`, five mechanism-motivated profiles, not COSMIC's) and indels; a clonal allele-specific copy-number profile (`random_copy_number`: gains, losses, copy-neutral LOH, amplifications, whole-genome doubling); early mutations at every copy of their allele, late and subclonal ones at one; microsatellites that slip; purity; tumour reads whose depth follows copy number and whose molecules carry what their clone carries; the normal's source; cohorts of catalogues by the Poisson model of signatures (`simulate_catalogues`) |

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
| structural (`check_genomics_simulate_structural`) | Every base of every rearranged chromosome lifted to the reference is the base it was taken from (0 of 1 204 900 wrong), the inserted bases exactly the novel insertion's; the copy number of every reference base by brute force over the lifted bases equals the reported step function, by haplotype (1+0, 0+0, 1+3, 2+0 after LOH); lengths change by what the variants add and remove; each kind makes its junctions with the breakends of the theory note (§1.4), LOH none; a trisomy and a monosomy; lifted bases differ from the reference at 0.7 · 10⁻³ under SNVs at 10⁻³; random variants as many as asked, valid, clear of N; pairs: fragments 450 ± 45, FR, mates their chromosome's bases but for their errors, junction crossings as predicted, shards equal, GC bias toward its optimum, no fragment touching N |
| tumours (`check_genomics_simulate_tumour`) | A tree against the sum rule refused; the clones' own fractions; every mutation's expected VAF ρφm/(2(1 − ρ) + ρN); early mutations at every copy of their allele, none on a lost one; each clone's substitutions following its signature (cosine ≥ 0.98) in their reference trinucleotides; tumour depth by segment following 2(1 − ρ) + ρN within 3%; the reads carrying the substitutions against Σ f · depth; a shard of tumour reads equal to the same reads of the whole sample; cohorts with the model's Poisson means |

## 5. References

- Dohm, J. C., Lottaz, C., Borodina, T. and Himmelbauer, H. (2008), Substantial biases in ultra-short read data sets from high-throughput DNA sequencing, *Nucleic Acids Research* 36, e105.
- Danecek, P. et al. (2011), The variant call format and VCFtools, *Bioinformatics* 27, 2156–2158.
- Thompson, J. D., Plewniak, F. and Poch, O. (1999), A comprehensive comparison of multiple sequence alignment programs, *Nucleic Acids Research* 27, 2682–2690.
- Weirather, J. L. et al. (2017), Comprehensive comparison of Pacific Biosciences and Oxford Nanopore Technologies and their applications to transcriptome analysis, *F1000Research* 6, 100.

## 6. Open Questions

1. **Error profiles from data.** The short-read model is linear in position; profiles learnt from real runs (ART's, or a sequencer's own quality calibration) would be more faithful.
2. **Paired reads** are simulated (`pairs.hpp`), FR only; mate-pair (RF) libraries, adapter read-through and PCR duplicates are not.
