# variant: Variants of a Diploid

*Status: draft (2026-10-24).* The calling of single-nucleotide variants and short indels of a diploid from its mapped reads: pileups, genotype likelihoods and their posterior (Li 2011), filters on depth and strand bias; VCF 4.2 written and read; the normalisation of indels (Tan, Abecasis and Kang 2015); the evaluation of calls against a truth set. Serves use (the calls of the atlas's Reads view and the genotypes of the evolution domain's population genetics) and learning (the theory note [mapping-and-variants.md](../theory/mapping-and-variants.md), §5–8).

## 1. Concepts

- **A site is a small inference.** At each position the reads give bases with their error probabilities; a diploid has one of three genotypes there (or six, with two alternate alleles). The likelihood of each genotype is a product over the reads, the prior comes from the heterozygosity θ, and the posterior gives the genotype, its quality GQ and the site's quality QUAL = −10 log₁₀ P(0/0 | D).
- **An indel has many names.** In a repeat, one insertion or deletion can be written at several positions with more or less context. Two variants are equivalent when they make the same haplotype; the left-aligned, parsimonious representation is the canonical member of each class (theory §7). Reads are normalised the same way: the gaps of each alignment are moved left without changing its score, so that every read carrying an indel reports it at one place.
- **Reads that cannot tell.** A read ending inside a repeat looks the same with an indel or without it; only reads spanning the indel's whole interval of ambiguity count for the reference allele.
- **Deterministic in parallel.** The genome is called in fixed windows on calx's dynamic loop, each window's calls in its own slot.

## 2. Content

| Header | Content |
|---|---|
| `normalize.hpp` | `normalize(position, ref, alts, reference)` → `Normalized{position, ref, alts, changed, ref_matches}`: Tan et al.'s Algorithm 1 (trim the right end, extend left on an empty allele, trim the left end), any number of alleles, VCF's padding after the event at the first base; `ref_agrees`; `apply(reference, position, ref, alt)` (the haplotype); `left_align_gaps(alignment, query, target)`: every gap of a path moved left while the score stays, idempotent |
| `genotype.hpp` | `BaseCall{base, quality}`; `log10_likelihoods(bases, alleles)` (Li's eq. 2 for two alleles, its n-allele form for more; VCF's order `genotype_index(j, k)` = k(k + 1)/2 + j); `log10_likelihoods(ref_count, alt_count, error)` (the counting form for indels); `genotype_prior(alleles, θ)` (1 − 3θ/2, θ, θ/2; θ² between two alternates); `posterior(likelihoods, prior)` → `GenotypeCall{best, allele_a, allele_b, log10_posterior, gq, qual, pl}`; `fisher_exact(a, b, c, d)` (two-sided, from log Γ), `strand_bias` (FS) |
| `pileup.hpp` | `PileupOptions{min_base_quality, min_mapping_quality, default_quality, left_align}`; `ReadView{bases, quality}` and `read_views(reads)`; `pileup(reference, contig, begin, end, mappings, reads, options)` → `Pileup{offset, bases, indels, spans}`: `Observation{base, quality, mapq, reverse, read}` per position in one flat array, `IndelObservation{position, length, bases, reverse, mapq, read}` at their anchors, `ReadSpan`s; `at(p)`, `depth(p)`, `for_each_spanning(b, e, f)`, `spanning(b, e)` |
| `caller.hpp` | `CallerOptions{pileup, theta, theta_indel, indel_error, min_alt, min_alt_fraction, second_allele_fraction, min_qual, min_depth, max_depth, max_strand_bias, max_strand_bias_indel, window}`; `Call{contig, position, ref, alts, qual, filters, depth, allele_depths, genotype, gq, pl, strand_bias, indel}`; `call_region(reference, contig, begin, end, mappings, reads, options)`; `call_variants(pool, contigs, mappings, reads, options)` (windows in parallel, calls in order); `caller_header(names, lengths, sample)` and `vcf_record(call, names)` (INFO DP, FS; FORMAT GT:GQ:DP:AD; FILTER PASS, LowDepth, HighDepth, StrandBias) |
| `vcf.hpp` | `VcfGenotype{alleles, phased}`, `parse_genotype`, `genotype_string`; `VcfRecord{chrom, position, id, ref, alts, qual, filters, info, format, samples, line}` with `info_value`, `sample_value`, `genotype(sample)`, `passed`; `VcfHeader{fileformat, meta, contigs, samples}`; `read_vcf(text, problems)` → `Vcf{header, records}`, never throwing on content; `vcf_line`, `write_vcf` |
| `evaluate.hpp` | `VariantType` (SNV, insertion, deletion, other), `type_of`; `TruthVariant{contig, position, ref, alt, dosage}`, `truth_variant`, `truth_set` (from `simulate::Variant`s); `compare(calls, truth, contigs, pass_only)` → `Comparison{snv, insertion, deletion, other, all, false_positives, false_negatives}` of `Counts{tp, fp, fn, genotype_match}` with precision, recall, F₁ and genotype concordance, both sides normalised, multi-allelic calls split |

## 3. Dependencies

`genomics::mapping` (the `Mapping`s it piles up), `genomics::align` (CIGARs and alignments), `genomics::sequence` (codes, the Phred scale, `Problem`); `calx::special_functions` (log Γ for Fisher's test), `calx::scheduler` (the pool, the dynamic loop). The checks also use `genomics::simulate`, `calx::random` and `calx::verification`.

## 4. Checks

`node domains/genomics/tests/run-tests.mjs`, or one program at a time from the build (`check_genomics_variant_<topic>`).

| Check | What it verifies |
|---|---|
| `normalize` | A CA deleted from CACACA, written in ten ways (the problem of Tan et al.'s Fig. 1), one haplotype and one normal form (2, GCA, G); an insertion at the repeat's end moved to its start; substitutions trimmed; a multi-allelic indel normalised as one; VCF's padding at the first base; a REF that is not the reference's reported; 400 random indels in random tandem repeats, all 40 649 representations of their haplotypes normalising to one form, which makes the haplotype and is the shortest representation of least position (brute force); gaps of alignments left-aligned (9M2D5M → 5M2D9M), scores kept, idempotent, no gap movable further |
| `genotype` | Li's eq. 2 against its closed forms (three reference bases at Q20 and one alternate at Q30; mixed qualities, a third allele left out); the counting form equal to eq. 2 at equal qualities; the three-allele form reducing to eq. 2, L(1/2) by hand; VCF's order of genotypes; the prior; QUAL, GQ, the posterior and PL against direct computation; 1000 reads in log space; Fisher's test on Fisher's tea-tasting table (p = 34/70) and on 400 tables against exact enumeration in integers; FS |
| `vcf` | The VCF 4.2 specification's example (§1.1): 5 records, 3 samples, 17 meta lines, POS, ALT ('.' and two), QUAL, FILTER, INFO flags, genotypes with their phasing, dropped trailing fields; written back to the same text; GT values (missing, haploid, partly missing, triploid, alleles above 9) and malformed ones; calls written and read back; malformed lines reported with their line numbers |
| `caller` | A 500 kb diploid with 509 SNVs and 198 indels at 30× of 150-base reads with 1% substitutions, mapped by the mapper: SNVs precision above 0.99 (measured 1.0000), recall above 0.98 (0.9961), genotypes above 0.98 (1.0000); insertions and deletions precision and recall above 0.9 (1.0000; 0.988 and 0.991), genotypes above 0.9 (0.977 and 1.000); all variants F₁ 0.997; the calls the same on one and eight threads; the pileup's mean depth 30.0; the reads of 188 of 196 truth indels reporting them at their normalised place; the calls written as VCF and read back |

## 5. Performance

Measured on the development machine (i7-11370H, WSL, Release, `-march=native`), 2026-10-24, with other builds running: the 500 kb diploid at 30× (100 000 reads of 150 bases) is mapped in 0.6–0.9 s and called in 0.15–0.19 s on eight threads; the pileups and calls cost a fifth of the mapping.

## 6. References

- Li, H. (2011), A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data, *Bioinformatics* 27, 2987–2993.
- Tan, A., Abecasis, G. R. and Kang, H. M. (2015), Unified representation of genetic variants, *Bioinformatics* 31, 2202–2204.
- Danecek, P., Auton, A., Abecasis, G. et al. (2011), The variant call format and VCFtools, *Bioinformatics* 27, 2156–2158; *The Variant Call Format Specification*, VCFv4.2 (samtools/hts-specs).
- DePristo, M. A., Banks, E., Poplin, R. et al. (2011), A framework for variation discovery and genotyping using next-generation DNA sequencing data, *Nature Genetics* 43, 491–498.
- Van der Auwera, G. A., Carneiro, M. O., Hartl, C. et al. (2013), From FastQ data to high-confidence variant calls: the Genome Analysis Toolkit best practices pipeline, *Current Protocols in Bioinformatics* 43, 11.10.1–11.10.33.
- Fisher, R. A. (1935), *The Design of Experiments*, Oliver and Boyd, ch. 2.
- Agresti, A. (1992), A survey of exact inference for contingency tables, *Statistical Science* 7, 131–153.
- Krusche, P., Trigg, L., Boutros, P. C. et al. (2019), Best practices for benchmarking germline small-variant calls in human genomes, *Nature Biotechnology* 37, 555–560.

## 7. Open Questions

1. **Read likelihoods by the pair HMM.** The design note's variant module scores reads against candidate haplotypes with the pair HMM (`align/pairhmm.hpp`), as GATK's HaplotypeCaller does; here the pileup's per-base model and a fixed indel error. Local reassembly or haplotype realignment would weigh an indel and the SNVs near it together.
2. **Base quality recalibration** is assumed done: the qualities are taken as calibrated, which the simulator's are.
3. **Several samples.** Joint calling of a cohort (Li 2011 estimates allele frequencies by EM and the site frequency spectrum from many samples, and uses them as the prior) would make the genotypes the evolution domain reads from many samples at once; here one sample.
4. **A hypergeometric law in calx** (pmf, cdf, tails) would serve Fisher's test here and enrichment tests elsewhere; Fisher's test is built from `log_choose`.
5. **Haplotype-level evaluation** (vcfeval, hap.py) would match complex variants written differently in the two sets; here alleles are matched after normalisation.
