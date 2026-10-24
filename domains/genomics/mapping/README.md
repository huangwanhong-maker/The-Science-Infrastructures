# mapping: Reads to a Genome

*Status: draft (2026-10-24).* A minimap2-like mapper of short and long reads to a genome: minimizer anchors (or BWA-MEM's super-maximal exact matches), chaining by dynamic programming with Li's gap cost, primary and secondary chains with their mapping quality, base-level alignment of the gaps between seeds and of the read's ends, soft clipping; SAM written and read; the accuracy of mappings against the simulator's truth. Serves use (the reads of the variant caller and of the atlas's Reads view) and learning (the theory note [mapping-and-variants.md](../theory/mapping-and-variants.md), §1–4).

## 1. Concepts

- **Seed, chain, align.** A read is placed in three steps, each cheaper than the next is costly (Li 2018).
  - *Seeds.* The read's minimizers, looked up in the genome's minimizer index, give exact matches of k bases. Alternatively, its super-maximal exact matches in the FMD-index give exact matches of any length (Li 2013).
  - *Chains.* A dynamic programme over the seeds, sorted along the genome, finds the chains of seeds that are colinear on both sequences. It scores each seed by the bases it adds and each gap between seeds by a concave cost.
  - *Alignment.* Only the best chains are aligned base by base, and only between their seeds and beyond their ends.
- **The strand mapped.** A seed on the reverse strand is turned into the coordinates of the read's reverse complement. Every chain and every alignment is then a forward one, of the read or of its reverse complement, as SAM writes it.
- **Mapping quality is a posterior made cheap.** A read whose best placement is barely better than its second is likely to be wrong.
  - Li's formula turns the ratio of the two chain scores into a Phred-scaled quality, refined here by the ratio of the alignment scores.
  - Seeded by SMEMs, the quality is BWA-MEM's, from the difference of the alignment scores.
  - Either way, a read inside an exact repeat gets 0.
- **Certified shortcuts.** A gap between seeds whose two sides have one length is matched without gaps when no gapped path can score more (a bound, not a guess). The other gaps go to the wavefront algorithm, or to bands certified optimal.
- **Deterministic in parallel.** Reads are mapped on calx's dynamic loop, each into its own slot. The mappings are the same on any number of threads, and the index is the same whatever its pool.

## 2. Content

| Header | Content |
|---|---|
| `chain.hpp` | `Seed{contig, reverse, target, query, span}`, ordered by its ends, and `seeds_of(anchors, read_length, k)` (reverse anchors to the reverse complement's coordinates). `ChainOptions{max_gap, bandwidth, max_predecessors, min_anchors, min_score, mask_level, secondary_ratio, best_n}`. `gap_cost(l, w̄)`: Li's γ in exact integers. `chain_step`, and `chain_dp(seeds, options)` → `ChainDp{f, predecessor}`: eqs. 1–2 on the seeds' ends, with the predecessor limit h and the maximum gap G. `chain(seeds, options)` → `Chain{contig, reverse, score, seeds, spans, primary, parent, sub_score, mapq}`: backtracking, §2.1.3. `select_chains(chains, options)`: fragments dropped, then primaries and secondaries by overlap on the query, §2.1.4. `mapping_quality(f₁, f₂, m)` = 40 (1 − f₂/f₁) min(1, m/10) ln f₁ |
| `mapper.hpp` | `Seeding{minimizers, smems}`. `MapperOptions{seeding, index, smem_min_length, smem_max_occ, chain, scoring, drop, clip_penalty, secondary, gap_wfa_limit}`, with three presets: `short_reads()` (k = 21, w = 11, minimap2's sr scoring), `short_reads_smem()` (SMEMs of at least 19 bases occurring at most 500 times) and `long_reads()` (k = 15, w = 10, map-ont). `Mapping{read, contig, target_begin, target_end, reverse, mapq, primary, supplementary, query_length, alignment, nm, chain_score, seeds}`, with `cigar()` (soft clips added), `clip_left()` and `clip_right()`. `Mapper(pool, contigs, options)` has `seeds(read)`; `map(read, id)` → the primary first, then supplementary parts, then secondaries; `map_all(pool, reads, grain)` in read order; and `align_chain(chain, strand)`. The gaps between seeds are aligned ungapped when certified, else by the wavefront algorithm or certified bands. The ends are aligned by X-drop and Z-drop, end to end when clipping gains less than the clipping penalty |
| `sam.hpp` | `sam_header(names, lengths, program, command_line)` (@HD, @SQ, @PG). `SamRecord{qname, flag, rname, position, mapq, cigar, rnext, pnext, tlen, seq, qual, tags}` and `SamTag`. `sam_record(mapping, name, read, quality, names)`: FLAG 0x10, 0x100, 0x800; POS 1-based; SEQ and QUAL on the strand mapped; NM, AS, tp, s1. `sam_unmapped`, `sam_line`, `write_sam`. `read_sam(text, problems)` → `Sam{header, references, records}`, never throwing on content. `mapping_of(record, names, mapping)` |
| `evaluate.hpp` | `Truth{contig, begin, end, reverse}` and `truth_of(read)`. `correct(mapping, truth)`: the right contig and strand, overlapping the true interval. `primary_of(mappings)`. `accuracy(mappings, truth, thresholds)` and `accuracy_of(mappings, reads)` → `Accuracy{reads, mapped, correct, unmapped, curve}`, with `AccuracyPoint{mapq, mapped, wrong, fraction_mapped, error_rate}` (Li 2018, Fig. 1). `Throughput{reads, bases, seconds}` |

## 3. Dependencies

- **Domain:** `genomics::index` (the minimizer index and its anchors, the FMD-index and its SMEMs); `genomics::align` (scorings, CIGARs, `rescore`, the wavefront algorithm, banded alignment, X-drop and Z-drop); `genomics::sequence` (codes, reverse complements).
- **Foundations:** `calx::scheduler` (the pool, the dynamic loop).
- **Checks:** they also use `genomics::simulate`, `calx::random` and `calx::verification`.

## 4. Checks

`node domains/genomics/tests/run-tests.mjs`, or one program at a time from the build (`check_genomics_mapping_<topic>`).

| Check | What it verifies |
|---|---|
| `chain` | The integer gap cost: symmetric, never above Li's γ(l) = 0.01 w̄ \|l\| + 0.5 log₂\|l\|, and within 1.5 of it for \|l\| ≤ 5000. f(i) of the DP equal to the best chain over every subset of seeds (300 random sets of 4–12 seeds), and the predecessor limit h = 1. Backtracked chains disjoint, increasing on both sequences, of the minimum seeds and score. Two equal copies one primary and one secondary of mapping quality 0; a unique placement 60; Li's formula at (100, 80, 5) = 18. Fragments dropped; two parts of a split read both primary; a weak secondary dropped but kept as f₂. Reverse seeds exact forward matches of the reverse complement |
| `mapper` | A simulated genome of 750 kb on two contigs: diverged repeats, tandem arrays, and three exact copies of 2.5 kb, one reverse-complemented. 20 000 short reads (150 b, 1% substitutions), seeded by minimizers and by SMEMs: those outside repeats correct above 0.995 (1.0000 for both); all of them above 0.98 (0.9920); every read inside an exact copy of MAPQ ≤ 3 (195 of 195, mean 0.00); MAPQ ≥ 20 wrong less than once in a thousand (1.5 · 10⁻⁴ and 2.0 · 10⁻⁴ at MAPQ ≥ 20, for 98.2% and 98.8% of the reads). 300 long reads (8% errors): correct above 0.97 (1.0000), and under 2% of their bases clipped (0.06%). Every CIGAR consistent with its read and interval, rescored to its AS, with NM its edit distance; the same mappings on one and eight threads |
| `sam` | The SAM specification's example (§1.1): six records; POS; CIGARs read back as written, with query lengths SEQ's; FLAGs; NM recomputed from the reference (3 for r001, and 1 for its mate, as its NM:i:1 says). 650 mappings written and read back the same (contig, POS, strand, MAPQ, flags, CIGAR with clips, AS, NM, tp), with SEQ and QUAL on the strand mapped, and unmapped reads. Malformed lines reported with their line numbers |
| `ecoli` | Escherichia coli K-12 MG1655 (NC_000913.3; skipped without the data cache). 40 000 short reads, seeded both ways: correct above 0.98 (0.9870; the rest lie in the rRNA operons and insertion sequences, at MAPQ 0), and wrong at MAPQ ≥ 20 less than once in a thousand (5.1 · 10⁻⁵ and 2.5 · 10⁻⁵). 400 long reads of 8 kb: correct above 0.98 (1.0000). The same mappings on one and eight threads; the time of each index and of the mapping |

## 5. Performance

Measured on the development machine (i7-11370H, WSL, Release, `-march=native`), 2026-10-24, with other builds running (load 3 to 4 of 8), so the ranges are of several runs.

On a 2 Mb simulated genome, one thread, a 150-base read takes 13 µs: 6 µs of seeding, 2 µs of chaining, the rest alignment. That is 75 000 reads a second. An 8 kb read at 8% errors takes 1.2 ms: 6.6 Mb/s, of which seeding is a third. The ungapped certificate for gaps between seeds halved the time of a short read, from 39 µs, and cut that of a long read by a third.

| Input (in the checks) | Index | 1 thread | 8 threads |
|---|---|---|---|
| Escherichia coli, 150-base reads, minimizers (k = 21, w = 11) | 0.11–0.15 s on 8 threads | 25 000–32 000 reads/s | 100 000–206 000 reads/s (15–31 Mb/s) |
| Escherichia coli, 150-base reads, SMEMs | 2.2–2.9 s (FMD-index, one thread) | 16 000–22 000 reads/s | 95 000–176 000 reads/s (14–26 Mb/s) |
| Escherichia coli, 8 kb reads at 8% errors (k = 15, w = 10) | 0.12–0.18 s | 460–640 reads/s (3.7–5.2 Mb/s) | 2 200–3 000 reads/s (17–24 Mb/s) |
| Simulated 750 kb, 150-base reads, minimizers | 0.02–0.05 s | 27 000–41 000 reads/s | 115 000–183 000 reads/s |
| Simulated 750 kb, 6 kb reads | — | 650–1 020 reads/s | 1 700–5 000 reads/s (10–31 Mb/s) |

## 6. References

- Li, H. (2018), Minimap2: pairwise alignment for nucleotide sequences, *Bioinformatics* 34, 3094–3100.
- Li, H. (2013), Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM, arXiv:1303.3997.
- Li, H. (2012), Exploring single-sample SNP and INDEL calling with whole-genome de novo assembly, *Bioinformatics* 28, 1838–1844 (SMEMs and the FMD-index).
- Li, H., Handsaker, B., Wysoker, A. et al. (2009), The Sequence Alignment/Map format and SAMtools, *Bioinformatics* 25, 2078–2079.
- The SAM/BAM Format Specification Working Group, *Sequence Alignment/Map Format Specification*, version 1.6.
- Li, H., Ruan, J. and Durbin, R. (2008), Mapping short DNA sequencing reads and calling variants using mapping quality scores, *Genome Research* 18, 1851–1858.
- Marco-Sola, S., Moure, J. C., Moreto, M. and Espinosa, A. (2021), Fast gap-affine pairwise alignment using the wavefront algorithm, *Bioinformatics* 37, 456–463.
- Roberts, M., Hayes, W., Hunt, B. R., Mount, S. M. and Yorke, J. A. (2004), Reducing storage requirements for biological sequence comparison, *Bioinformatics* 20, 3363–3369.
- Zhang, Z., Schwartz, S., Wagner, L. and Miller, W. (2000), A greedy algorithm for aligning DNA sequences, *Journal of Computational Biology* 7, 203–214.

## 7. Open Questions

1. **Re-seeding of SMEMs.** BWA-MEM re-seeds long SMEMs to find the shorter matches that a repeat hides, and samples the occurrences of repetitive ones; neither is written. A read whose SMEMs all occur more than 500 times is unmapped.
2. **Two-piece affine gaps and the Z-drop inside chains.** minimap2 aligns long gaps with a second, flatter piece of the gap cost, and breaks an alignment at an inversion by the Z-drop between anchors. Here there is one affine piece, and the Z-drop at the ends only.
3. **Paired-end reads.** Mates, their insert sizes and the rescue of an unmapped mate are not handled; every read is single-end.
4. **Spliced alignment** (minimap2's splice mode) for transcripts belongs with the expression module.
5. **Faster seeding.** `sequence::minimizers` runs at about 18 ns a base (147 µs for 8 kb), on a `std::deque`; a ring buffer of w entries would make seeding several times faster. The gaps between seeds are short and many, and an inter-sequence kernel on calx lanes (align's planned `inter.hpp`) could align eight or sixteen at a time.
