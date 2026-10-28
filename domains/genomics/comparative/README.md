# comparative: Whole-Genome Alignment by Maximal Unique Matches

*Status: draft (2026-10-24).* Two genomes aligned whole by MUMmer's method (Delcher et al. 1999): their maximal unique matches on both strands, the heaviest chain of matches that keeps the order of both genomes, the gaps between its anchors, the collinear blocks forward and reversed, and the dot plot. It serves use, to compare strains, species and assemblies, and learning, as the meeting of suffix arrays with increasing subsequences. The theory is [theory/signals-classification-and-comparison.md](../theory/signals-classification-and-comparison.md), §4.

## 1. Concepts

- **A MUM occurs once in each genome and extends at neither end.** Found in one pass over the suffix array and LCP array of the two genomes' runs free of N, each run with a separator of its own: an LCP peak of exactly two suffixes, one from each genome, whose preceding bases differ. The suffix array and LCP array are morphe's (SA-IS, Kasai), not written again.
- **MUMs never nest**, so their starts in each genome are distinct, and a chain increasing in both genomes' starts increases in their ends.
- **The chain is a heaviest increasing subsequence.** Sorted by their starts in A, the MUMs whose starts in B also increase are an alignment that keeps both orders; the heaviest, by length, is found by morphe's `heaviest_chain` in O(m log m). The mapper's chaining (`mapping/chain.hpp`) does not fit: its DP looks back over a bounded number of predecessors, within a gap cut-off and a band, as reads need; two genomes need the global chain, whose gaps may be as long as a rearrangement.
- **The gaps say what changed**: one base in each genome is a single-base substitution, equal lengths a stretch of substitutions, unequal lengths an indel. A reversed segment is a chain of the reverse strand. On a 200 kb genome against a copy with a 30 kb inversion and 232 SNVs and indels, the chains give three blocks, the reversed one the inversion to the base, every isolated SNV as a single-base gap and every isolated indel as a gap of its length.
- **Circular genomes begin where their records begin.** The mitochondria of human and chimpanzee align along one forward chain over 0.964 of each, every MUM on one diagonal modulo the chimpanzee's length; their records begin 582 bases apart on the circle, and with the chimpanzee's rotated by that offset the chain spans 0.999 of both.

## 2. Content

| Header | Content |
|---|---|
| [`mums.hpp`](include/genomics/comparative/mums.hpp) | `Match{a, b, length, reverse}` with `a_end`, `b_end`; `maximal_unique_matches(A, B, min_length, reverse)` (morphe's `concatenate`, SA-IS and LCP array) |
| [`whole_genome.hpp`](include/genomics/comparative/whole_genome.hpp) | `GapKind{snp, substitutions, indel}`; `Gap{a, b, length_a, length_b}` with `kind()`, `indel()`; `Chain{reverse, anchors, gaps, weight}`; `strand_b`; `chain_matches(matches, length_b)` (morphe's `heaviest_chain`); `Block{reverse, a_begin, a_end, b_begin, b_end, anchors, matched}`, `blocks(chain, max_gap)`; `GenomeOptions{min_length, max_gap, both_strands}`; `align_genomes(A, B, options)` → `GenomeAlignment{forward, reverse, forward_chain, reverse_chain, blocks}` with `coverage_a`, `coverage_b`; `DotSegment`, `dot_plot(alignment)` |

## 3. Dependencies

`genomics::sequence` (codes, reverse complements); morphe `formal` (concatenation, SA-IS, the LCP array) and `combinatorics` (heaviest increasing subsequences). The checks use `genomics::simulate` and calx `random`.

## 4. Checks

| Check | What it verifies |
|---|---|
| `mums` (5) | On 400 random pairs sharing segments on both strands, with repeats and bases N, the MUMs of each strand equal those found by brute force (6270 in all); every MUM an exact match on its strand, free of N; no MUM's interval contains another's in A or in B; a 5 kb genome against itself one forward MUM, itself, and no reversed one |
| `genomes` (12) | The chain against the quadratic recurrence on 40 related pairs with moved segments. 200 kb against a copy with an inversion of [60 000, 90 000) and then SNVs and indels (`simulate`, with the truth of its liftover): the bases of the chained anchors where the truth puts them but for 31 of 199 760 within 8 bases of an anchor's end (bases matching by chance across a breakpoint or beside an insertion, and indels placed elsewhere in a repeat); three blocks, forward, reversed and forward, the reversed one [60 000, 90 000) (within 10 asserted, exact measured); 152 isolated SNVs each a single-base gap at its place, the 178 single-base gaps all at SNVs, 38 isolated indels each a gap of its length; the dot plot one segment a MUM at 45°, the chained marked; coverage 0.999. With the cache, human (NC_012920.1) against chimpanzee (NC_001643.1): MUMs of 15 and of 20, one forward chain spanning at least 0.95 of each (0.964), every forward MUM on its diagonal modulo the chimpanzee's length, no reversed MUM; rotated by the chain's offset (−582), at least 0.99 (0.999) |

## 5. References

- Delcher, A. L., Kasif, S., Fleischmann, R. D., Peterson, J., White, O. and Salzberg, S. L. (1999), Alignment of whole genomes, *Nucleic Acids Research* 27, 2369–2376.
- Jacobson, G. and Vo, K.-P. (1992), Heaviest increasing/common subsequence problems, *CPM 1992*, LNCS 644 (morphe's `increasing.hpp`).
- Nong, G., Zhang, S. and Chan, W. H. (2011) and Kasai, T. et al. (2001), as morphe's `suffix_array.hpp` cites them.

## 6. Open Questions

1. **Local chains.** One heaviest chain a strand recovers one collinear backbone a strand. Translocations, several inversions and duplications need clusters of nearby anchors chained locally and then joined, as MUMmer 3's nucmer does; the mapper's chaining with genome-sized gap limits is a starting point.
2. **Circular genomes.** The rotation that closes the circle is the offset of the main diagonal modulo |B|; an option could find it from a first pass and report coordinates modulo the circle.
3. **The gaps' bases.** A gap is classified by its lengths; aligning its bases (`align/`, Gotoh or the wavefront algorithm) would count its substitutions and place its indels, as MUMmer's alignment of the gaps does.
4. **Maximal exact matches that repeat.** MUMs drop matches that occur twice in a genome; for repetitive genomes, maximal matches unique in one genome only (MUMmer's options) or all maximal exact matches would anchor more of the alignment.
