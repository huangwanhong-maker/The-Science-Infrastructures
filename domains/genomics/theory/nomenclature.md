# The Description of Sequence Variants

*Status: draft (2026-11-01).* The theory behind `genomics/nomenclature`: the HGVS Nomenclature (version 21.1) as a formal language; what a description means; the normal forms of HGVS, VCF and SPDI as choices from one equivalence class; the projection of a genomic change through a transcript to its protein; and the places where the standard is ambiguous, underspecified or inconsistent. Learning and research.

## 1. Concepts

A genome that departs from its reference is described relative to it. The HGVS Nomenclature, maintained by the HGVS Variant Nomenclature Committee under HUGO (den Dunnen et al. 2016; Hart et al. 2024), is the international standard for sequence variants; VCF (samtools/hts-specs) and NCBI's SPDI (Holmes et al. 2020) are two data models for the same objects; ISCN (the karyotype module) describes chromosomes. Three questions organise the theory: what the well-formed descriptions are (a language), what each denotes (a set of haplotypes), and which description of a denotation the standard chooses (a normal form).

## 2. The Language

**Definition 1 (alphabet and grammar).** Let Σ be the bytes. The grammar G of `grammar.hpp` is a context-free grammar over Σ with start symbol `hgvs`; L(G) ⊂ Σ* is its language. Its rules fall into a head (descriptions, references, counts, separators, the protein family) and four instances of one nucleotide template, for the coordinate types {g, m, o}, c, n and r, which differ in their positions (plain numbers, pter, qter, cen for g., m., o.; UTR signs and intronic offsets for c. and r.; offsets alone for n.) and in their alphabet (IUPAC capitals for DNA, lower case for RNA; standards.md).

The structure of a description, from the outside in:

$$
\text{hgvs} = d\,(\texttt{\^{}}\,d)^* \mid [d(;d)^*](;[d(;d)^*])^*,\qquad d = (\text{ref}:)^?\,t\,.\,\ell,
$$

where ℓ is a list of alternatives `^` of groups joined by `(;)` (unknown phase), a group is alleles in trans `[…];[…]` or a chain, an allele is chains separated by `;` (cis), `(;)` or `,` (transcripts of one allele), a chain is items joined by `::`, and an item is an edit, a predicted choice `(e^e)`, or a mixture `=/e`. An edit is a location and a change: substitution `>`, `del`, `dup`, `ins`, `delins`, `inv`, `=`, a repeated sequence or copy number `[n]`, `sup`, methylation `|gom`; for proteins also `fs`, `ext` at either end and `Met1?`.

**Proposition 1 (unambiguity on the standard's examples).** Every one of the 702 descriptions the standard writes as examples (hgvs-examples.tsv) has exactly one derivation in G.

*Proof.* By computation: the checks count the derivations of each example over the compiled grammar with a memoised recursion on spans (ε-rules included, saturated at 2), independent of Earley's parser; each count is 1. ∎

*Remark.* Unambiguity of a context-free grammar is undecidable in general; the grammar is built to make the usual sources impossible: a predicted list `[(a;b)]` requires a separator, so that `[(a)]` is only a predicted item; a range of residues inserted (`insGly4_Ser6`) requires two positions, so that `insTer12` is only a stop at a length; three-letter and one-letter sequences cannot be mixed; an identifier has two characters at least, so `c.260` is not a reference. 600 random descriptions of every kind confirm it beyond the examples.

**Definition 2 (conditions of meaning).** L_HGVS ⊂ L(G) is cut out by conditions a context-free grammar does not express, checked on the value (`check`): a substitution names one position and two different bases; a range names two different positions, 5′ to 3′ (except across the origin of `o.` and `m.`); an insertion lies between adjacent positions; an inversion covers more than one; variants in cis that touch are one delins, and two changes in one codon one base apart are one delins; `=` stands on an allele only when it holds no variant; a frameshift starts at its first changed residue; a change of the initiating Met1 is not a substitution; a delins is minimal. Some conditions need the reference (the 3′ rule; positions beyond a transcript): they are conditions on the pair (description, reference), §3.

## 3. The Semantics

**Definition 3.** A *change* of a reference r ∈ {A, C, G, T}ⁿ is a triple (p, D, I) with D = r[p, p + |D|); it makes the haplotype h = r[0, p)·I·r[p + |D|, n). Two changes are *equivalent* when they make the same haplotype. A change is *minimal* when D and I share neither their first nor their last base. An exact g. description denotes a change (`change_of`): its location replaced by what its edit says; an uncertain one (`(a_b)`, `?`, `N[15]`) denotes the set of changes its readings allow; an allele in cis the composition of its changes, in trans a pair of haplotypes.

Let lcp and lcs be the lengths of the longest common prefix and suffix of r and h, and m = |h|.

**Theorem 1 (the class of a change).** Let (p, D, I) be a change with h ≠ r.
1. If lcp + lcs < min(n, m), the minimal changes making h are exactly one: (lcp, r[lcp, n − lcs), h[lcp, m − lcs)).
2. Otherwise h arises from a pure deletion (n > m) or insertion (n < m) of L = |n − m| bases, and its minimal changes are the placements p ∈ [p_min, p_max] with p_min = m − lcs, p_max = lcp (deletions) or p_min = n − lcs, p_max = lcp (insertions); they form an interval.

*Proof.* A change (q, D′, I′) makes h exactly when r[0, q) = h[0, q) and r[q + |D′|, n) = h[m − (n − q − |D′|), m), i.e. q ≤ lcp and n − q − |D′| ≤ lcs. A minimal one cannot extend its window by a common base on either side, so its window is maximal under these constraints: q is as large as possible and its suffix as long as possible. In case 1 the constraints q ≤ lcp and q + |D′| ≥ n − lcs have the unique tightest solution q = lcp, q + |D′| = n − lcs, and both sides are non-empty because lcp + lcs < min(n, m). In case 2 one side may be empty: for a deletion |I′| = 0 and |D′| = L, so q ≤ lcp and q ≥ n − L − lcs = m − lcs, every such q giving a minimal change (its sides are D′ and ε); the interval is non-empty because lcp + lcs ≥ m. The insertion is symmetric with |D′| = 0. ∎

The sliding computation of `placements` (a deletion of r[p, p + L) moves right while r[p] = r[p + L]) reaches the same interval, since each step preserves the two equalities above.

**Theorem 2 (three normal forms).** For a pure indel, (i) VCF's normal form (Tan, Abecasis and Kang 2015, Algorithm 1) is the placement p_min, written with the base before it; (ii) the HGVS Nomenclature's 3′ rule ("the most 3′ position possible of the reference sequence is arbitrarily assigned to have been changed", general.md) is p_max; (iii) SPDI's contextual allele, by VOCA ("Correcting overprecision", Holmes et al. 2020), deletes the hull r[p_min, p_max + |D|) and inserts the haplotype over it.

*Proof.* (i) Tan et al. prove their form left-aligned (no equivalent representation has a smaller position) and parsimonious, i.e. the least placement, with VCF's padding base (variant/normalize.hpp, mapping-and-variants.md §7). (ii) The rightmost placement is the definition. (iii) VOCA's shrinking trims the common ends; its blossom repeats the bud periodically from the change's position, and the largest interval matching it is the union of all placements: a deletion of a period-L bud can slide exactly while r[i] = r[i + L], which is the blossom's agreement with r one period on; an insertion likewise. ∎

The checks verify both theorems by enumerating every minimal representation of 600 random changes in tandem repeats.

**Corollary (the HGVS description of an indel).** At p_max, an insertion of u is a duplication of r[p_max − L, p_max) exactly when p_max − p_min ≥ L: if u equals the L bases before p_max, inserting it before them gives the same haplotype, so p_max − L is a placement; conversely a span of L placements contains a full period. This is the standard's "a copy … inserted directly 3′ of the original copy" (DNA/duplication.md).

## 4. Transcripts

**Definition 4.** A transcript is a sequence of exons E₁, …, E_k (genomic intervals, in the transcript's order, on one strand), with a CDS [b, e) in transcript coordinates. The c. number of a transcript coordinate t is t − b + 1 for b ≤ t < e, −(b − t) for t < b, *(t − e + 1) for t ≥ e; there is no c.0 (numbering.md). An intronic base at distances d₅ from the upstream exon's last base and d₃ from the downstream exon's first is numbered from the upstream exon with +d₅ when d₅ ≤ d₃, from the downstream with −d₃ otherwise, so that the central base of an odd intron takes '+' (numbering.md: "the central nucleotide is numbered relative to the upstream exon followed by a +").

**Proposition 2 (strand symmetry).** Let ρ be reverse complementation of the genome (position g ↦ n − 1 − g) and T̄ the transcript with its exons mirrored and its strand reversed. Then the c. number of g under T equals that of n − 1 − g under T̄, and the projection of a change c under T equals that of ρ(c) under T̄.

*Proof.* The transcript coordinates of corresponding bases are equal by construction; the numbering depends on them alone. The projection describes the change on the transcript's strand, which is the same sequence in both cases. ∎ The checks verify both on 420 bases and 400 changes.

**The 3′ rule on a transcript** is the 3′ rule in the transcript's direction (general.md: it "applies to ALL descriptions … of a given variant"), so a minus-strand transcript takes the genomic p_min: the standard's g.32343183del is c.5697del. Placing the change on the genome, the exon–exon exception of numbering.md (c.3921del, not c.3922del) is not an exception: the intron separates the exons, and the change cannot slide across it. On a transcript-only reference it is imposed by cutting the sequence at the end of the exon holding the 5′-most placement.

**The protein.** The variant transcript is translated with the record's genetic code from the start codon to the first stop; P (reference, with its stop) and Q (observed) are compared as protein/*.md prescribes (`protein.hpp`): an extension when the first difference is at P's stop, a frameshift when the change's length is not divisible by three and the first difference precedes the stop (its new stop counted from the first changed residue as 1), in-frame changes by the class of Theorem 1 on the protein (the most C-terminal residue changed), a duplication by its corollary.

## 5. Structural Variants

A breakend is a base with the side of it a junction is on; an adjacency joins two breakends with any inserted bases, read across the junction; the adjacency read backwards swaps its breakends and reverse-complements the inserted bases, an involution whose orbits are the adjacencies themselves. VCF's four bracket forms are the four pairs of strands (§5.4); a symbolic <DEL>, <DUP:TANDEM> and <INV> are one, one and two adjacencies. An HGVS description of a rearrangement is the delins of the derivative chromosome, named by the chromosome whose centromere it keeps (DNA/complex.md: "variant descriptions are always in the forward orientation … determined by the chromosomal origin of the intact centromere").

## 6. Where the Standard Is Ambiguous, Underspecified or Inconsistent

The examples are from the standard's sources at commit 78c03ac7e4 (version 21.1); each was found by the checks.

1. **Its own syntax is not exact.** syntax.yaml's frameshift form `aa_position "fs" ["Ter" position]` leaves out the new residue its example has (p.Arg97ProfsTer23); its C-terminal extension `"Ter" aa_position aa "ext" …` repeats the residue (aa_position already holds it). grammar.md states that an exact EBNF is "expected in a future release".
2. **Repeat or duplication.** "When more than one additional copies are inserted directly 3′ of the original copy" a repeat is written (DNA/duplication.md); an inserted u = v^k (k ≥ 2) is one copy of u or k copies of v. Mutalyzer 3 takes the primitive root v (4_12CAG[5]); so does `describe`.
3. **Which copy an inserted range names.** g.17_18ins5_16 (DNA/insertion.md) names one of two occurrences of the inserted sequence; Mutalyzer 3 names the other (17_18ins1_12). The nearest is taken here, which agrees with the standard's example.
4. **Population knowledge decides the form.** p.Ala2[9] is "preferred over p.Ala11del when the repeat is variable in the population" (protein/repeated.md): the description depends on facts outside the sequences.
5. **The protein principle and its examples.** Protein descriptions "should be derived from comparing the variant protein sequence with the reference protein sequence" (general.md), yet how many residues a truncating delins replaces, p.(Pro578_Lys579delinsLeuTer) for c.1733_1735delinsTTT, is not determined by the two proteins: any number of residues after the first change could be said deleted, since Q stops. The standard's examples count the codons the DNA change touches; so do `predict` and VEP (HBB c.20_22delinsTTT: p.(Glu7_Glu8delinsValTer)), while Mutalyzer 3 deletes the rest of the protein (p.(Glu7_His147delinsVal)), the form the standard rejects ("not … p.(Met3_Ile3418delinsGly)", protein/insertion.md).
6. **The extension's count.** "Ter17 … at position 17 of the added sequence" (protein/extension.md) does not say whether the residue replacing the stop is position 1 of the added sequence, while the frameshift's rule counts its first changed residue as codon 1. Mutalyzer 3 counts the added sequence from the residue after the one replacing the stop (HBB c.444A>T, TAA → TAT: Tyr at 148, the new stop at 169, p.(Ter148TyrextTer21)); `predict` follows it, so the two counts differ by one in their conventions.
7. **Unknown start.** p.(Met1?) for a lost start (protein/substitution.md); Mutalyzer 3 writes p.?.
8. **Transcript flanks.** Positions beyond a transcript are "not allowed" with the transcript as reference, but NC_000023.10(NM_004006.1):c.-128354C>T is valid: c. numbering is extrapolated beyond the transcript without the standard defining it.
9. **Introns on transcripts.** numbering.md forbids intronic positions on a transcript reference, while DNA/repeated.md writes NM_000492.3:c.1210-12_1210-6T[7].
10. **RNA repeats.** RNA/repeated.md calls r.-125_-123cug[4] redundant, while the same page's proposal allows only that form (r.123_191cag[23]); its composite form r.456_465[4]466_489[9]490_499[3] is marked invalid in the sentence that introduces it.
11. **Alleles inside positions.** DNA/alleles.md forbids c.2376[G>C];[G>C]; RNA/alleles.md recommends r.76[a>c];[=].
12. **Translocations.** The t(3;14) example writes the inserted bases before the chromosome 14 piece in der(3), while ISCN puts the piece first; it writes one inverted range ascending (g.29745314_qterinv) and one descending (g.36969141_pterinv); its t(2;11)(q31.1;q22.3) gives chromosome 2's breakpoint as 17450009 in HGVS and 174500098 in ISCN (2q31.1 lies near 174.5 Mb). Here inverted ranges are written 5′ to 3′ with inv, and pieces in junction order.
13. **Lapses of letters and markup.** An RNA linker with a DNA 't' (guaugauuuuuuaataa, RNA/adjoined_transcript.md); r.(100_150)delN[15] with a DNA N; an insertion without a position (g.[…;insG], DNA/complex.md); an insertion between non-adjacent residues (p.[Ser68_Ala74insSerGln], protein/alleles.md); r.124_500delins[AB053210.2:r.1289-365_1289-73inv] marked invalid in the sentence that recommends it (RNA/inversion.md); Pyl missing from the amino-acid table (standards.md).

## 7. References

den Dunnen et al. 2016, Human Mutation 37, 564–569; Hart et al. 2024, Genome Medicine 16, 149; HGVS Nomenclature 21.1 (https://hgvs-nomenclature.org); Holmes et al. 2020, Bioinformatics 36, 1902–1907; Tan, Abecasis and Kang 2015, Bioinformatics 31, 2202–2204; The Variant Call Format Specification v4.3 (samtools/hts-specs); McLaren et al. 2016, Genome Biology 17, 122; Earley 1970, Communications of the ACM 13, 94–102; Vis et al. 2015, Bioinformatics 31, 3751–3757 (the description extractor).
