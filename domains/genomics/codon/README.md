# codon: Codon Usage and Expression

*Status: draft (2026-11-01).* The usage of synonymous codons by a genome's genes and what it says of their expression: synonymous families under any of NCBI's genetic codes, third positions, the effective number of codons (Wright 1990; Sun, Yang and Xia 2013), RCBS and DCBS, the codon adaptation index (Sharp and Li 1987), Fop, MILC, Karlin and Mrázek's predicted highly expressed genes, the tRNA adaptation index from the genome's tRNA genes with the S value of translational selection (dos Reis, Savva and Wernisch 2004) and stAI's penalties learnt from the genome alone (Sabi and Tuller 2014); and the agreement of every index with measured protein abundance (PaxDb). Serves use (the Genome Workbench paints the indices on a genome's map), learning and research (the theory note [codon-usage.md](../theory/codon-usage.md)).

## 1. Concepts

- **A gene's synonymous choices are a sample from a distribution its expression has shaped.** Highly expressed genes of fast-growing microbes use the codons their abundant tRNAs read (Ikemura 1981); the indices measure how far a gene's choices lean the way a reference set of highly expressed genes leans (CAI, Fop, MILC, Karlin and Mrázek's E), the way the tRNA pool leans (tAI), or away from uniform (N_c, DCBS).
- **The genetic code is a parameter.** Families, degeneracy classes, the stops and the single-codon amino acids follow the code a gene is read with: under the vertebrate mitochondrial code (2) methionine and tryptophan have two codons and CAI counts them.
- **What a codon means decides whether it is counted.** Genes read with a frameshift, recoded (selenocysteine), partial or pseudogenes are left out; a stop completed by polyadenylation keeps its gene.
- **Prediction is tested on measurement.** On *E. coli* K-12, CAI ranks PaxDb's protein abundances with Spearman ρ = 0.58 and tAI with 0.56.

## 2. Content

| Header | Content |
|---|---|
| `genes.hpp` | `CodingGene{locus_tag, gene, product, protein_id, start, stop, reverse, wraps, table, cds, feature}`, `coding_genes(record, GeneOptions{default_table, min_codons})` → `CodingGenes{accession, length, circular, genes, skipped, problems}`; `codon_counts`, `pooled_counts`; reference sets: `genes_by_product(genes, regex)`, `RIBOSOMAL_PROTEINS`, `genes_by_name`, `ecoli_chaperones()`, `ecoli_translation_factors()` (Karlin and Mrázek's Table 2); `TrnaGene`, `trna_genes(record)` (from /anticodon or RefSeq's notes, initiators marked), `amino_acid_letter` |
| `usage.hpp` | `Families`, `synonymous_families(code, split_by_box)`; `third_positions` (GC3, GC3s); `enc_wright(counts, code, cap)` (reusing `sequence::effective_number_of_codons`), `enc_sun(counts, code, split, pseudocount)`; `rcbs`, `dcbs`; `rscu` |
| `adaptation.hpp` | `relative_adaptiveness(reference, code, absent = 0.5)` → `Adaptiveness{w, counted}`, `cai`, `adaptiveness_profile`; `optimal_by_adaptiveness`, `optimal_by_contrast(high, low, code, α)` (χ²), `fop`; `milc(counts, reference, code, pseudocount)`; `codon_bias` (B(F\|G)), `ExpressionStandards`, `predicted_expression` → `PredictedExpression{b_*, e_rp, e_ch, e_tf, e, phx}` |
| `trna.hpp` | `Isoacceptor`, `trna_pool(trnas)`; `Penalties` (dos Reis et al.'s Table 2), `reads(isoacceptor, code, s)` (wobble rules, cognate, lysidine), `trna_adaptiveness(pool, code, s)` → `TrnaAdaptiveness{W, w, counted, w_mean, unread}`, `tai`; `wright_expected_nc`, `NcCurve`, `fit_nc_curve` (eqs. 8–9, by simulation or from points), `enc_wright_defined`; `s_value(genes, w, code, curve, permutations, seed, pool)` → `SValue{s, tai, psi, nc, gc3s, p_value}`; `optimise_penalties(genes, target, pool, code, StaiOptions, threads)` → `StaiResult` (stAI's hill climbing) |
| `abundance.hpp` | `read_paxdb(text)` → `AbundanceTable`; `agreement(genes, index, table)` → `Agreement{matched, spearman, pearson_log, …}` |

## 3. Dependencies

`genomics::sequence` (codes, genetic codes, GenBank, N_c), `genomics::annotation` (features' extents across the origin); `tyche::descriptive` (Pearson, Spearman's ρ); `calx::optimisation` (Nelder–Mead), `calx::random`, `calx::scheduler`, `calx::special_functions` (χ²).

## 4. Checks

`node domains/genomics/tests/run-tests.mjs`, or one program at a time (`check_genomics_codon_<topic>`). Fixtures: `tests/fixtures/codon-ecoli-genes.tsv`, `codon-weights.tsv`, made by `make_codon_fixtures.py` with Biopython 1.88 and codon-bias 0.5.0.

| Check | What it verifies |
|---|---|
| `usage` (23) | The families of the standard code (61 codons, Wright's classes 2 × 1, 9 × 2, 1 × 3, 5 × 4, 3 × 6), split by boxes into 23 (Sun et al.), and of the mitochondrial code (60; 22 split); GC3 and GC3s by hand; Wright's N_c of one codon per amino acid (20), by hand on two amino acids (55.4); Sun's N_c exactly the number of sense codons for uniform usage at any count (61; 60 under code 2) and tending to the number of families (23); RCBS and DCBS by hand; *E. coli* K-12's 4294 counted genes, the same as Biopython's parse; Wright's N_c of 3110 genes against codon-bias (6 × 10⁻¹⁴); Sun's N_c of the pooled genome against codon-bias (2 × 10⁻⁷ relative; codon-bias weights by n + m, we by n); the human mitochondrion's 13 genes under code 2 |
| `adaptation` (18) | Sharp and Li's w = X/X_max and the 0.5 of absent codons, CAI by hand, Met and Trp counted under code 2; Fop by hand; the χ² of the optimal-codon contrast by hand; MILC by hand and against its own usage (−C); B(F\|G) by hand; *E. coli*'s 56 ribosomal proteins as the reference; w of 61 codons and the CAI of 4294 genes against Biopython (5 × 10⁻¹⁶, below 10⁻¹²); w of the human mitochondrion's 60 codons against Biopython; Karlin and Mrázek on *E. coli*: 275 of 3763 genes ≥ 100 codons PHX (7.3%; the paper 306 of 3898, 8%), largest E 2.98 (2.66), all 36 long ribosomal proteins PHX |
| `trna` (24) | The wobble rules by hand (I:U, I:C, I:A, C:G, G:U, lysidine; tRNA-Met kept off AUA), w by eq. 2 with the geometric mean, tAI without methionine; Wright's f and dos Reis et al.'s f₁; *E. coli*'s 86 tRNA genes (4 initiators, 1 Sec), 81 elongators in 40 isoacceptors; tAI's W of 59 codons against codon-bias's (proportional; AUA by hand: codon-bias reads every C₃₄ as lysidine); eq. 9's fit recovering a planted curve; *E. coli*'s genes simulated without selection: N_c 62.6 at GC3s ½ (uniform usage); the S value 0.645, beyond all 999 permutations, the same on one and four threads; stAI raising ρ(tAI, DCBS) from 0.459 to 0.476 over 35 starts, the same on one and four threads |
| `abundance` (8) | PaxDb's format; 4090 of 4294 genes matched to PaxDb 6.0; Spearman ρ with protein abundance: CAI 0.579 (above 0.4), tAI 0.557 (above 0.25), Fop 0.518, MELP 0.477, DCBS 0.308, E 0.263, N_c −0.378; the ribosomal proteins all over ten times the median abundance |

## 5. Performance

Measured in WSL on the development machine (i7-11370H), Release: *E. coli*'s 4294 genes counted and indexed in well under a second; the S value with 999 permutations 0.14 s on four threads; stAI's hill climbing (35 starts, 9198 evaluations) 2.2 s on four threads.

## 6. References

- Ikemura, T. (1981), Correlation between the abundance of *Escherichia coli* transfer RNAs and the occurrence of the respective codons in its protein genes, *J. Mol. Biol.* 146, 1–21; (1985) *Mol. Biol. Evol.* 2, 13–34.
- Sharp, P. M., Tuohy, T. M. F. and Mosurski, K. R. (1986), Codon usage in yeast, *Nucleic Acids Res.* 14, 5125–5143.
- Sharp, P. M. and Li, W.-H. (1987), The codon adaptation index, *Nucleic Acids Res.* 15, 1281–1295.
- Wright, F. (1990), The 'effective number of codons' used in a gene, *Gene* 87, 23–29.
- Karlin, S. and Mrázek, J. (2000), Predicted highly expressed genes of diverse prokaryotic genomes, *J. Bacteriol.* 182, 5238–5250.
- dos Reis, M., Savva, R. and Wernisch, L. (2004), Solving the riddle of codon usage preferences: a test for translational selection, *Nucleic Acids Res.* 32, 5036–5044; the R package tAI 0.21.
- Supek, F. and Vlahoviček, K. (2005), Comparison of codon usage measures and their applicability in prediction of microbial gene expressivity, *BMC Bioinformatics* 6, 182; Elek, A., Kuzman, M. and Vlahoviček, K. (2019), coRdon, Bioconductor.
- Roymondal, U., Das, S. and Sahoo, S. (2009), Predicting gene expression level from relative codon usage bias, *DNA Res.* 16, 13–30.
- Sun, X., Yang, Q. and Xia, X. (2013), An improved implementation of effective number of codons (N_c), *Mol. Biol. Evol.* 30, 191–196.
- Sabi, R. and Tuller, T. (2014), Modelling the efficiency of codon–tRNA interactions based on codon usage bias, *DNA Res.* 21, 511–526.
- Huang, Q., Szklarczyk, D., Oehninger, J. and von Mering, C. (2026), PaxDb v6.0, *Nucleic Acids Res.* 54, D427–D439.

## 7. Open Questions

1. **Mitochondrial decoding.** tAI's rules do not describe superwobbling U₃₄ (one tRNA for a four-fold box); a mitochondrial table of readings would make tAI meaningful there.
2. **dos Reis et al.'s Table 3.** Our simulation of genes without selection gives N_c ≈ 62.6 at GC3s ½, where their refit gives 59; their simulation and CodonW's handling of F̄ ≤ 0 are not described well enough to find why. The S value uses their suggested constants.
3. **CodonW.** Not built here; its outputs (Nc, Fop with Ikemura's optimal codons, CBI) would be a further reference.
4. **The codon context.** Codon-pair bias, ramp regions near the start, and the folding energy of the 5′ end, which also set translation, are not modelled.
