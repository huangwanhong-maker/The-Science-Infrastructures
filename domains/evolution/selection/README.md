# selection: Signatures of Selection in Coding Sequences

*Status: draft. Module of the evolution domain.* What selection leaves in protein-coding DNA: the ratio of non-synonymous to synonymous divergence between two sequences (dN/dS, by Nei and Gojobori's counting method), and the comparison of that ratio between differences fixed between species and polymorphisms within them (the McDonald–Kreitman test). The theory is [theory/selection.md](../theory/selection.md). Serves use and learning.

## 1. Concepts

- **Synonymous change is the neutral yardstick.** A change that keeps the amino acid is invisible to selection on the protein, so its rate is the mutation rate; a non-synonymous rate below it (ω < 1) is purifying selection, above it (ω > 1) positive selection.
- **Sites are fractions.** A codon position is counted as synonymous in the proportion of its three possible changes that are synonymous, so a sequence has S synonymous and N non-synonymous sites, S + N = 3 per codon.
- **Fixed and polymorphic differences are the same mutations seen at two times.** Under neutrality, the replacement-to-synonymous ratio is the same for both; adaptive fixation raises it among fixed differences, weakly deleterious variation among polymorphisms.

## 2. Content

| Header | Content |
|---|---|
| `dnds.hpp` | `CodonTable(code)`: for each of the 64 codons of an NCBI genetic code, its amino acid and synonymous sites (NG86's s, n = 3 − s); for each pair, synonymous and non-synonymous differences averaged over the pathways that avoid stop codons. `nei_gojobori(x, y[, code or table])` → `NeiGojobori {S, N, Sd, Nd, pS, pN, dS, dN, omega, saturated_S, saturated_N, codons, skipped}`, with the Jukes–Cantor correction |
| `mcdonald_kreitman.hpp` | `mk_table(species1, species2[, code])` → `MKTable {Dn, Ds, Pn, Ps, codons, complex_codons, skipped}` from aligned coding sequences; `mcdonald_kreitman(table)` → `MKTest` (neutrality index, α = 1 − NI, direction of selection, Fisher's exact two-sided P, G and Williams-corrected G with their χ² P); `fisher_exact(a, b, c, d)` |

Codons with an ambiguous base, a gap or a stop codon are left out pair by pair (dN/dS) or column by column (MK). Sequences are read in frame from their first base.

## 3. Dependencies

- **Domains:** genomics (`sequence`: codes, the NCBI genetic codes, codon indices).
- **Foundations:** calx — `special_functions` (log Γ for the hypergeometric probabilities, the χ² tail), `verification` (checks).

## 4. Checks

| Check | What it verifies |
|---|---|
| `dnds` | TTA's 2/3 and 7/3 sites (NG86, p. 419); a fourfold-degenerate codon's 1; ATG and TGG none, TGG 1/3 under the vertebrate mitochondrial code; s + n = 3 and symmetric differences for every codon; GTT/GTA (1, 0) and TTT/GTA (0.5, 1.5), NG86's examples; TCA/TGG with the stop pathway left out (1, 1); a twelve-codon pair worked by hand, S = 67/6, N = 149/6, S_d = 3.5, N_d = 2.5, and its d_S, d_N, ω in closed form; synonymous change alone gives ω = 0; gaps and stops left out; saturation reported |
| `mcdonald_kreitman` | The Adh table (7, 17, 2, 42) of McDonald and Kreitman 1991: Williams-corrected G = 7.4327 with P = 0.0064, their published P = 0.006; uncorrected G = 7.9071; Fisher's two-sided P = 0.007327 (as R's `fisher.test`); NI = 0.1156, α = 0.8844; Fisher's tea-tasting table (P = 34/70); a table read from sequences built by hand, with a column of two variable sites and columns left out |

## 5. References

- Nei, M. and Gojobori, T. (1986), Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions, *Mol. Biol. Evol.* 3, 418–426.
- Nei, M. and Kumar, S. (2000), *Molecular Evolution and Phylogenetics*, Oxford, ch. 4.
- Yang, Z. (2007), PAML 4, *Mol. Biol. Evol.* 24, 1586–1591; Yang, Z. and Nielsen, R. (2000), *Mol. Biol. Evol.* 17, 32–43.
- McDonald, J. H. and Kreitman, M. (1991), Adaptive protein evolution at the *Adh* locus in *Drosophila*, *Nature* 351, 652–654.
- Rand, D. M. and Kann, L. M. (1996), *Mol. Biol. Evol.* 13, 735–748; Smith, N. G. C. and Eyre-Walker, A. (2002), *Nature* 415, 1022–1024; Stoletzki, N. and Eyre-Walker, A. (2011), *Mol. Biol. Evol.* 28, 63–69.
- Sokal, R. R. and Rohlf, F. J. (1995), *Biometry*, 3rd ed., Freeman; Williams, D. A. (1976), *Biometrika* 63, 33–37; Fisher, R. A. (1935), *The Design of Experiments*, Oliver and Boyd.
- McDonald, J. H. (2014), *Handbook of Biological Statistics*, 3rd ed., Sparky House (biostathandbook.com), "Fisher's exact test".

## 6. Open Questions

1. **A published worked example of NG86's counts.** The paper gives d_S and d_N of globins (its Table 2) but not their sequences' S, N, S_d and N_d; the checks use its codon examples and a pair worked by hand. Nei and Kumar (2000, ch. 4) has a worked example that could be added when a copy is at hand.
2. **Likelihood methods.** Codon models (Goldman and Yang 1994; Muse and Gaut 1994) estimate ω by maximum likelihood, and site and branch-site models test for positive selection; they need the phylogeny module's likelihood over 61 states.
3. **A hypergeometric distribution in calx.** Fisher's exact test forms the hypergeometric probabilities from calx's `log_choose`; a `hypergeometric` probability and tail in `calx::special_functions` would serve other tests (enrichment in the genomics domain's annotation).
4. **Polymorphism in more than two species**, and the asymptotic MK test (Messer and Petrov 2013), which uses the frequencies of polymorphisms rather than their counts.
