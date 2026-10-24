# Selection: Codon Sites, dN/dS and the McDonald–Kreitman Test

*Status: draft. Serves learning and use.* How selection on proteins shows in the DNA that codes for them: the geometry of the genetic code on the space of codons, the counting of synonymous and non-synonymous sites and differences (Nei and Gojobori 1986), the ratio ω = dN/dS, and the McDonald–Kreitman test, which compares that ratio between divergence and polymorphism. The code is [`selection/`](../selection/).

## 1. The Genetic Code on the Space of Codons

**Definition 1.1.** Let 𝔅 = {A, C, G, T} and 𝒞 = 𝔅³, the 64 codons. A *genetic code* is a map γ: 𝒞 → 𝒜 ∪ {∗} to the twenty amino acids and the stop signal; its *sense codons* are 𝒞_s = γ⁻¹(𝒜). The NCBI tables (`genomics::sequence::genetic_code`) are 27 such maps.

**Definition 1.2 (the Hamming graph).** H(3, 4) is the graph on 𝒞 whose edges join codons differing at exactly one position; each codon has 9 neighbours, 3 at each position. An edge {c, c′} is *synonymous* when γ(c) = γ(c′), and the synonymous edges are those inside the fibres of γ.

A single-nucleotide mutation is a step along an edge of H(3, 4). Selection on the protein sees only whether the step leaves the fibre: the fibres of the standard code are mostly unions of codons differing at the third position (the wobble), and so most synonymous edges are third-position edges.

## 2. Synonymous Sites

**Definition 2.1 (Nei and Gojobori 1986).** For a codon c and a position i ∈ {1, 2, 3}, let f_i(c) be the fraction of the three changes at position i that are synonymous. The codon has s(c) = Σ_i f_i(c) synonymous sites and n(c) = 3 − s(c) non-synonymous ones. A sequence of codons c₁, …, c_r has S = Σ_j s(c_j) and N = 3r − S.

So s(c) is the number of synonymous neighbours of c in H(3, 4) divided by 3: if mutation were uniform over the nine single-nucleotide changes, s(c)/3 would be the probability that a mutation of c is synonymous. A change to a stop codon counts as non-synonymous, as in the paper (NG86, p. 419: TTA, leucine, has s = 1/3 + 0 + 1/3 = 2/3 and n = 7/3, though two of its second-position neighbours, TAA and TGA, are stops).

*Examples (standard code).* A fourfold-degenerate codon (GCN, alanine) has s = 1, all from its third position. ATG (methionine) and TGG (tryptophan), alone in their fibres, have s = 0. Under the vertebrate mitochondrial code, where TGA is tryptophan too, TGG has s = 1/3. The checks hold the code to each.

*Remark (what "site" means).* A site here is a fraction of a nucleotide position, a measure of opportunity: S is the expected number of positions at which a uniformly random mutation would be synonymous. It is the counting analogue of the rates of a codon model, and the methods that weight changes by transition bias and codon frequencies (Li, Wu and Luo 1985; Yang and Nielsen 2000) refine the same idea.

## 3. Differences and Pathways

Two sense codons c, c′ at Hamming distance k differ by k mutations, which may have occurred in any of the k! orders. Each order is a geodesic of H(3, 4) from c to c′, and along each, every step is synonymous or not.

**Definition 3.1.** Let Γ(c, c′) be the geodesics from c to c′ whose intermediate codons are sense codons. The synonymous and non-synonymous differences of the pair are the means over Γ(c, c′) of the numbers of synonymous and non-synonymous steps,

s_d(c, c′) = (1/|Γ|) Σ_{g ∈ Γ} #syn(g),  n_d(c, c′) = k − s_d(c, c′).

For a pair of aligned sequences, S_d = Σ_j s_d(c_j, c′_j) and N_d = Σ_j n_d(c_j, c′_j).

*Example (NG86).* TTT (Phe) and GTA (Val) differ at positions 1 and 3. The geodesic TTT → GTT → GTA has one non-synonymous step (Phe → Val) and one synonymous (Val → Val); TTT → TTA → GTA has two non-synonymous ones (Phe → Leu → Val). So s_d = 1/2 and n_d = 3/2.

*Remark (stops).* The paper does not say how it treated geodesics through stop codons. They are excluded here, as in PAML (Yang 2007) and MATLAB's `dnds`, since a lineage cannot pass through a nonsense codon; TCA (Ser) and TGG (Trp), whose geodesic through TGA is excluded, then have s_d = n_d = 1 rather than 1/2 and 3/2. When every geodesic passes a stop, PAML's convention is followed.

The weighting is uniform over geodesics, which is a simplification: a geodesic through a radical amino-acid change is less likely to have been taken, and weighting by a codon model's transition probabilities (as the likelihood methods do) removes the arbitrariness.

## 4. Distances and ω

**Definition 4.1.** p_S = S_d/S and p_N = N_d/N are the proportions of synonymous and non-synonymous difference. Treating each class of sites as an independent JC69 process (theory/phylogeny.md, §2) gives

d_S = −(3/4) ln(1 − 4p_S/3),  d_N = −(3/4) ln(1 − 4p_N/3),  ω = d_N/d_S

(NG86, eqs. 1–3). Code: `dnds.hpp`, `nei_gojobori`.

*Interpretation.* Synonymous substitutions accumulate at roughly the mutation rate, so d_S measures time; d_N measures time discounted by the probability that a non-synonymous mutation is fixed. Under neutrality ω = 1, purifying selection gives ω < 1 (typical genes have ω between 0.05 and 0.3), and positive selection on a whole gene gives ω > 1, which is rare averaged over all its codons; most positive selection is found by models that let ω vary among sites or branches (Nielsen and Yang 1998).

*Caveat.* The counting method assumes equal rates of transition and transversion, and uniform codon usage. With a transition bias, the third-position transitions that are synonymous happen more often than the count of sites expects, so S is underestimated, d_S overestimated and ω underestimated (Yang and Nielsen 2000). With saturation (p_S ≥ 3/4) the correction is undefined; the code reports it and gives an infinite d_S rather than a number.

## 5. The McDonald–Kreitman Test

**Proposition 5.1 (the neutral expectation; McDonald and Kreitman 1991).** Suppose synonymous mutations are neutral, and a fraction f of replacement mutations is neutral and the rest strongly deleterious (never polymorphic at appreciable frequency, never fixed). Let μ_s and μ_n be the synonymous and replacement mutation rates per locus. Then the expected numbers of fixed differences between two species and of polymorphic sites within a sample satisfy

E[D_n]/E[D_s] = E[P_n]/E[P_s] = f μ_n / μ_s.

*Proof (sketch).* Under neutrality, the number of substitutions on a lineage is the mutation rate times its length (Kimura's result that the substitution rate equals the mutation rate), and the number of segregating sites in a sample is the mutation rate times the total length of the sample's genealogy. Both are a mutation rate times a length that is the same for the two classes of mutation, since they share the genealogy and the divergence time: the classes differ only in their rates, f μ_n and μ_s. So the ratio of the classes is the same in both rows. ∎

So under neutrality the 2 × 2 table

|  | replacement | synonymous |
|---|---|---|
| fixed | D_n | D_s |
| polymorphic | P_n | P_s |

has, in expectation, independent rows and columns, and a test of independence is a test of neutrality that needs no model of demography: whatever the population's history, it changes the lengths, which cancel.

*Departures.* Adaptive substitution adds replacement differences to the fixed row only, raising D_n/D_s: the neutrality index NI = (P_n/P_s)/(D_n/D_s) falls below 1, and α = 1 − NI estimates the proportion of replacement substitutions driven by positive selection (Smith and Eyre-Walker 2002). Slightly deleterious replacement mutations segregate but rarely fix, raising P_n/P_s and biasing α downwards, the more so when rare polymorphisms are counted (Eyre-Walker 2002, Messer and Petrov 2013). The direction of selection DoS = D_n/(D_n + D_s) − P_n/(P_n + P_s) (Stoletzki and Eyre-Walker 2011) is the same comparison as a difference of proportions, defined when a cell is zero.

**The tests.** Conditionally on the margins, under independence the top-left cell is hypergeometric: P(D_n = a) = C(R₁, a) C(R₂, C₁ − a) / C(n, C₁), with R the row sums, C the column sums and n the total. *Fisher's exact test* sums the probabilities of the tables with these margins that are no more probable than the observed one (two-sided); conditioning on the margins removes the nuisance parameters, the row and column proportions, exactly, since the margins are sufficient for them under independence (Fisher 1935; Lehmann and Romano 2005, Testing Statistical Hypotheses, ch. 4). The *G-test* compares G = 2 Σ O ln(O/E) with χ²₁, and Williams' correction divides G by q = 1 + (n Σ 1/R_i − 1)(n Σ 1/C_j − 1)/(6n) to bring its distribution nearer χ²₁ in small samples (Williams 1976). Code: `mcdonald_kreitman.hpp`, with the hypergeometric probabilities from calx's log Γ.

*The Adh data.* McDonald and Kreitman found, at the alcohol dehydrogenase locus of *Drosophila melanogaster*, *D. simulans* and *D. yakuba*, D_n = 7, D_s = 17, P_n = 2, P_s = 42: replacement changes were 29 % of the fixed differences but 5 % of the polymorphisms. NI = (2/42)/(7/17) = 0.116, α = 0.88; the G-test with Williams' correction gives G = 7.43 and P = 0.0064, which they reported as P = 0.006, and Fisher's exact test P = 0.0073. The checks reproduce these numbers.

*Remark (classification of sites in the code).* `mk_table` counts nucleotide sites: a site varying within either species is a polymorphism, one fixed for different bases in the two species a fixed difference. A site is synonymous when, in every codon observed at that position of the alignment, changing it among its observed bases keeps the amino acid. For a codon column with one variable site this is McDonald and Kreitman's rule; for columns with several (reported as `complex_codons`), the rule is conservative, calling a change a replacement if any context makes it one, where a pathway method would split it.

## 6. References

- Nei, M. and Gojobori, T. (1986), *Mol. Biol. Evol.* 3, 418–426.
- Nei, M. and Kumar, S. (2000), *Molecular Evolution and Phylogenetics*, Oxford University Press.
- Li, W.-H., Wu, C.-I. and Luo, C.-C. (1985), *Mol. Biol. Evol.* 2, 150–174.
- Yang, Z. and Nielsen, R. (2000), *Mol. Biol. Evol.* 17, 32–43; Nielsen, R. and Yang, Z. (1998), *Genetics* 148, 929–936; Yang, Z. (2007), PAML 4, *Mol. Biol. Evol.* 24, 1586–1591.
- McDonald, J. H. and Kreitman, M. (1991), *Nature* 351, 652–654.
- Smith, N. G. C. and Eyre-Walker, A. (2002), *Nature* 415, 1022–1024; Eyre-Walker, A. (2002), *Genetics* 162, 2017–2024; Messer, P. W. and Petrov, D. A. (2013), *PNAS* 110, 8615–8620; Stoletzki, N. and Eyre-Walker, A. (2011), *Mol. Biol. Evol.* 28, 63–69; Rand, D. M. and Kann, L. M. (1996), *Mol. Biol. Evol.* 13, 735–748.
- Fisher, R. A. (1935), *The Design of Experiments*, Oliver and Boyd; Williams, D. A. (1976), *Biometrika* 63, 33–37; Sokal, R. R. and Rohlf, F. J. (1995), *Biometry*, 3rd ed., Freeman; Lehmann, E. L. and Romano, J. P. (2005), *Testing Statistical Hypotheses*, 3rd ed., Springer.
