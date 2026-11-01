# Codon Usage and Expression

*Status: draft (2026-11-01). The theory of the `codon` module ([README](../codon/README.md)): synonymous families as the fibres of the genetic code, the indices of codon bias and what they estimate, the tRNA adaptation index and the test of translational selection, and what measured abundance says of them.*

## 1. The Genetic Code as a Map and Its Fibres

**Definition 1.1.** A genetic code is a map c : 𝒞 → 𝒜 ∪ {∗} from the 64 codons 𝒞 = {A, C, G, T}³ to the amino acids and the stop. Its *synonymous families* are the fibres c⁻¹(a), a ∈ 𝒜; the *sense codons* are c⁻¹(𝒜). The *degeneracy* of a is |c⁻¹(a)|.

NCBI's 27 tables (`sequence/genetic_code.hpp`) are such maps; they differ in a few fibres. In the standard code the degeneracies are 1 (Met, Trp), 2 (nine amino acids), 3 (Ile), 4 (five) and 6 (Leu, Ser, Arg), Wright's classes; in the vertebrate mitochondrial code (table 2) AGA and AGG are stops, AUA joins AUG as methionine and UGA joins UGG as tryptophan, so every amino acid has at least two codons.

**Definition 1.2.** A *box* is a set of four codons sharing their first two bases. Splitting each fibre by boxes gives the *codon families* of Sun, Yang and Xia (2013): Leu = TTR ⊔ CTN, Ser = TCN ⊔ AGY, Arg = CGN ⊔ AGR in the standard code (`synonymous_families(code, true)`: 23 families).

*Rationale.* Within a box the third position alone varies, and one tRNA with a wobbling first anticodon base can read several of its codons; across boxes the codons are read by unrelated tRNAs. A six-fold amino acid is two decoding problems.

## 2. Indices of Bias

Let n = (n_x)_{x∈𝒞} be a gene's codon counts and, for a family F, n_F = Σ_{x∈F} n_x and p_x = n_x/n_F.

**Definition 2.1 (RSCU and CAI; Sharp and Li 1987, eqs. 1–4).** Over a reference set with counts X, RSCU_x = X_x/(X_F/|F|) and w_x = X_x / max_{y∈F} X_y, an absent codon taking X = ½. The codon adaptation index of a gene is the geometric mean CAI = exp((1/L) Σ_x n_x ln w_x) over the L codons of families with more than one codon.

**Proposition 2.2.** CAI ∈ (0, 1], and CAI = 1 exactly when the gene uses only the reference set's most frequent codon of each amino acid.

*Proof.* w_x ∈ (0, 1] with equality iff x maximises X on its family; a geometric mean of numbers in (0, 1] lies in (0, 1] and equals 1 iff all are 1. ∎

**Definition 2.3 (Wright's N_c; Wright 1990).** For an amino acid observed n_F ≥ 2 times, F̂ = (n_F Σ p_x² − 1)/(n_F − 1), the unbiased estimate of the homozygosity Σ π_x² of its codon distribution π. With F̄_k the mean F̂ over the amino acids of degeneracy k, N_c = Σ_k K_k/F̄_k, K_k the number of amino acids of degeneracy k.

**Proposition 2.4.** If every amino acid of a gene uses one codon, N_c = |𝒜| (20 in the standard code). If π is uniform on every family, E F̂ = 1/|F|, and N_c with F̄_k replaced by its expectation is the number of sense codons (61).

*Proof.* One codon: Σ p² = 1, F̂ = 1, N_c = Σ_k K_k. Uniform: E[n Σ p²] = 1 + (n − 1)/|F| for multinomial counts, so E F̂ = 1/|F| and Σ_k K_k k = 61. ∎

*Caveat.* F̂ is unbiased, 1/F̄ is not, and F̂ ≤ 0 is possible (an amino acid seen n ≤ |F| times in n different codons). For isoleucine, alone in its class, that makes 1/F̄₃ undefined; Wright capped N_c at 61, and we keep the cap (`trna.hpp`'s remark explains why dos Reis et al.'s uncapped values break the S value).

**Definition 2.5 (Sun, Yang and Xia 2013).** On codon families of m codons, F_CF = Σ_i ((n_i + 1)/(n + m))², and F̄_m = Σ_j n_j F_CF,j / Σ_j n_j; N_c = Σ_m K_m/F̄_m.

**Proposition 2.6.** For uniform usage (n_i = n/m) Sun's N_c equals the number of sense codons for every n; for one codon per family it tends to the number of families as the counts grow.

*Proof.* Uniform: (n/m + 1)/(n + m) = 1/m, so F_CF = 1/m and N_c = Σ_m K_m m. One codon: F_CF = ((n + 1)² + (m − 1))/(n + m)² → 1. ∎ (Both are checked in `check_usage`.)

**Definition 2.7 (MILC; Supek and Vlahoviček 2005).** With expected within-family frequencies g from a reference set, M_a = 2 Σ_{x∈a} n_x ln(n_x/(n_a g_x)), the G statistic of the gene's family a, and MILC = Σ_a M_a/L − C, C = Σ_a (|a| − 1)/L − ½ over the families present. Under the null the G statistic has expectation |a| − 1 (its χ² limit), so C removes the length bias that afflicts scaled χ². We follow the authors' implementation (coRdon), whose factor 2 and sign of ½ differ from the paper's printed equations.

**Definition 2.8 (Karlin and Mrázek 2000, eqs. 1–3).** B(F|G) = Σ_a p_a(F) Σ_{x∈a} |f_x − g_x|, the L¹ distance of within-family frequencies averaged over F's amino acids; E_S(g) = B(g|C)/B(g|S); E(g) = B(g|C)/(½B(g|RP) + ¼B(g|CH) + ¼B(g|TF)); a gene is *predicted highly expressed* when two of E_RP, E_CH, E_TF exceed 1.05 and E ≥ 1.

**Definition 2.9 (RCBS and DCBS; Roymondal et al. 2009; Sabi and Tuller 2014).** With f the gene's codon frequencies and f₁, f₂, f₃ its base frequencies at the three positions, d_x = f_x/(f₁f₂f₃) − 1; RCBS = (Π_i (1 + d_{x_i}))^{1/L} − 1, DCBS = (1/L) Σ_i max(f/(f₁f₂f₃), f₁f₂f₃/f). They need no reference set and control for compositional bias position by position.

## 3. The tRNA Adaptation Index

**Definition 3.1.** An isoacceptor is a pair (amino acid, anticodon N₃₄N₃₅N₃₆) with its gene copy number tGCN. It *reads* codon x = comp(N₃₆)comp(N₃₅)z with penalty s when (N₃₄, z) is one of the pairings of dos Reis et al.'s Table 1 — I:U, I:C, I:A for a genomic A₃₄ deaminated to inosine; G:C, G:U; U:A, U:G; C:G; and lysidine's L:A for the C₃₄ of a bacterial tRNA-Ile — and c(x) is its amino acid (a *cognate* read).

**Definition 3.2 (dos Reis, Savva and Wernisch 2004, eqs. 1–3).** W_x = Σ_j (1 − s_xj) tGCN_j over the isoacceptors reading x; w_x = W_x/max W, and the geometric mean of the non-zero w where W_x = 0; tAI_g = exp((1/L) Σ_x n_x ln w_x) over sense codons other than methionine's.

*Remark.* The cognate rule is ours; dos Reis et al.'s formulas in codon order realise it for the standard code, and their R package sets W(AUA) to 1 − s_L:A whatever the genome's count of tRNA-Ile2. codon-bias 0.5 applies the lysidine pairing to every C₃₄, which overcounts NNA codons; the checks compare with it under its eukaryotic rules.

**Definition 3.3 (the S value; eqs. 5–8).** With x_g a gene's GC3s and N_c,g its effective number, ψ_g = f₁(x_g) − N_c,g, f₁(x) = a + x + b/(x² + (c − x)²) the expected N_c of genes without selection (Table 3's suggested a = −6, b = 34, c = 1.025). S = cor(tAI, ψ) over the genes; under the null of no co-adaptation, permuting w over the codons leaves S's distribution centred near 0, and the permutation p-value is (b + 1)/(B + 1).

**Proposition 3.4.** If the w are permuted uniformly at random over the counted codons, the observed S is exchangeable with the permuted ones under the null hypothesis that the assignment of w to codons is unrelated to ψ, so P(p ≤ α) ≤ α.

*Proof.* The observed assignment is one of the equally likely permutations under the null; the rank of the observed statistic among B + 1 exchangeable values is uniform, and (b + 1)/(B + 1) is a valid p-value (Phipson and Smyth 2010). ∎

**stAI (Sabi and Tuller 2014).** Without expression data, the wobble penalties are chosen to maximise Spearman's ρ(tAI, DCBS) over the genes by hill climbing from 35 starts; on *E. coli* the climb raises ρ from 0.459 to 0.476.

## 4. What the Indices Predict

On *Escherichia coli* K-12 (4294 coding genes of NC_000913.3, 4090 matched to PaxDb 6.0's integrated proteome), Spearman's ρ with protein abundance is: CAI (ribosomal proteins as reference) 0.579, tAI 0.557, Fop 0.518, MELP 0.477, DCBS 0.308, Karlin and Mrázek's E 0.263, Wright's N_c −0.378. The S value is 0.645 (p = 0.001 against 999 permutations). Karlin and Mrázek's criterion predicts 275 of 3763 genes of at least 100 codons highly expressed (their 2000 count, on the annotation of the time: 306 of 3898).

*Caveat.* These are associations of a sequence statistic with an abundance that transcription, translation initiation and decay also set; ρ ≈ 0.6 is what codon usage alone can carry in a fast-growing bacterium, and much less in organisms with weak translational selection (dos Reis et al. find S near 0 in humans).

## 5. Place in the Field

The indices here are the classical ones (Ikemura 1981; Sharp and Li 1987; Wright 1990; Karlin and Mrázek 2000; dos Reis et al. 2004; Supek and Vlahoviček 2005; Sun et al. 2013; Sabi and Tuller 2014). Later work models translation as a process: elongation rates from ribosome profiling (the normalised translational efficiency of Pechmann and Frydman 2013), the 5′ ramp of slow codons (Tuller et al. 2010), mRNA folding at the start codon, and codon-pair bias (Coleman et al. 2008); and deep models predict protein abundance from sequence with the codons as one input among many. Those are planned, not written.

## 6. References

See the module's [README](../codon/README.md), §6.
