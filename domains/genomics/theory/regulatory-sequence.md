# Regulatory Sequence and Expression

*Status: draft (2026-11-01). The theory of the `regulation` module ([README](../regulation/README.md)): the σ⁷⁰ promoter as a two-box model with a variable spacer and its EM, the thermodynamic model of occupancy and its regulation factors, and sequence-to-expression models learnt from data, placed in the landscape of the field.*

## 1. The σ⁷⁰ Promoter as a Mixture

**The biology.** Bacterial RNA polymerase holoenzyme with σ⁷⁰ binds two hexamers upstream of the start site +1, the −35 box (consensus TTGACA) and the −10 box (TATAAT), 15 to 19 bases apart, 17 most often and best (Hawley and McClure 1983; Harley and Reynolds 1987). σ⁷⁰'s region 4.2 reads the −35 box, region 2.4 the −10 box, which also melts to open the DNA; an extended −10 (TGn) and an AT-rich UP element, read by the α subunits' C-terminal domains, add strength (Ross et al. 1993).

**Definition 1.1.** For windows x of the region −U … +D about a start site, a *two-box model* is θ = (P₁₀, P₃₅, π_r, π_s, γ, q): position probability matrices of widths w₁₀, w₃₅, a law π_r of the −10 box's offset r from +1, a law π_s of the spacer s, the probability γ that a window holds a promoter, and a background composition q. Under θ a window is background with probability 1 − γ; otherwise z = (r, s) is drawn from π_r ⊗ π_s and the boxes' bases from P₁₀ and P₃₅, the rest from q.

**Proposition 1.2.** The likelihood of a window is Pr(x | θ) = Pr₀(x) Λ(x), Pr₀ the background's and Λ(x) = (1 − γ) + γ Σ_z π_r(r) π_s(s) R(x, z), with R the product of the boxes' ratios P/q at their positions.

*Proof.* Conditioning on z, the bases outside the boxes have the background's probability and those inside P's; dividing by Pr₀(x) leaves the ratios inside the boxes. Summing over z and the background alternative gives Λ. ∎

**Proposition 1.3 (EM).** With the posterior τ(z | x) = γ π π R(x, z)/Λ(x), the update P₁₀[k, b] ∝ Σ_x Σ_z τ(z | x) [x_{a(z)+k} = b] + β q_b (and likewise P₃₅), π ∝ Σ τ + 1, γ = mean of Σ_z τ, does not decrease the log posterior Σ_x log Λ(x) + log prior, the prior being Dirichlet with these pseudocounts.

*Proof.* This is EM for a finite mixture with a conjugate prior (Dempster, Laird and Rubin 1977, Theorem 1; the M-step maximises the expected complete-data log posterior in closed form, each matrix row and each law being a multinomial). ∎

*Remark.* EM finds a local maximum. MEME's sites, which start it, may sit a few bases off the true box (a phase shift; Lawrence et al. 1993 shift their sampler for the same reason); `learn_sigma70` starts EM from MEME's sites moved by δ ∈ {−2, …, 2} for each box and keeps the model of largest posterior.

**Scanning and a threshold.** The score of a −10 position is log₂ R₁₀ + max_s (log₂ R₃₅ + log₂ π_s/max π_s); random sequence of composition q gives the null law of the score, whose upper quantile at a chosen rate of false hits is the threshold (`promoter_threshold`).

## 2. Occupancy and Regulation

**Definition 2.1 (Bintu et al. 2005).** A promoter architecture has binding sites, each for one species X with n_X copies and specific energy Δε (relative to the N_NS non-specific sites), overlapping sites excluding each other, and bound pairs interacting with factors ω. The weight of a state σ (a set of bound sites without overlaps) is Π_{i∈σ} (n_{X(i)}/N_NS) e^{−Δε_i} Π_{i<j∈σ} ω_ij.

**Proposition 2.2.** If polymerase P has one site, overlapped by no other species' site but those that exclude it, then p_bound = 1/(1 + (N_NS/(P F_reg)) e^{Δε_pd}) with F_reg = Z_reg^{P-bound}/Z_reg^{P-free}, the ratio of the regulators' partition functions with polymerase bound and free.

*Proof.* Factor Z by the polymerase site's state: Z = Z_free + p Z_bound, p = (P/N_NS)e^{−Δε_pd}, where Z_bound sums the regulators' states compatible with the bound polymerase, weighted by their interaction factors with it. Then p_bound = p Z_bound/Z = 1/(1 + Z_free/(p Z_bound)). ∎

**Corollary 2.3 (Bintu et al., eqs. 8, 9, 10 and Table 1).** For one activator site (weight a) interacting with polymerase by f: F_reg = (1 + af)/(1 + a). For a repressor that excludes polymerase: F_reg = 1/(1 + r). For two repressors that both exclude polymerase but not each other, interacting by ω: F_reg = 1/(1 + r₁ + r₂ + ω r₁ r₂). A helper site h interacting with the activator by ω acts as an activator of weight a(1 + hω)/(1 + h).

*Proof.* Compute Z_bound and Z_free in each case: (1 + af)/(1 + a); 1/(1 + r) (only the empty regulator state is compatible); 1/(1 + r₁ + r₂ + ωr₁r₂); for the helper, Z_bound = 1 + h + af(1 + hω), Z_free = 1 + h + a(1 + hω), and dividing both by 1 + h gives the activation form with a′ = a(1 + hω)/(1 + h). ∎

*Remark.* Sequence enters through Δε: a site's log-odds score under a weight matrix approximates −Δε/k_BT up to scale and offset (Berg and von Hippel 1987), the bridge between §1 and §2. Kinney et al. (2010) inferred the lac promoter's energy matrices for RNAP and CRP from sort-seq data and fitted exactly such a model, with an interaction energy between the two.

## 3. Learning Expression from Sequence

**Definition 3.1.** Given sequences x_i with measured activities y_i, a model is a map f_θ from sequences to activities chosen to minimise held-out error. A *linear k-mer model* has f(x) = β₀ + Σ_w β_w c_w(x), c_w the count of word w (or gapped word; Ghandi et al. 2014), with β by the elastic net (Friedman, Hastie and Tibshirani 2010). A *convolutional network* has a first layer of filters W_o ∈ ℝ^{4×w} applied to the one-hot sequence, h_o(t) = b_o + Σ_{b,k} W_o[b, k] x[b, t + k]: on a one-hot sequence h_o(t) is a position weight matrix's score of the window at t.

**Proposition 3.2.** A linear k-mer model with k ≥ w represents every additive position-free weight matrix model of width w summed over positions (Σ_t PWM score of window t), and nothing that depends on where a word lies; a convolutional network with pooling of width m sees position at resolution m.

*Proof.* Σ_t Σ_k W[x_{t+k}, k] is a sum over the w-mers occurring in x of a function of the w-mer; with k = w this is Σ_w β_w c_w with β_w the matrix's score of w. Counts forget positions; pooling of width m keeps the window of m in which a filter fired. ∎

**Held-out evaluation.** The DREAM Challenge's test set mixes native promoters, designed sequences and their single-nucleotide and motif variants (Rafi et al. 2025), so random splits put near copies in both training and test and reward memory; `grouped_split` keeps sequences sharing the variable region's bases 17–36 and 77–96 on one side.

## 4. Place in the Field

The landscape, from classical to deep, and where this module stands:

- **Weight matrices and their thermodynamics.** PWMs (Stormo 2000), their energy reading (Berg and von Hippel 1987), occupancy models (Shea and Ackers 1985; Bintu et al. 2005), sort-seq energy matrices (Kinney et al. 2010), and genome-scale σ⁷⁰ models fitted to massively parallel assays (Urtecho et al. 2019; LaFleur, Hossain and Salis 2022, a 346-parameter free-energy model). Here: a two-box σ⁷⁰ model learnt by EM from Cappable-seq start sites (Ettwiller et al. 2016) and the exact occupancy of any architecture; no fit to an *E. coli* reporter assay (Urtecho et al.'s GEO deposit could not be fetched on 2026-11-01).
- **k-mer models.** gkm-SVM (Ghandi et al. 2014) and its successors; here the gapped k-mer features with a penalised linear model, not the kernel SVM.
- **Convolutional networks.** DeepBind (Alipanahi et al. 2015), Basset (Kelley et al. 2016), and the yeast models of de Boer et al. (2020) and Vaishnav et al. (2022) trained on tens of millions of random promoters; Xpresso (Agarwal and Shendure 2020) for native promoters. Here: a small CNN of DeepBind's shape, trained on thousands of sequences in C++ on tyche's `learning`.
- **Long-context models.** Enformer (Avsec et al. 2021; 200 kb, attention), Borzoi (Linder et al. 2025; RNA coverage) and AlphaGenome (Avsec et al. 2025; one megabase, many modalities) predict tracks along whole loci; the DREAM Challenge (Rafi et al. 2025) found that design choices (data handling, augmentation, architecture details) matter as much as size. None of this is attempted here: the module's networks are for learning what a model of regulatory sequence is, checked to recover planted codes, and evaluated honestly on a stated subset.

## 5. References

See the module's [README](../regulation/README.md), §6.
