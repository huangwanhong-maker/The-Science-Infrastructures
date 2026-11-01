# Association: From Genotype to Trait and to Expression

*Status: draft. Serves learning, use and research.* The models behind the module [`association/`](../association/) and the statistics it takes from tyche (`regression`'s SuSiE and variance components, `testing`'s q-values and permutation minima, `regression`'s elastic net):

- the linear model of a trait on many variants, and why one variant at a time is fitted first;
- the projection off covariates, which makes a scan cost one inner product a site;
- the permutation law of a minimum, and the false discovery rate over genes;
- fine-mapping as Bayesian variable selection, and colocalisation as a joint one;
- expression predicted from genotypes, and its association with a trait from summary statistics;
- heritability as a variance component, and polygenic prediction as its posterior mean.

The module's README lists the checks that hold each statement to numbers; §8 says where the code implements it.

**Conventions.** n samples, m sites; g_ij ∈ {0, 1, 2} the dosage of the derived (alternate) allele, p_j its frequency, z_ij = (g_ij − 2p_j)/√(2p_j(1 − p_j)) the standardised dosage; C the n × k covariate matrix (an intercept among its columns), M = I − C(CᵀC)⁺Cᵀ the projection off its columns, r = rank C; ‖·‖ the Euclidean norm; ỹ = My.

## 1. The Model of a Trait

**Definition 1.1 (the additive model).** y = Cγ + Σ_j β_j g_j + e, e ~ N(0, σ_e² I) independent of the genotypes. The genetic value is G = Σ_j β_j g_j, the (narrow-sense, SNP) heritability h² = Var G/(Var G + σ_e²).

Under the infinitesimal limit of many small effects, β_j standardised ~ N(0, h²/m) independently, G has covariance h² A with A = ZZᵀ/m the genetic relationship matrix (Yang et al. 2011, eq. 3): this is §6's variance-component model. Under sparsity, a few β_j are non-zero: §4's fine-mapping model. The simulator (`cohort.hpp`) draws a sparse direct part on standardised dosages (the α = −1 model) and a mediated part through gene expression (§5), and scales them to the h² asked.

*Remark (why one variant at a time).* With m ≫ n the joint model is not identifiable by least squares; the marginal regression of y on each g_j estimates β_j plus the effects of the variants in linkage disequilibrium with it, Σ_k r_jk β_k (for standardised genotypes). Scans find regions; fine-mapping (§4) resolves them.

## 2. Scans and the Projection off Covariates

**Proposition 2.1 (Frisch–Waugh–Lovell).** In y = Cγ + βg + e fitted by least squares, β̂ = g̃ᵀỹ/g̃ᵀg̃ and the residual sum of squares is ỹᵀỹ − β̂²g̃ᵀg̃, with g̃ = Mg. So t = β̂/se on n − r − 1 degrees of freedom is the full regression's.

*Proof.* The normal equations' block for γ gives Cγ̂ = (I − M)(y − β̂g); substituting into the block for β gives gᵀM(y − β̂g) = 0, so β̂ = gᵀMy/gᵀMg = g̃ᵀỹ/g̃ᵀg̃ as M is a symmetric idempotent. The residual is M(y − β̂g) = ỹ − β̂g̃, whose squared norm is ỹᵀỹ − 2β̂g̃ᵀỹ + β̂²g̃ᵀg̃ = ỹᵀỹ − β̂²g̃ᵀg̃. ∎

So ỹ is formed once, and each site costs O(nr) to project and O(n) to test. The partial correlation r = g̃ᵀỹ/(‖g̃‖‖ỹ‖) gives t = r√(ν/(1 − r²)), ν = n − r − 1, which is Matrix eQTL's identity (Shabalin 2012).

**Proposition 2.2 (the score test of a logistic model).** For logit μ = Cγ + βg, the score at β = 0 and the null estimate γ̂₀ is U = gᵀ(y − μ̂₀), and its variance with γ estimated is V = gᵀWg − gᵀWC(CᵀWC)⁻¹CᵀWg, W = diag μ̂₀(1 − μ̂₀); U²/V → χ²₁ under the null.

*Proof.* The efficient information for β is the Schur complement of the γ block in the Fisher information [[CᵀWC, CᵀWg], [gᵀWC, gᵀWg]] (Cox and Hinkley 1974, §9.3); the score test refers U² to it. ∎

**Proposition 2.3 (genomic control).** If every statistic is a χ²₁ multiplied by a common λ (an inflation by cryptic relatedness or stratification), the median of the observed statistics is λ F⁻¹_{χ²₁}(½), and dividing by λ̂ = median/0.4549 restores the null (Devlin and Roeder 1999).

*Caveat.* A polygenic trait inflates the median without confounding: every site in LD with some causal variant has a non-central statistic. λ_GC is then not a measure of confounding; LD score regression's intercept is (planned).

## 3. The Least of Correlated Tests

Let a gene's window hold p variants and T = max_j r_j² its strongest association. Under the null, permuting the expression (residualised on C) leaves the genotypes' correlations intact and makes every permutation of the samples equally likely, so the permutation law of T is its exact conditional null law, and (b + 1)/(B + 1) is a valid p-value (Phipson and Smyth 2010).

**Proposition 3.1.** If the p variants were independent, the least of their p-values would be Beta(1, p): P(min ≤ x) = 1 − (1 − x)^p.

*Proof.* Each p-value is uniform; the minimum exceeds x when all do. The density p(1 − x)^{p−1} is Beta(1, p)'s. ∎

With linkage the minimum behaves as the least of p_eff < p independent tests, and FastQTL's beta law Beta(a, b), fitted by maximum likelihood to the permutation minima, absorbs both b ≈ p_eff and a small departure a ≠ 1; its distribution function at the observed minimum is the adjusted p. The beta family is an exponential family in (a − 1, b − 1), so its log-likelihood is concave and Newton's method from the moments converges (tyche `fit_beta`).

The genes' adjusted p-values are then a family of tests, whose false discovery rate is controlled by Storey's q-values: q(p) = π̂₀ m p/#{p_i ≤ p}, made monotone (Storey and Tibshirani 2003).

## 4. Fine-Mapping

**Definition 4.1 (the single-effect regression).** y = x_γ β + e with γ ~ Mult(1, π) and β ~ N(0, σ₀²). For each j, with b̂_j = x_jᵀy/x_jᵀx_j and s_j² = σ²/x_jᵀx_j, the Bayes factor against β = 0 is BF_j = N(b̂_j; 0, σ₀² + s_j²)/N(b̂_j; 0, s_j²).

**Proposition 4.2.** The posterior of γ is α_j ∝ π_j BF_j, and given γ = j, β ~ N(σ₁ⱼ² b̂_j/s_j², σ₁ⱼ²) with σ₁ⱼ² = (1/s_j² + 1/σ₀²)⁻¹.

*Proof.* The likelihood of y given γ = j and β depends on β through b̂_j alone (it is sufficient for β in that model), b̂_j | β ~ N(β, s_j²); integrating β against its prior gives the marginal N(b̂_j; 0, σ₀² + s_j²) up to a factor common to all j, and the conjugate normal update gives the conditional posterior. ∎

This is Wakefield's approximate Bayes factor when b̂ and its variance come from any regression (`finemap_abf`), and the exact posterior when one variant is causal. SuSiE sums L such effects, b = Σ_l b_l, and fits q(b) = Π_l q_l(b_l) by coordinate ascent: each q_l is the exact posterior of a single-effect regression on the expected residual y − X Σ_{l′≠l} b̄_{l′} (Wang et al. 2020, Proposition 1), so each step raises the evidence lower bound. A credible set of effect l is the smallest set of variants holding α_l's mass 0.95; its purity, the least |r| between two of its members, tells a localised effect from a diffuse one.

**Proposition 4.3 (summary statistics suffice).** For standardised columns and trait, the marginal z-scores and the LD matrix R determine XᵀX = (n − 1)R and Xᵀy = (n − 1)r with r_j = z_j/√(z_j² + n − 2).

*Proof.* For one standardised column, t² = (n − 2)r²/(1 − r²) for the sample correlation r; solve for r. Then Xᵀy/(n − 1) is the vector of correlations and XᵀX/(n − 1) the correlation matrix. ∎

So SuSiE-RSS on (z, R, n) is SuSiE on the sufficient statistics (Zou et al. 2022, the PVE-adjusted z-scores); with in-sample R it reproduces the individual-level fit, which the checks see.

## 5. Colocalisation

Two traits measured at the same p variants (a disease and a gene's expression). If each has at most one causal variant, the configurations are pairs (j, k) ∈ ({0} ∪ [p])², and the five hypotheses H₀ (0, 0), H₁ (j, 0), H₂ (0, k), H₃ (j, k) with j ≠ k, H₄ (j, j) partition them.

**Proposition 5.1 (Giambartolomei et al. 2014).** With independent studies, per-variant priors p₁, p₂, p₁₂ and Bayes factors BF¹, BF², the posterior odds of H₃ and H₄ against H₀ are p₁p₂(Σ_j BF¹_j Σ_k BF²_k − Σ_j BF¹_j BF²_j) and p₁₂ Σ_j BF¹_j BF²_j.

*Proof.* Independent studies factor the likelihood of a configuration, so its Bayes factor against (0, 0) is BF¹_j BF²_k (with BF_0 = 1). Summing over the configurations of each hypothesis and weighting by their priors gives the odds; the distinct pairs are all pairs less the diagonal. ∎

coloc.susie (Wallace 2021) applies this to each pair of SuSiE effects, one per trait, relaxing the one-variant assumption.

## 6. Expression Predicted from Genotypes

PrediXcan (Gamazon et al. 2015) learns, in a panel where expression is measured, weights w with Ê = Σ_l w_l g_l from the cis window by the elastic net, and tests in a separate study the association of the trait with Ê. If expression mediates the trait, y = aE + …, then y and Ê covary through the cis-genetic part of E, which is what the test sees: a TWAS tests the association of the genetically determined part of expression.

**Proposition 6.1 (S-PrediXcan from summary statistics; Barbeira et al. 2018, eq. 2).** For the regression of y on Ê = Σ_l w_l g_l, with σ_l the sd of g_l and σ_g² = wᵀΓw the variance of Ê (Γ the covariance of the dosages),

z_g ≈ Σ_l w_l (σ_l/σ_g) z_l,

where z_l = β̂_l/se_l are the marginal statistics of the scan.

*Proof (sketch).* The slope of y on Ê is Cov(Ê, y)/Var Ê = Σ_l w_l Cov(g_l, y)/σ_g², and Cov(g_l, y) = β̂_l σ_l²; the slope's standard error is σ_y/(√n σ_g) and that of β̂_l is σ_y/(√n σ_l) when each variant explains little of y, so z_g = Σ_l w_l β̂_l σ_l²/σ_g · √n/σ_y = Σ_l w_l (σ_l/σ_g) z_l. ∎

With Γ = diag σ · R · diag σ this is FUSION's z = w̃ᵀz/√(w̃ᵀRw̃), w̃ = σ ∘ w (Gusev et al. 2016). The approximation is exact up to the share of variance each variant explains, which is why the checks find S-PrediXcan within a few per cent of the individual-level z.

*Remark (honest accuracy).* A model's cross-validated R² is optimistic when λ is chosen on the same folds; nested cross-validation (PredictDB's) chooses λ inside each training fold, so the held-out R² estimates what the model does on new samples, and it is bounded above by the gene's cis-heritability, which the checks hold it to.

## 7. Heritability and Polygenic Prediction

**Proposition 7.1 (the variance-component model).** Under §1's infinitesimal limit, y ~ N(Cγ, σ_g² A + σ_e² I), and h² = σ_g²/(σ_g² + σ_e²). REML maximises the likelihood of the error contrasts, free of γ (Patterson and Thompson 1971); with A = U S Uᵀ the covariance is diagonal in the rotated data, so one decomposition gives the exact restricted likelihood at every h² (tyche `reml_spectral`).

**Proposition 7.2 (Haseman–Elston).** For standardised y and i ≠ j, E[y_i y_j] = h² A_ij, so the least-squares slope of y_i y_j on A_ij estimates h² (a method of moments).

*Proof.* E[y_i y_j] = Cov(y_i, y_j) = σ_g² A_ij + σ_e² δ_ij with σ_g² + σ_e² = 1. ∎

**Proposition 7.3 (the liability scale; Lee et al. 2011, eq. 23).** If a case–control trait is the indicator that a Gaussian liability exceeds its 1 − K quantile, the heritability of the 0/1 trait in a sample with case fraction P is h²_o = h²_l z² P(1 − P)/(K²(1 − K)²), z the normal density at the threshold.

LDpred-inf (Vilhjálmsson et al. 2015) is the posterior mean of the standardised effects under the infinitesimal prior β ~ N(0, (h²/m) I) given the marginal estimates β̃ = Xᵀy/(n − 1) and the LD D: since β̃ | β ~ N(Dβ, D/n) approximately, the posterior mean is (D + (m/(nh²)) I)⁻¹ β̃, which with in-sample D is ridge regression with penalty m/(nh²), the best linear unbiased predictor. Clumping and thresholding instead keeps one variant per LD clump below a p-value threshold, a crude but robust approximation of sparsity.

## 8. Where the Code Is

| Section | Code |
|---|---|
| §1 the simulation | `cohort.hpp`: `simulate_study` (the coalescent's founders, Li–Stephens mosaics, the trait, the expression) |
| §2 scans | `scan.hpp`: `Residualizer`, `linear_scan`, `logistic_scan`, `genomic_control`, `qq_points`; `ld.hpp`: `principal_components`, `ld_correlation`, `ld_scores`, `clump` |
| §3 permutations | `eqtl.hpp`: `cis_eqtl`, `expression_components`; tyche `testing` |
| §4 fine-mapping | `finemap.hpp`: `log_abf`, `finemap_abf`, `finemap_susie`, `finemap_rss`, `against_truth`; tyche `regression/susie.hpp` |
| §5 colocalisation | `coloc.hpp`: `coloc_abf`, `coloc_from_lbf`, `coloc_susie` |
| §6 prediction | `predict.hpp`: `train_expression_models`, `predict_expression`, `associate`, `s_predixcan`; tyche `regression/penalized.hpp` |
| §7 heritability, scores | `heritability.hpp`: `genetic_relationship`, `snp_heritability`, `liability_scale`; `prs.hpp`: `clump_and_threshold`, `ldpred_inf`, `polygenic_score`, `score_r2`; tyche `regression/variance_components.hpp` |
