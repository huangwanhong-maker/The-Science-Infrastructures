# haplotypes

*Status: draft. A module of the evolution domain.* Haplotypes read against a reference panel: Li and Stephens's copying model, which sees a haplotype as an imperfect mosaic of the panel's; its posterior over the haplotype copied, which imputes the alleles a target was not typed at; its pair model of a diploid's unphased genotypes, which phases them; and the haplotype scans for selection, extended haplotype homozygosity (EHH) and the integrated haplotype score (iHS). Serves use (tested imputation, phasing and scans), learning (every algorithm held to brute force and to the process that generates its data) and research. The theory is [theory/haplotypes.md](../theory/haplotypes.md).

## 1. Concepts

- **A haplotype is a mosaic of others.** Going along the chromosome, the target copies one panel haplotype at a time, switching at the rate of recombination and miscopying at the rate of mutation (Li and Stephens 2003). That is a hidden Markov model whose states are the panel's K haplotypes, and its jumps are to a haplotype drawn uniformly: the transition matrix is (1 − r)I + rP, P the uniform projection, so every step of forward, backward and Viterbi costs O(K) rather than O(K²).
- **Untyped sites are missing data.** A site the target was not typed at is emitted with probability 1, so the chain runs over every site of the panel and imputation is the posterior predictive there. Adding or removing such a site leaves the rest of the model unchanged, because the jumps compose as a semigroup.
- **A diploid is two copying processes.** Its genotypes are the sums of two haplotypes, each a mosaic of the panel. Over the K² ordered pairs of haplotypes the transition is a Kronecker product, and the step is again cheap: O(K²), from row and column sums.
- **Phase is a property of a path.** The pair model is symmetric in its two copies, so the marginal posterior that one copy carries an allele at a heterozygous site is exactly ½. Phasing therefore takes the joint maximum over the copying paths and the haplotypes.
- **Selection lengthens haplotypes.** An allele that rose fast carries the chromosome it arose on, and recombination has had little time to break it up. EHH measures that length about a core site by allele, iHH integrates it over the genetic map, and iHS compares the two alleles, standardised against sites of the same frequency (Sabeti et al. 2002; Voight et al. 2006).
- **Truth by simulation.** The panels and targets are the coalescent's samples (`population/coalescent.hpp`), and the sweeps are those of `population/sweep.hpp`, so imputation, phasing and the scans are measured against a known truth.

## 2. Content

| Header | Content |
|---|---|
| `copying.hpp` | `CopyingModel` (ρ, θ or ε), `li_stephens_theta`, `miscopying`; `Chain` (jumps r_j = 1 − exp(−ρ d_j/K), ε) and `copying_chain(panel, model, map)`; `forward`, `backward`, `log_likelihood` (O(K) memory), `posterior` and `viterbi` in O(N K); `path_log_probability`; `forward_dense`, the O(N K²) reference; `sample_mosaic`, a target drawn from the model |
| `phasing.hpp` | The pair model: `pair_log_likelihood` in O(N K²) and `pair_log_likelihood_dense` in O(N K⁴); `for_each_pair_posterior` and `pair_posterior`; `phase` (the joint maximum over the copying paths and haplotypes, with its alleles) and `phased_log_probability`; `SwitchErrors`, `switch_errors` |
| `imputation.hpp` | `impute_haplotype` (the dosages at untyped sites from the copying posterior), `impute_genotypes` (genotype probabilities and dosages of an unphased diploid from the pair posterior), `dosage_r2` |
| `ehh.hpp` | `carriers`, `ehh_at`, `ehh` (`EhhCurve`, both sides, to a cut-off), by refining the partition of the carriers site by site; `ihh` (the trapezoid rule over the genetic map, cut-off 0.05, truncation flagged); `ihs_scan` (`IhsSite`: frequency, iHH of each allele, ln(iHH_A/iHH_D)), in parallel; `IhsStandardisation`, `fit_ihs_standardisation`, `standardise` |

The selective sweeps the scans are checked on are simulated by `population/sweep.hpp`: a Wright–Fisher population whose standing variation is a coalescent sample, carried forward with selection and recombination.

## 3. Dependencies

- **Domains:** genomics (`hmm`: its `Path` and log helpers; its homogeneous HMM checks the copying chain on an evenly spaced map).
- **Within evolution:** `population` (`Haplotypes`, the coalescent and the sweep for the checks).
- **Foundations:** calx: `random` (streams), `scheduler` (the scan's dynamic loop); `verification` for the checks.

The copying chain is written here rather than on `genomics::hmm`, whose `Model` is homogeneous with emissions fixed by the state: the copying chain's jumps change with the map and its emissions with the panel's alleles at each site, and its uniform jumps make a step O(K), where the general recursion is O(K²).

## 4. Checks

Each check program runs in seconds: `node domains/evolution/tests/run-tests.mjs` builds and runs them all.

| Check | What it verifies |
|---|---|
| copying | θ̃(2) = 1, θ̃(3) = 2/3, ε = 1/6 at K = 2; on 36 random panels (K = 2 … 4, N = 4 … 6, missing alleles) the forward, backward, O(K)-memory and dense likelihoods equal log Σ over all K^N paths, the posterior the enumerated marginals and Viterbi's score the best path's, its path attaining it (all within 1e-12); a target no path can emit; the O(N K) forward against the O(N K²) one on coalescent haplotypes; the likelihood, posterior and Viterbi score against genomics::hmm on an evenly spaced map, with the site folded into the symbol; calibration on 400 mosaics drawn from the model, E[γ_j(Z_j)] = E[Σ_k γ_j(k)²]; a panel of 1000 haplotypes in about 120 ms |
| phasing | Switch errors by hand; on 16 random panels the pair likelihood (fast and dense) and posterior against all K^{2N} pairs of paths, and phasing's score against every pair of paths and assignment of heterozygous alleles, its paths and haplotypes attaining it and adding up to the genotypes; O(K²) against O(K⁴) on coalescent sites; the posterior symmetric in the two copies; 100 coalescent diploids phased on panels of 100 with a switch error rate of 0.014 (random phase 0.49; asserted below 0.1); diploids drawn from the model phased at 0.010 (asserted below 0.05) |
| imputation | On the coalescent's samples (panels of 200, a quarter of the sites typed) haplotypes imputed with concordance 0.994 against the panel's major allele's 0.880, at common sites 0.993 against 0.776, squared error a twentieth of the frequency's, mean r² 0.966 at common sites; unphased diploids (panels of 100) at 0.974 against Hardy–Weinberg's mode 0.778; typed sites kept; genotype laws summing to one; the posterior predictive unbiased under the model, E[x_j − d_j] = 0 for haplotypes and diploids drawn from it |
| ehh | A six-haplotype example by hand: EHH of both alleles, iHH 10/3 and 8/3 (7/3 with the cut-off at 0.5), the unstandardised iHS ln 0.8, truncation, the map's scale; the partition walk against EHH by pairs of carriers on coalescent haplotypes with missing alleles (exact); EHH 1 at the core and non-increasing; a scan the same on 1 and 8 threads; iHS of mean 0 and variance 1 in each bin by construction; partial sweeps to 0.7 (2N = 1000, s = 0.05): 30 of 39 selected sites at iHS < −2 (mean −2.7), against 4.6% of neutral sites at \|iHS\| > 2 |

## 5. Performance

Measured on the development machine (i7-11370H, WSL):

- forward–backward and Viterbi over a panel of 1000 haplotypes at 4571 sites: about 120 ms;
- a diploid phased on a panel of 100 haplotypes at about 230 sites (10⁴ pairs of states): about 10 ms;
- 100 iHS scans of 100 haplotypes at about 300 sites each, on 8 threads: under 2 s.

## 6. Open Questions

1. **Panels of thousands.** The pair model is O(N K²), so it suits panels of up to a few hundred. Larger ones are phased by iterating haploid copying against the other samples' current estimates, and imputed from a subset of the panel chosen by long matches to the target. The positional Burrows–Wheeler transform (Durbin 2014) finds those matches; it is being written in morphe's `formal`, and its use here would come after.
2. **Estimating ρ and θ.** Both are parameters here. Li and Stephens estimate ρ, and its hotspots, by maximising their approximate likelihood; an EM of the copying chain (expected switches and miscopyings from the posterior) would estimate both from a panel.
3. **Other scans.** XP-EHH (Sabeti et al. 2007), which compares a site's EHH between populations, and nSL, which measures distance in sites rather than on the map, would reuse the same walk.
4. **Gaps.** A long gap between sites is integrated over like any interval. Implementations of iHS often penalise long gaps or skip them; a rule for gaps would come with real data.

## 7. References

- Durbin, R. (2014), Efficient haplotype matching and storage using the positional Burrows–Wheeler transform (PBWT), *Bioinformatics* 30, 1266–1272.
- Li, N. and Stephens, M. (2003), Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data, *Genetics* 165, 2213–2233.
- Marchini, J., Howie, B., Myers, S., McVean, G. and Donnelly, P. (2007), A new multipoint method for genome-wide association studies by imputation of genotypes, *Nature Genetics* 39, 906–913.
- Maynard Smith, J. and Haigh, J. (1974), The hitch-hiking effect of a favourable gene, *Genetical Research* 23, 23–35.
- Rabiner, L. R. (1989), A tutorial on hidden Markov models and selected applications in speech recognition, *Proceedings of the IEEE* 77, 257–286.
- Sabeti, P. C. et al. (2002), Detecting recent positive selection in the human genome from haplotype structure, *Nature* 419, 832–837.
- Voight, B. F., Kudaravalli, S., Wen, X. and Pritchard, J. K. (2006), A map of recent positive selection in the human genome, *PLoS Biology* 4, e72.
