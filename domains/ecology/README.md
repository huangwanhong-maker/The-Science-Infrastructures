# ecology

*Status: planned. Field, height 5, on geography, climate and evolution.* Organisms among each other and in their environment: populations and their growth, the interactions that join species into communities, behaviour, biodiversity and its measures, the distributions of species and the biogeography that explains them, the flows of energy and matter through ecosystems, and the fisheries that harvest them. Behavioural ecology, biodiversity, biogeography and fisheries are its modules. Placed in the order with the other domains of life at the user's request (2026-10-31, [ADR 0061](../../development/decisions/0061-the-domains-of-life.md)); its work is tier 3.7 of the [README's order](../../README.md#46-the-order-of-implementation).

## 1. Purpose

- **A population is a dynamical system, often a structured one.** For one species: logistic growth, delays and harvest. For a population structured by age or stage: a non-negative projection matrix A. When A is primitive, its dominant eigenvalue λ is the asymptotic rate of growth, its right eigenvector the stable structure and its left eigenvector the reproductive values (Perron–Frobenius; Caswell 2001).
- **A community is a system of interacting populations.** Its models are Lotka–Volterra and consumer–resource equations; the Lotka–Volterra equations for n species are equivalent to the replicator equation for n + 1 strategies (Hofbauer and Sigmund 1998). A large random community whose self-regulation is scaled to one is stable, almost surely, only when σ√(nC) < 1, for n species, connectance C and interaction strength σ (May 1972). Food webs follow simple rules of niche (Williams and Martinez 2000).
- **Diversity is a family of means.** Hill's numbers ^qD = (Σᵢ pᵢ^q)^{1/(1−q)} are the effective numbers of species of order q: richness at q = 0, the exponential of Shannon's entropy as q → 1, the inverse of Simpson's index at q = 2 (Hill 1973). Rarefaction and estimators correct for what a sample misses (Chao 1984).
- **Distributions follow climate and history.** Where a species lives is a point process over the land whose intensity depends on climate and terrain; MaxEnt is such a model, a Poisson process fitted by penalised likelihood (Renner and Warton 2013). Island biogeography balances immigration against extinction (MacArthur and Wilson 1967); neutral theory asks what follows when species are ecologically equivalent (Hubbell 2001).
- **An ecosystem is a network of stocks and flows.** Energy enters by photosynthesis and passes up trophic levels, with losses at each (Lindeman 1942); carbon and nitrogen cycle through compartments.

## 2. Scope

| Module | Scope | Status |
|---|---|---|
| `populations` | Single populations: logistic, delayed and harvested growth; matrix population models (Leslie, Lefkovitch), their sensitivities and elasticities; demographic stochasticity and extinction; spread by Fisher–KPP | planned |
| `communities` | Interacting populations: Lotka–Volterra and consumer–resource models, the replicator equation, food webs and the niche model, May's stability of random communities, assembly | planned |
| `behaviour` | Behavioural ecology: optimal foraging and the marginal value theorem (Charnov 1976), movement as random walks read from tracks, evolutionary games of behaviour | planned |
| `biodiversity` | Hill numbers, rarefaction and extrapolation, Chao's estimators, beta diversity and its partitions, species–area relations, neutral theory and Ewens's sampling formula | planned |
| `distributions` | Occurrences (GBIF, OBIS) against climate and terrain: species distribution models as point processes, ranges, biomes on climate's classes, island biogeography | planned |
| `ecosystems` | Stocks and flows: trophic levels, primary production from remote sensing (geography's NDVI), nutrient cycles as compartment models | planned |
| `fisheries` | Stock assessment: surplus production (Schaefer), stock and recruitment (Beverton–Holt, Ricker), yield per recruit, maximum sustainable yield | planned |

## 3. Dependencies

| Source | What ecology uses from it |
|---|---|
| geography | Terrain, land cover, remote sensing (NDVI), zonal statistics, regions |
| climate | Climate classes and normals (Köppen–Geiger), and through it meteorology and oceanography |
| evolution | Phylogeny (phylogenetic diversity), population genetics |
| physis | `dynamics`, `chemistry` (reaction–diffusion, for spread), `statistical_mechanics` (birth–death processes, master equations) |
| tyche | Regression and generalised linear models, spatial statistics, point processes (`processes`, planned) |
| calx | Linear algebra (dominant eigenvalues), integrators, random numbers |

Physiology will be declared when `ecosystems` needs photosynthesis, or the metabolic theory of ecology (Brown et al. 2004).

## 4. Used By

No work uses it yet. The Earth Atlas's vegetation (Sentinel-2's NDVI) and Köppen–Geiger classes are what `distributions` and `ecosystems` will read first.

## 5. Open Questions

- Whether ecology should stand on climate, which makes it the highest domain in the order (height 5, through air–sea), or on meteorology and oceanography directly until climate's bridge is complete. Climate's classes are needed early, so climate is declared.
- Microbial communities: microbiology's `communities`, or here once ecology declares microbiology.
- Eco-evolutionary dynamics, where ecological and evolutionary time meet: a module here, since ecology already stands on evolution.
- GBIF's occurrences carry CC0, CC BY or CC BY-NC record by record; the data library must keep each record's licence.

## 6. References

- Begon, M. and Townsend, C. R. (2021), *Ecology: From Individuals to Ecosystems*, 5th ed., Wiley.
- Kot, M. (2001), *Elements of Mathematical Ecology*, Cambridge University Press.
- May, R. M. (1973), *Stability and Complexity in Model Ecosystems*, Princeton University Press.
- Caswell, H. (2001), *Matrix Population Models*, 2nd ed., Sinauer.
- Hofbauer, J. and Sigmund, K. (1998), *Evolutionary Games and Population Dynamics*, Cambridge University Press.
- May, R. M. (1972), Will a large complex system be stable?, *Nature* 238, 413–414.
- Williams, R. J. and Martinez, N. D. (2000), Simple rules yield complex food webs, *Nature* 404, 180–183.
- Hill, M. O. (1973), Diversity and evenness: a unifying notation and its consequences, *Ecology* 54, 427–432.
- Chao, A. (1984), Nonparametric estimation of the number of classes in a population, *Scandinavian Journal of Statistics* 11, 265–270.
- MacArthur, R. H. and Wilson, E. O. (1967), *The Theory of Island Biogeography*, Princeton University Press.
- Hubbell, S. P. (2001), *The Unified Neutral Theory of Biodiversity and Biogeography*, Princeton University Press.
- Phillips, S. J., Anderson, R. P. and Schapire, R. E. (2006), Maximum entropy modeling of species geographic distributions, *Ecological Modelling* 190, 231–259.
- Renner, I. W. and Warton, D. I. (2013), Equivalence of MAXENT and Poisson point process models for species distribution modeling in ecology, *Biometrics* 69, 274–281.
- Lindeman, R. L. (1942), The trophic-dynamic aspect of ecology, *Ecology* 23, 399–417.
- Charnov, E. L. (1976), Optimal foraging, the marginal value theorem, *Theoretical Population Biology* 9, 129–136.
- Schaefer, M. B. (1954), Some aspects of the dynamics of populations important to the management of the commercial marine fisheries, *Bulletin of the Inter-American Tropical Tuna Commission* 1, 27–56.
- Beverton, R. J. H. and Holt, S. J. (1957), *On the Dynamics of Exploited Fish Populations*, Fishery Investigations, Series II, 19, HMSO.
- Ricker, W. E. (1954), Stock and recruitment, *Journal of the Fisheries Research Board of Canada* 11, 559–623.
- Brown, J. H., Gillooly, J. F., Allen, A. P., Savage, V. M. and West, G. B. (2004), Toward a metabolic theory of ecology, *Ecology* 85, 1771–1789.
