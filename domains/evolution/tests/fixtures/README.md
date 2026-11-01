# evolution/tests/fixtures: Reference Data

Reference values the evolution domain's checks read (`EVOLUTION_FIXTURES`), each made by the script beside it.

| File | Content | Source | Retrieved |
|---|---|---|---|
| `coloc.txt` | Three seeded scenarios of two traits' estimates and variances at 60 variants in AR(1) linkage (a shared causal variant, two distinct ones, one trait null) and R's `coloc.abf` posteriors PP.H0–PP.H4 and SNP.PP.H4 (first trait quantitative with sdY = 1, second case–control; default priors) | `make_coloc_fixtures.py`, running R 4.4.3 with coloc 5.2.3 (Giambartolomei et al. 2014; CRAN, GPL-3, used as a reference program, not distributed) | 2026-11-01 |
