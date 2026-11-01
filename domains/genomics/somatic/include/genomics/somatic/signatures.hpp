// Mutational signatures: their extraction de novo from a cohort's catalogues by non-negative
// matrix factorisation with the stability of each rank, the matching of signatures to a reference
// set, the refitting of a tumour's catalogue to known signatures by non-negative least squares,
// and the bootstrap of its exposures. The factorisation, the stability analysis and NNLS are
// tyche's (multivariate/nmf.hpp, regression/nnls.hpp), the assignment calx's; this file says what
// they mean for catalogues.
//
// - The model (Alexandrov et al. 2013, Cell Reports 3, 246–259; Nature 500, 415–421). A cohort's
//   catalogues V (96 classes × n tumours) are counts with means W H: W's columns the signatures,
//   distributions over the classes, H the exposures, the mutations each signature contributed to
//   each tumour. The Poisson likelihood of V is the Kullback–Leibler objective of Lee and Seung's
//   updates, the one SigProfilerExtractor uses by default (Islam et al. 2022, Cell Genomics 2,
//   100179). The factorisation is not identifiable in general: W H = (W D P)(P⁻¹ D⁻¹ H) for any
//   permutation P and positive diagonal D, which normalising W's columns to sum 1 removes up to
//   the permutation; beyond that, rotations W Q, Q⁻¹ H that keep both non-negative exist unless the
//   factors are sufficiently sparse (theory/the-somatic-genome.md, §5). The stability of a rank is
//   how reproducible its signatures are over resampled catalogues (tyche's nmf_stability).
// - The rank: the largest stable one (average stability ≥ 0.80 and every cluster's ≥ 0.20,
//   SigProfilerExtractor's thresholds for genomes); SigProfilerExtractor chooses among the stable
//   ranks by a test of their reconstruction errors, not written here.
// - Matching (match_signatures): the assignment of extracted to reference signatures maximising
//   the total cosine similarity (calx's Hungarian method), each pair's cosine reported; the
//   measure of Islam et al. 2022's comparisons with COSMIC.
// - Refitting (refit): exposures e ≥ 0 minimising ‖S e − c‖₂ for a catalogue c and signatures S
//   (columns summing to 1), by tyche's NNLS; the cosine of S e and c says how well the signatures
//   explain the catalogue (SigProfilerAssignment's measure; its forward and backward selection of
//   signatures is not written).
// - The bootstrap (bootstrap_exposures): B catalogues drawn as Multinomial(Σc, c/Σc), calx's exact
//   sampler, keyed by the seed and the replicate, each refitted; the exposures' means, standard
//   deviations and 2.5% and 97.5% quantiles (tyche's type-7 quantiles). It measures the sampling
//   noise of a catalogue of finite size, not the model's error (Huang, Wojtowicz and Przytycka 2018,
//   Bioinformatics 34, 330–337, bootstrap the same way).
//
// Complexity: the factorisations' (O(m n k) per update); NNLS O(m k²) per fit.

#pragma once

#include <genomics/sequence/contexts.hpp>

#include <tyche/descriptive/summary.hpp>
#include <tyche/multivariate/nmf.hpp>
#include <tyche/regression/nnls.hpp>

#include <calx/linalg/rectangular.hpp>
#include <calx/optimisation/assignment.hpp>
#include <calx/random/multinomial.hpp>
#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace genomics::somatic {

using calx::linalg::Rectangular;

// The matrix of profiles (each of m classes) as m × n columns.
inline Rectangular profile_matrix(const std::vector<std::vector<double>>& profiles) { return Rectangular::from_columns(profiles); }

inline std::vector<std::vector<double>> columns_of(const Rectangular& M) {
  std::vector<std::vector<double>> out(M.cols, std::vector<double>(M.rows));
  for (std::size_t j = 0; j < M.cols; j++)
    for (std::size_t i = 0; i < M.rows; i++) out[j][i] = M(i, j);
  return out;
}

struct Extraction {
  std::vector<tyche::multivariate::RankStability> ranks;  // rank_min, …, rank_max
  std::size_t chosen = 0;                                 // the index in `ranks` of the rank chosen
  bool stable = false;                                    // whether any rank was stable
  const tyche::multivariate::RankStability& solution() const { return ranks.at(chosen); }
};

// De novo extraction from catalogues (one per tumour, each over the same classes).
inline Extraction extract_signatures(const std::vector<std::vector<double>>& catalogues, std::size_t rank_min, std::size_t rank_max,
                                     const tyche::multivariate::StabilityOptions& o = {}) {
  if (catalogues.empty()) throw std::invalid_argument("genomics::somatic::extract_signatures: at least one catalogue");
  Extraction e;
  e.ranks = tyche::multivariate::nmf_rank_stability(profile_matrix(catalogues), rank_min, rank_max, o);
  for (std::size_t k = 0; k < e.ranks.size(); k++)
    if (e.ranks[k].stable) {
      e.chosen = k;
      e.stable = true;
    }
  return e;
}

struct SignatureMatch {
  std::vector<std::size_t> reference;  // for each extracted signature, its reference
  std::vector<double> cosine;          // and their cosine similarity
  double mean_cosine = 0;
};

// The assignment of each extracted signature to a distinct reference one of greatest total cosine.
inline SignatureMatch match_signatures(const std::vector<std::vector<double>>& extracted, const std::vector<std::vector<double>>& reference) {
  if (extracted.size() > reference.size()) throw std::invalid_argument("genomics::somatic::match_signatures: no more extracted signatures than references");
  Rectangular cost(extracted.size(), reference.size());
  for (std::size_t i = 0; i < extracted.size(); i++)
    for (std::size_t j = 0; j < reference.size(); j++) cost(i, j) = -tyche::multivariate::cosine_similarity(extracted[i], reference[j]);
  const auto a = calx::optimisation::assignment(cost);
  SignatureMatch m;
  m.reference = a.column;
  for (std::size_t i = 0; i < extracted.size(); i++) {
    m.cosine.push_back(-cost(i, a.column[i]));
    m.mean_cosine += m.cosine.back();
  }
  if (!extracted.empty()) m.mean_cosine /= static_cast<double>(extracted.size());
  return m;
}

struct Refit {
  std::vector<double> exposures;     // mutations attributed to each signature
  std::vector<double> reconstruction;
  double cosine = 0;                 // of the reconstruction and the catalogue
  double residual = 0;               // ‖S e − c‖₂
};

// The exposures of a catalogue to signatures (each normalised to sum 1 here).
inline Refit refit(std::span<const double> catalogue, const std::vector<std::vector<double>>& signatures) {
  if (signatures.empty()) throw std::invalid_argument("genomics::somatic::refit: at least one signature");
  std::vector<std::vector<double>> s = signatures;
  for (auto& x : s) {
    if (x.size() != catalogue.size()) throw std::invalid_argument("genomics::somatic::refit: signatures over the catalogue's classes");
    double t = 0;
    for (double v : x) t += v;
    if (t > 0)
      for (double& v : x) v /= t;
  }
  const Rectangular S = Rectangular::from_columns(s);
  const auto r = tyche::regression::nnls(S, catalogue);
  Refit f;
  f.exposures = r.x;
  f.residual = r.residual_norm;
  f.reconstruction.assign(catalogue.size(), 0);
  for (std::size_t i = 0; i < catalogue.size(); i++)
    for (std::size_t j = 0; j < s.size(); j++) f.reconstruction[i] += S(i, j) * r.x[j];
  f.cosine = tyche::multivariate::cosine_similarity(f.reconstruction, std::vector<double>(catalogue.begin(), catalogue.end()));
  return f;
}

struct ExposureInterval {
  double estimate = 0, mean = 0, sd = 0, low = 0, high = 0;  // the refit's, then the bootstrap's
};

struct Bootstrap {
  std::vector<ExposureInterval> exposures;
  std::size_t replicates = 0;
};

inline Bootstrap bootstrap_exposures(std::span<const double> catalogue, const std::vector<std::vector<double>>& signatures, std::size_t replicates, std::uint64_t seed) {
  const Refit base = refit(catalogue, signatures);
  Bootstrap b;
  b.replicates = replicates;
  const std::size_t k = signatures.size();
  double total = 0;
  for (double x : catalogue) total += x;
  const auto n = static_cast<std::int64_t>(std::llround(total));
  std::vector<std::vector<double>> draws(k);
  for (std::size_t r = 0; r < replicates && n > 0; r++) {
    calx::random::Stream rng(seed, 0x424f4f54u /* "BOOT" */, static_cast<std::uint32_t>(r), 0);
    const auto counts = calx::random::multinomial(rng, n, catalogue);
    std::vector<double> c(counts.begin(), counts.end());
    const Refit f = refit(c, signatures);
    for (std::size_t j = 0; j < k; j++) draws[j].push_back(f.exposures[j]);
  }
  b.exposures.resize(k);
  for (std::size_t j = 0; j < k; j++) {
    ExposureInterval& e = b.exposures[j];
    e.estimate = base.exposures[j];
    if (draws[j].empty()) {
      e.mean = e.low = e.high = e.estimate;
      continue;
    }
    e.mean = tyche::descriptive::mean(draws[j]);
    e.sd = draws[j].size() > 1 ? tyche::descriptive::standard_deviation(draws[j]) : 0;
    e.low = tyche::descriptive::quantile(draws[j], 0.025);
    e.high = tyche::descriptive::quantile(draws[j], 0.975);
  }
  return b;
}

}  // namespace genomics::somatic
