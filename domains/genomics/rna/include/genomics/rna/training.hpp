// Training a structure grammar from sequences alone: expectation maximisation by morphe's
// inside–outside, untied or with the rule probabilities tied along the orbits of a group acting
// on the grammar (the equivariant EM of the user's question Q1), the held-out likelihood, and
// the likelihood-ratio test of a tie (theory/rna-structure.md §5).
//
// EM. The E-step is morphe's expected_counts, the expected number C_r of uses of each rule in the
// derivations of each sequence, summed over the sequences (Lari and Young 1990, The estimation of
// stochastic context-free grammars using the Inside-Outside algorithm, Computer Speech and
// Language 4, 35–56); the M-step sets θ_r = C_r / Σ_{r′ of the same left side} C_r′. Each iteration
// never decreases the likelihood (Dempster, Laird and Rubin 1977, Maximum likelihood from
// incomplete data via the EM algorithm, J. R. Statist. Soc. B 39, 1–38).
//
// The tie. A group G acting on the grammar's symbols and fixing its nonterminals permutes the
// rules of each left side (morphe::formal::rule_orbits); the tied model asks θ to be constant on
// each orbit O. Its M-step is
//
//   θ_r = C_O / (|O| Σ_{O′ of the same left side} C_O′)   for r ∈ O,   C_O = Σ_{r ∈ O} C_r,
//
// which is the exact maximiser of the expected complete log-likelihood Σ_r C_r log θ_r on the
// tied parameters (theory, Theorem 5.2), so the tied EM also never decreases the likelihood.
// A tied θ makes P(g·x) = P(x) for every g ∈ G and every sequence x (theory, Proposition 5.4).
// Each iteration costs one inside–outside per sequence, O(Σ_k n_k³ |R_binary|).
//
// The test. With ℓ̂ the maximised log-likelihoods of the untied and the tied models on the same
// sequences, Λ = 2(ℓ̂_untied − ℓ̂_tied) is asymptotically χ² with as many degrees of freedom as
// the tie removes, Σ_A (#rules of A − #orbits of A) (Wilks 1938, The large-sample distribution of
// the likelihood ratio for testing composite hypotheses, Ann. Math. Statist. 9, 60–62). The
// asymptotics assume interior parameters and maximised likelihoods; EM finds local maxima, and
// the untied fit is started from the tied one so that Λ ≥ 0.
//
// Determinism: the E-step runs per sequence on calx's thread pool, each sequence's counts in its
// own slot, summed in order, so the result is the same on any number of threads.

#pragma once

#include <morphe/formal/equivariant.hpp>

#include <calx/scheduler/parallel_for.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/special_functions/gamma.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "scfg.hpp"

namespace genomics::rna {

// A partition of the rules into orbits, each within one left side.
struct Tie {
  std::string name = "untied";
  std::vector<std::uint32_t> orbit_of;  // per rule
  std::vector<std::size_t> size;        // per orbit
  std::vector<std::size_t> lhs;         // per orbit
  std::size_t orbits() const { return size.size(); }
};

// Every rule its own orbit.
inline Tie untied(const StructureGrammar& sg) {
  Tie t;
  for (std::uint32_t r = 0; r < sg.rules(); r++) {
    t.orbit_of.push_back(r);
    t.size.push_back(1);
    t.lhs.push_back(sg.info[r].lhs);
  }
  return t;
}

// The orbits of the rules under a group acting on the grammar (it must fix every nonterminal, so
// that each orbit keeps one left side); std::invalid_argument when the grammar is not invariant.
inline Tie orbit_tie(const StructureGrammar& sg, const mf::GrammarAction& action, std::string name) {
  if (!action.fixes_nonterminals()) throw std::invalid_argument("genomics::rna::orbit_tie: the action must fix every nonterminal");
  const auto orb = mf::rule_orbits(sg.grammar, action);
  Tie t;
  t.name = std::move(name);
  t.orbit_of = orb.orbit_of;
  t.size = orb.size;
  t.lhs.assign(orb.size.size(), SIZE_MAX);
  for (std::uint32_t r = 0; r < sg.rules(); r++) {
    auto& l = t.lhs[orb.orbit_of[r]];
    if (l != SIZE_MAX && l != sg.info[r].lhs) throw std::logic_error("genomics::rna::orbit_tie: an orbit with two left sides");
    l = sg.info[r].lhs;
  }
  return t;
}

// The free parameters of a model: Σ_A (#orbits of A − 1).
inline std::size_t free_parameters(const Tie& t) {
  std::vector<std::size_t> per;
  for (std::size_t o = 0; o < t.orbits(); o++) {
    if (t.lhs[o] >= per.size()) per.resize(t.lhs[o] + 1, 0);
    per[t.lhs[o]]++;
  }
  std::size_t f = 0;
  for (std::size_t a : per) f += a ? a - 1 : 0;
  return f;
}

// Whether θ is constant on every orbit, to a tolerance.
inline bool is_tied(const Tie& t, const Parameters& theta, double tolerance = 1e-12) {
  std::vector<double> first(t.orbits(), -1);
  for (std::size_t r = 0; r < theta.size(); r++) {
    double& f = first[t.orbit_of[r]];
    if (f < 0) f = theta[r];
    else if (std::fabs(f - theta[r]) > tolerance) return false;
  }
  return true;
}

// The orbit averages of θ: a tied θ with the same sums on each left side.
inline Parameters project(const Tie& t, const Parameters& theta) {
  std::vector<double> mean(t.orbits(), 0);
  for (std::size_t r = 0; r < theta.size(); r++) mean[t.orbit_of[r]] += theta[r];
  for (std::size_t o = 0; o < t.orbits(); o++) mean[o] /= static_cast<double>(t.size[o]);
  Parameters out(theta.size());
  for (std::size_t r = 0; r < theta.size(); r++) out[r] = mean[t.orbit_of[r]];
  return out;
}

// ---- E-step, M-step, Q ---------------------------------------------------------------------------------
struct Expectation {
  std::vector<double> counts;  // Σ over the sequences of the expected rule counts
  double log_likelihood = 0;   // Σ log P(x)
};

inline Expectation expectation(const StructureGrammar& sg, const Parameters& theta, const std::vector<Codes>& data,
                               calx::scheduler::ThreadPool* pool = nullptr) {
  std::vector<std::vector<double>> part(data.size());
  std::vector<double> logp(data.size(), 0);
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t i = begin; i < end; i++) part[i] = expected_rule_counts(sg, theta, data[i], &logp[i]);
  };
  if (pool) calx::scheduler::parallel_for(*pool, data.size(), 1, body);
  else body(0, data.size(), 0);
  Expectation e;
  e.counts.assign(sg.rules(), 0);
  for (std::size_t i = 0; i < data.size(); i++) {
    for (std::size_t r = 0; r < sg.rules(); r++) e.counts[r] += part[i][r];
    e.log_likelihood += logp[i];
  }
  return e;
}

// The maximiser of Σ_r C_r log θ_r over the parameters tied by t; a left side whose counts are
// all zero keeps the previous θ.
inline Parameters maximisation(const StructureGrammar& sg, const Tie& t, const std::vector<double>& counts, const Parameters& previous) {
  std::vector<double> orbit(t.orbits(), 0);
  for (std::size_t r = 0; r < counts.size(); r++) orbit[t.orbit_of[r]] += counts[r];
  std::vector<double> side(sg.nonterminals(), 0);
  for (std::size_t o = 0; o < t.orbits(); o++) side[t.lhs[o]] += orbit[o];
  Parameters theta(counts.size());
  for (std::size_t r = 0; r < counts.size(); r++) {
    const std::size_t o = t.orbit_of[r], a = sg.info[r].lhs;
    theta[r] = side[a] > 0 ? orbit[o] / (static_cast<double>(t.size[o]) * side[a]) : previous[r];
  }
  return theta;
}

// Q(θ) = Σ_r C_r log θ_r, the expected complete log-likelihood up to a constant (−∞ when a
// counted rule has θ_r = 0).
inline double expected_complete_log_likelihood(const std::vector<double>& counts, const Parameters& theta) {
  double q = 0;
  for (std::size_t r = 0; r < counts.size(); r++) {
    if (counts[r] == 0) continue;
    if (!(theta[r] > 0)) return -std::numeric_limits<double>::infinity();
    q += counts[r] * std::log(theta[r]);
  }
  return q;
}

// ---- the iteration --------------------------------------------------------------------------------------
struct TrainOptions {
  int max_iterations = 40;
  double tolerance = 1e-5;  // stop when the log-likelihood rises by less than this per base
};

struct Training {
  Parameters theta;
  std::vector<double> log_likelihood;  // before each M-step: ℓ(θ_0), ℓ(θ_1), …, and ℓ(theta) last
  int iterations = 0;
  bool converged = false;
  bool monotone = true;  // no iteration lowered ℓ beyond rounding (1e-9 relative)
};

// EM from theta0 (projected onto the tie first, so that every iterate is tied).
inline Training train(const StructureGrammar& sg, const Tie& t, const Parameters& theta0, const std::vector<Codes>& data,
                      const TrainOptions& o = {}, calx::scheduler::ThreadPool* pool = nullptr) {
  validate(sg, theta0);
  std::size_t bases = 0;
  for (const auto& x : data) bases += x.size();
  Training out;
  out.theta = project(t, theta0);
  const double tol = o.tolerance * static_cast<double>(std::max<std::size_t>(bases, 1));
  for (int it = 0; it < o.max_iterations; it++) {
    const auto e = expectation(sg, out.theta, data, pool);
    const double l = e.log_likelihood;
    if (!out.log_likelihood.empty()) {
      const double prev = out.log_likelihood.back();
      if (l < prev - 1e-9 * std::fabs(prev)) out.monotone = false;
      if (l - prev < tol) {
        out.log_likelihood.push_back(l);
        out.converged = true;
        break;
      }
    }
    out.log_likelihood.push_back(l);
    out.theta = maximisation(sg, t, e.counts, out.theta);
    out.iterations = it + 1;
  }
  if (!out.converged) out.log_likelihood.push_back(expectation(sg, out.theta, data, pool).log_likelihood);
  return out;
}

// Σ log P(x) over the sequences.
inline double log_likelihood(const StructureGrammar& sg, const Parameters& theta, const std::vector<Codes>& data,
                             calx::scheduler::ThreadPool* pool = nullptr) {
  std::vector<double> l(data.size(), 0);
  const auto body = [&](std::size_t begin, std::size_t end, std::size_t) {
    for (std::size_t i = begin; i < end; i++) l[i] = log_likelihood(sg, theta, data[i]);
  };
  if (pool) calx::scheduler::parallel_for(*pool, data.size(), 1, body);
  else body(0, data.size(), 0);
  double s = 0;
  for (double v : l) s += v;
  return s;
}

// ---- the likelihood-ratio test -----------------------------------------------------------------------------
struct RatioTest {
  double statistic = 0;  // Λ = 2(ℓ_full − ℓ_tied)
  std::size_t df = 0;
  double p_value = 1;    // P(χ²_df ≥ Λ)
};

inline RatioTest likelihood_ratio(double full, std::size_t full_free, double tied, std::size_t tied_free) {
  if (full_free <= tied_free) throw std::invalid_argument("genomics::rna::likelihood_ratio: the tied model must have fewer parameters");
  RatioTest t;
  t.statistic = std::max(0.0, 2 * (full - tied));
  t.df = full_free - tied_free;
  t.p_value = calx::special_functions::chi_square_q(t.statistic, static_cast<double>(t.df));
  return t;
}

}  // namespace genomics::rna
