// Checks of secondary structures and of the one recursion for additive weights: dot-bracket read
// and written back; the rules a structure breaks; the number of structures against their
// enumeration and, for minimum loop 0, against an independent count (the paired positions chosen,
// then matched without crossing in a Catalan number of ways); Nussinov's maximal pairing and the
// maximum expected accuracy structure against the best of every structure; the expected
// accuracy of the MEA optimum never decreasing in γ; the base-pair distance a metric.

#include <genomics/rna/pairing.hpp>
#include <genomics/rna/structure.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::rna;
using morphe::arithmetic::Integer;

namespace {

Codes random_rna(calx::random::Stream& r, std::size_t n) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(4));
  return x;
}

Integer binomial(std::size_t n, std::size_t k) {
  Integer b(1);
  for (std::size_t i = 0; i < k; i++) b = b * Integer(static_cast<long long>(n - i)) / Integer(static_cast<long long>(i + 1));
  return b;
}
Integer catalan(std::size_t k) { return binomial(2 * k, k) / Integer(static_cast<long long>(k + 1)); }

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(2026, 101, 0, 0);

  // ---- dot-bracket and validity -------------------------------------------------------------------
  {
    bool round = true;
    std::size_t seen = 0;
    for (std::size_t n = 0; n <= 9; n++)
      enumerate_structures(n, 0, [](std::size_t, std::size_t) { return true; }, [&](const Structure& y) {
        round = round && parse_dot_bracket(dot_bracket(y)) == y && violation(y).empty();
        seen++;
      });
    s.require("dot-bracket reads back every structure of up to nine bases", round, std::to_string(seen) + " structures");
    bool refused = true;
    for (const char* bad : {"(()", "())", ")(", "(.x.)"}) {
      try {
        parse_dot_bracket(bad);
        refused = false;
      } catch (const std::invalid_argument&) {}
    }
    s.require("unbalanced brackets and unknown characters are refused", refused);
    Structure crossing(6);
    crossing.partner = {3, 4, -1, 0, 1, -1};  // (0, 3) and (1, 4) cross
    const Codes x{sequence::G, sequence::C, sequence::A, sequence::C, sequence::G, sequence::A};
    Structure hairpin = parse_dot_bracket("(...)."), wobble_pair = parse_dot_bracket("(..)..");
    const Codes gu{sequence::G, sequence::A, sequence::A, sequence::T, sequence::A, sequence::A};
    s.require("crossing pairs, short loops and non-canonical pairs are named",
              violation(crossing) == "two pairs cross" && violation(hairpin, 4) == "a loop shorter than the minimum" &&
                  violation(hairpin, 3, &x) == "a pair that is not canonical" && violation(wobble_pair, 2, &gu).empty() &&
                  violation(wobble_pair, 2, &gu, false) == "a pair that is not canonical",
              violation(crossing) + "; " + violation(hairpin, 4) + "; " + violation(hairpin, 3, &x));
  }

  // ---- the number of structures ---------------------------------------------------------------------
  {
    bool agree = true, motzkin = true;
    std::string row;
    for (std::size_t m = 0; m <= 3; m++)
      for (std::size_t n = 0; n <= 14; n++) {
        long long count = 0;
        enumerate_structures(n, m, [](std::size_t, std::size_t) { return true; }, [&](const Structure&) { count++; });
        agree = agree && structure_count(n, m) == Integer(count);
        if (m == 0) {
          Integer independent(0);
          for (std::size_t k = 0; 2 * k <= n; k++) independent = independent + binomial(n, 2 * k) * catalan(k);
          motzkin = motzkin && independent == structure_count(n, 0);
        }
        if (m == 3 && n == 14) row = "S_3(14) = " + structure_count(n, m).str();
      }
    s.require("S_m(n) by its recursion equals the enumeration (m ≤ 3, n ≤ 14)", agree, row);
    s.require("with minimum loop 0 it is Σ_k C(n, 2k) Catalan(k), the paired positions chosen and matched without crossing", motzkin);
  }

  // ---- Nussinov against the best of every structure ---------------------------------------------------
  {
    bool best = true, valid = true;
    std::size_t tested = 0;
    for (int trial = 0; trial < 300; trial++) {
      const std::size_t n = 1 + rng.below(14);
      const Codes x = random_rna(rng, n);
      NussinovOptions o;
      o.min_loop = rng.below(4);
      o.wobble = rng.below(2) == 1;
      const auto r = nussinov(x, o);
      std::size_t most = 0;
      enumerate_structures(n, o.min_loop, [&](std::size_t i, std::size_t j) { return canonical(x[i], x[j], o.wobble); },
                           [&](const Structure& y) { most = std::max(most, y.pair_count()); });
      best = best && static_cast<std::size_t>(r.score) == most;
      valid = valid && violation(r.structure, o.min_loop, &x, o.wobble).empty() && r.structure.pair_count() == most;
      tested++;
    }
    s.require("Nussinov's count is the most canonical pairs of any structure (300 sequences up to 14 bases, loops 0–3, with and without G·U)", best,
              std::to_string(tested) + " sequences");
    s.require("its traceback is a structure of canonical pairs with that many pairs and loops of the minimum or more", valid);
  }

  // ---- MEA against the best of every structure, and its monotonicity in γ ------------------------------------
  {
    bool best = true, monotone = true, accuracy = true;
    double worst = 0;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t n = 2 + rng.below(11);
      // Random pair probabilities with every row summing to at most 1.
      PairProbabilities bp;
      bp.n = n;
      bp.p.assign(n * n, 0);
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = i + 1; j < n; j++)
          if (rng.below(3)) bp.p[i * n + j] = bp.p[j * n + i] = rng.uniform();
      double top = 0;
      for (std::size_t i = 0; i < n; i++) {
        double r = 0;
        for (std::size_t j = 0; j < n; j++) r += bp.p[i * n + j];
        top = std::max(top, r);
      }
      for (auto& v : bp.p) v /= top * (1 + rng.uniform());
      const std::size_t min_loop = rng.below(3);
      double previous_pairs = -1;
      for (double gamma : {0.25, 0.5, 1.0, 2.0, 4.0, 8.0}) {
        const auto r = maximum_expected_accuracy(bp, gamma, min_loop);
        double most = -1;
        enumerate_structures(n, min_loop, [&](std::size_t i, std::size_t j) { return bp.at(i, j) > 0; },
                             [&](const Structure& y) { most = std::max(most, expected_accuracy(bp, y, gamma)); });
        worst = std::max(worst, std::fabs(r.score - most));
        best = best && std::fabs(r.score - most) <= 1e-12 * std::max(1.0, most);
        accuracy = accuracy && std::fabs(expected_accuracy(bp, r.structure, gamma) - r.score) <= 1e-12 * std::max(1.0, r.score);
        double paired = 0;
        for (const auto& [i, j] : r.structure.pairs()) paired += 2 * bp.at(i, j);
        monotone = monotone && paired >= previous_pairs - 1e-12;
        previous_pairs = paired;
      }
    }
    s.require("the MEA structure attains the largest expected accuracy of every structure (200 instances, six γ)", best);
    s.check("its largest deviation from the enumeration", worst, 1e-12);
    s.require("its score is the expected accuracy of the structure it returns", accuracy);
    s.require("the expected number of correct paired positions Σ 2 p_ij of the optimum never decreases as γ grows", monotone);
  }

  // ---- the base-pair distance ------------------------------------------------------------------------------
  {
    std::vector<Structure> all;
    enumerate_structures(8, 1, [](std::size_t, std::size_t) { return true; }, [&](const Structure& y) { all.push_back(y); });
    bool metric = true;
    for (int t = 0; t < 2000; t++) {
      const auto& a = all[rng.below(all.size())];
      const auto& b = all[rng.below(all.size())];
      const auto& c = all[rng.below(all.size())];
      const std::size_t ab = base_pair_distance(a, b);
      metric = metric && ab == base_pair_distance(b, a) && (ab == 0) == (a == b) && base_pair_distance(a, c) <= ab + base_pair_distance(b, c);
    }
    s.require("the base-pair distance is a metric (symmetric, zero only between equal structures, the triangle inequality)", metric,
              std::to_string(all.size()) + " structures of 8 bases");
  }

  return s.done();
}
