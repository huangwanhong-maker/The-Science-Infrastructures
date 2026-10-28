// Checks of the structure grammars against their definitions: which symmetries each admits and
// the pair sets' closures; the number of derivations of a word (morphe's counting semiring) equal
// to the number of structures the grammar should derive, so that each is derived once; on short
// sequences under random probabilities, the inside probability, the pair probabilities, the
// expected rule counts, the Viterbi structure and the MEA structure against the enumeration of
// every structure; Knudsen and Hein's grammar as a sub-model of the reversible one; derivations
// refused where the grammar has none; sampling; and the invariance P(g·x) = P(x) of a tied model.

#include <genomics/rna/scfg.hpp>
#include <genomics/rna/symmetry.hpp>
#include <genomics/rna/training.hpp>

#include <morphe/formal/equivariant.hpp>

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

// Random rule probabilities, each left side normalised, every rule at least 0.05 of the mean.
Parameters random_parameters(const StructureGrammar& sg, calx::random::Stream& r) {
  Parameters t(sg.rules());
  std::vector<double> sum(sg.nonterminals(), 0);
  for (std::size_t k = 0; k < t.size(); k++) sum[sg.info[k].lhs] += t[k] = 0.05 + r.uniform();
  for (std::size_t k = 0; k < t.size(); k++) t[k] /= sum[sg.info[k].lhs];
  return t;
}

std::vector<mf::Symbol> word(const Codes& x) { return {x.begin(), x.end()}; }

struct Variant {
  Shape shape;
  bool lonely;
  PairSet pairs;
};

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(2026, 102, 0, 0);
  const std::vector<Variant> variants{{Shape::knudsen_hein, true, ALL_PAIRS}, {Shape::reversible, true, ALL_PAIRS}, {Shape::knudsen_hein, false, ALL_PAIRS},
                                      {Shape::reversible, false, ALL_PAIRS},  {Shape::reversible, false, canonical_pairs()}};

  // ---- the grammars and their symmetries ---------------------------------------------------------------
  {
    const auto& kh = knudsen_hein();
    const auto& khr = knudsen_hein_reversible();
    const auto& canon = structure_grammar(Shape::reversible, false, canonical_pairs());
    const auto sym = pair_symmetries();  // ρ, κ, κρ, the strand group
    s.require("Knudsen–Hein's grammar: 39 rules, 66 in binary form; the reversible one 41 and 70; canonical pairs without lonely ones 27 and 42",
              kh.rules() == 39 && kh.binary.binary.size() == 66 && khr.rules() == 41 && khr.binary.binary.size() == 70 && canon.rules() == 27 &&
                  canon.binary.binary.size() == 42,
              std::to_string(kh.rules()) + ", " + std::to_string(kh.binary.binary.size()) + "; " + std::to_string(khr.rules()) + ", " +
                  std::to_string(khr.binary.binary.size()) + "; " + std::to_string(canon.rules()) + ", " + std::to_string(canon.binary.binary.size()));
    s.require("Knudsen–Hein's grammar admits the complement κ but not reversal ρ, κρ or the strand group (S → L S reversed is no rule)",
              admits(kh, sym[1]) && !admits(kh, sym[0]) && !admits(kh, sym[2]) && !admits(kh, sym[3]));
    s.require("the reversible grammar admits all four", admits(khr, sym[0]) && admits(khr, sym[1]) && admits(khr, sym[2]) && admits(khr, sym[3]));
    s.require("with the six canonical pairs it admits ρ alone: κ and κρ carry G·U to C·A and A·C",
              admits(canon, sym[0]) && !admits(canon, sym[1]) && !admits(canon, sym[2]) && !admits(canon, sym[3]));
    const std::string closed = pair_names(closure(sym[1], canonical_pairs()));
    s.require("the closure of the canonical pairs under κ, κρ and the strand group adds A·C and C·A; under ρ it adds nothing",
              closed == "AC AU CA CG GC GU UA UG" && closure(sym[2], canonical_pairs()) == closure(sym[1], canonical_pairs()) &&
                  closure(sym[3], canonical_pairs()) == closure(sym[1], canonical_pairs()) && closure(sym[0], canonical_pairs()) == canonical_pairs(),
              closed);
    bool burnside = true;
    std::string orbits;
    for (const auto& g : sym) {
      const auto orb = mf::rule_orbits(khr.grammar, grammar_action(g));
      burnside = burnside && orb.burnside == orb.representatives.size();
      orbits += (orbits.empty() ? "" : "; ") + g.name + " " + std::to_string(orb.representatives.size());
    }
    s.require("the rule orbits of the reversible grammar number as Burnside's lemma counts them", burnside, orbits);
  }

  // ---- each structure derived once ---------------------------------------------------------------------------
  {
    bool once = true, closed_form = true;
    std::string row;
    for (const auto& v : variants) {
      const auto& sg = structure_grammar(v.shape, v.lonely, v.pairs);
      const std::vector<Integer> unit(sg.rules(), Integer(1));
      for (std::size_t n = 1; n <= 12; n++) {
        const Codes x = random_rna(rng, n);
        const auto in = mf::inside<mf::Counting>(sg.binary, unit, word(x));
        long long structures = 0;
        enumerate_structures(n, 2, [](std::size_t, std::size_t) { return true; }, [&](const Structure& y) { structures += derivable(sg, x, y); });
        once = once && in.at(0, 0, n) == Integer(structures);
        if (v.lonely && v.pairs == ALL_PAIRS) closed_form = closed_form && in.at(0, 0, n) == structure_count(n, 2);
        if (n == 12 && v.shape == Shape::reversible && !v.lonely && v.pairs == ALL_PAIRS) row = "no lonely pairs, 12 bases: " + in.at(0, 0, n).str();
      }
    }
    s.require("the number of derivations of a word is the number of structures each grammar should derive (five grammars, 1–12 bases)", once, row);
    s.require("with lonely pairs allowed and every pair a rule it is S_2(n), the structures with loops of two or more", closed_form);
  }

  // ---- probabilities against the enumeration -------------------------------------------------------------------
  {
    bool inside = true, pairs = true, counts = true, viterbi_ok = true, mea_ok = true;
    double worst_inside = 0, worst_pairs = 0;
    std::size_t words = 0;
    for (const auto& v : variants) {
      const auto& sg = structure_grammar(v.shape, v.lonely, v.pairs);
      for (int trial = 0; trial < 12; trial++) {
        const Parameters theta = random_parameters(sg, rng);
        const std::size_t n = 4 + rng.below(7);
        const Codes x = random_rna(rng, n);
        double z = 0, best = -1e300;
        std::vector<double> pp(n * n, 0), expect(sg.rules(), 0);
        std::vector<std::pair<double, Structure>> all;
        enumerate_structures(n, 2, [](std::size_t, std::size_t) { return true; }, [&](const Structure& y) {
          if (!derivable(sg, x, y)) return;
          const double lp = log_joint(sg, theta, x, y), p = std::exp(lp);
          z += p;
          best = std::max(best, lp);
          for (const auto& [i, j] : y.pairs()) pp[i * n + j] += p;
          const auto c = rule_counts(sg, x, y);
          for (std::size_t r = 0; r < c.size(); r++) expect[r] += p * c[r];
          all.push_back({lp, y});
        });
        words++;
        const double l = log_likelihood(sg, theta, x);
        worst_inside = std::max(worst_inside, std::fabs(std::exp(l) - z) / z);
        inside = inside && std::fabs(std::exp(l) - z) <= 1e-12 * z;
        const auto bp = pair_probabilities(sg, theta, x);
        for (std::size_t i = 0; i < n; i++)
          for (std::size_t j = i + 1; j < n; j++) {
            worst_pairs = std::max(worst_pairs, std::fabs(bp.at(i, j) - pp[i * n + j] / z));
            pairs = pairs && std::fabs(bp.at(i, j) - pp[i * n + j] / z) <= 1e-12;
          }
        double lp = 0;
        const auto c = expected_rule_counts(sg, theta, x, &lp);
        for (std::size_t r = 0; r < c.size(); r++) counts = counts && std::fabs(c[r] - expect[r] / z) <= 1e-10;
        const auto f = viterbi(sg, theta, x);
        viterbi_ok = viterbi_ok && std::fabs(f.log_probability - best) <= 1e-10 && std::fabs(log_joint(sg, theta, x, f.structure) - best) <= 1e-10;
        // MEA over every structure with loops of two or more, on the grammar's pair probabilities.
        for (double gamma : {0.5, 1.0, 3.0}) {
          const auto m = mea_fold(sg, theta, x, gamma);
          double most = -1;
          enumerate_structures(n, 2, [&](std::size_t i, std::size_t j) { return bp.at(i, j) > 0; },
                               [&](const Structure& y) { most = std::max(most, expected_accuracy(bp, y, gamma)); });
          mea_ok = mea_ok && std::fabs(m.score - most) <= 1e-12 * std::max(1.0, most);
        }
      }
    }
    s.require("the inside probability is the sum of P(x, y) over the structures (five grammars, 60 words of 4–10 bases)", inside, std::to_string(words) + " words");
    s.check("its largest relative deviation", worst_inside, 1e-12);
    s.require("the pair probabilities from inside and outside are the enumerated marginals", pairs);
    s.check("their largest deviation", worst_pairs, 1e-12);
    s.require("the expected rule counts are the enumerated ones", counts);
    s.require("the Viterbi structure attains the most probable structure, and its derivation has that probability", viterbi_ok);
    s.require("the MEA structure on the grammar's pair probabilities attains the best expected accuracy (γ 0.5, 1, 3)", mea_ok);
  }

  // ---- Knudsen–Hein inside the reversible grammar -------------------------------------------------------------
  {
    bool same = true, joint = true;
    double worst = 0;
    for (bool lonely : {true, false})
      for (int trial = 0; trial < 20; trial++) {
        const auto& a = structure_grammar(Shape::knudsen_hein, lonely);
        const Parameters theta = random_parameters(a, rng);
        const Parameters rev = as_reversible(a, theta);
        const auto& b = structure_grammar(Shape::reversible, lonely);
        validate(b, rev);
        const Codes x = random_rna(rng, 6 + rng.below(30));
        const double la = log_likelihood(a, theta, x), lb = log_likelihood(b, rev, x);
        worst = std::max(worst, std::fabs(la - lb));
        same = same && std::fabs(la - lb) <= 1e-12 * std::fabs(la);
        if (x.size() <= 10)
          enumerate_structures(x.size(), 2, [](std::size_t, std::size_t) { return true; }, [&](const Structure& y) {
            if (derivable(a, x, y)) joint = joint && std::fabs(log_joint(a, theta, x, y) - log_joint(b, rev, x, y)) <= 1e-12;
          });
      }
    s.require("Knudsen–Hein's parameters carried to the reversible grammar give every sequence the same probability", same);
    s.check("the largest difference of log-likelihoods (nats)", worst, 1e-10);
    s.require("and every structure the same probability", joint);
  }

  // ---- derivations refused --------------------------------------------------------------------------------------
  {
    const auto& canon = structure_grammar(Shape::reversible, false, canonical_pairs());
    const auto& kh = knudsen_hein();
    const Codes x = sequence::encode_dna("GGAAACCAAAA");
    const auto refused = [&](const StructureGrammar& g, const char* db) {
      try {
        derivation(g, x, parse_dot_bracket(db));
        return false;
      } catch (const std::invalid_argument&) {
        return true;
      }
    };
    s.require("derivations are refused for a lonely pair without lonely pairs, a loop of one base, and a pair outside the pair set",
              refused(canon, "(.....)....") && !refused(kh, "(.....)....") && refused(kh, "((.))......") && refused(canon, "((...).)...") &&
                  !refused(canon, "((...))....") && refused(canon, "(.((...)))."));
  }

  // ---- sampling ----------------------------------------------------------------------------------------------------
  {
    const auto& sg = structure_grammar(Shape::reversible, false, canonical_pairs());
    const Parameters theta = starting_parameters(sg, {0.85, 0.6, 1});
    bool valid = true;
    std::size_t first_is_single = 0;
    const std::size_t draws = 3000;
    for (std::size_t k = 0; k < draws; k++) {
      const auto smp = sample(sg, theta, rng);
      valid = valid && derivable(sg, smp.sequence, smp.structure) && std::isfinite(log_joint(sg, theta, smp.sequence, smp.structure));
      first_is_single += derivation(sg, smp.sequence, smp.structure).front() == sg.S_L;
    }
    const double f = static_cast<double>(first_is_single) / static_cast<double>(draws), p = theta[sg.S_L];
    const double z = (f - p) / std::sqrt(p * (1 - p) / static_cast<double>(draws));
    s.require("sampled sequences carry structures the grammar derives, of finite probability", valid);
    s.check("the first rule is S → L as often as its probability (standard scores)", z, 4.0);
  }

  // ---- a tied model's invariance -------------------------------------------------------------------------------------
  {
    const auto& khr = knudsen_hein_reversible();
    bool invariant = true, asymmetric = false;
    double worst = 0;
    std::string names;
    for (const auto& g : pair_symmetries()) {
      const Tie t = orbit_tie(khr, g);
      for (int trial = 0; trial < 6; trial++) {
        const Parameters free = random_parameters(khr, rng), tied = project(t, free);
        validate(khr, tied);
        const Codes x = random_rna(rng, 8 + rng.below(30));
        const double l = log_likelihood(khr, tied, x);
        for (std::uint32_t e = 1; e < g.letters->order(); e++) {
          const double d = std::fabs(log_likelihood(khr, tied, act(g, e, x)) - l);
          worst = std::max(worst, d);
          invariant = invariant && d <= 1e-12 * std::fabs(l);
          asymmetric = asymmetric || std::fabs(log_likelihood(khr, free, act(g, e, x)) - log_likelihood(khr, free, x)) > 1e-6;
        }
      }
      names += (names.empty() ? "" : "; ") + g.name;
    }
    s.require("θ constant on the rule orbits gives P(g·x) = P(x) for every element of ρ, κ, κρ and the strand group", invariant, names);
    s.check("the largest difference of log-likelihoods (nats)", worst, 1e-10);
    s.require("untied, the probabilities of x and g·x differ", asymmetric);
  }

  return s.done();
}
