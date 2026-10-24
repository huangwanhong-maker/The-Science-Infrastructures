// Checks of pairhmm.hpp: forward and backward give the same P(x, y); P(x, y) is morphe's generic
// sum over paths in the probability semiring with the model's weights, and the sum of the
// probabilities of every enumerated alignment; Viterbi is the generic sum in the tropical semiring
// of the log weights, and its path has the probability it claims; the posteriors lie in [0, 1],
// sum to at most 1 along a row and a column, and sum to the expected number of aligned pairs
// computed by enumeration; the maximum-expected-accuracy alignment has at least the expected pairs
// of the Viterbi alignment; long sequences neither underflow nor disagree; and the target
// frequencies of BLOSUM62 at its λ form a distribution.

#include <genomics/align/pairhmm.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/align/statistics.hpp>
#include <genomics/simulate/family.hpp>
#include <genomics/simulate/genome.hpp>

#include <morphe/formal/alignment.hpp>
#include <morphe/formal/semiring.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;
namespace mf = morphe::formal;

namespace {

Codes random_codes(calx::random::Stream& r, std::size_t n) {
  Codes x(n);
  for (auto& c : x) c = static_cast<Code>(r.below(4));
  return x;
}

// The model's transitions as alignment.hpp's table [from: begin, M, I, D][to: end, M, I, D], in
// probability (log = false) or log probability.
std::array<std::array<double, 4>, 4> transitions(const PairHmm& h, bool log) {
  const double mm = 1 - 2 * h.delta - h.tau, xm = 1 - h.epsilon - h.tau;
  std::array<std::array<double, 4>, 4> t = {{{h.tau, mm, h.delta, h.delta}, {h.tau, mm, h.delta, h.delta}, {h.tau, xm, h.epsilon, 0}, {h.tau, xm, 0, h.epsilon}}};
  if (log)
    for (auto& row : t)
      for (auto& v : row) v = v > 0 ? std::log(v) : LOG_ZERO;
  return t;
}

// The expected aligned pairs of an alignment under a posterior (m × n).
double pairs_of(const Alignment& a, const std::vector<double>& post, std::size_t n) {
  double total = 0;
  std::size_t i = 0, j = 0;
  for (const auto& o : a.cigar)
    for (std::uint32_t k = 0; k < o.length; k++) {
      if (o.op == 'M') total += post[i++ * n + j++];
      else (o.op == 'I' ? i : j)++;
    }
  return total;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x70686d6d, 0, 0);
  const PairHmm h = PairHmm::dna(0.9, 0.05, 0.4, 0.02);

  {
    double fb = 0, generic = 0, enumerated = 0, viterbi_gap = 0, path_gap = 0, expected = 0;
    bool bounds = true, rows = true, mea = true;
    for (int rep = 0; rep < 80; rep++) {
      const Codes x = random_codes(r, r.below(7)), y = random_codes(r, r.below(7));
      const PairTables f = forward(h, x, y), b = backward(h, x, y);
      fb = std::max(fb, std::fabs(f.log_probability - b.log_probability));
      // The probability semiring.
      const auto prob = mf::three_state<mf::Arithmetic<double>>(
          transitions(h, false), [&](Code a, Code c) { return h.p[a * 4u + c]; }, [&](Code a) { return h.q[a]; }, [&](Code c) { return h.q[c]; });
      const double P = mf::path_sum(x, y, prob, mf::Mode::global);
      generic = std::max(generic, std::fabs(std::log(P) - f.log_probability));
      double total = 0, pairs = 0;
      mf::enumerate_paths(x.size(), y.size(), mf::Mode::global, [&](const mf::EditPath& p) {
        const double w = mf::path_weight(p, x, y, prob);
        total += w;
        pairs += w * static_cast<double>(std::count(p.ops.begin(), p.ops.end(), 'M'));
      });
      enumerated = std::max(enumerated, std::fabs(std::log(total) - f.log_probability));
      // The tropical semiring of the log weights.
      const auto logw = mf::three_state<mf::MaxPlus<double>>(
          transitions(h, true), [&](Code a, Code c) { return h.log_p(a, c); }, [&](Code a) { return h.log_q(a); }, [&](Code c) { return h.log_q(c); });
      const PairViterbi v = viterbi(h, x, y);
      viterbi_gap = std::max(viterbi_gap, std::fabs(mf::path_sum(x, y, logw, mf::Mode::global) - v.log_probability));
      mf::EditPath vp;
      for (const auto& o : v.alignment.cigar) vp.ops.append(o.length, o.op);
      path_gap = std::max(path_gap, std::fabs(mf::path_weight(vp, x, y, logw) - v.log_probability));
      // Posteriors.
      const auto post = posterior(h, x, y);
      double sum = 0;
      for (double q : post) {
        bounds = bounds && q >= 0 && q <= 1 + 1e-12;
        sum += q;
      }
      for (std::size_t i = 0; i < x.size(); i++) {
        double row = 0;
        for (std::size_t j = 0; j < y.size(); j++) row += post[i * y.size() + j];
        rows = rows && row <= 1 + 1e-12;
      }
      for (std::size_t j = 0; j < y.size(); j++) {
        double col = 0;
        for (std::size_t i = 0; i < x.size(); i++) col += post[i * y.size() + j];
        rows = rows && col <= 1 + 1e-12;
      }
      expected = std::max(expected, std::fabs(sum - pairs / total));
      const ExpectedAccuracy e = mea_alignment(post, x.size(), y.size());
      mea = mea && e.expected_pairs + 1e-12 >= pairs_of(v.alignment, post, y.size()) && std::fabs(e.expected_pairs - pairs_of(e.alignment, post, y.size())) < 1e-12;
    }
    s.check("pair HMM: forward and backward give the same log P(x, y)", fb, 1e-12);
    s.check("pair HMM: log P(x, y) is the generic sum over paths in the probability semiring (+, ×)", generic, 1e-12);
    s.check("pair HMM: log P(x, y) is the log of the sum over every enumerated alignment", enumerated, 1e-12);
    s.check("pair HMM: Viterbi is the generic sum in the tropical semiring of log weights", viterbi_gap, 1e-12);
    s.check("pair HMM: the Viterbi path has the log probability it claims", path_gap, 1e-12);
    s.require("pair HMM: posteriors lie in [0, 1]", bounds);
    s.require("pair HMM: posteriors sum to at most 1 along each row and each column", rows);
    s.check("pair HMM: the posteriors sum to the expected number of aligned pairs, by enumeration", expected, 1e-10);
    s.require("pair HMM: the MEA alignment's expected pairs are at least the Viterbi alignment's, and are what it reports", mea);
  }
  {
    // Long homologous sequences: P(x, y) is far below the smallest double, and log space keeps it.
    const Codes root = simulate::random_sequence(1000, 0.5, 5);
    const auto f = simulate::evolve_family(root, 2, 0.05, 0.02, 5);
    const PairTables fw = forward(h, f.sequences[0], f.sequences[1]), bw = backward(h, f.sequences[0], f.sequences[1]);
    const double rel = std::fabs(fw.log_probability - bw.log_probability) / std::fabs(fw.log_probability);
    s.check("pair HMM: on 1 kb homologous sequences, forward = backward in log space (relative)", rel, 1e-12);
    s.require("pair HMM: there P(x, y) underflows double (log P < −745) but its log is finite", std::isfinite(fw.log_probability) && fw.log_probability < -745,
              "log P = " + std::to_string(fw.log_probability));
    const auto post = posterior(h, f.sequences[0], f.sequences[1]);
    const ExpectedAccuracy e = mea_alignment(post, f.sequences[0].size(), f.sequences[1].size());
    s.require("pair HMM: the MEA alignment of related sequences has high expected accuracy (> 0.8)", e.accuracy > 0.8, std::to_string(e.accuracy));
  }
  {
    const auto p = robinson_frequencies();
    const KarlinAltschul k = karlin_altschul(matrix("BLOSUM62"), p, p);
    const Matrix B = matrix("BLOSUM62");
    double total = 0;
    for (int a = 0; a < 24; a++)
      for (int b = 0; b < 24; b++) total += p[static_cast<std::size_t>(a)] * p[static_cast<std::size_t>(b)] * std::exp(k.lambda * B(static_cast<Code>(a), static_cast<Code>(b)));
    const PairHmm protein = PairHmm::from_matrix(B, p, k.lambda, 0.02, 0.4, 0.01);
    double ps = 0;
    for (double v : protein.p) ps += v;
    s.check("pair HMM from BLOSUM62: q_a q_b e^{λ s_ab} at Karlin–Altschul's λ sums to 1 before normalising", total - 1, 1e-12);
    s.check("pair HMM from BLOSUM62: its pair emissions sum to 1", ps - 1, 1e-12);
    const Codes x = sequence::encode_protein("HEAGAWGHEE"), y = sequence::encode_protein("PAWHEAE");
    const PairTables fw = forward(protein, x, y), bw = backward(protein, x, y);
    s.check("pair HMM from BLOSUM62: forward = backward on HEAGAWGHEE and PAWHEAE", fw.log_probability - bw.log_probability, 1e-12);
  }
  return s.done();
}
