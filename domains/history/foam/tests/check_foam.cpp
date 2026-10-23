// Checks of the spin foam of a history at G = ℤ₂ and of the other complexes (theory note §6–7).
//
// On records small enough to know: the refined complex of three ties between two meetings is a
// prism, a sphere; the heat bath agrees with exact enumeration, prior and posterior; the spins
// from the character expansion, P(f ∈ S) = t(⟨h⟩ − t)/(1 − t²), equal those of the foams
// enumerated one by one, and every foam is a ℤ/2 2-cycle of K̂ (the Gauss law); a Wilson loop
// around a circle no cell fills vanishes; the posterior of a sign between two attested ones is
// tanh 2β, and of a sign closing a triangle tanh β_Δ; switching a being changes no cell's mean.
// SU(2) kinematics: dimensions of invariants against Clebsch–Gordan and the Catalan numbers.
// The other complexes: the shells filled kill H₂; Dowker's two complexes have the same
// homology; Atkin's Q-vector of a small relation; persistence in time ends at the Betti numbers
// over ℤ/2. And on every canon: the attested signs never move, and the spins lie in [0, 1].

#include <history/foam/model.hpp>
#include <history/foam/models.hpp>
#include <history/foam/refined.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

using namespace history;
namespace mt = morphe::topology;

namespace {

const std::string HEAD = R"({"kind": "canon", "id": "test", "title": "test", "codebook": 1}
{"kind": "source", "id": "s", "cite": "s"}
)";
std::string being(const std::string& id) { return R"({"kind": "being", "id": ")" + id + R"(", "name": ")" + id + R"(", "type": "person"})" + "\n"; }
std::string occasion(const std::string& id, const std::string& date, const std::string& present) {
  return R"({"kind": "occasion", "id": ")" + id + R"(", "date": ")" + date + R"(", "type": "meeting", "present": )" + present + R"(, "sources": ["s"], "confidence": 1})" + "\n";
}
std::string tie(const std::string& a, const std::string& b, const std::string& attested, const std::string& colour = "letters") {
  return R"({"kind": "tie", "id": ")" + a + "~" + b + R"(", "a": ")" + a + R"(", "b": ")" + b + R"(", "type": "friendship", "colour": ")" + colour + R"(", "attested": )" + attested +
         R"(, "sources": ["s"], "confidence": 1})" + "\n";
}

std::vector<std::size_t> ranks(const mt::CWComplex& X) {
  std::vector<std::size_t> b;
  for (const auto& h : mt::homology(X, false)) b.push_back(h.rank);
  return b;
}
std::string text(const std::vector<std::size_t>& v) {
  std::string s;
  for (std::size_t x : v) s += (s.empty() ? "" : " ") + std::to_string(x);
  return s;
}
double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
  double m = 0;
  for (std::size_t i = 0; i < a.size(); i++) m = std::max(m, std::fabs(a[i] - b[i]));
  return m;
}

bool signed_foam_of(const HistoryModel& M) {
  return std::any_of(M.cells.begin(), M.cells.end(), [](const HistoryModel::Cell& c) { return c.beta < 0; });
}

// A model with every coupling set by hand.
HistoryModel with_couplings(const HistoryComplex& K, const RefinedComplex& R, double lens, double balance) {
  Couplings C;
  C.lens["all"] = lens;
  C.balance = balance;
  return history_model(K, R, C);
}

}  // namespace

int main() {
  calx::verification::Suite s;
  const std::string all = R"([["a", ""], ["b", ""], ["c", ""]])";

  // The prism: three ties between two meetings of three beings.
  {
    const HistoryComplex K = construct(read_record_text(HEAD + being("a") + being("b") + being("c") + occasion("p", "1800", all) + occasion("q", "1801", all) +
                                                        tie("a", "b", R"([["p", 1], ["q", 1]])") + tie("b", "c", R"([["p", 1], ["q", -1]])") + tie("a", "c", R"([["p", -1], ["q", 1]])")));
    const RefinedComplex R = refine(K);
    s.require("the refined prism: 6 presences, 3 + 6 edges, 3 lenses + 2 triangles", R.cw.vertices == 6 && R.cw.edges.size() == 9 && R.temporal == 3 && R.cw.faces.size() == 5 && R.triangles.size() == 2);
    s.require("the refined prism is a sphere (Betti 1 0 1)", ranks(R.cw) == std::vector<std::size_t>{1, 0, 1}, text(ranks(R.cw)));

    const Couplings C = estimate_couplings(K, R);
    s.require("couplings: letters from 3 lenses attested at both ends, 2 flips; balance from 2 triangles, both frustrated", C.estimates[0].of == "letters" && C.estimates[0].n == 3 && C.estimates[0].k == 2 && C.estimates.back().of == "balance" && C.estimates.back().n == 2 && C.estimates.back().k == 2);

    const HistoryModel M = with_couplings(K, R, 0.7, 0.4);
    const Ensemble prior = exact(M, false), post = exact(M, true);
    Sampling S;
    S.sweeps = 60000;
    const Ensemble mc_prior = sample(M, false, S), mc_post = sample(M, true, S);
    s.check("heat bath against enumeration, prior: ⟨h_f⟩", max_difference(mc_prior.cell, prior.cell), 0.02);
    s.check("heat bath against enumeration, posterior: ⟨u_r⟩", max_difference(mc_post.rung, post.rung), 0.02);
    bool fixed = true;
    for (std::size_t r = 0; r < M.rung.size(); r++)
      if (M.attested[r]) fixed = fixed && mc_post.rung[r] == M.attested[r];
    s.require("the posterior never moves an attested sign", fixed);

    std::size_t foams = 0;
    const auto by_foams = foam_spins_by_enumeration(M, &foams);
    const FoamSpace space = foam_space(R);
    s.require("the space of foams of the prism: dimension 1, every cell in its support", space.dimension == 1 && std::all_of(space.support.begin(), space.support.end(), [](bool b) { return b; }));
    s.check("spins: the character expansion's identity against the foams enumerated", max_difference(spins(M, prior), by_foams), 1e-12);
    s.require("the foams of the prism: the empty one and the whole sphere (2 = 2^b₂ over ℤ/2)", foams == 2);
    // Every foam is a ℤ/2 2-cycle of K̂: the sum of all five cells has every edge twice.
    std::vector<int> parity(R.cw.edges.size(), 0);
    for (const auto& f : R.cw.faces)
      for (const auto& step : f.boundary) parity[step.edge] ^= 1;
    s.require("the Gauss law: the whole prism is a 2-cycle over ℤ/2, temporal edges included", std::all_of(parity.begin(), parity.end(), [](int p) { return p == 0; }));
  }

  // A Wilson loop around a circle of encounters that nothing fills.
  {
    const HistoryComplex K = construct(read_record_text(HEAD + being("a") + being("b") + being("c") + occasion("o1", "1801", R"([["a", ""], ["b", ""]])") +
                                                        occasion("o2", "1802", R"([["b", ""], ["c", ""]])") + occasion("o3", "1803", R"([["c", ""], ["a", ""]])") +
                                                        tie("a", "b", R"([["o1", 1]])") + tie("b", "c", R"([["o2", 1]])") + tie("c", "a", R"([["o3", 1]])")));
    const RefinedComplex R = refine(K);
    const HistoryModel M = with_couplings(K, R, 0.9, 0.9);
    Sampling S;
    S.sweeps = 20000;
    const Ensemble E = sample(M, false, S, {{0, 1, 2}});
    s.require("a circle no cell fills: H₁ of K̂ is ℤ", ranks(R.cw) == std::vector<std::size_t>{1, 1}, text(ranks(R.cw)));
    s.check("its Wilson loop vanishes", E.loop[0], 4 * E.loop_error[0] + 0.02);
  }

  // Closed forms of the posterior.
  {
    const std::string ab = R"([["a", ""], ["b", ""]])";
    const HistoryComplex K = construct(read_record_text(HEAD + being("a") + being("b") + occasion("p", "1800", ab) + occasion("q", "1801", ab) + occasion("r", "1802", ab) +
                                                        tie("a", "b", R"([["p", 1], ["r", 1]])")));
    const RefinedComplex R = refine(K);
    const double beta = 0.6;
    const Ensemble E = exact(with_couplings(K, R, beta, 0), true);
    s.check("a sign between two attested +1, coupled by β to each: ⟨u⟩ = tanh 2β", E.rung[1] - std::tanh(2 * beta), 1e-12);
    const HistoryComplex T = construct(read_record_text(HEAD + being("a") + being("b") + being("c") + occasion("p", "1800", all) +
                                                        tie("a", "b", R"([["p", 1]])") + tie("b", "c", R"([["p", 1]])") + tie("a", "c", R"([["p", 1]])")));
    RefinedComplex RT = refine(T);
    HistoryModel MT = with_couplings(T, RT, 0, 0.8);
    MT.attested[2] = 0;  // the third sign unattested
    const Ensemble ET = exact(MT, true);
    s.check("a sign closing a triangle of two attested +1 under balance β_Δ: ⟨u⟩ = tanh β_Δ", ET.rung[2] - std::tanh(0.8), 1e-12);
  }

  // Switching a being changes no cell's mean, prior or posterior (exactly).
  {
    const std::string rec = HEAD + being("a") + being("b") + being("c") + occasion("p", "1800", all) + occasion("q", "1801", all) + occasion("w", "1802", all) +
                            tie("a", "b", R"([["p", 1], ["w", -1]])") + tie("b", "c", R"([["p", 1], ["q", -1], ["w", 1]])") + tie("a", "c", R"([["q", -1], ["w", 1]])");
    Record r = read_record_text(rec), flipped = r;
    for (auto& t : flipped.ties)
      if (t.a == 1 || t.b == 1)
        for (auto& a : t.attested) a.sign = -a.sign;
    const HistoryComplex K = construct(r), KF = construct(flipped);
    const RefinedComplex R = refine(K), RF = refine(KF);
    const Couplings C = estimate_couplings(K, R), CF = estimate_couplings(KF, RF);
    s.require("switching: the couplings estimated are the same", std::fabs(C.balance - CF.balance) < 1e-15 && std::fabs(C.of("letters") - CF.of("letters")) < 1e-15);
    const Ensemble E = exact(history_model(K, R, C), true), EF = exact(history_model(KF, RF, CF), true);
    s.check("switching b: every cell's posterior mean unchanged", max_difference(E.cell, EF.cell), 1e-12);
  }

  // SU(2) kinematics: dim Inv(V_{j₁} ⊗ … ⊗ V_{jₙ}).
  {
    auto half = [](std::size_t n) { return std::vector<int>(n, 1); };
    s.require("dim Inv: (½)⊗2 = 1, (½)⊗4 = 2, (½)⊗6 = 5, (½)⊗8 = 14, the Catalan numbers",
              invariant_dimension(half(2)) == 1 && invariant_dimension(half(4)) == 2 && invariant_dimension(half(6)) == 5 && invariant_dimension(half(8)) == 14);
    s.require("dim Inv: an odd number of spins ½ has none", invariant_dimension(half(3)) == 0 && invariant_dimension(half(5)) == 0);
    s.require("dim Inv: (½, ½, 1) = 1, (1, 1, 1) = 1, (½, 1) = 0, (½, 1, 2) = 0 (the triangle fails)",
              invariant_dimension({1, 1, 2}) == 1 && invariant_dimension({2, 2, 2}) == 1 && invariant_dimension({1, 2}) == 0 && invariant_dimension({1, 2, 4}) == 0);
  }

  // The other complexes, on small records.
  {
    const HistoryComplex K = construct(read_record_text(HEAD + being("a") + being("b") + being("c") + occasion("p", "1800", all) + occasion("q", "1801", all) +
                                                        tie("a", "b", R"([["p", 1], ["q", 1]])") + tie("b", "c", R"([["p", 1], ["q", 1]])") + tie("a", "c", R"([["p", 1], ["q", 1]])")));
    const auto S = shells(K);
    const mt::CWComplex K3 = shells_filled(K, S);
    s.require("K³: the theta-sphere's shell filled is a ball (Betti 1 0 0 0)", ranks(K3) == std::vector<std::size_t>{1, 0, 0, 0}, text(ranks(K3)));
  }
  {
    // Four beings and five occasions, each being at no more than four.
    const HistoryComplex K = construct(read_record_text(HEAD + being("a") + being("b") + being("c") + being("d") +
                                                        occasion("o1", "1801", R"([["a", ""], ["b", ""]])") + occasion("o2", "1802", R"([["b", ""], ["c", ""]])") +
                                                        occasion("o3", "1803", R"([["c", ""], ["d", ""]])") + occasion("o4", "1804", R"([["d", ""], ["a", ""]])") +
                                                        occasion("o5", "1805", R"([["a", ""], ["b", ""], ["c", ""]])")));
    const auto DB = dowker_beings(K), DO = dowker_occasions(K);
    std::vector<std::size_t> hb, ho;
    for (const auto& h : morphe::algebra::homology(DB.chain_complex(), false)) hb.push_back(h.rank);
    for (const auto& h : morphe::algebra::homology(DO.chain_complex(), false)) ho.push_back(h.rank);
    hb.resize(4, 0);
    ho.resize(4, 0);
    s.require("Dowker's theorem: the complexes on the beings and on the occasions have the same homology", hb == ho && hb[1] == 1, text(hb) + " | " + text(ho));
    s.require("Atkin's Q-vector: Q₀ = 1; Q₁ = 3 (o1 and o2 1-near through o5, o3 and o4 alone); Q₂ = 1", text(q_analysis(K)) == "1 3 1", text(q_analysis(K)));
  }

  // Every canon.
  for (const auto& entry : std::filesystem::directory_iterator(HISTORY_CANONS)) {
    if (entry.path().extension() != ".jsonl") continue;
    const std::string name = entry.path().filename().string();
    const HistoryComplex K = construct(read_record_file(entry.path().string()));
    const RefinedComplex R = refine(K);
    const HistoryModel M = history_model(K, R, estimate_couplings(K, R));
    Sampling S;
    S.burn = 200;
    S.sweeps = 20000;
    const Ensemble post = sample(M, true, S), prior = sample(M, false, S);
    const FoamSpace space = foam_space(R);
    const auto b2 = morphe::algebra::betti(R.cw.chain_complex(), 2);
    s.require("canon " + name + ": the foams are the 2-cycles of K̂ over ℤ/2 (their dimension is b₂)", b2.size() > 2 && space.dimension == b2[2], std::to_string(space.dimension) + " foams");
    bool fixed = true;
    for (std::size_t r = 0; r < M.rung.size(); r++)
      if (M.attested[r]) fixed = fixed && post.rung[r] == M.attested[r];
    s.require("canon " + name + ": the posterior keeps every attested sign", fixed, std::to_string(M.rung.size()) + " rungs, " + std::to_string(M.cells.size()) + " cells");
    // Griffiths–Kelly–Sherman: with every coupling positive, ⟨h_f⟩ ≥ tanh β_f for every cell; a
    // mean below it would betray a sampler that does not mix.
    if (!signed_foam_of(M)) {
      double worst = 0;
      for (std::size_t c = 0; c < M.cells.size(); c++) worst = std::min(worst, prior.cell[c] - std::tanh(M.cells[c].beta));
      s.check("canon " + name + ": the GKS inequality ⟨h_f⟩ ≥ tanh β_f on every cell, within the sampling error", std::min(0.0, worst + 0.03), 0.0);
    }
    const auto j = spins(M, prior, &space);
    const double lo = *std::min_element(j.begin(), j.end()), hi = *std::max_element(j.begin(), j.end());
    const bool signed_foam = std::any_of(M.cells.begin(), M.cells.end(), [](const HistoryModel::Cell& c) { return c.beta < 0; });
    s.require("canon " + name + ": the spins of the prior lie in [0, 1] within the sampling error (or the foam is signed)", signed_foam || (lo > -0.05 && hi < 1.05), std::to_string(lo) + " .. " + std::to_string(hi));
    // Corollary 6.2: with every coupling positive, no cell exceeds tanh β_f / (1 + tanh β_f).
    if (!signed_foam) {
      double over = 0;
      for (std::size_t c = 0; c < M.cells.size(); c++) {
        const double t = std::tanh(M.cells[c].beta);
        over = std::max(over, j[c] - t / (1 + t));
      }
      s.check("canon " + name + ": no spin exceeds t_f/(1 + t_f) (Corollary 6.2), within the sampling error", std::max(0.0, over - 0.02), 0.0);
    }
    const auto D = persistence_in_time(K);
    const auto b = mt::betti(K.cw, 2);
    bool ends = true;
    for (std::size_t d = 0; d < b.size(); d++) ends = ends && D.betti(d, K.order.size()) == b[d];
    s.require("canon " + name + ": persistence in time ends at the Betti numbers over ℤ/2", ends, text(b));
  }
  return s.done();
}
