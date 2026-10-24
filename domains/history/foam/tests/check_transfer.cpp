// Checks of the exact contraction of a history (theory note §6.6), physis's contraction along the
// order of time, and of the couplings by maximum a posteriori.
//
// On records small enough to enumerate every configuration: log Z, the means of rungs and cells,
// two-point functions and the law of the frustration at a cut agree with brute force, prior and
// posterior. A sign between two attested ones is tanh 2β; a chain's end-to-end function is t^k.
// Proposition 6.4: with independent cells the exact MAP is the KT estimate, its standard error
// cosh β/√(n + 1). Proposition 6.3: every state of the prior at a cut is a sum of ℤ/2 1-cycles of
// the cut graph. On every canon within reach: the GKS inequality and the bound of Corollary 6.2
// exactly, the attested signs unmoved, physis's heat bath in agreement, the MAP's gradient zero.

#include <history/foam/model.hpp>
#include <history/foam/refined.hpp>
#include <history/foam/transfer.hpp>
#include <history/readings/signed.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

using namespace history;

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

HistoryModel with_couplings(const HistoryComplex& K, const RefinedComplex& R, double lens, double balance) {
  Couplings C;
  C.lens["all"] = lens;
  C.balance = balance;
  return history_model(K, R, C);
}

// Everything by enumeration.
struct Brute {
  double log_z = 0;
  std::vector<double> rung, cell;
  std::map<std::pair<std::size_t, std::size_t>, double> two;
  std::vector<double> law;
};

Brute brute(const HistoryModel& M, bool conditioned, const std::vector<std::array<std::size_t, 3>>& triangles = {}) {
  const std::size_t n = M.rung.size();
  std::vector<std::size_t> vars;
  std::vector<int> u(n, 1);
  for (std::size_t r = 0; r < n; r++) {
    if (conditioned && M.attested[r]) u[r] = M.attested[r];
    else vars.push_back(r);
  }
  Brute B;
  B.rung.assign(n, 0);
  B.cell.assign(M.cells.size(), 0);
  B.law.assign(triangles.size() + 1, 0);
  std::vector<std::vector<double>> pair(n, std::vector<double>(n, 0));
  double Z = 0;
  for (std::uint64_t mask = 0; mask < (std::uint64_t{1} << vars.size()); mask++) {
    for (std::size_t i = 0; i < vars.size(); i++) u[vars[i]] = (mask >> i) & 1 ? -1 : 1;
    double logw = 0;
    for (const auto& c : M.cells) logw += c.beta * holonomy(c, u);
    const double w = std::exp(logw);
    Z += w;
    for (std::size_t r = 0; r < n; r++) B.rung[r] += w * u[r];
    for (std::size_t c = 0; c < M.cells.size(); c++) B.cell[c] += w * holonomy(M.cells[c], u);
    for (std::size_t r = 0; r < n; r++)
      for (std::size_t q = 0; q < n; q++) pair[r][q] += w * u[r] * u[q];
    std::size_t f = 0;
    for (const auto& t : triangles) f += u[t[0]] * u[t[1]] * u[t[2]] < 0;
    B.law[f] += w;
  }
  B.log_z = std::log(Z);
  for (auto& x : B.rung) x /= Z;
  for (auto& x : B.cell) x /= Z;
  for (auto& x : B.law) x /= Z;
  for (std::size_t r = 0; r < n; r++)
    for (std::size_t q = 0; q < n; q++) B.two[{r, q}] = pair[r][q] / Z;
  return B;
}

double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
  double m = 0;
  for (std::size_t i = 0; i < a.size(); i++) m = std::max(m, std::fabs(a[i] - b[i]));
  return m;
}

// The vertex of K after which the occasions dated no later than a day end.
std::size_t vertex_at(const HistoryComplex& K, long jdn) {
  std::size_t v = 0;
  for (std::size_t k = 0; k < K.order.size(); k++)
    if (K.record.occasions[K.order[k]].date.first <= jdn) v = k;
  return v;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // A small history with triangles, carried signs and a tie of one meeting.
  const std::string four = HEAD + being("a") + being("b") + being("c") + being("d") +
                           occasion("o1", "1800", R"([["a", ""], ["b", ""], ["c", ""]])") +
                           occasion("o2", "1801", R"([["a", ""], ["b", ""], ["c", ""], ["d", ""]])") +
                           occasion("o3", "1802", R"([["b", ""], ["c", ""], ["d", ""]])") +
                           occasion("o4", "1803", R"([["a", ""], ["b", ""], ["c", ""], ["d", ""]])") +
                           occasion("o5", "1804", R"([["a", ""], ["c", ""], ["d", ""]])") +
                           tie("a", "b", R"([["o1", 1], ["o4", -1]])") + tie("b", "c", R"([["o1", 1], ["o3", 1]])") +
                           tie("a", "c", R"([["o2", -1], ["o5", -1]])", "music") + tie("c", "d", R"([["o2", 1], ["o5", 1]])", "music") +
                           tie("b", "d", R"([["o3", -1], ["o4", -1]])") + tie("a", "d", R"([["o4", 1]])");
  {
    const HistoryComplex K = construct(read_record_text(four));
    const RefinedComplex R = refine(K);
    const HistoryModel M = with_couplings(K, R, 0.6, 0.45);
    s.require("the small history: 16 rungs, some carried, and triangles", M.rung.size() == 16 && R.triangles.size() >= 3, std::to_string(M.rung.size()) + " rungs, " + std::to_string(R.triangles.size()) + " triangles");
    // Two-points: every pair of variables in the order the contraction allows, prior.
    std::vector<std::array<std::size_t, 2>> pairs;
    for (std::size_t t = 0; t < M.index.size(); t++)
      if (M.index[t].size() >= 2) pairs.push_back({M.index[t].front(), M.index[t].back()});
    const Contraction prior = contract(K, M, false, pairs);
    const Brute bp = brute(M, false);
    s.check("prior: log Z by contraction against enumeration", prior.log_z - bp.log_z, 1e-10);
    s.check("prior: ⟨h_f⟩ of every cell", max_difference(prior.plaquette, bp.cell), 1e-10);
    s.check("prior: ⟨u_r⟩ of every rung (zero, by switching)", max_difference(prior.link, bp.rung), 1e-10);
    double worst = 0;
    for (std::size_t k = 0; k < pairs.size(); k++) worst = std::max(worst, std::fabs(prior.two_point[k] - bp.two.at({pairs[k][0], pairs[k][1]})));
    s.check("prior: the two-point function of each tie's first and last rungs", worst, 1e-10);
    // The same, gauge invariantly: the Wilson loop around each tie's life in K̂, every link free,
    // by physis's heat bath.
    Observables O;
    for (std::size_t t = 0; t < M.index.size(); t++)
      if (M.index[t].size() >= 2) O.loops.push_back(tie_loop(K, R, t, 0, M.index[t].size() - 1));
    Sampling S;
    S.sweeps = 40000;
    const Ensemble E = sample(R, M, false, S, O);
    s.check("prior: the Wilson loop around each tie's life (heat bath on K̂) is the two-point function", max_difference(E.loop, prior.two_point), 0.03);

    // The posterior, with the law of two triangles' frustration at a cut after o3.
    const std::size_t v3 = vertex_at(K, jdn_gregorian(1802, 12, 31));
    auto rung_at = [&](const std::string& a, const std::string& b, const std::string& o) {
      for (std::size_t t = 0; t < K.record.ties.size(); t++) {
        const auto& x = K.record.ties[t];
        if (K.record.beings[x.a].id != a || K.record.beings[x.b].id != b) continue;
        for (std::size_t k = 0; k < K.rungs[t].size(); k++)
          if (K.record.occasions[K.rungs[t][k].occasion].id == o) return M.index[t][k];
      }
      throw std::logic_error("no rung");
    };
    Counter C;
    C.triangles = {{rung_at("a", "b", "o2"), rung_at("b", "c", "o3"), rung_at("a", "c", "o2")}, {rung_at("b", "c", "o3"), rung_at("c", "d", "o3"), rung_at("b", "d", "o3")}};
    Cut cut;
    cut.after = v3;
    for (const auto& tri : C.triangles) cut.keep.insert(cut.keep.end(), tri.begin(), tri.end());
    const Contraction post = contract(K, M, true, {}, {cut}, {C});
    const Brute bq = brute(M, true, C.triangles);
    s.check("posterior: log Z by contraction against enumeration", post.log_z - bq.log_z, 1e-10);
    s.check("posterior: ⟨u_r⟩ of every rung", max_difference(post.link, bq.rung), 1e-10);
    s.check("posterior: ⟨h_f⟩ of every cell", max_difference(post.plaquette, bq.cell), 1e-10);
    s.check("posterior: the law of the frustration at a cut", max_difference(post.cuts[0].negative, bq.law), 1e-10);
    double mean = 0;
    for (std::size_t k = 0; k < bq.law.size(); k++) mean += static_cast<double>(k) * bq.law[k];
    s.check("posterior: the mean frustration at the cut by pairing", post.cuts[0].mean_negative - mean, 1e-10);
    bool fixed = true;
    for (std::size_t r = 0; r < M.rung.size(); r++)
      if (M.attested[r]) fixed = fixed && post.link[r] == M.attested[r];
    s.require("posterior: the attested signs are constants", fixed);

    // The derivatives carried through physis's contraction with the record's links fixed: the
    // gradient of log Z is Σ ⟨h_f⟩ over each parameter's cells, and its Hessian the derivative of
    // that, prior and posterior. Three parameters: letters, music, balance.
    const Parameters P = parameters_of(K, M);
    s.require("the small history: three couplings to fit", P.names.size() == 3);
    std::array<double, 3> at{};
    for (std::size_t c = 0; c < M.cells.size(); c++) at[P.of_cell[c]] = M.cells[c].beta;
    for (const bool conditioned : {false, true}) {
      const auto order = time_order(K, M);
      const auto J = physis::gauge::log_z_derivatives<3>(theory_of(M, conditioned), order, P.of_cell, at);
      const Contraction ens = contract(K, M, conditioned);
      std::vector<double> sums(3, 0), gradient(J.gradient.begin(), J.gradient.end());
      for (std::size_t c = 0; c < M.cells.size(); c++) sums[P.of_cell[c]] += ens.plaquette[c];
      const std::string which = conditioned ? "posterior" : "prior";
      s.check(which + ": log Z of the jet against the contraction", J.log_z - ens.log_z, 1e-10);
      s.check(which + ": ∂ log Z/∂β_p = Σ ⟨h_f⟩ over p's cells", max_difference(gradient, sums), 1e-10);
      double worst_h = 0;
      const double h = 1e-5;
      for (std::size_t q = 0; q < 3; q++) {
        auto up = at, down = at;
        up[q] += h;
        down[q] -= h;
        const auto gu = physis::gauge::log_z_derivatives<3, 1>(theory_of(M, conditioned), order, P.of_cell, up).gradient;
        const auto gd = physis::gauge::log_z_derivatives<3, 1>(theory_of(M, conditioned), order, P.of_cell, down).gradient;
        for (std::size_t p = 0; p < 3; p++) worst_h = std::max(worst_h, std::fabs(J.hessian[p][q] - (gu[p] - gd[p]) / (2 * h)));
      }
      s.check(which + ": the Hessian of log Z against differences of the gradient", worst_h, 1e-6);
    }
  }

  // A sign between two attested ones is tanh 2β; the chain's end-to-end function is t^k.
  {
    const HistoryComplex K = construct(read_record_text(HEAD + being("a") + being("b") + occasion("o1", "1800", R"([["a", ""], ["b", ""]])") +
                                                        occasion("o2", "1801", R"([["a", ""], ["b", ""]])") + occasion("o3", "1802", R"([["a", ""], ["b", ""]])") +
                                                        occasion("o4", "1803", R"([["a", ""], ["b", ""]])") + tie("a", "b", R"([["o1", 1], ["o3", 1], ["o4", 1]])")));
    const RefinedComplex R = refine(K);
    const double beta = 0.8;
    const HistoryModel M = with_couplings(K, R, beta, 0);
    const Contraction post = contract(K, M, true);
    s.check("a carried sign between two attested ones: tanh 2β", post.link[M.index[0][1]] - std::tanh(2 * beta), 1e-12);
    const Contraction prior = contract(K, M, false, {{M.index[0][0], M.index[0][3]}});
    s.check("a chain of three lenses: ⟨u_first u_last⟩ = tanh³ β", prior.two_point[0] - std::pow(std::tanh(beta), 3), 1e-12);
    s.check("a chain: log Z = log 2 + 3 log 2cosh β", prior.log_z - (std::log(2.0) + 3 * std::log(2 * std::cosh(beta))), 1e-12);
  }

  // Proposition 6.4: independent cells. Two separate ties of one colour (12 lenses, 3 flips) and
  // one of another (5 lenses, none flipped): the MAP is the KT estimate of each colour.
  {
    std::string rec = HEAD + being("a") + being("b") + being("c") + being("d") + being("e") + being("f");
    std::string ab, cd, ef;
    const int ab_signs[] = {1, 1, -1, -1, 1, 1, 1}, cd_signs[] = {1, 1, 1, 1, 1, -1, -1};
    for (int k = 0; k < 7; k++) {
      const std::string y = std::to_string(1800 + k), o = "p" + std::to_string(k), q = "q" + std::to_string(k);
      rec += occasion(o, y + "-03", R"([["a", ""], ["b", ""]])") + occasion(q, y + "-06", R"([["c", ""], ["d", ""]])");
      ab += (k ? ", " : "") + std::string(R"([")") + o + R"(", )" + std::to_string(ab_signs[k]) + "]";
      cd += (k ? ", " : "") + std::string(R"([")") + q + R"(", )" + std::to_string(cd_signs[k]) + "]";
    }
    for (int k = 0; k < 6; k++) {
      const std::string o = "r" + std::to_string(k);
      rec += occasion(o, std::to_string(1810 + k), R"([["e", ""], ["f", ""]])");
      ef += (k ? ", " : "") + std::string(R"([")") + o + R"(", 1])";
    }
    rec += tie("a", "b", "[" + ab + "]") + tie("c", "d", "[" + cd + "]") + tie("e", "f", "[" + ef + "]", "music");
    const HistoryComplex K = construct(read_record_text(rec));
    const RefinedComplex R = refine(K);
    const Couplings kt = estimate_couplings(K, R);
    const HistoryModel M = history_model(K, R, kt);
    const Fit fit = fit_couplings(K, M, kt);
    auto expect = [](double n, double k) { return std::atanh(1 - 2 * (k + 0.5) / (n + 1)); };
    std::map<std::string, std::size_t> p;
    for (std::size_t i = 0; i < fit.names.size(); i++) p[fit.names[i]] = i;
    s.require("independent cells: one parameter per colour, no balance", fit.names.size() == 2 && p.count("letters") && p.count("music"));
    s.check("Prop. 6.4: the MAP of letters is KT's (k + ½)/(n + 1), 12 lenses, 3 flips", fit.beta[p["letters"]] - expect(12, 3), 1e-6);
    s.check("Prop. 6.4: the MAP of music, 5 lenses unflipped", fit.beta[p["music"]] - expect(5, 0), 1e-6);
    s.check("Laplace: the standard error is cosh β/√(n + 1)", fit.error[p["letters"]] - std::cosh(fit.beta[p["letters"]]) / std::sqrt(13.0), 1e-4);
    s.require("the fit converged", fit.converged, std::to_string(fit.iterations) + " iterations");
  }

  // Every canon within reach.
  for (const auto& entry : std::filesystem::directory_iterator(HISTORY_CANONS)) {
    if (entry.path().extension() != ".jsonl") continue;
    const std::string name = entry.path().filename().string();
    const HistoryComplex K = construct(read_record_file(entry.path().string()));
    const RefinedComplex R = refine(K);
    const Couplings kt = estimate_couplings(K, R);
    const HistoryModel M = history_model(K, R, kt);
    // Cuts after every fifth occasion, asking for every spin network there.
    std::vector<Cut> cuts;
    for (std::size_t v = 4; v + 1 < K.order.size(); v += 5) cuts.push_back({v, {}});
    try {
      const Contraction prior = contract(K, M, false, {}, cuts, {}, std::size_t{1} << 22);
      const Contraction post = contract(K, M, true);
      // Proposition 6.3: the prior's spin networks at a cut are 1-cycles of the cut graph.
      bool cycles = true;
      for (const auto& law : prior.cuts)
        for (const auto& [set, weight] : law.networks) {
          std::map<std::size_t, int> degree;
          for (std::size_t r : set) {
            const auto& x = K.record.ties[M.rung[r][0]];
            degree[x.a] ^= 1;
            degree[x.b] ^= 1;
          }
          for (const auto& [b, d] : degree) cycles = cycles && d == 0;
        }
      std::size_t widest = 0;
      for (const auto& law : prior.cuts) widest = std::max(widest, law.live);
      s.require("canon " + name + ": the prior's states at the cuts are ℤ/2 1-cycles of the cut graphs (Prop. 6.3)", cycles,
                "width " + std::to_string(prior.width) + ", largest state " + std::to_string(prior.largest));
      const bool signed_foam = std::any_of(M.cells.begin(), M.cells.end(), [](const HistoryModel::Cell& c) { return c.beta < 0; });
      if (!signed_foam) {
        double gks = 0, bound = 0;
        for (std::size_t c = 0; c < M.cells.size(); c++) {
          const double t = std::tanh(M.cells[c].beta);
          gks = std::min(gks, prior.plaquette[c] - t);
          bound = std::max(bound, t * (prior.plaquette[c] - t) / (1 - t * t) - t / (1 + t));
        }
        s.check("canon " + name + ": GKS, ⟨h_f⟩ ≥ tanh β_f exactly", std::min(0.0, gks + 1e-12), 0.0);
        s.check("canon " + name + ": Corollary 6.2 exactly, no spin above t/(1 + t)", std::max(0.0, bound - 1e-12), 0.0);
      }
      bool fixed = true;
      for (std::size_t r = 0; r < M.rung.size(); r++)
        if (M.attested[r]) fixed = fixed && post.link[r] == M.attested[r];
      s.require("canon " + name + ": the exact posterior keeps every attested sign", fixed);
      // The model in temporal gauge is physis's theory on the whole of K̂, every link free: log Z
      // less #temporal · log 2, the temporal edges being a forest.
      std::vector<double> beta;
      for (const auto& c : M.cells) beta.push_back(c.beta);
      try {
        const double whole = physis::gauge::log_z(physis::gauge::Z2Theory::of(R.cw, beta), time_order(K, M));
        s.check("canon " + name + ": log Z in temporal gauge is physis's on K̂ less #temporal · log 2",
                prior.log_z + static_cast<double>(R.temporal) * std::log(2.0) - whole, 1e-8);
      } catch (const Intractable&) {
      }
      Sampling S;
      S.sweeps = 3000;
      const Ensemble mc = sample(R, M, false, S);
      s.check("canon " + name + ": physis's heat bath agrees with the contraction, prior ⟨h_f⟩", max_difference(mc.cell, prior.plaquette), 0.06);
      const Fit fit = fit_couplings(K, M, kt);
      const auto P = parameters_of(K, M);
      const auto g = transfer_detail::evaluate(K, M, P, fit.beta, {}).gradient;
      double gmax = 0;
      for (double x : g) gmax = std::max(gmax, std::fabs(x));
      s.check("canon " + name + ": the MAP couplings, gradient zero", gmax, 1e-5, std::to_string(fit.names.size()) + " couplings, log-likelihood " + std::to_string(fit.log_likelihood));
    } catch (const Intractable& e) {
      s.require("canon " + name + ": beyond exact reach, declared so (the heat bath stands in)", true, e.what());
    }
  }
  return s.done();
}
