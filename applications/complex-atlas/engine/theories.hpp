// The Complex Atlas's engine: field theories that generate spacetime (physis decision 0006),
// from their actions through their Feynman expansions to their large-N solutions.
//
// Theories in invariant form (scalar, vector, matrix, tensor and tensorial group field
// theories) are physis field_theory's; theories in coloured form (the simplicial tensor model,
// Boulatov's and Ooguri's group field theories) are enumerated here over their coloured
// graphs, whose isomorphism classes come from a canonical form: a connected coloured graph is
// rigid once a vertex is fixed, so relabelling breadth first from each black vertex and keeping
// the least result names its class (O(p² (D+1)) per graph).

#pragma once

#include <algorithm>
#include <climits>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include <physis/field_theory/invariant.hpp>
#include <physis/field_theory/large_n.hpp>
#include <physis/quantum_gravity/gft.hpp>
#include <physis/quantum_gravity/tensor_models.hpp>

#include "generation.hpp"

namespace atlas {

namespace pf = physis::field_theory;
namespace pq = physis::quantum_gravity;
using morphe::arithmetic::Integer;
using morphe::arithmetic::Rational;

struct TheoryInfo {
  const char* id;
  const char* title;
  const char* family;
  const char* description;
  const char* action;
};

inline const std::vector<TheoryInfo>& theory_infos() {
  static const std::vector<TheoryInfo> infos = {
      {"scalar", "Complex scalar", "QFT", "A complex scalar with a (φ̄φ)^k interaction, zero-dimensional: its Feynman graphs are plain graphs, matchings of φ with φ̄.",
       "S[φ] = φ̄φ + λ (φ̄φ)^k"},
      {"vector", "Vector model", "QFT", "A U(N)-invariant vector with (φ̄·φ)^k: each closed index line gives N; at large N the cactus graphs dominate and a gap equation sums them.",
       "S[φ] = φ̄·φ + λ (φ̄·φ)^k"},
      {"matrix", "Matrix model", "QFT", "A complex N × N matrix with tr (MM†)^k: its Feynman graphs are ribbon graphs, surfaces of genus g weighted N^{2−2g}; the planar ones dominate, and at the critical coupling they become continuous 2D quantum gravity.",
       "S[M] = tr MM† + λ tr (MM†)^k"},
      {"tensor-melonic", "Quartic melonic tensor model", "QFT", "A tensor of rank d with the quartic melonic bubble of every colour: Feynman graphs are (d+1)-coloured graphs, pseudomanifolds weighted by their Gurau degree; the melonic ones dominate and resum into branched polymers.",
       "S[φ] = φ̄·φ + λ N^{1−d} Σ_c I_{melon,c}(φ, φ̄)"},
      {"tensor-wheel", "Sextic wheel tensor model", "QFT", "A rank-3 tensor with the wheel K₃,₃, a sextic bubble that is not melonic.", "S[φ] = φ̄·φ + λ I_{K₃,₃}(φ, φ̄)"},
      {"tensor-necklace", "Quartic necklace tensor model", "QFT", "A rank-3 tensor with the quartic necklace.", "S[φ] = φ̄·φ + λ I_{necklace}(φ, φ̄)"},
      {"coloured-simplicial", "Coloured simplicial tensor model", "QFT", "D + 1 tensors whose interaction is a D-simplex: every Feynman graph is a (D+1)-coloured graph and triangulates a D-dimensional pseudomanifold (Gurau 2011).",
       "S[φ] = Σ_c φ̄_c·φ_c + λ φ_0 ⋯ φ_D + λ̄ φ̄_0 ⋯ φ̄_D"},
      {"gft-boulatov", "Boulatov's group field theory", "GFT", "Fields on G³ with the gauge projection and the tetrahedral interaction: Feynman graphs are 4-coloured, their 2-complexes spin foams of three-dimensional gravity (Boulatov 1992); amplitudes here for a finite group.",
       "S[φ] = Σ_c ∫ φ̄_c P φ_c + λ ∫ φ_0 φ_1 φ_2 φ_3 + c.c.,  φ_c : G³ → ℂ"},
      {"gft-ooguri", "Ooguri's group field theory", "GFT", "Fields on G⁴ with the 4-simplex interaction: Feynman graphs are 5-coloured, their 2-complexes the spin foams of four-dimensional BF theory (Ooguri 1992); amplitudes for a finite group.",
       "S[φ] = Σ_c ∫ φ̄_c P φ_c + λ ∫ φ_0 ⋯ φ_4 + c.c.,  φ_c : G⁴ → ℂ"},
      {"gft-melonic", "Quartic melonic group field theory", "GFT", "A tensorial group field theory: fields on G^d with the gauge projection and quartic melonic bubbles; its condensed phase is the setting of GFT cosmology.",
       "S[φ] = ∫ φ̄ P φ + λ Σ_c ∫ I_{melon,c}(φ, φ̄),  φ : G^d → ℂ"},
  };
  return infos;
}

inline bool coloured_theory(const std::string& id) { return id == "coloured-simplicial" || id == "gft-boulatov" || id == "gft-ooguri"; }
inline bool gft_theory(const std::string& id) { return id.rfind("gft-", 0) == 0; }

// The invariant theory behind an id (not for the coloured ones).
inline pf::InvariantTheory invariant_theory(const std::string& id, const Json& p) {
  const auto k = static_cast<std::size_t>(integer(p, "k", 2, 1, 4));
  const auto rank = static_cast<std::size_t>(integer(p, "rank", 3, 3, 5));
  if (id == "scalar") return pf::theory_by_name("scalar", 0, k);
  if (id == "vector") return pf::theory_by_name("vector", 1, k);
  if (id == "matrix") return pf::theory_by_name("matrix", 2, k);
  if (id == "tensor-melonic" || id == "gft-melonic") return pf::theory_by_name(p["colours"].text("all") == "one" ? "melonic-1" : "melonic", rank);
  if (id == "tensor-wheel") return pf::theory_by_name("wheel", 3);
  if (id == "tensor-necklace") return pf::theory_by_name("necklace", 3);
  throw std::invalid_argument("no theory \"" + id + "\"");
}

inline std::size_t coloured_dimension(const std::string& id, const Json& p) {
  if (id == "gft-boulatov") return 3;
  if (id == "gft-ooguri") return 4;
  return static_cast<std::size_t>(integer(p, "D", 3, 2, 4));
}

inline morphe::algebra::FiniteGroup gauge_group(const Json& p) {
  const std::string g = p["group"].text("Z2");
  if (g == "trivial") return morphe::algebra::cyclic(1);
  if (g == "Z2") return morphe::algebra::cyclic(2);
  if (g == "Z3") return morphe::algebra::cyclic(3);
  if (g == "S3") return morphe::algebra::symmetric(3);
  if (g == "Q8") return morphe::algebra::quaternion();
  throw std::invalid_argument("the group is trivial, Z2, Z3, S3 or Q8");
}

inline Json params_of(const std::string& id) {
  Json ps = Json::array();
  if (id == "scalar" || id == "vector" || id == "matrix") ps.push(param("k", "int", 1, 4, 2));
  if (id == "tensor-melonic" || id == "gft-melonic") ps.push(param("rank", "int", 3, 5, 3));
  if (id == "coloured-simplicial") ps.push(param("D", "int", 2, 4, 3));
  if (gft_theory(id)) ps.push(param("group", "text", 0, 0, "Z2"));
  return ps;
}

inline Json theories() {
  Json out = Json::array();
  for (const auto& t : theory_infos()) {
    const std::string id = t.id;
    Json bubbles = Json::array();
    std::size_t rank = 0;
    if (coloured_theory(id)) {
      // The simplicial interaction has no bubble: every colour is a propagator, each vertex a
      // D-simplex (the page draws it).
      rank = coloured_dimension(id, Json::object());
    } else {
      const auto T = invariant_theory(id, Json::object());
      rank = T.rank;
      for (const auto& I : T.interactions) bubbles.push(Json::object({{"name", I.name}, {"coloured_graph", coloured_graph_json(I.bubble)}, {"scaling", static_cast<double>(I.scaling)}}));
    }
    out.push(Json::object({{"id", t.id}, {"title", t.title}, {"family", t.family}, {"rank", static_cast<double>(rank)}, {"description", t.description},
                           {"action", t.action}, {"params", params_of(id)}, {"bubbles", bubbles}}));
  }
  return Json::object({{"theories", out}});
}

// ---- the expansion ----------------------------------------------------------------------------------------------------

// (−1)^n P: a coefficient of (−λ)^n as one of λ^n.
inline pf::Polynomial signed_order(pf::Polynomial P, std::size_t n) {
  if (n % 2)
    for (auto& [f, c] : P) c = -c;
  return P;
}
inline pf::Polynomial times(const pf::Polynomial& a, const pf::Polynomial& b) {
  pf::Polynomial c;
  for (const auto& [i, x] : a)
    for (const auto& [j, y] : b) c[i + j] += x * y;
  return c;
}
// Z = exp F as series in one variable: n z_n = Σ_{k=1}^{n} k w_k z_{n−k}, z_0 = 1.
inline std::vector<pf::Polynomial> exponentiate(const std::vector<pf::Polynomial>& w) {
  std::vector<pf::Polynomial> z(w.size());
  z[0] = {{0, Rational(1)}};
  for (std::size_t n = 1; n < w.size(); n++) {
    pf::Polynomial acc;
    for (std::size_t k = 1; k <= n; k++)
      for (const auto& [f, c] : times(w[k], z[n - k])) acc[f] += c * Rational(static_cast<long long>(k));
    for (auto& [f, c] : acc) c = c / Rational(static_cast<long long>(n));
    z[n] = acc;
  }
  return z;
}

inline Json polynomial_json(const pf::Polynomial& P) {
  Json terms = Json::array();
  for (auto it = P.rbegin(); it != P.rend(); ++it)
    if (!it->second.is_zero()) terms.push(Json::object({{"power", static_cast<double>(it->first)}, {"coeff", it->second.str()}}));
  return terms;
}

// A canonical form of a connected coloured graph: the least relabelling, breadth first from a
// black root, over all roots.
inline std::vector<std::vector<std::size_t>> canonical_sigma(const mt::ColouredGraph& G) {
  const std::size_t p = G.p, C = G.colours;
  std::vector<std::vector<std::size_t>> inv(C, std::vector<std::size_t>(p));
  for (std::size_t c = 0; c < C; c++)
    for (std::size_t b = 0; b < p; b++) inv[c][G.sigma[c][b]] = b;
  std::vector<std::vector<std::size_t>> best;
  for (std::size_t root = 0; root < p; root++) {
    std::vector<std::size_t> nb(p, SIZE_MAX), nw(p, SIZE_MAX), order;
    std::size_t blacks = 0, whites = 0;
    nb[root] = blacks++;
    order.push_back(root);
    for (std::size_t q = 0; q < order.size(); q++) {
      const std::size_t b = order[q];
      for (std::size_t c = 0; c < C; c++) {
        const std::size_t w = G.sigma[c][b];
        if (nw[w] == SIZE_MAX) {
          nw[w] = whites++;
          for (std::size_t d = 0; d < C; d++) {
            const std::size_t b2 = inv[d][w];
            if (nb[b2] == SIZE_MAX) nb[b2] = blacks++, order.push_back(b2);
          }
        }
      }
    }
    if (blacks != p) return {};  // not connected
    std::vector<std::vector<std::size_t>> s(C, std::vector<std::size_t>(p));
    for (std::size_t c = 0; c < C; c++)
      for (std::size_t b = 0; b < p; b++) s[c][nb[b]] = nw[G.sigma[c][b]];
    if (best.empty() || s < best) best = s;
  }
  return best;
}

inline Json expand(const Json& request) {
  const std::string id = request["theory"].text();
  const Json& p = request["params"];
  const auto order = static_cast<std::size_t>(integer(request, "order", 2, 0, 6));
  const auto limit = static_cast<std::size_t>(integer(request, "limit", 200, 1, 2000));
  Json graphs = Json::array(), z = Json::array(), fe = Json::array();
  std::size_t total = 0;
  bool truncated = false;
  long long lead_power = LLONG_MIN;
  std::string lead_coeff;
  if (coloured_theory(id)) {
    const std::size_t D = coloured_dimension(id, p);
    if ((D == 3 && order > 4) || (D == 4 && order > 3) || (D == 2 && order > 5)) throw std::invalid_argument("the coloured expansion is enumerated up to order 5, 4 or 3 for D = 2, 3, 4");
    const bool gft = gft_theory(id);
    const auto group = gft ? gauge_group(p) : morphe::algebra::cyclic(1);
    std::vector<pf::Polynomial> Fs(order + 1);
    // Free energy at every order up to the requested one: Σ over connected classes of N^F / |Aut|.
    for (std::size_t n = 1; n <= order; n++) {
      std::map<std::vector<std::vector<std::size_t>>, std::size_t> classes;
      std::map<std::vector<std::vector<std::size_t>>, mt::ColouredGraph> rep;
      pq::for_each_coloured_graph(D, n, [&](const mt::ColouredGraph& Gm) {
        auto key = canonical_sigma(Gm);
        if (key.empty()) return;
        if (!classes[key]++) rep.emplace(key, Gm);
      });
      pf::Polynomial F;
      const Integer nf = morphe::arithmetic::factorial(n);
      for (const auto& [key, labelled] : classes) {
        const auto& Gm = rep.at(key);
        const long long faces = static_cast<long long>(Gm.face_count());
        F[faces] += Rational(Integer(static_cast<long long>(labelled)), nf);
        if (n != order) continue;
        total += labelled;
        if (graphs.size() >= limit) { truncated = true; continue; }
        Json g = Json::object({{"id", static_cast<double>(total - 1)}, {"coloured_graph", coloured_graph_json(Gm)}, {"counts", Json::array({static_cast<double>(n)})},
                               {"connected", true}, {"n_power", static_cast<double>(faces)}, {"n_power_scaled", static_cast<double>(pq::amplitude_exponent(Gm))},
                               {"degree", static_cast<double>(Gm.degree())}, {"genus", Json()},
                               {"automorphisms", morphe::arithmetic::exact_div(nf, Integer(static_cast<long long>(labelled))).to_double()},
                               {"weight", Rational(Integer(static_cast<long long>(labelled)), nf).str()}, {"amplitude", Json()}});
        if (gft) g.set("amplitude", pq::boulatov_amplitude(Gm, group).to_double());
        graphs.push(g);
      }
      Fs[n] = F;
      fe.push(Json::object({{"order", static_cast<double>(n)}, {"terms", polynomial_json(F)}, {"variable", "λλ̄"}}));
      if (n == order && !F.empty()) lead_power = F.rbegin()->first, lead_coeff = F.rbegin()->second.str();
    }
    const auto Zs = exponentiate(Fs);
    for (std::size_t n = 1; n <= order; n++) z.push(Json::object({{"order", static_cast<double>(n)}, {"terms", polynomial_json(Zs[n])}, {"variable", "λλ̄"}}));
  } else {
    const auto T = invariant_theory(id, p);
    const bool gft = gft_theory(id);
    const auto group = gft ? gauge_group(p) : morphe::algebra::cyclic(1);
    // Every distribution of the order among the interactions.
    const std::size_t I = T.interactions.size();
    std::vector<std::size_t> counts(I, 0);
    std::function<void(std::size_t, std::size_t)> each = [&](std::size_t i, std::size_t left) {
      if (i + 1 == I) {
        counts[i] = left;
        std::size_t P = 0;
        for (std::size_t t = 0; t < I; t++) P += counts[t] * T.interactions[t].bubble.p;
        if (P > 8) throw std::invalid_argument("more than 8 fields of each kind at this order: choose a lower order");
        const auto E = pf::expand(T, counts, true);
        Integer denominator = 1;
        for (auto c : counts) denominator *= morphe::arithmetic::factorial(c);
        for (const auto& cl : E.classes) {
          total += cl.labelled;
          if (graphs.size() >= limit) { truncated = true; continue; }
          Json c = Json::array();
          for (auto x : counts) c.push(static_cast<double>(x));
          Json g = Json::object({{"id", static_cast<double>(total - 1)}, {"coloured_graph", coloured_graph_json(cl.graph)}, {"counts", c},
                                 {"connected", cl.connected}, {"n_power", static_cast<double>(cl.faces)}, {"n_power_scaled", static_cast<double>(cl.scaled)},
                                 {"degree", T.rank >= 3 ? Json(static_cast<double>(cl.degree)) : Json()}, {"genus", T.rank == 2 ? Json(static_cast<double>(cl.degree)) : Json()},
                                 {"automorphisms", cl.automorphisms.to_double()},
                                 {"weight", Rational(Integer(static_cast<long long>(cl.labelled)), denominator).str()}, {"amplitude", Json()}});
          if (gft) g.set("amplitude", pq::gft_amplitude(pq::tensorial_two_complex(E, cl.sigma0, T.rank), group).to_double());
          graphs.push(g);
        }
        return;
      }
      for (std::size_t c = 0; c <= left; c++) counts[i] = c, each(i + 1, left - c);
    };
    if (order > 0) each(0, order);
    if (I > 1) {
      // Several bubbles, every coupling λ: Z_n and F_n summed over the distributions of the order.
      for (std::size_t n = 0; n <= order; n++) {
        pf::Polynomial Zn, Fn;
        std::function<void(std::size_t, std::size_t)> dist = [&](std::size_t i, std::size_t left) {
          if (i + 1 == I) {
            counts[i] = left;
            const auto E = pf::expand(T, counts, false);
            Integer denominator = 1;
            for (auto c : counts) denominator *= morphe::arithmetic::factorial(c);
            for (const auto& [f, k] : E.all) Zn[f] += Rational(k, denominator);
            for (const auto& [f, k] : E.connected) Fn[f] += Rational(k, denominator);
            return;
          }
          for (std::size_t c = 0; c <= left; c++) counts[i] = c, dist(i + 1, left - c);
        };
        dist(0, n);
        if (n == 0) continue;
        Zn = signed_order(Zn, n), Fn = signed_order(Fn, n);
        z.push(Json::object({{"order", static_cast<double>(n)}, {"terms", polynomial_json(Zn)}}));
        fe.push(Json::object({{"order", static_cast<double>(n)}, {"terms", polynomial_json(Fn)}}));
        if (n == order && !Fn.empty()) lead_power = Fn.rbegin()->first, lead_coeff = Fn.rbegin()->second.str();
      }
    }
    if (I == 1) {
      const auto S = pf::series(T, order);
      // The series in λ itself: (−1)^n times the coefficients of (−λ)^n.
      for (std::size_t n = 1; n <= order; n++) {
        z.push(Json::object({{"order", static_cast<double>(n)}, {"terms", polynomial_json(signed_order(S.z[n], n))}}));
        fe.push(Json::object({{"order", static_cast<double>(n)}, {"terms", polynomial_json(signed_order(S.free_energy[n], n))}}));
      }
      if (order >= 1) {
        const auto Fo = signed_order(S.free_energy[order], order);
        if (!Fo.empty()) lead_power = Fo.rbegin()->first, lead_coeff = Fo.rbegin()->second.str();
      }
    }
  }
  // The leading graphs first.
  std::vector<Json> sorted(graphs.as_array().begin(), graphs.as_array().end());
  std::stable_sort(sorted.begin(), sorted.end(), [](const Json& a, const Json& b) { return a["n_power_scaled"].as_number() > b["n_power_scaled"].as_number(); });
  Json g2 = Json::array();
  for (auto& x : sorted) g2.push(x);
  return Json::object({{"theory", id}, {"order", static_cast<double>(order)}, {"graphs", g2}, {"total", static_cast<double>(total)}, {"truncated", truncated},
                       {"series", Json::object({{"Z", z}, {"F", fe}})},
                       {"leading", lead_power == LLONG_MIN ? Json() : Json::object({{"power", static_cast<double>(lead_power)}, {"coeff", lead_coeff}})}});
}

// ---- the large-N or condensed phase --------------------------------------------------------------------------------------

inline Json large_n(const Json& request) {
  const std::string id = request["theory"].text();
  const Json& p = request["params"];
  // c_k x_c^k, formed exactly before rounding, so that series growing like x_c^{−k} stay finite.
  auto scaled = [](const auto& c, const Rational& xc) {
    std::vector<double> d;
    Rational power = 1;
    for (const auto& x : c) d.push_back((Rational(x) * power).to_double()), power *= xc;
    return d;
  };
  if (id == "matrix") {
    // The Hermitian quartic model V = M²/2 + g M⁴ (Brézin, Itzykson, Parisi and Zuber 1978).
    const double gc = pf::bipz_critical;
    double g = number(p, "g", 0.02);
    if (g <= gc) g = gc * 0.999;
    const double a2 = pf::bipz_a2(g), a = std::sqrt(a2);
    std::vector<double> x, rho, gs, e0;
    for (int i = 0; i <= 240; i++) {
      const double l = -2 * a + 4 * a * i / 240.0;
      x.push_back(l), rho.push_back(pf::bipz_density(l, g));
    }
    for (int i = 0; i <= 200; i++) {
      const double gg = gc * (1 - 1e-6) + (0.1 - gc) * std::pow(i / 200.0, 1.5);
      gs.push_back(gg), e0.push_back(pf::bipz_free_energy(gg));
    }
    const auto c = pf::bipz_series(400);
    const double alpha = pf::exponent_from_series(scaled(c, Rational(-1, 48)), 1.0);
    Json coeffs = Json::array();
    for (std::size_t k = 0; k < 12 && k < c.size(); k++) coeffs.push(c[k].str());
    return Json::object({{"kind", "matrix"}, {"model", "Hermitian quartic, V = M²/2 + g M⁴"}, {"g", g}, {"a2", a2},
                         {"density", Json::object({{"x", Json(x)}, {"rho", Json(rho)}})}, {"free_energy", Json::object({{"g", Json(gs)}, {"E0", Json(e0)}})},
                         {"critical", Json::object({{"g_c", gc}, {"gamma", -0.5}, {"gamma_estimate", 3 - alpha}})}, {"series", coeffs}});
  }
  if (id == "vector") {
    const double mu = 1, lc = pf::vector_critical(mu);
    std::vector<double> l, gap;
    for (int i = 0; i <= 200; i++) {
      const double x = lc * (1 - 1e-9) + (1.0 - lc) * i / 200.0;
      l.push_back(x), gap.push_back(pf::vector_two_point(x, mu));
    }
    return Json::object({{"kind", "vector"}, {"coupling", Json(l)}, {"gap", Json(gap)}, {"critical", Json::object({{"lambda_c", lc}})}});
  }
  if (id == "scalar") return Json::object({{"kind", "none"}, {"note", "a zero-dimensional complex scalar has no large-N limit"}});
  // Melonic: G = 1 + x G^m. The quartic melonic tensor model has m = 2 (x = −2Dλ, D colours);
  // the coloured simplicial models m = D + 1 (x = λλ̄).
  int m = 2;
  if (coloured_theory(id)) m = static_cast<int>(coloured_dimension(id, p)) + 1;
  const double xc = pf::melonic_critical(m);
  std::vector<double> xs, G;
  for (int i = 0; i <= 200; i++) {
    const double x = -xc + 2 * xc * i / 200.0;
    xs.push_back(x), G.push_back(pf::melonic_two_point(std::min(x, xc), m));
  }
  const auto c = pf::melonic_coefficients(m, 400);
  const Rational xr(morphe::arithmetic::pow(Integer(m - 1), static_cast<unsigned long long>(m - 1)), morphe::arithmetic::pow(Integer(m), static_cast<unsigned long long>(m)));
  const double alpha = pf::exponent_from_series(scaled(c, xr), 1.0);
  Json coeffs = Json::array();
  for (std::size_t k = 0; k < 12; k++) coeffs.push(c[k].str());
  Json out = Json::object({{"kind", gft_theory(id) ? "condensate" : "melonic"}, {"m", m}, {"lambda", Json(xs)}, {"G", Json(G)},
                           {"critical", Json::object({{"lambda_c", xc}, {"gamma", 0.5}, {"gamma_estimate", 2 - alpha}})}, {"series", coeffs},
                           {"equation", "G = 1 + x G^" + std::to_string(m)}});
  return out;
}

}  // namespace atlas
