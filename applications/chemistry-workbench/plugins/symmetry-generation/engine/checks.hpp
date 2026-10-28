// The manuscript's claims checked: for each claim, named by an id that the verification note's
// sections follow (theory/formal-verification.md), the computations that confirm or refute it.

#pragma once

#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "fields.hpp"
#include "grammars.hpp"
#include "selection.hpp"

namespace chem {

struct Check {
  std::string name;
  bool ok = false;
  std::string detail;
};

using Distribution = std::vector<Rational>;
inline Rational total_variation(const Distribution& p, const Distribution& q) {
  Rational s = 0;
  for (std::size_t i = 0; i < p.size(); i++) s = s + (p[i] > q[i] ? p[i] - q[i] : q[i] - p[i]);
  return s / Rational(2);
}
// The law p K of a distribution through a stochastic matrix K (rows: inputs).
inline Distribution through(const Distribution& p, const std::vector<Distribution>& K) {
  Distribution out(K[0].size(), Rational(0));
  for (std::size_t i = 0; i < p.size(); i++)
    for (std::size_t j = 0; j < out.size(); j++) out[j] = out[j] + p[i] * K[i][j];
  return out;
}
inline Distribution random_law(std::mt19937& rng, std::size_t n) {
  Distribution p(n);
  long long total = 0;
  std::vector<long long> w(n);
  for (auto& x : w) total += x = 1 + rng() % 9;
  for (std::size_t i = 0; i < n; i++) p[i] = Rational(Integer(w[i]), Integer(total));
  return p;
}

inline std::string yes(bool b) { return b ? "yes" : "no"; }
inline std::string num(double x) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.6g", x);
  return buf;
}

inline Json checks() {
  std::map<std::string, std::vector<Check>> R;
  const auto add = [&](const std::string& id, std::string name, bool ok, std::string detail = "") { R[id].push_back({std::move(name), ok, std::move(detail)}); };

  // ---- Formalism I ----
  const Json ch = stereocentre();
  add("prop-9.1", "two terminal classes, one reflection orbit, one constitutional class",
      ch["terminalClasses"].as_number() == 2 && ch["orbits"].as_number() == 1 && ch["constitutionalClasses"].as_number() == 1,
      num(ch["terminalClasses"].as_number()) + " classes, " + num(ch["orbits"].as_number()) + " orbit");
  const Json& geo = ch["geometry"];
  add("prop-9.1", "τ(X₊) = +1, τ(X₋) = −1 = τ(pX₊)", geo["tauPlus"].as_number() == 1 && geo["tauMinus"].as_number() == -1);
  add("prop-9.1", "X₊ and X₋ are congruent only by an improper isometry", geo["congruences"].as_number() == 1 && geo["congruenceDet"].as_number() == -1,
      num(geo["congruences"].as_number()) + " congruence, det " + num(geo["congruenceDet"].as_number()));
  add("prop-9.1", "the configuration's point group is C1 (it is chiral)", geo["pointGroup"].as_string() == "C1", geo["pointGroup"].as_string());
  add("tau-cip", "with the order (H, F, Cl, Br), τ = +1 is the CIP descriptor R", geo["cipPlus"].as_string() == "R" && geo["cipMinus"].as_string() == "S",
      "X₊ is " + geo["cipPlus"].as_string() + ", X₋ is " + geo["cipMinus"].as_string());
  add("tau-cip", "an odd relabelling of the ligands (F ↔ Cl) flips τ as a reflection does", geo["tauAfterSwap"].as_number() == -1);
  add("prop-4.1", "the stereocentre's rules are closed under p and their derivations commute with it", ch["closed"].as_bool() && ch["equivariant"].as_bool());
  add("prop-4.1", "the rules {p₊} alone are not closed, and equivariance fails (the hypothesis is needed)", !ch["halfClosed"].as_bool() && !ch["halfEquivariant"].as_bool());
  const Json gly = glycine(3);
  add("prop-4.1", "the coupling rule commutes with S₃ relabelling at every derivable graph", gly["labelled"]["equivariant"].as_bool());
  add("glycine", "labelled units: 6 terminal backbones, each by 2 derivations (12 in all), one S₃ orbit",
      gly["labelled"]["terminalClasses"].as_number() == 6 && gly["labelled"]["terminalDerivations"].as_number() == 12 && gly["labelled"]["orbits"].as_number() == 1,
      num(gly["labelled"]["terminalClasses"].as_number()) + " backbones, " + num(gly["labelled"]["terminalDerivations"].as_number()) + " derivations");
  add("glycine", "identical units: one backbone class reached by 12 derivations",
      gly["identical"]["terminalClasses"].as_number() == 1 && gly["identical"]["terminalDerivations"].as_number() == 12);
  add("glycine", "two waters created; 3 Gly → Gly₃ + 2 H₂O balances C, H, N, O", gly["labelled"]["water"].as_number() == 2 && gly["balanced"].as_bool());
  const Json pr = proton_transfer();
  add("eq-stoich", "CH₃COOH + H₂O → CH₃COO⁻ + H₃O⁺: AN = 0 and zᵀN = 0", pr["balanced"].as_bool());
  add("def-4.2", "the edit keeps every atom and changes two oxygens' charges (a relabelling of kept vertices)",
      pr["before"].dump() == "{\"C\":2,\"H\":6,\"O\":3,\"charge\":0}" && pr["after"]["O"].as_number() == 3 && pr["after"]["charge"].as_number() == 0,
      "before " + pr["before"].dump() + ", after " + pr["after"].dump());
  add("counting", "matches = |Aut L| × occurrences: 2 matches of the edit, |Aut L| = 2, one product",
      pr["matches"].as_number() == 2 && pr["automorphismsOfL"].as_number() == 2 && pr["productClasses"].as_number() == 1);
  const Json eth = ethanol();
  bool ethanol_ok = false;
  for (const auto& c : eth["space"]["classes"].as_array())
    if (c["terminal"].as_bool()) ethanol_ok = c["derivations"].as_number() == 1 && c["graph"]["vertices"].size() == 3;
  add("grammar-syntax", "A0 ⇒ CH₃–Tail₁ ⇒ CH₃–CH₂–OH: one terminal graph, one derivation", ethanol_ok);

  // ---- dynamics ----
  FrankParameters P;
  P.precursor = 24;
  const Json fr = frank(P);
  add("eq-stoich", "Frank's network conserves a + ℓ + r + 2d (the left kernel of N)", fr["conservation"].dump() == "[[1,1,1,2]]" || fr["conservation"].dump() == "[[-1,-1,-1,-2]]",
      fr["conservation"].dump());
  add("prop-4.2", "the Markov generator commutes with the mirror L ↔ R, exactly", fr["chain"]["commutes"].as_bool() && fr["mirrorSymmetric"].as_bool(),
      num(fr["chain"]["states"].as_number()) + " states");
  add("prop-4.2", "the law of outcomes from (A, L, R) = (8, 1, 1) is exactly mirror-symmetric (rational arithmetic)", frank_law_exactly_symmetric(8));
  add("prop-4.2", "from (24, 1, 1, 0) the outcome law is symmetric while each outcome is asymmetric",
      fr["chain"]["asymmetry"].as_number() < 1e-12 && fr["chain"]["extinct"].as_number() < 1, "P(both extinct) = " + num(fr["chain"]["extinct"].as_number()));
  add("frank-m", "ṁ = (k_i s/2) m (1 − m²) along the integrated solution", fr["ode"]["identityError"].as_number() < 1e-9, "relative error " + num(fr["ode"]["identityError"].as_number()));
  add("frank-m", "a + ℓ + r + 2d is conserved by the integration", fr["ode"]["conservationError"].as_number() < 1e-9);
  {
    // The minor enantiomer decays at about k_i δ∞, slowly when the excess was small: from m₀ = 0.1
    // the limit is reached well within t = 400.
    FrankParameters Q;
    Q.t_end = 400, Q.l0 = 0.011, Q.r0 = 0.009;
    const Json f2 = frank(Q);
    add("frank-m", "from m₀ = 0.1 the excess tends to 1 (the minor enantiomer is consumed)", f2["ode"]["mFinal"].as_number() > 0.999, "m(400) = " + num(f2["ode"]["mFinal"].as_number()));
    Q.l0 = Q.r0 = 0.01;
    const Json f3 = frank(Q);
    double worst = 0;
    for (const auto& x : f3["ode"]["m"].as_array()) worst = std::max(worst, std::fabs(x.as_number()));
    add("lemma-4.1", "from a racemic start the deterministic solution stays racemic (m ≡ 0)", worst == 0, "max |m| = " + num(worst));
    Q.l0 = 0.0101;
    const Json f4 = frank(Q);
    bool sign = true;
    for (std::size_t i = 1; i < f4["ode"]["m"].size(); i++) sign = sign && f4["ode"]["m"][i].as_number() > 0;
    add("lemma-4.1", "the sign of ℓ − r never changes (δ̇ = k a δ), so the stabilizer is constant in finite time", sign);
  }
  {
    const double a = 0.7, b = 2;
    const double m = std::sqrt(a / b);
    add("normal-form", "ṁ = am − bm³: f′(±√(a/b)) = −2a and f′(0) = a", std::fabs((a - 3 * b * m * m) + 2 * a) < 1e-12);
  }

  // ---- Formalism II ----
  const Json cx = complexes();
  const auto homology_ranks = [](const Json& c) {
    std::string s;
    for (const auto& h : c["homology"].as_array()) s += (s.empty() ? "" : ",") + num(h["rank"].as_number());
    return s;
  };
  add("app-L", "the six-cycle: Betti (1, 1)", homology_ranks(cx["ring"]) == "1,1", homology_ranks(cx["ring"]));
  add("app-L", "the filled six-cycle: Betti (1, 0, 0), ∂₁∂₂ = 0", homology_ranks(cx["filled"]) == "1,0,0" && cx["filled"]["boundaryOfBoundary"].as_bool());
  add("prop-9.2", "the dipole: (V, E, F) = (2, 4, 6), χ = 4", cx["dipole"]["counts"].dump() == "[2,4,6]" && cx["dipole"]["euler"].as_number() == 4);
  add("prop-9.2", "∂₁∂₂ = 0 and homology ℤ, 0, ℤ³ (torsion-free)",
      cx["dipole"]["boundaryOfBoundary"].as_bool() && cx["dipole"]["homology"][1]["group"].as_string() == "0" && cx["dipole"]["homology"][2]["group"].as_string() == "ℤ^3",
      cx["dipole"]["homology"].dump());
  const Json w2 = wick(2, "Z2");
  add("prop-9.2", "of the 24 gluings of two vertices, 4 give six faces (the dipole), 4 four, 16 three", w2["gluingFaces"].dump() == "{\"3\":16,\"4\":4,\"6\":4}",
      w2["gluingFaces"].dump());
  add("eq-complex-expansion", "⟨V²⟩ over ℤ₂: Σ over the 105 pairings of |G|^F equals the Gaussian moment", w2["sum"].as_string() == w2["moment"].as_string(), w2["sum"].as_string());
  const Json w1 = wick(1, "Z3");
  add("eq-complex-expansion", "⟨V⟩ over ℤ₃: Σ over pairings of |G|^F equals the Gaussian moment", w1["sum"].as_string() == w1["moment"].as_string(), w1["sum"].as_string());
  {
    bool ok = true;
    for (const auto& c : w2["classes"].as_array()) ok = ok && c["pairings"].as_number() * c["automorphisms"].as_number() == 32;
    add("eq-complex-expansion", "each class Γ holds |Aut v|² 2!/|Aut Γ| = 32/|Aut Γ| pairings", ok, num(w2["classes"].size()) + " classes");
  }
  const Json fz = field("Z2", Json());
  add("example-5.1", "V_tet[1] = 64 over ℤ₂", fz["value"].as_number() == 64);
  {
    std::vector<double> odd(8);
    for (std::size_t p = 0; p < 8; p++) odd[p] = pq::coordinates(p, 2, 3)[0] ? -1 : 1;
    const Json fo = field("Z2", Json(odd));
    add("example-5.1", "V_tet[(−1)^a] = 0 over ℤ₂", fo["value"].as_number() == 0);
  }
  {
    std::mt19937 rng(5);
    std::vector<double> psi(8);
    for (auto& x : psi) x = static_cast<double>(static_cast<int>(rng() % 9) - 4);
    const Json fr2 = field("Z2", Json(psi));
    std::string inv;
    for (const auto& t : fr2["translations"].as_array()) inv += t["invariant"].as_bool() ? "1" : "0";
    add("prop-C.1", "V_tet is invariant under the diagonal translation (1,1,1)", fr2["translations"][7]["invariant"].as_bool());
    add("symmetry-GxG", "and under exactly the four (g, h, g) of ℤ₂³: (000), (010), (101), (111)", inv == "10100101", inv);
    std::string modes;
    for (const auto& m : fz["modes"].as_array()) modes += m["V"].as_number() == 64 ? "1" : "0";
    add("flat-eta", "V_tet along the characters (−1)^{s·x}: 0 along (−1)^a, 64 along (−1)^b", fz["modes"][4]["V"].as_number() == 0 && fz["modes"][2]["V"].as_number() == 64,
        "s = 000…111: " + modes);
    bool killed = true, nonzero = false;
    for (const auto& m : fr2["modes"].as_array())
      if (m["odd"].as_bool()) killed = killed && m["etaProjected"].as_number() == 0, nonzero = nonzero || m["eta"].as_number() != 0;
    add("gauge-eta", "the gauge projection annihilates every odd order parameter over ℤ₂ (the diagonal translation is the gauge)", killed && nonzero);
    add("eta-incomplete", "(−1)^b is odd under U_p and has η = 0: η = 0 does not imply symmetry",
        [&] {
          std::vector<double> b(8);
          for (std::size_t p = 0; p < 8; p++) b[p] = pq::coordinates(p, 2, 3)[1] ? -1 : 1;
          const Json fb = field("Z2", Json(b));
          return fb["modes"][4]["eta"].as_number() == 0 && fb["modes"][2]["eta"].as_number() == 1;
        }());
    add("regulator", "V_tet takes negative values, so e^{−λV} is not integrable for any λ ≠ 0 without a regulator", fz["negative"]["V"].as_number() < 0,
        "V = " + num(fz["negative"]["V"].as_number()) + " on " + fz["negative"]["psi"].dump());
  }
  {
    const double alpha = -1.3, beta = 0.4;
    const double eta = std::sqrt(-alpha / (2 * beta));
    const double V = alpha * eta * eta + beta * eta * eta * eta * eta;
    add("prop-5.1", "η± = ±√(−α/2β), V(η±) = −α²/4β, V″(η±) = −4α", std::fabs(V + alpha * alpha / (4 * beta)) < 1e-12 &&
                                                                         std::fabs(2 * alpha + 12 * beta * eta * eta + 4 * alpha) < 1e-12);
  }

  // ---- probability kernels (§5.4, §6, §10, §11, App. F) ----
  {
    std::mt19937 rng(17);
    bool contraction = true, mixture = true;
    for (int trial = 0; trial < 200; trial++) {
      const auto p = random_law(rng, 5), q = random_law(rng, 5);
      std::vector<Distribution> K(5);
      for (auto& row : K) row = random_law(rng, 4);
      contraction = contraction && total_variation(through(p, K), through(q, K)) <= total_variation(p, q);
      // Mixtures: d(Σκ E_h, Σκ Ẽ_h) ≤ Σ κ d(E_h, Ẽ_h).
      const auto kappa = random_law(rng, 3);
      std::vector<Distribution> E(3), F(3);
      for (auto& e : E) e = random_law(rng, 4);
      for (auto& e : F) e = random_law(rng, 4);
      Rational bound = 0;
      for (int h = 0; h < 3; h++) bound = bound + kappa[h] * total_variation(E[h], F[h]);
      mixture = mixture && total_variation(through(kappa, E), through(kappa, F)) <= bound;
    }
    add("prop-10.1", "d_TV(PK, QK) ≤ d_TV(P, Q) for 200 random laws and kernels (exact)", contraction);
    add("prop-10.1", "the mixture bound, for 200 random cases (exact)", mixture);
    // Pushforward covariance: ℤ₂ acting on complexes {0,1,2,3} by 0↔1, 2↔3 and on states {0,1,2}
    // by 0↔1; an invariant ensemble and an equivariant readout give an invariant output law.
    const std::vector<std::size_t> gk{1, 0, 3, 2}, gx{1, 0, 2};
    Distribution pk{Rational(1, 6), Rational(1, 6), Rational(1, 3), Rational(1, 3)};
    std::vector<Distribution> Rk{{Rational(1, 2), Rational(1, 4), Rational(1, 4)}, {}, {Rational(1, 5), Rational(1, 5), Rational(3, 5)}, {}};
    for (int k : {0, 2}) {
      Rk[gk[k]].assign(3, Rational(0));
      for (int x = 0; x < 3; x++) Rk[gk[k]][gx[x]] = Rk[k][x];
    }
    const auto Q = through(pk, Rk);
    add("prop-5.2", "an invariant ensemble and an equivariant readout give an invariant law Q(gA) = Q(A)", Q[0] == Q[1], "Q = (" + Q[0].str() + ", " + Q[1].str() + ", " + Q[2].str() + ")");
    // The constant readout erases every distinction.
    std::vector<Distribution> constant(4, Distribution{Rational(1), Rational(0), Rational(0)});
    add("prop-5.2", "a constant readout is equivariant yet maps both broken sectors to one state",
        through(Distribution{1, 0, 0, 0}, constant) == through(Distribution{0, 1, 0, 0}, constant));
    // Unrestricted readout: R(·|k) = P* reproduces any target from any ensemble.
    const Distribution target{Rational(1, 7), Rational(2, 7), Rational(4, 7)};
    add("prop-11.1", "the constant kernel R(·|k) = P* gives P* from any normalized ensemble",
        through(pk, std::vector<Distribution>(4, target)) == target && through(random_law(rng, 4), std::vector<Distribution>(4, target)) == target);
    // Factorization (Prop. 6.1): a kernel constant on the fibres of T factors; Example 6.1's δ.
    const Rational delta(3, 10);
    const Distribution plus{(Rational(1) + delta) / Rational(2), (Rational(1) - delta) / Rational(2)},
        minus{(Rational(1) - delta) / Rational(2), (Rational(1) + delta) / Rational(2)};
    add("prop-6.1", "Example 6.1: d_TV between the two responses is δ, though T_r(x₊) = T_r(x₋)", total_variation(plus, minus) == delta, "δ = " + delta.str());
    // Contrastive models (App. F): equal under C₀, |q − p| apart under C*.
    const Rational pp(1, 3), qq(3, 4);
    add("prop-8.1", "the two-state models agree in C₀ and differ by |q − p| in C*", total_variation({Rational(1) - pp, pp}, {Rational(1) - pp, pp}) == Rational(0) &&
                                                                                       total_variation({Rational(1) - qq, qq}, {Rational(1) - pp, pp}) == qq - pp);
  }

  // ---- selection rules on molecules (§7.2, App. E, App. I) ----
  {
    std::map<std::string, Json> pg;
    for (const auto& m : catalogue()) pg[m.id] = point_group_json(m);
    add("prop-7.1", "(E)-1,2-dichloroethene: C2h, no fixed vector (no permanent dipole)",
        pg["dichloroethene-e"]["schoenflies"].as_string() == "C2h" && pg["dichloroethene-e"]["dipole"].as_number() == 0);
    add("prop-7.1", "(Z)-1,2-dichloroethene: C2v, one fixed direction (a dipole is allowed)",
        pg["dichloroethene-z"]["schoenflies"].as_string() == "C2v" && pg["dichloroethene-z"]["dipole"].as_number() == 1);
    add("prop-7.1", "CHBrClF: C1, chiral, every vector allowed", pg["chbrclf-plus"]["schoenflies"].as_string() == "C1" && pg["chbrclf-plus"]["chiral"].as_bool());
    add("prop-7.1", "benzene D6h, methane Td, SF₆ Oh: no dipole",
        pg["benzene"]["schoenflies"].as_string() == "D6h" && pg["methane"]["schoenflies"].as_string() == "Td" && pg["sulfur-hexafluoride"]["schoenflies"].as_string() == "Oh" &&
            pg["benzene"]["dipole"].as_number() + pg["methane"]["dipole"].as_number() + pg["sulfur-hexafluoride"]["dipole"].as_number() == 0);
    add("prop-7.1", "hydrogen peroxide (skew): C2, chiral, a dipole along the axis",
        pg["hydrogen-peroxide"]["schoenflies"].as_string() == "C2" && pg["hydrogen-peroxide"]["chiral"].as_bool() && pg["hydrogen-peroxide"]["dipole"].as_number() == 1);
  }

  Json out = Json::object();
  for (const auto& [id, checks] : R) {
    Json cj = Json::array();
    for (const auto& c : checks) cj.push(Json::object({{"name", c.name}, {"ok", c.ok}, {"detail", c.detail}}));
    out.set(id, cj);
  }
  return out;
}

}  // namespace chem
