// Checks of the shared problem's header (economy.hpp), against closed forms and an independent recomputation:
//
//   1. the encoding: decode(encode(p)) = p for the plan as built and for random plans; every point of the cube
//      decodes to a valid plan;
//   2. the distances: zero on the diagonal, symmetric, the triangle inequality on every triple, and no greater with
//      every link made than with none;
//   3. Parseval: the week's lines of the spending sum to its mean square, so the shares of the mean, daily and other
//      lines sum to one, and the spread is 1/(1 + CV²);
//   4. the value against a recomputation written apart (a second loop over hours and frontages, with competition and
//      catalysis written out), for random plans under random flows;
//   5. closed forms: every frontage vacant gives V = 0 and access at the cap; one open frontage alone, without
//      catalysis, captures P π φ s m an hour; two frontages of one kind on one segment share 1 − (1 − p)² of the
//      passers (each p(1 − p/2));
//   6. the value is linear in the flows when nothing else changes (twice the passers, twice the margin).
//
//   check_economy <results/problem-inputs.json> [results/problem.json]

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>

#include <calx/random/stream.hpp>

#include "economy.hpp"

namespace {

int failures = 0, checks = 0;
void check(bool ok, const std::string& what) {
  checks++;
  if (!ok) failures++, std::printf("FAIL  %s\n", what.c_str());
}
bool near(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol * std::max(1.0, std::max(std::fabs(a), std::fabs(b))); }

// The value written apart from Problem::evaluate.
double value_apart(const economy::Problem& P, const economy::Plan& plan, const std::vector<double>& passers) {
  using namespace economy;
  const auto Dist = P.distances(plan.links);
  const std::size_t n = P.size();
  double margin = 0;
  for (std::size_t h = 0; h < week; h++) {
    std::vector<double> p(n, 0.0);
    for (std::size_t f = 0; f < n; f++) {
      const auto k = plan.kind[f];
      if (!trades(k) || !plan.open(f, h)) continue;
      double kappa = 0;
      for (std::size_t g = 0; g < n; g++) {
        if (g == f || !plan.open(g, h)) continue;
        const double d = Dist->at(f, g);
        if (d <= 3 * P.P.coupling) kappa += std::exp(-d / P.P.coupling) * P.P.affinity[static_cast<std::size_t>(plan.kind[g])][static_cast<std::size_t>(k)];
      }
      const double t = static_cast<double>(h % 24);
      const std::size_t b = t >= 18 ? 1 : (t >= 7 ? 0 : 2);
      p[f] = P.P.propensity[static_cast<std::size_t>(k)] * P.P.period[static_cast<std::size_t>(k)][b] * std::min(2.0, std::max(0.5, 1 + P.P.catalysis * kappa));
    }
    for (std::size_t f = 0; f < n; f++) {
      if (p[f] <= 0) continue;
      double none = 1, sum = 0;
      for (std::size_t g = 0; g < n; g++)
        if (p[g] > 0 && plan.kind[g] == plan.kind[f] && P.frontages[g].segment == P.frontages[f].segment) none *= 1 - p[g], sum += p[g];
      const auto k = static_cast<std::size_t>(plan.kind[f]);
      margin += passers[f * week + h] * p[f] * (1 - none) / sum * P.P.spend[k] * P.P.margin[k];
    }
  }
  double cost = 0;
  for (std::size_t f = 0; f < n; f++) {
    const auto k = static_cast<std::size_t>(plan.kind[f]);
    if (plan.kind[f] == vacant) continue;
    double hours = 0, days = 0;
    for (std::size_t h = 0; h < week; h++) hours += plan.open(f, h);
    for (std::size_t d = 0; d < 7; d++) {
      bool any = false;
      for (std::size_t t = 0; t < 24; t++) any = any || plan.open(f, 24 * d + t);
      days += any;
    }
    cost += P.P.staff[k] * P.P.wage * (hours + days) + P.P.outgoings[k] * P.frontages[f].trading_area;
    if (plan.kind[f] != P.frontages[f].kind0) cost += P.P.refit_weekly(P.frontages[f].trading_area);
  }
  for (std::size_t q = 0; q < P.links.size(); q++)
    if (plan.links[q]) cost += P.P.link_weekly(P.links[q].length);
  return margin - cost;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace economy;
  if (argc < 2) {
    std::cerr << "usage: check_economy <results/problem-inputs.json> [results/problem.json]\n";
    return 2;
  }
  const auto read = [](const char* path) {
    std::ifstream in(path, std::ios::binary);
    std::stringstream b;
    b << in.rdbuf();
    return Json::parse(b.str());
  };
  const Json in = read(argv[1]);
  Json pj;
  if (argc > 2) pj = read(argv[2]);
  const Problem P(in, argc > 2 ? &pj : nullptr);
  const std::size_t n = P.size(), L = P.links.size();
  calx::random::Stream s(20261105, 1, 0, 0);
  Kernels K{&P};

  // 1. The encoding.
  {
    const Plan a = P.as_built();
    check(P.decode(P.encode(a)) == a, "decode(encode(as built)) = as built");
    for (int t = 0; t < 200; t++) {
      const Plan p = K.sample(s);
      if (!(P.decode(P.encode(p)) == p)) { check(false, "decode(encode(p)) = p for a random plan"); break; }
      std::vector<double> u(P.dimension());
      for (auto& x : u) x = s.uniform();
      try {
        P.validate(P.decode(u));
      } catch (const std::exception& e) {
        check(false, std::string("a point of the cube decodes to a valid plan: ") + e.what());
      }
    }
    check(true, "encoding round trips");
  }
  // 2. Distances.
  {
    const auto D0 = P.distances(hv::Bits(L, 0));
    const auto D1 = P.distances(hv::Bits(L, 1));
    bool diag = true, sym = true, tri = true, mono = true;
    for (std::size_t f = 0; f < n; f++) {
      diag = diag && D0->at(f, f) == 0;
      for (std::size_t g = 0; g < n; g++) {
        sym = sym && near(D0->at(f, g), D0->at(g, f), 1e-12);
        mono = mono && D1->at(f, g) <= D0->at(f, g) + 1e-6;
        for (std::size_t q = 0; q < n; q++) tri = tri && D0->at(f, g) <= D0->at(f, q) + D0->at(q, g) + 1e-6;
      }
    }
    check(diag, "d(f, f) = 0");
    check(sym, "d symmetric");
    check(tri, "the triangle inequality on every triple");
    check(mono, "every link made: no distance grows");
    std::printf("distances: median %.1f m with no link, %zu frontages with neighbours within 3 coupling lengths\n", D0->at(0, n / 2),
                static_cast<std::size_t>(std::count_if(D0->near.begin(), D0->near.end(), [](const auto& v) { return !v.empty(); })));
  }
  // Random flows: a smooth day at each frontage, scaled per frontage.
  std::vector<double> flows(n * week);
  for (std::size_t f = 0; f < n; f++) {
    const double scale = 200 + 3000 * s.uniform();
    for (std::size_t h = 0; h < week; h++) {
      const double t = static_cast<double>(h % 24);
      flows[f * week + h] = scale * (0.1 + std::exp(-0.5 * std::pow((t - 13) / 4, 2))) * (is_weekend(h) ? 0.7 : 1.0);
    }
  }
  // 3 and 4. Parseval, and the value written apart.
  {
    double worst = 0;
    bool parseval = true, shares = true;
    for (int t = 0; t < 30; t++) {
      Plan p = t == 0 ? P.as_built() : K.sample(s);
      const auto o = P.evaluate(p, flows);
      const double v = value_apart(P, p, flows);
      worst = std::max(worst, std::fabs(o.value - v) / std::max(1.0, std::fabs(v)));
      double sum = 0, ms = 0, mean = 0;
      for (double x : o.power) sum += x;
      for (double x : o.spending_by_hour) ms += x * x / week, mean += x / week;
      parseval = parseval && near(sum, ms, 1e-9);
      shares = shares && near(o.spread + o.daily_share + o.weekly_share, 1, 1e-12) && (ms == 0 || near(o.spread, mean * mean / ms, 1e-9));
    }
    check(parseval, "Parseval: the week's lines sum to the mean square");
    check(shares, "the shares sum to one; the spread is 1/(1 + CV^2)");
    check(worst < 1e-9, "the value equals its recomputation written apart");
    std::printf("value against its recomputation: largest relative difference %.2e\n", worst);
  }
  // 5. Closed forms.
  {
    Plan none = P.as_built();
    for (auto& k : none.kind) k = vacant;
    const auto o = P.evaluate(none, flows);
    check(o.value == 0 && o.spending == 0, "every frontage vacant: V = 0");
    check(near(o.access, P.P.access_cap), "every frontage vacant: access at the cap");
    // One food frontage alone, open all day every day, without catalysis (it has no open neighbour).
    std::size_t f = 0;
    while (f < n && !P.frontages[f].fronts) f++;
    Plan one = none;
    one.kind[f] = food;
    for (std::size_t q = 0; q < 4; q++) one.hours[4 * f + q] = q % 2 == 0 ? 0 : 24;
    const auto o1 = P.evaluate(one, flows);
    double expect = 0;
    for (std::size_t h = 0; h < week; h++) expect += flows[f * week + h] * P.P.propensity[food] * P.P.period[food][band(h)] * P.P.spend[food] * P.P.margin[food];
    check(near(o1.margin, expect, 1e-12), "one frontage alone captures P pi phi s m");
    // Two of one kind on one segment.
    std::size_t g = n;
    for (auto q : P.on_segment[f]) { g = q; break; }
    if (g < n) {
      Parameters Q = P.P;
      Plan two = one;
      two.kind[g] = food;
      for (std::size_t q = 0; q < 4; q++) two.hours[4 * g + q] = q % 2 == 0 ? 0 : 24;
      // With equal flows at both and no catalysis, each takes p(1 - p/2) of the passers.
      std::vector<double> equal(n * week, 0.0);
      for (std::size_t h = 0; h < week; h++) equal[f * week + h] = equal[g * week + h] = 1000;
      Problem copy(in);
      for (auto& row : copy.P.affinity) row.fill(0.0);
      const auto o2 = copy.evaluate(two, equal);
      double ex2 = 0;
      for (std::size_t h = 0; h < week; h++) {
        const double pp = Q.propensity[food] * Q.period[food][band(h)];
        ex2 += 2 * 1000 * pp * (1 - pp / 2) * Q.spend[food] * Q.margin[food];
      }
      check(near(o2.margin, ex2, 1e-12), "two of one kind on one segment share 1 - (1 - p)^2");
    }
    // 6. Linear in the flows.
    std::vector<double> twice(flows);
    for (auto& x : twice) x *= 2;
    const Plan p = P.as_built();
    const auto a = P.evaluate(p, flows), b = P.evaluate(p, twice);
    check(near(b.margin, 2 * a.margin, 1e-12), "twice the passers, twice the margin");
    std::printf("as built under the check's flows: V = %.0f A$/week, spending %.0f, spread %.3f, mix %.3f nats, access %.1f m\n", a.value, a.spending, a.spread, a.mix, a.access);
  }
  std::printf("%d checks, %d failed\n", checks, failures);
  return failures ? 1 : 0;
}
