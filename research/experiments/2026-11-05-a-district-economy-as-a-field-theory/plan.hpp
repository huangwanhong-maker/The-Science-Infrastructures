// The plan as a field coupled to the people's (README §5): exact local conditionals of the value, and a heat bath on
// them, annealed to its zero-temperature limit (iterated conditional modes).
//
// The value of a plan σ under the people's field is (economy.hpp, model.hpp)
//
//   V(σ) = Σ_h Σ_f A_f(h; σ) [P_f(h; σ)]₊ − Σ_f C_f(σ_f) − Σ_links,
//   A_f(h) = [open] p_f c_f s_k m_k,   P(h) = P⁰(h) + R(h) (p(h) − p⁰(h)),
//
// a sum over hours of terms that couple a frontage only to its catalysis neighbours (within 3ξ), its competitors on
// its segment and, through the people's linear response R, to everyone at the same hour. A change at one frontage f
// changes, at each hour, p on N(f) = {f} ∪ near(f) and A on the segments of N(f), S(f); its exact change is
//
//   ΔV_h = Σ_g (A'_g [P_g + Σ_{g'∈N(f)} R_gg'(h) Δp_g']₊ − A_g [P_g]₊),   A' = A off S(f),
//
// n |N(f)| products an hour. Where no passer count is held at zero ([·]₊), this is Σ_{g∈N(f)} Ā_g Δp_g plus local
// terms, with Ā(h) = R(h)ᵀ A(h) the adjoint of the response: the value's linear response to each frontage's
// attraction, the connected correlation of the value with the coupling (README §5.1). And since the hours enter only through sums over hours,
// f's conditional value is, for each kind k, a constant plus the sums of an hourly margin over its weekday and its
// weekend interval, less their staffing: every one of a frontage's 6 × 600 × 600 states is valued exactly from 6 × 168
// local terms. The heat bath draws f's state from exp(V/T) restricted to f, the encoding's states as they are (a
// vacant frontage's 600² hour genes all alike), and the links are drawn by exact full evaluations.
//
// Provisional code of the study's own: the incremental value is a second implementation of the problem's value
// (economy.hpp), held equal to it by check_field.

#pragma once

#include <chrono>
#include <numeric>

#include <calx/random/stream.hpp>

#include "model.hpp"

namespace field {

struct LocalWork {
  double seconds = 0;     // spent on local conditionals and their updates
  std::size_t full = 0;   // full evaluations (rebuilds of the state)
  std::size_t states = 0; // frontage states valued exactly
};

class PlanState {
 public:
  const Model* model = nullptr;
  const Problem* P = nullptr;
  Plan plan;
  std::shared_ptr<const Flows> L;
  std::shared_ptr<const economy::Distances> D;
  std::size_t n = 0;
  std::vector<double> kappa, p, c, A, Pass, Abar;  // n × 168; Pass the passers before they are held at zero
  double passers(std::size_t i) const { return std::max(0.0, Pass[i]); }
  std::vector<char> open;
  std::vector<double> cost;  // per frontage
  double links = 0, value = 0;
  std::vector<std::size_t> segment_of;                 // frontage ↦ its segment's index
  std::vector<std::vector<std::size_t>> members;       // segment ↦ its frontages

  PlanState(const Model& m, Plan x) : model(&m), P(m.problem), plan(std::move(x)), n(m.problem->size()) {
    segment_of.assign(n, 0);
    for (std::size_t s = 0; s < P->segments.size(); s++) {
      members.push_back(P->segments[s].second);
      for (auto f : P->segments[s].second) segment_of[f] = s;
    }
    rebuild();
  }

  double frontage_cost(std::size_t f, std::int64_t k, double weekday_span, double weekend_span) const {
    if (k == economy::vacant) return 0;
    const auto& Q = P->P;
    const auto kk = static_cast<std::size_t>(k);
    const double hours = 5 * weekday_span + 2 * weekend_span, days = 5.0 * (weekday_span > 0) + 2.0 * (weekend_span > 0);
    double c0 = Q.staff[kk] * Q.wage * (hours + days) + Q.outgoings[kk] * P->frontages[f].trading_area;
    if (k != P->frontages[f].kind0) c0 += Q.refit_weekly(P->frontages[f].trading_area);
    return c0;
  }

  void rebuild() {
    L = model->flows(plan.links);
    D = P->distances(plan.links);
    const auto& Q = P->P;
    kappa.assign(n * week, 0.0), p.assign(n * week, 0.0), c.assign(n * week, 0.0), A.assign(n * week, 0.0), Abar.assign(n * week, 0.0);
    open.assign(n * week, 0);
    for (std::size_t f = 0; f < n; f++)
      for (std::size_t h = 0; h < week; h++) open[f * week + h] = plan.open(f, h);
    for (std::size_t h = 0; h < week; h++) {
      for (std::size_t f = 0; f < n; f++) {
        double kap = 0;
        const auto k = static_cast<std::size_t>(std::max<std::int64_t>(0, plan.kind[f]));
        for (const auto& [g, w] : D->near[f])
          if (open[g * week + h]) kap += w * Q.affinity[static_cast<std::size_t>(plan.kind[g])][k];
        kappa[f * week + h] = kap;
        if (open[f * week + h] && economy::trades(plan.kind[f])) p[f * week + h] = chance(plan.kind[f], h, kap);
      }
      for (std::size_t f = 0; f < n; f++) {
        c[f * week + h] = share(f, h, [&](std::size_t g) { return p[g * week + h]; }, [&](std::size_t g) { return plan.kind[g]; });
        A[f * week + h] = worth(plan.kind[f], p[f * week + h], c[f * week + h]);
      }
    }
    Pass = L->P0;
    for (std::size_t h = 0; h < week; h++)
      for (std::size_t f = 0; f < n; f++) {
        const double* r = &model->R[(h * n + f) * n];
        double s = 0;
        for (std::size_t g = 0; g < n; g++) s += r[g] * (p[g * week + h] - model->p0[g * week + h]);
        Pass[f * week + h] = Pass[f * week + h] + s;
      }
    for (std::size_t h = 0; h < week; h++)
      for (std::size_t f = 0; f < n; f++) {
        const double a = A[f * week + h];
        if (a == 0) continue;
        const double* r = &model->R[(h * n + f) * n];
        for (std::size_t g = 0; g < n; g++) Abar[g * week + h] += a * r[g];
      }
    cost.assign(n, 0.0);
    for (std::size_t f = 0; f < n; f++) cost[f] = frontage_cost(f, plan.kind[f], plan.span(f, false), plan.span(f, true));
    links = 0;
    for (std::size_t j = 0; j < P->links.size(); j++)
      if (plan.links[j]) links += Q.link_weekly(P->links[j].length);
    value = -links;
    for (std::size_t i = 0; i < n * week; i++) value += A[i] * passers(i);
    for (double x : cost) value -= x;
  }

  double chance(std::int64_t k, std::size_t h, double kap) const {
    const auto& Q = P->P;
    const auto kk = static_cast<std::size_t>(k);
    return Q.propensity[kk] * Q.period[kk][economy::band(h)] * std::clamp(1 + Q.catalysis * kap, Q.catalysis_floor, Q.catalysis_ceiling);
  }
  double worth(std::int64_t k, double pf, double cf) const {
    if (!(pf > 0)) return 0;
    const auto kk = static_cast<std::size_t>(k);
    return pf * cf * P->P.spend[kk] * P->P.margin[kk];
  }
  // Competition: f's share of its segment's same-kind stops (economy.hpp).
  template <class Pf, class Kf>
  double share(std::size_t f, std::size_t, Pf pof, Kf kof) const {
    const double pf = pof(f);
    if (!(pf > 0)) return 0;
    double none = 1 - pf, sum = pf;
    for (auto g : P->on_segment[f]) {
      const double pg = pof(g);
      if (pg > 0 && kof(g) == kof(f)) none *= 1 - pg, sum += pg;
    }
    return (1 - none) / sum;
  }

  // The exact change of hour h's value when f is put in configuration (is_open, k) at h, the rest as it is. On
  // return, the touched frontages' new p and A are in the scratch arrays (for apply()).
  struct Scratch {
    std::vector<double> pn, an, kn, dP;
    std::vector<char> touched;
    std::vector<std::size_t> changed, affected;
  };
  Scratch scratch() const {
    Scratch s;
    s.pn.assign(n, 0.0), s.an.assign(n, 0.0), s.kn.assign(n, 0.0), s.dP.assign(n, 0.0), s.touched.assign(n, 0);
    return s;
  }

  double delta_hour(std::size_t f, std::size_t h, bool is_open, std::int64_t k, Scratch& s) const {
    const auto& Q = P->P;
    const std::size_t base = h;  // index helper
    (void)base;
    s.changed.clear(), s.affected.clear();
    const auto at = [&](std::size_t g) { return g * week + h; };
    const bool was_open = open[at(f)];
    const std::int64_t k_was = plan.kind[f];
    // New p on N(f).
    const auto touch = [&](std::size_t g) {
      if (!s.touched[g]) s.touched[g] = 1, s.pn[g] = p[at(g)], s.kn[g] = kappa[at(g)], s.an[g] = A[at(g)], s.affected.push_back(g);
    };
    touch(f);
    {
      double kap = 0;
      if (is_open)
        for (const auto& [g, w] : D->near[f])
          if (open[at(g)]) kap += w * Q.affinity[static_cast<std::size_t>(plan.kind[g])][static_cast<std::size_t>(k)];
      s.kn[f] = is_open ? kap : kappa[at(f)];
      s.pn[f] = is_open && economy::trades(k) ? chance(k, h, kap) : 0.0;
    }
    // Every catalysis neighbour's κ follows f (a closed one's too, so that κ stays exact for later changes).
    for (const auto& [g, w] : D->near[f]) {
      const auto kg = static_cast<std::size_t>(std::max<std::int64_t>(0, plan.kind[g]));
      const double dk = (is_open ? w * Q.affinity[static_cast<std::size_t>(k)][kg] : 0.0) - (was_open ? w * Q.affinity[static_cast<std::size_t>(k_was)][kg] : 0.0);
      if (dk == 0) continue;
      touch(g);
      s.kn[g] = kappa[at(g)] + dk;
      if (open[at(g)] && economy::trades(plan.kind[g])) s.pn[g] = chance(plan.kind[g], h, s.kn[g]);
    }
    for (auto g : s.affected)
      if (s.pn[g] != p[at(g)]) s.changed.push_back(g);
    // Competition on the segments of the changed frontages (and f's, whose kind may change).
    const auto kind_of = [&](std::size_t g) { return g == f ? (is_open ? k : k_was) : plan.kind[g]; };
    const auto p_of = [&](std::size_t g) { return s.touched[g] ? s.pn[g] : p[at(g)]; };
    std::vector<std::size_t> segs{segment_of[f]};
    for (auto g : s.changed)
      if (std::find(segs.begin(), segs.end(), segment_of[g]) == segs.end()) segs.push_back(segment_of[g]);
    for (auto sg : segs)
      for (auto g : members[sg]) touch(g);
    // The new A on the affected set; the passers' change everywhere; the change of the hour's value.
    const double* Rh = &model->R[h * n * n];
    for (auto g : s.affected) s.an[g] = worth(kind_of(g), p_of(g), share(g, h, p_of, kind_of));
    double dv = 0;
    for (std::size_t g = 0; g < n; g++) {
      double dP = 0;
      for (auto q : s.changed) dP += Rh[g * n + q] * (s.pn[q] - p[at(q)]);
      s.dP[g] = dP;
      const double a_new = s.touched[g] ? s.an[g] : A[at(g)];
      if (a_new == 0 && A[at(g)] == 0) continue;
      dv += a_new * std::max(0.0, Pass[at(g)] + dP) - A[at(g)] * std::max(0.0, Pass[at(g)]);
    }
    return dv;
  }

  void clear(Scratch& s) const {
    for (auto g : s.affected) s.touched[g] = 0;
    s.affected.clear(), s.changed.clear();
  }

  // Puts f in configuration (is_open, k) at hour h, updating the cached arrays (Pass and Ā exactly).
  void apply_hour(std::size_t f, std::size_t h, bool is_open, std::int64_t k, Scratch& s) {
    const double dv = delta_hour(f, h, is_open, k, s);
    const auto at = [&](std::size_t g) { return g * week + h; };
    if (!s.changed.empty())
      for (std::size_t g = 0; g < n; g++) Pass[at(g)] += s.dP[g];
    const double* Rh = &model->R[h * n * n];
    for (auto g : s.affected) {
      const double da = s.an[g] - A[at(g)];
      if (da != 0)
        for (std::size_t q = 0; q < n; q++) Abar[at(q)] += da * Rh[g * n + q];
      p[at(g)] = s.pn[g], A[at(g)] = s.an[g], kappa[at(g)] = s.kn[g];
    }
    open[at(f)] = is_open;
    value += dv;
    clear(s);
  }

  // One frontage's conditional: every state's value relative to the current one, and a draw from exp(V/T) (T > 0)
  // or its mode (T = 0). Returns the change of value.
  double update(std::size_t f, double T, calx::random::Stream& z, LocalWork& work) {
    const auto T0 = std::chrono::steady_clock::now();
    const auto& Q = P->P;
    Scratch s = scratch();
    // ΔV_h of f closed, and the margins of f open as each kind.
    std::vector<double> closed(week);
    std::array<std::vector<double>, economy::kinds> margin;
    for (std::size_t h = 0; h < week; h++) {
      closed[h] = delta_hour(f, h, false, plan.kind[f], s);
      clear(s);
    }
    for (std::size_t k = 1; k < economy::kinds; k++) {
      margin[k].assign(week, 0.0);
      for (std::size_t h = 0; h < week; h++) {
        margin[k][h] = delta_hour(f, h, true, static_cast<std::int64_t>(k), s) - closed[h];
        clear(s);
      }
    }
    double base = -cost[f];  // relative to the current state: + Σ closed − (the new state's cost − the current's)
    for (double x : closed) base += x;
    // Each kind's best and its log-partition over the hours.
    struct Choice { double value = -1e300; std::int64_t s = 0, l = 0; };
    struct KindTable { std::vector<double> wd, we; double fixed = 0; };
    std::array<KindTable, economy::kinds> tab;
    const auto interval_values = [&](std::size_t k, bool weekend) {
      std::vector<double> out(600);
      std::array<double, 24> M{};
      for (std::size_t d = weekend ? 5 : 0; d < (weekend ? 7u : 5u); d++)
        for (std::size_t t = 0; t < 24; t++) M[t] += margin[k][24 * d + t];
      const double days = weekend ? 2 : 5;
      const double per = Q.staff[k] * Q.wage * days;
      for (std::size_t st = 0; st < 24; st++) {
        double run = 0;
        for (std::size_t len = 0; len <= 24; len++) {
          if (len > 0 && st + len - 1 < 24) run += M[st + len - 1];
          const double span = len > 0 ? static_cast<double>(std::min<std::size_t>(24, st + len) - st) : 0.0;
          out[st * 25 + len] = len > 0 ? run - per * (span + 1) : 0.0;
        }
      }
      return out;
    };
    for (std::size_t k = 1; k < economy::kinds; k++) {
      tab[k].wd = interval_values(k, false);
      tab[k].we = interval_values(k, true);
      tab[k].fixed = Q.outgoings[k] * P->frontages[f].trading_area + (static_cast<std::int64_t>(k) != P->frontages[f].kind0 ? Q.refit_weekly(P->frontages[f].trading_area) : 0.0);
    }
    const auto lse = [](const std::vector<double>& v, double T) {
      double mx = -1e300;
      for (double x : v) mx = std::max(mx, x);
      double s = 0;
      for (double x : v) s += std::exp((x - mx) / T);
      return mx + T * std::log(s);
    };
    std::int64_t k_new = plan.kind[f];
    std::int64_t ws = plan.hours[4 * f], wl = plan.hours[4 * f + 1], es = plan.hours[4 * f + 2], el = plan.hours[4 * f + 3];
    if (T > 0) {
      std::array<double, economy::kinds> F{};
      F[0] = base + 2 * T * std::log(600.0);
      for (std::size_t k = 1; k < economy::kinds; k++) F[k] = base - tab[k].fixed + lse(tab[k].wd, T) + lse(tab[k].we, T);
      double mx = *std::max_element(F.begin(), F.end()), Z = 0;
      for (double x : F) Z += std::exp((x - mx) / T);
      double u = z.uniform() * Z;
      std::size_t pick = economy::kinds - 1;
      for (std::size_t k = 0; k < economy::kinds; k++) {
        u -= std::exp((F[k] - mx) / T);
        if (u <= 0) { pick = k; break; }
      }
      k_new = static_cast<std::int64_t>(pick);
      const auto draw = [&](const std::vector<double>& v) {
        double m2 = *std::max_element(v.begin(), v.end()), z2 = 0;
        for (double x : v) z2 += std::exp((x - m2) / T);
        double w = z.uniform() * z2;
        for (std::size_t i = 0; i < v.size(); i++) {
          w -= std::exp((v[i] - m2) / T);
          if (w <= 0) return i;
        }
        return v.size() - 1;
      };
      if (pick == 0) {
        ws = static_cast<std::int64_t>(z.below(24)), wl = static_cast<std::int64_t>(z.below(25)), es = static_cast<std::int64_t>(z.below(24)), el = static_cast<std::int64_t>(z.below(25));
      } else {
        const auto a = draw(tab[pick].wd), b = draw(tab[pick].we);
        ws = static_cast<std::int64_t>(a / 25), wl = static_cast<std::int64_t>(a % 25), es = static_cast<std::int64_t>(b / 25), el = static_cast<std::int64_t>(b % 25);
      }
    } else {
      // The mode; the current state on a tie.
      const double current = 0;  // the current state's relative value is 0 by construction
      double best = current + 1e-9;
      bool moved = false;
      if (base > best && plan.kind[f] != economy::vacant) best = base, k_new = economy::vacant, moved = true;
      for (std::size_t k = 1; k < economy::kinds; k++) {
        const auto iw = static_cast<std::size_t>(std::max_element(tab[k].wd.begin(), tab[k].wd.end()) - tab[k].wd.begin());
        const auto ie = static_cast<std::size_t>(std::max_element(tab[k].we.begin(), tab[k].we.end()) - tab[k].we.begin());
        const double v = base - tab[k].fixed + tab[k].wd[iw] + tab[k].we[ie];
        if (v > best) best = v, k_new = static_cast<std::int64_t>(k), ws = static_cast<std::int64_t>(iw / 25), wl = static_cast<std::int64_t>(iw % 25),
                      es = static_cast<std::int64_t>(ie / 25), el = static_cast<std::int64_t>(ie % 25), moved = true;
      }
      if (!moved) {
        work.states += 6 * 600 * 2;
        work.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
        return 0;
      }
    }
    // Apply the new state hour by hour.
    const double before = value;
    Plan next = plan;
    next.kind[f] = k_new, next.hours[4 * f] = ws, next.hours[4 * f + 1] = wl, next.hours[4 * f + 2] = es, next.hours[4 * f + 3] = el;
    for (std::size_t h = 0; h < week; h++) {
      const bool o = next.open(f, h);
      if (o == static_cast<bool>(open[f * week + h]) && (!o || k_new == plan.kind[f])) continue;
      apply_hour(f, h, o, o ? k_new : plan.kind[f], s);
    }
    // The kind may change where f stays closed: it then matters to no hour; record it.
    plan = next;
    const double c_new = frontage_cost(f, k_new, plan.span(f, false), plan.span(f, true));
    value += cost[f] - c_new;
    cost[f] = c_new;
    work.states += 6 * 600 * 2;
    work.seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - T0).count();
    return value - before;
  }
};

}  // namespace field
