// Kepler's third law rediscovered from JPL's ephemeris DE440s, by genetic programming typed by physical dimension
// (heuresis programs) and, for comparison, by the same search untyped.
//
// Part 1, the data. For the barycentres of the planetary systems (Mercury to Pluto; the Earth as the Earth–Moon
// barycentre) and for the Moon about the Earth, read through the astronomy domain's SPK reader:
//  - the state of each planet in Jacobi coordinates, relative to the barycentre of the Sun and the planets inside
//    its orbit, with μ = G(M_sun + Σ_inside m_i + m_planet) (Murray and Dermott 1999, Solar System Dynamics, §6.2);
//    and, for comparison, heliocentric, with μ = G(M_sun + m_planet). The Moon is geocentric, μ = G(M_E + M_M). GM
//    are DE440's (Park et al. 2021, AJ 161, 105; the astronomy domain's bodies.hpp);
//  - the osculating semi-major axis a = 1/(2/r − v²/μ) (Murray and Dermott 1999, eq. 2.134) at every sample of the
//    span, and its time mean;
//  - the sidereal period, from the ephemeris alone: the ecliptic longitude (the J2000 ecliptic, mean obliquity
//    84381.406″, IAU 2006: Capitaine, Wallace and Chapront 2003, A&A 412, 567; IERS Conventions 2010, eq. 5.40),
//    unwrapped, gives 2π over the slope of its least-squares line where the span holds at least 20 revolutions,
//    and otherwise the mean time for it to advance by 2π from up to 400 starts spread over the span, each return
//    found by Brent's method (calx solvers). The returns are biased where the apsides turn fast (the Moon's in
//    8.85 years), the line where there are few revolutions; each is used where it is not.
// Part 2, the search. Cases (a, GM, e, T) in two systems of units: the astronomical unit (IAU 2012 Resolution B2,
// 149 597 870 700 m) and the Julian year (365.25 d of 86 400 s), Kepler's own, in which the Sun's GM is nearly 4π²;
// and SI. The mean osculating eccentricity e (dimensionless) is offered as a variable the law does not need: without
// it, a and GM admit no dimensionless group, so every typed term of the sort of a time would be C √(a³/GM) (the Π
// theorem) and the typed search would only choose C. Two problems:
//  - "two centres": the nine planets and the Moon; the variables a (a length), GM (L³ T⁻²) and e; the target T;
//  - "one centre": the nine planets alone; the variables a and e, and an ephemeral constant k of the sort of a
//    gravitational parameter (L³ T⁻²) drawn in [1, 100] au³ yr⁻² (its range converted with its unit), so that a
//    typed term of the sort of a time can be formed at all.
// Each is searched by GP typed by dimension and by the same GP with every symbol dimensionless (the untyped k keeps
// the range [1, 100] in either system, since nothing tells it that it is a GM), 30 seeds each; the error is the
// RMS relative error of T. A term is the law when it has the law's form (two centres: equal to 2π √(a³/GM) at probe
// points to 10⁻⁹; one centre: proportional to a^{3/2} to 10⁻⁹), its error is below 1 % on the data, and it stays
// below 1 % when the units change (constants converted by their declared dimension): a law, not a coincidence of
// the units. Every term met is kept by size; the non-dominated set of (size, error) over all seeds is the front of
// formulas (heuresis multiobjective's dominance).
//
// Usage: kepler <path to de440s.bsp> <results directory>

#include <astronomy/ephemerides/bodies.hpp>
#include <astronomy/ephemerides/spk.hpp>
#include <calx/solvers/roots.hpp>
#include <heuresis/multiobjective/dominance.hpp>
#include <heuresis/programs/programs.hpp>
#include <heuresis/search/search.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <string>
#include <thread>
#include <vector>

namespace ae = astronomy::ephemerides;
namespace hp = heuresis::programs;
namespace hs = heuresis::search;
namespace hm = heuresis::multiobjective;
namespace dim = physis::units::dim;
using physis::units::Dimension;

constexpr double PI = std::numbers::pi;
constexpr double AU = 149597870700.0;                          // m, IAU 2012 Resolution B2
constexpr double DAY = 86400.0;                                // s
constexpr double YEAR = 365.25 * DAY;                          // s, the Julian year
constexpr double OBLIQUITY = 84381.406 / 3600.0 * PI / 180.0;  // rad, IAU 2006 mean obliquity at J2000.0

// ---- part 1: orbits from the ephemeris --------------------------------------------------------------------------

struct Target {
  std::string name;
  int target;
  int centre;                 // the Moon's: 399; a planet's: 0 (Jacobi) or 10 (heliocentric)
  std::vector<int> interior;  // Jacobi: the Sun and the planets' barycentres inside the orbit
  double mu;                  // m³ s⁻²
  double step_days;
};

struct Orbit {
  std::string name, frame;
  int target = 0;
  double mu = 0;
  double a_mean = 0, a_sd = 0, a_min = 0, a_max = 0;  // m
  double e_mean = 0, e_sd = 0;                         // the osculating eccentricity
  double period_fit = 0, period_returns = 0, period_returns_sd = 0, period = 0;  // s
  std::string method;
  int starts = 0;
  double revolutions = 0;
};

double wrap(double x) {
  while (x > PI) x -= 2 * PI;
  while (x <= -PI) x += 2 * PI;
  return x;
}

double gm_of(int id) {
  const auto g = [](const char* key) { return ae::body(key)->gm; };
  switch (id) {
    case 10: return g("sun");
    case 1: return g("mercury");
    case 2: return g("venus");
    case 3: return g("earth") + g("moon");
    case 4: return g("mars");
    case 5: return g("jupiter");
    case 6: return g("saturn");
    case 7: return g("uranus");
    case 8: return g("neptune");
    case 9: return g("pluto");
  }
  throw std::invalid_argument("no GM for body " + std::to_string(id));
}

// The target's state in its frame: relative to its centre, or to the barycentre of the bodies inside its orbit.
ae::State relative(const ae::Spk& spk, const Target& b, double et) {
  if (b.interior.empty()) return spk.state(b.target, b.centre, et);
  ae::State p = spk.barycentric(b.target, et), c;
  double m = 0;
  for (int i : b.interior) {
    const auto s = spk.barycentric(i, et);
    const double g = gm_of(i);
    for (int k = 0; k < 3; k++) c.position[k] += g * s.position[k], c.velocity[k] += g * s.velocity[k];
    m += g;
  }
  for (int k = 0; k < 3; k++) p.position[k] -= c.position[k] / m, p.velocity[k] -= c.velocity[k] / m;
  return p;
}

double longitude(const ae::State& s) {
  const double ye = std::cos(OBLIQUITY) * s.position[1] + std::sin(OBLIQUITY) * s.position[2];
  return std::atan2(ye, s.position[0]);
}

Orbit measure(const ae::Spk& spk, const Target& b, double et0, double et1) {
  Orbit o;
  o.name = b.name, o.target = b.target, o.mu = b.mu;
  o.frame = b.target == 301 ? "geocentric" : b.interior.empty() ? "heliocentric" : "jacobi";
  const double h = b.step_days * DAY;
  const std::size_t n = std::size_t((et1 - et0) / h) + 1;
  std::vector<double> t(n), raw(n), lam(n);
  double sa = 0, sa2 = 0, se = 0, se2 = 0;
  o.a_min = std::numeric_limits<double>::infinity(), o.a_max = -o.a_min;
  for (std::size_t i = 0; i < n; i++) {
    t[i] = et0 + double(i) * h;
    const auto s = relative(spk, b, t[i]);
    const double r = std::hypot(s.position[0], s.position[1], s.position[2]);
    const double v2 = s.velocity[0] * s.velocity[0] + s.velocity[1] * s.velocity[1] + s.velocity[2] * s.velocity[2];
    const double a = 1 / (2 / r - v2 / b.mu);
    sa += a, sa2 += a * a, o.a_min = std::min(o.a_min, a), o.a_max = std::max(o.a_max, a);
    // the eccentricity vector ((v² − μ/r) r − (r·v) v)/μ (Murray and Dermott 1999, eq. 2.133)
    const double rv = s.position[0] * s.velocity[0] + s.position[1] * s.velocity[1] + s.position[2] * s.velocity[2];
    double ev[3];
    for (int k = 0; k < 3; k++) ev[k] = ((v2 - b.mu / r) * s.position[k] - rv * s.velocity[k]) / b.mu;
    const double ecc = std::hypot(ev[0], ev[1], ev[2]);
    se += ecc, se2 += ecc * ecc;
    raw[i] = longitude(s);
    lam[i] = i == 0 ? raw[0] : lam[i - 1] + wrap(raw[i] - raw[i - 1]);
    if (i > 0 && !(lam[i] > lam[i - 1])) throw std::runtime_error(b.name + ": the longitude does not advance at a sample; refine the step");
  }
  o.a_mean = sa / double(n);
  o.a_sd = std::sqrt(std::max(sa2 / double(n) - o.a_mean * o.a_mean, 0.0));
  o.e_mean = se / double(n);
  o.e_sd = std::sqrt(std::max(se2 / double(n) - o.e_mean * o.e_mean, 0.0));
  o.revolutions = (lam[n - 1] - lam[0]) / (2 * PI);
  // The least-squares line through the unwrapped longitude.
  {
    double tm = 0, lm = 0;
    for (std::size_t i = 0; i < n; i++) tm += t[i], lm += lam[i];
    tm /= double(n), lm /= double(n);
    double stl = 0, stt = 0;
    for (std::size_t i = 0; i < n; i++) stl += (t[i] - tm) * (lam[i] - lm), stt += (t[i] - tm) * (t[i] - tm);
    o.period_fit = 2 * PI / (stl / stt);
  }
  // The returns: from starts whose next revolution ends within the span.
  std::size_t last = 0;
  for (std::size_t i = 0; i < n; i++)
    if (lam[n - 1] >= lam[i] + 2 * PI + 1e-6) last = i;
  const std::size_t count = std::min<std::size_t>(400, last + 1);
  std::vector<double> periods;
  std::size_t j = 1;
  for (std::size_t k = 0; k < count; k++) {
    const std::size_t i = count == 1 ? 0 : std::size_t(std::llround(double(k) * double(last) / double(count - 1)));
    const double goal = lam[i] + 2 * PI;
    if (j <= i) j = i + 1;
    while (j < n && lam[j] < goal) j++;
    if (j >= n) break;
    const double base = lam[j - 1], braw = raw[j - 1];
    auto f = [&](double tt) { return base + wrap(longitude(relative(spk, b, tt)) - braw) - goal; };
    periods.push_back(calx::solvers::brent(f, t[j - 1], t[j], 1e-3).x - t[i]);
  }
  double s = 0, s2 = 0;
  for (double p : periods) s += p, s2 += p * p;
  o.starts = int(periods.size());
  o.period_returns = s / double(periods.size());
  o.period_returns_sd = periods.size() > 1 ? std::sqrt(std::max((s2 - s * s / double(periods.size())) / double(periods.size() - 1), 0.0)) : 0.0;
  const bool fit = o.revolutions >= 20;
  o.period = fit ? o.period_fit : o.period_returns;
  o.method = fit ? "line" : "returns";
  return o;
}

double kepler_period(const Orbit& o) { return 2 * PI * std::sqrt(o.a_mean * o.a_mean * o.a_mean / o.mu); }

// ---- part 2: the search -----------------------------------------------------------------------------------------

struct Case {
  double a, gm, e, T;  // in the variant's units (e dimensionless)
};

const Dimension GM_DIM = dim::length.pow(3) / dim::time.pow(2);

struct Units {
  std::string name;
  double length, time;  // the unit's size in SI: (AU, YEAR) or (1, 1)
};
const Units KEPLER{"au, yr", AU, YEAR}, SI{"SI", 1, 1};

struct Variant {
  std::string problem;  // "two centres" or "one centre"
  bool two_centres;
  Units units;
  std::vector<Case> cases;
};

double factor(const Dimension& d, const Units& u) { return std::pow(u.length, d.e[1]) * std::pow(u.time, d.e[2]); }

hp::Signature signature(bool typed, bool two_centres, const Units& u) {
  hp::Signature s(4);
  for (auto p : {hp::library::add(), hp::library::subtract(), hp::library::multiply(), hp::library::divide(), hp::library::square(), hp::library::cube(), hp::library::root()})
    s.add(p);
  const auto sort = [&](Dimension d) { return typed ? hp::of(d) : hp::dimensionless(); };
  s.add(hp::Variable{"a", sort(dim::length)});
  if (two_centres) s.add(hp::Variable{"GM", sort(GM_DIM)});
  s.add(hp::Variable{"e", hp::dimensionless()});
  s.add(hp::Constant{"2", hp::dimensionless(), 2});
  s.add(hp::Constant{"pi", hp::dimensionless(), PI});
  s.add(hp::Constant{"c", hp::dimensionless(), 0, true, -1, 1});
  if (!two_centres) {
    // k in [1, 100] au³ yr⁻²: converted to the variant's units when typed, a bare number when not.
    const double f = typed ? factor(GM_DIM, KEPLER) / factor(GM_DIM, u) : 1.0;
    s.add(hp::Constant{"k", sort(GM_DIM), 0, true, 1 * f, 100 * f});
  }
  return s;
}

std::vector<double> inputs_of(const Case& c, bool two_centres) { return two_centres ? std::vector<double>{c.a, c.gm, c.e} : std::vector<double>{c.a, c.e}; }

// The RMS relative error of T and the per-case errors; a failure where a value is not finite.
hs::Evaluation score(const hp::Signature& sig, const hp::Tree& t, const Variant& v) {
  std::vector<double> e(v.cases.size());
  double s = 0;
  for (std::size_t i = 0; i < v.cases.size(); i++) {
    const double y = hp::interpret(sig, t, inputs_of(v.cases[i], v.two_centres));
    if (!std::isfinite(y)) return hs::failure("not finite on case " + std::to_string(i));
    const double r = (y - v.cases[i].T) / v.cases[i].T;
    e[i] = std::fabs(r), s += r * r;
  }
  return hs::evaluation({s / double(v.cases.size())}, {}, std::move(e));
}

// The RMS relative error of a term after a change of units from `from` to `to`: the inputs converted by their
// physical dimension, the term's constants by their declared sort (an untyped constant is a bare number).
double error_in(const hp::Signature& sig, hp::Tree t, const Variant& v, const Units& to) {
  for (auto& n : t)
    if (sig.is_constant(n.symbol)) n.value *= factor(n.sort.dimension, v.units) / factor(n.sort.dimension, to);
  const double fa = factor(dim::length, v.units) / factor(dim::length, to), fg = factor(GM_DIM, v.units) / factor(GM_DIM, to), ft = factor(dim::time, v.units) / factor(dim::time, to);
  double s = 0;
  for (const auto& c : v.cases) {
    const double y = hp::interpret(sig, t, v.two_centres ? std::vector<double>{c.a * fa, c.gm * fg, c.e} : std::vector<double>{c.a * fa, c.e});
    if (!std::isfinite(y)) return std::numeric_limits<double>::infinity();
    const double r = (y - c.T * ft) / (c.T * ft);
    s += r * r;
  }
  return std::sqrt(s / double(v.cases.size()));
}

const Units& other(const Units& u) { return u.name == SI.name ? KEPLER : SI; }

// The law's form at probe points (in the variant's units), its error on the data, and its error in the other units.
bool is_law(const hp::Signature& sig, const hp::Tree& t, const Variant& v, double rms) {
  if (!(rms <= 0.01)) return false;
  static const double pa[] = {0.003, 0.4, 1.7, 9.0, 31.0, 77.0};
  static const double pg[] = {1.3e-4, 39.0, 0.02, 3.0, 400.0, 1.0};
  static const double pe[] = {0.0, 0.055, 0.21, 0.5, 0.01, 0.8};
  const double fa = AU / v.units.length, fg = factor(GM_DIM, KEPLER) / factor(GM_DIM, v.units);
  if (v.two_centres) {
    for (int i = 0; i < 6; i++) {
      const double a = pa[i] * fa, g = pg[i] * fg;
      const double y = hp::interpret(sig, t, std::vector<double>{a, g, pe[i]});
      if (!(std::fabs(y / (2 * PI * std::sqrt(a * a * a / g)) - 1) <= 1e-9)) return false;
    }
  } else {
    const double y0 = hp::interpret(sig, t, std::vector<double>{fa, 0.1});
    if (!(std::isfinite(y0) && y0 > 0)) return false;
    for (int i = 0; i < 6; i++) {
      const double y = hp::interpret(sig, t, std::vector<double>{pa[i] * fa, pe[i]});
      if (!(std::fabs(y / (y0 * std::pow(pa[i], 1.5)) - 1) <= 1e-9)) return false;
    }
  }
  return error_in(sig, t, v, other(v.units)) <= 0.01;
}

// ---- formulas written for the page (matplotlib's mathtext) ------------------------------------------------------

std::string number(double v) {
  char buf[64];
  if (std::fabs(v - 2 * PI) < 1e-12) return "2\\pi";
  if (std::fabs(v - PI) < 1e-12) return "\\pi";
  if (std::fabs(v - 4 * PI * PI) < 1e-10) return "4\\pi^2";
  if (std::fabs(v - std::round(v)) < 1e-12 && std::fabs(v) < 1e6) std::snprintf(buf, sizeof buf, "%.0f", v);
  else std::snprintf(buf, sizeof buf, "%.4g", v);
  std::string s = buf;
  const auto e = s.find('e');
  if (e != std::string::npos) s = s.substr(0, e) + "\\times10^{" + std::to_string(std::stoi(s.substr(e + 1))) + "}";
  return s;
}

struct Printer {
  const hp::Signature& sig;
  const hp::Tree& t;
  bool constant_only(std::size_t i) const {
    const std::size_t e = hp::subtree_end(sig, t, i);
    for (std::size_t k = i; k < e; k++)
      if (sig.is_variable(t[k].symbol)) return false;
    return true;
  }
  double fold(std::size_t i) const {
    const std::size_t e = hp::subtree_end(sig, t, i);
    const hp::Tree sub(t.begin() + std::ptrdiff_t(i), t.begin() + std::ptrdiff_t(e));
    return hp::interpret(sig, sub, std::vector<double>(4, 0.0));
  }
  // The text and its precedence: 0 a sum, 1 a product or quotient, 2 an atom or a power.
  std::pair<std::string, int> at(std::size_t i) const {
    const auto s = t[i].symbol;
    if (sig.is_variable(s)) return {sig.name(s), 2};
    if (constant_only(i)) {
      const double v = fold(i);
      return {number(v), v < 0 ? 0 : 2};
    }
    const std::string name = sig.name(s);
    const std::size_t a = i + 1;
    if (sig.arity(s) == 1) {
      const auto [x, p] = at(a);
      if (name == "sqrt") return {"\\sqrt{" + x + "}", 2};
      const bool bare = p == 2 && x.find("\\sqrt") == std::string::npos && x.find('^') == std::string::npos;
      return {(bare ? x : "(" + x + ")") + (name == "sq" ? "^2" : "^3"), 2};
    }
    const std::size_t b = hp::subtree_end(sig, t, a);
    const auto [x, px] = at(a);
    const auto [y, py] = at(b);
    if (name == "+") return {x + " + " + y, 0};
    if (name == "-") return {x + " - " + (py == 0 ? "(" + y + ")" : y), 0};
    if (name == "/") return {"\\frac{" + x + "}{" + y + "}", 1};
    return {(px == 0 ? "(" + x + ")" : x) + "\\," + (py == 0 ? "(" + y + ")" : y), 1};
  }
  std::string text() const { return at(0).first; }
};

// ---- runs ---------------------------------------------------------------------------------------------------------

struct Kept {
  double error = std::numeric_limits<double>::infinity();
  hp::Tree tree;
};

struct RunResult {
  const Variant* variant = nullptr;
  bool typed = false;
  std::uint64_t seed = 0;
  double law_at = -1;  // the evaluation at which the law first appeared
  double best_error = 0, best_error_other = 0;
  std::size_t best_size = 0;
  bool best_is_law = false;
  std::string best_prefix, best_tex;
  std::vector<std::pair<double, double>> progress;  // (evaluations, RMS error)
  std::map<std::size_t, Kept> by_size;
  std::size_t failures = 0, evaluations = 0;
};

RunResult search_once(const Variant& v, bool typed, std::uint64_t seed, double budget) {
  const hp::Signature sig = signature(typed, v.two_centres, v.units);
  hp::GpOptions o;
  o.target = typed ? hp::of(dim::time) : hp::dimensionless();
  o.population = 500;
  o.min_depth = 2, o.max_depth = 5;
  o.limits = hp::Limits{10, 40, 3};
  o.crossover = 0.8, o.subtree = 0.1, o.point = 0.1;
  o.point_rate = 0.15;
  o.tournament = 7, o.elites = 2;
  hp::GeneticProgramming gp(sig, o, hs::Draws(seed));
  RunResult r;
  r.variant = &v, r.typed = typed, r.seed = seed;
  auto objective = [&](const hp::Tree& t) {
    r.evaluations++;
    auto e = score(sig, t, v);
    if (e.failed()) {
      r.failures++;
      return e;
    }
    const double rms = std::sqrt(e.objectives[0]);
    Kept& k = r.by_size[t.size()];
    if (rms < k.error) k.error = rms, k.tree = t;
    if (r.law_at < 0 && is_law(sig, t, v, rms)) r.law_at = double(r.evaluations);
    return e;
  };
  const auto out = hs::run(gp, objective, hs::Budget{budget});
  for (const auto& p : out.progress) r.progress.push_back({p.evaluations, std::sqrt(p.value)});
  if (out.best) {
    const hp::Tree& t = out.best->genotype;
    r.best_error = std::sqrt(out.best->evaluation.value());
    r.best_error_other = error_in(sig, t, v, other(v.units));
    r.best_is_law = is_law(sig, t, v, r.best_error);
    r.best_size = t.size();
    r.best_prefix = hp::to_string(sig, t);
    r.best_tex = Printer{sig, t}.text();
  }
  return r;
}

std::string json_string(const std::string& s) {
  std::string o = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') o += '\\';
    o += c;
  }
  return o + "\"";
}

std::string num(double v) {
  if (!std::isfinite(v)) return "null";
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "usage: kepler <de440s.bsp> <results directory>\n";
    return 2;
  }
  const std::filesystem::path out_dir = argv[2];
  std::filesystem::create_directories(out_dir);
  const auto spk = ae::Spk::read(argv[1]);
  const auto [c0, c1] = spk.coverage();
  const double et0 = c0 + DAY, et1 = c1 - DAY;
  std::fprintf(stderr, "DE440s: ET %.6e to %.6e s past J2000 (%.2f Julian years)\n", et0, et1, (et1 - et0) / YEAR);

  // The bodies: each planet in Jacobi coordinates and heliocentric; the Moon geocentric.
  const std::vector<std::pair<std::string, double>> planets = {{"Mercury", 0.5}, {"Venus", 1}, {"Earth", 1}, {"Mars", 2}, {"Jupiter", 10},
                                                               {"Saturn", 20}, {"Uranus", 30}, {"Neptune", 30}, {"Pluto", 30}};
  std::vector<Orbit> jacobi, helio;
  const auto t_start = std::chrono::steady_clock::now();
  for (int k = 0; k < int(planets.size()); k++) {
    const int id = k + 1;
    if (!spk.has(id)) {
      std::cerr << planets[k].first << ": not in the kernel, skipped\n";
      continue;
    }
    std::vector<int> inside{10};
    double m = gm_of(10);
    for (int i = 1; i < id; i++) inside.push_back(i), m += gm_of(i);
    jacobi.push_back(measure(spk, Target{planets[k].first, id, 0, inside, m + gm_of(id), planets[k].second}, et0, et1));
    helio.push_back(measure(spk, Target{planets[k].first, id, 10, {}, gm_of(10) + gm_of(id), planets[k].second}, et0, et1));
  }
  const Orbit moon = measure(spk, Target{"Moon", 301, 399, {}, ae::body("earth")->gm + ae::body("moon")->gm, 0.25}, et0, et1);
  jacobi.push_back(moon);
  std::fprintf(stderr, "orbits measured in %.1f s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count());
  {
    std::ofstream f(out_dir / "orbits.csv");
    f << "body,frame,mu_m3s2,a_mean_m,a_sd_m,a_min_m,a_max_m,e_mean,e_sd,period_s,method,period_line_s,period_returns_s,period_returns_sd_s,starts,revolutions,kepler_period_s,residual,"
         "helio_a_mean_m,helio_mu_m3s2,helio_residual\n";
    for (std::size_t i = 0; i < jacobi.size(); i++) {
      const auto& o = jacobi[i];
      f << o.name << ',' << o.frame << ',' << num(o.mu) << ',' << num(o.a_mean) << ',' << num(o.a_sd) << ',' << num(o.a_min) << ',' << num(o.a_max) << ',' << num(o.e_mean) << ',' << num(o.e_sd) << ',' << num(o.period) << ','
        << o.method << ',' << num(o.period_fit) << ',' << num(o.period_returns) << ',' << num(o.period_returns_sd) << ',' << o.starts << ',' << num(o.revolutions) << ','
        << num(kepler_period(o)) << ',' << num(o.period / kepler_period(o) - 1) << ',';
      if (i < helio.size()) {
        const auto& h = helio[i];
        f << num(h.a_mean) << ',' << num(h.mu) << ',' << num(h.period / kepler_period(h) - 1) << '\n';
      } else {
        f << ",,\n";
      }
      std::fprintf(stderr, "%-8s %-10s a %.9e m (sd %.1e)  T %.9e s by %-7s (line %.9e, returns %.9e)  T/T_K - 1 = %+.2e%s\n", o.name.c_str(), o.frame.c_str(), o.a_mean,
                   o.a_sd / o.a_mean, o.period, o.method.c_str(), o.period_fit, o.period_returns, o.period / kepler_period(o) - 1,
                   i < helio.size() ? ("  (heliocentric " + std::to_string(helio[i].period / kepler_period(helio[i]) - 1) + ")").c_str() : "");
    }
  }

  // The variants: two problems in two systems of units.
  std::vector<Variant> variants;
  for (const Units* u : {&KEPLER, &SI})
    for (bool two : {true, false}) {
      Variant v{two ? "two centres" : "one centre", two, *u, {}};
      for (const auto& o : jacobi) {
        if (!two && o.name == "Moon") continue;
        v.cases.push_back({o.a_mean / u->length, o.mu / factor(GM_DIM, *u), o.e_mean, o.period / u->time});
      }
      variants.push_back(std::move(v));
    }
  // The law's errors on the data, for reference: with each body's μ, and with the Sun's GM alone.
  double law_two = 0, law_one = 0;
  const double gm_sun = gm_of(10);
  int n_two = 0, n_one = 0;
  for (const auto& o : jacobi) {
    law_two += std::pow(o.period / kepler_period(o) - 1, 2), n_two++;
    if (o.name != "Moon") law_one += std::pow(o.period / (2 * PI * std::sqrt(o.a_mean * o.a_mean * o.a_mean / gm_sun)) - 1, 2), n_one++;
  }
  law_two = std::sqrt(law_two / n_two), law_one = std::sqrt(law_one / n_one);
  std::fprintf(stderr, "the law's RMS relative error: two centres %.3e (each body's μ); one centre %.3e (the Sun's GM, %.6f au³ yr⁻²)\n", law_two, law_one, gm_sun / factor(GM_DIM, KEPLER));

  // The runs: 4 variants × typed and untyped × 30 seeds, on the machine's threads.
  const double budget = 30000;
  const int seeds = 30;
  struct Job {
    const Variant* v;
    bool typed;
    std::uint64_t seed;
  };
  std::vector<Job> jobs;
  for (const auto& v : variants)
    for (bool typed : {true, false})
      for (int s = 0; s < seeds; s++) jobs.push_back({&v, typed, std::uint64_t(1000 + s)});
  std::vector<RunResult> results(jobs.size());
  std::atomic<std::size_t> next{0};
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<std::thread> pool;
  const unsigned threads = std::max(1u, std::min(8u, std::thread::hardware_concurrency()));
  for (unsigned w = 0; w < threads; w++)
    pool.emplace_back([&] {
      for (std::size_t i; (i = next++) < jobs.size();) results[i] = search_once(*jobs[i].v, jobs[i].typed, jobs[i].seed, budget);
    });
  for (auto& th : pool) th.join();
  std::fprintf(stderr, "%zu runs of %.0f evaluations in %.1f s\n", jobs.size(), budget, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());

  {
    std::ofstream f(out_dir / "gp_runs.json");
    f << "{\n  \"budget\": " << num(budget) << ",\n  \"population\": 500,\n  \"seeds\": " << seeds << ",\n  \"law_error\": {\"two centres\": " << num(law_two)
      << ", \"one centre\": " << num(law_one) << "},\n  \"gm_sun_au3_yr2\": " << num(gm_sun / factor(GM_DIM, KEPLER)) << ",\n  \"runs\": [\n";
    for (std::size_t i = 0; i < results.size(); i++) {
      const auto& r = results[i];
      f << "    {\"problem\": " << json_string(r.variant->problem) << ", \"units\": " << json_string(r.variant->units.name) << ", \"typed\": " << (r.typed ? "true" : "false")
        << ", \"seed\": " << r.seed << ", \"law_at\": " << num(r.law_at) << ", \"best_error\": " << num(r.best_error) << ", \"best_error_other_units\": " << num(r.best_error_other)
        << ", \"best_is_law\": " << (r.best_is_law ? "true" : "false") << ", \"best_size\": " << r.best_size << ", \"best\": " << json_string(r.best_prefix)
        << ", \"best_tex\": " << json_string(r.best_tex) << ", \"failures\": " << r.failures << ", \"evaluations\": " << r.evaluations << ", \"progress\": [";
      for (std::size_t k = 0; k < r.progress.size(); k++) f << (k ? ", " : "") << "[" << num(r.progress[k].first) << ", " << num(r.progress[k].second) << "]";
      f << "]}" << (i + 1 < results.size() ? "," : "") << "\n";
    }
    f << "  ]\n}\n";
  }

  // The fronts: for each variant and typing, the best term of each size over all seeds, and their non-dominated set.
  std::ofstream ff(out_dir / "gp_fronts.json");
  ff << "[\n";
  bool first_row = true;
  for (const auto& v : variants)
    for (bool typed : {true, false}) {
      const hp::Signature sig = signature(typed, v.two_centres, v.units);
      std::map<std::size_t, Kept> best;
      for (const auto& r : results) {
        if (r.variant != &v || r.typed != typed) continue;
        for (const auto& [size, k] : r.by_size)
          if (k.error < best[size].error) best[size] = k;
      }
      std::vector<hm::Point> points;
      std::vector<std::size_t> sizes;
      for (const auto& [size, k] : best) points.push_back({double(size), k.error}), sizes.push_back(size);
      for (std::size_t i : hm::non_dominated(points)) {
        const Kept& k = best[sizes[i]];
        const double e2 = error_in(sig, k.tree, v, other(v.units));
        const bool law = is_law(sig, k.tree, v, k.error);
        const std::string prefix = hp::to_string(sig, k.tree);
        // The constant the formula implies: for two centres C = T/√(a³/GM) at a = 1 au, GM the Sun's, e = 0.05 (2π for
        // the law); for one centre the GM it implies, 4π² a³/T² at a = 1 au, in au³ yr⁻².
        const double fa = AU / v.units.length, fg = factor(GM_DIM, KEPLER) / factor(GM_DIM, v.units), ft = YEAR / v.units.time;
        const double gs = gm_sun / factor(GM_DIM, KEPLER);
        double implied = 0;
        if (v.two_centres) implied = hp::interpret(sig, k.tree, std::vector<double>{fa, gs * fg, 0.05}) / ft / std::sqrt(1 / gs);
        else implied = 4 * PI * PI / std::pow(hp::interpret(sig, k.tree, std::vector<double>{fa, 0.05}) / ft, 2);
        std::fprintf(stderr, "%-11s %-6s %-7s size %2zu  error %.3e  other units %.3e %s %s\n", v.problem.c_str(), v.units.name.c_str(), typed ? "typed" : "untyped", sizes[i], k.error,
                     e2, law ? "LAW" : "   ", prefix.c_str());
        ff << (first_row ? "" : ",\n") << "  {\"problem\": " << json_string(v.problem) << ", \"units\": " << json_string(v.units.name) << ", \"typed\": " << (typed ? "true" : "false")
           << ", \"size\": " << sizes[i] << ", \"error\": " << num(k.error) << ", \"error_other_units\": " << num(e2) << ", \"law\": " << (law ? "true" : "false") << ", \"implied\": " << num(implied)
           << ", \"formula\": " << json_string(prefix) << ", \"tex\": " << json_string(Printer{sig, k.tree}.text()) << "}";
        first_row = false;
      }
    }
  ff << "\n]\n";

  for (const auto& v : variants)
    for (bool typed : {true, false}) {
      int found = 0, best_law = 0;
      double evals = 0;
      for (const auto& r : results)
        if (r.variant == &v && r.typed == typed) {
          if (r.law_at >= 0) found++, evals += r.law_at;
          best_law += r.best_is_law;
        }
      std::fprintf(stderr, "%-11s %-6s %-7s the law met on %2d of %d seeds (mean evaluations to it %.0f); the final best is the law on %d\n", v.problem.c_str(), v.units.name.c_str(),
                   typed ? "typed" : "untyped", found, seeds, found ? evals / found : 0.0, best_law);
    }
  return 0;
}
