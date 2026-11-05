// The first physical substrates: one example of each simulator of the fluids, chemical and molecular families
// (ADR 0076), written as JSON to <results>/ for figures.py. Every number comes from the domain's simulators and
// physis; nothing is fitted here.
//
//   substrates <results directory> [section ...]     sections: convection rotating faraday walkers tank capillary
//                                                    foam bubbles phillips oregonator spiral array programs
//                                                    electrodes strands seesaw pen (default: all)

#include <physical_computing/chemical/bz_array.hpp>
#include <physical_computing/chemical/electrochemical.hpp>
#include <physical_computing/chemical/oregonator.hpp>
#include <physical_computing/chemical/programs.hpp>
#include <physical_computing/fluids/bubble_logic.hpp>
#include <physical_computing/fluids/capillary.hpp>
#include <physical_computing/fluids/convection_cell.hpp>
#include <physical_computing/fluids/faraday.hpp>
#include <physical_computing/fluids/foam.hpp>
#include <physical_computing/fluids/phillips.hpp>
#include <physical_computing/fluids/rotating_tank.hpp>
#include <physical_computing/fluids/tank.hpp>
#include <physical_computing/fluids/walkers.hpp>
#include <physical_computing/molecular/pen.hpp>
#include <physical_computing/molecular/seesaw.hpp>
#include <physical_computing/molecular/strand_displacement.hpp>

#include <physis/chemistry/field_kinetics.hpp>
#include <physis/dynamics/continuation.hpp>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace physical_computing;
namespace pd = physis::dynamics;

// ---- a small JSON writer ----------------------------------------------------------------------------------------------
struct Json {
  std::ostringstream o;
  bool first = true;
  Json() { o.precision(9), o << "{"; }
  Json& key(const std::string& k) { o << (first ? "" : ",") << "\n\"" << k << "\":"; first = false; return *this; }
  Json& num(const std::string& k, double v) { key(k); write(v); return *this; }
  Json& str(const std::string& k, const std::string& v) { key(k) ; o << "\"" << v << "\""; return *this; }
  Json& arr(const std::string& k, const std::vector<double>& v) { key(k); list(v); return *this; }
  Json& arr2(const std::string& k, const std::vector<std::vector<double>>& v) {
    key(k), o << "[";
    for (std::size_t i = 0; i < v.size(); i++) o << (i ? "," : ""), list(v[i]);
    o << "]";
    return *this;
  }
  Json& raw(const std::string& k, const std::string& v) { key(k), o << v; return *this; }
  void write(double v) { if (std::isfinite(v)) o << v; else o << "null"; }
  void list(const std::vector<double>& v) {
    o << "[";
    for (std::size_t i = 0; i < v.size(); i++) o << (i ? "," : ""), write(v[i]);
    o << "]";
  }
  std::string text() { return o.str() + "\n}\n"; }
};

static std::filesystem::path out_dir;
static void save(const std::string& name, Json& j) {
  std::ofstream f(out_dir / (name + ".json"));
  f << j.text();
  std::printf("wrote %s.json\n", name.c_str());
}

// ---- fluids -----------------------------------------------------------------------------------------------------------
static void convection() {
  using physis::continuum::Plate;
  Json j;
  std::vector<double> a;
  for (int k = 0; k <= 55; k++) a.push_back(1.0 + 0.1 * k);
  j.arr("a", a);
  for (auto [name, bot, top] : {std::tuple{"rigid_rigid", Plate::rigid, Plate::rigid}, std::tuple{"rigid_free", Plate::rigid, Plate::free}, std::tuple{"free_free", Plate::free, Plate::free}}) {
    std::vector<double> ra;
    for (double x : a) ra.push_back(physis::continuum::marginal_rayleigh(x, bot, top));
    const auto c = physis::continuum::critical_rayleigh(bot, top);
    j.arr(std::string(name), ra).arr(std::string(name) + "_critical", {c.wavenumber, c.rayleigh});
  }
  // The cell's states against its heating, by Lorenz's truncation: the Nusselt number's range after transients.
  fluids::ConvectionCell cell({});
  const double dTc = cell.critical_temperature_difference();
  std::vector<double> r, lo, hi;
  for (int k = 0; k <= 160; k++) {
    const double rr = 0.4 + 0.22 * k;
    cell.set_heating(rr * dTc);
    const auto p = cell.lorenz();
    const pd::Lorenz f{p.sigma, p.r, p.b};
    pd::State x{1, 1, rr};
    pd::integrate(f, 0.0, 60.0, x, 0.002);
    double nmin = 1e300, nmax = -1e300;
    pd::integrate(f, 0.0, 40.0, x, 0.002, [&](double, const pd::State& y) {
      const double nu = physis::continuum::lorenz_nusselt(y[2], p.r);
      nmin = std::min(nmin, nu), nmax = std::max(nmax, nu);
    });
    r.push_back(rr), lo.push_back(nmin), hi.push_back(nmax);
  }
  j.arr("r", r).arr("nusselt_min", lo).arr("nusselt_max", hi).num("hopf_r", cell.hopf_r()).num("dTc", dTc).num("prandtl", cell.prandtl());
  // An orbit at r = 28 relative to the cell's onset.
  cell.set_heating(28 * dTc);
  const auto p = cell.lorenz();
  const pd::Lorenz f{p.sigma, p.r, p.b};
  pd::State x{1, 1, 28};
  pd::integrate(f, 0.0, 20.0, x, 0.002);
  std::vector<double> X, Z;
  std::size_t n = 0;
  pd::integrate(f, 0.0, 40.0, x, 0.002, [&](double, const pd::State& y) { if (n++ % 4 == 0) X.push_back(y[0]), Z.push_back(y[2]); });
  j.arr("orbit_x", X).arr("orbit_z", Z);
  save("convection", j);
}

static void rotating() {
  Json j;
  const auto w = fluids::Fluid::water();
  fluids::RotatingTankOptions base;
  auto at = [&](double theta, bool drag) {
    fluids::RotatingTankOptions q = base;
    q.ekman_drag = drag;
    q.temperature_difference = theta * q.rotation * q.rotation * q.gap * q.gap / (9.80665 * w.thermal_expansion * q.depth);
    return fluids::RotatingTank(q);
  };
  std::vector<double> ms;
  for (int m = 1; m <= 12; m++) ms.push_back(m);
  j.arr("m", ms);
  std::vector<double> thetas{0.3, 0.65, 1.0, 1.4};
  j.arr("thetas", thetas);
  std::vector<std::vector<double>> g_free, g_drag;
  for (double th : thetas) {
    std::vector<double> a, b;
    const auto t0 = at(th, false), t1 = at(th, true);
    for (double m : ms) a.push_back(t0.growth(m)), b.push_back(t1.growth(m));
    g_free.push_back(a), g_drag.push_back(b);
  }
  j.arr2("growth_inviscid", g_free).arr2("growth_ekman", g_drag);
  const fluids::RotatingTank t(base);
  std::vector<double> crit;
  for (double m : ms) crit.push_back(t.inviscid_critical_thermal_rossby(static_cast<int>(m)));
  j.arr("critical_theta_by_m", crit).num("theta", t.thermal_rossby()).num("taylor", t.taylor());
  // A nonlinear run in the wave regime, from small noise.
  fluids::RotatingTankOptions q = base;
  q.nx = 128, q.ny = 16;
  const fluids::RotatingTank tank(q);
  const auto [m, rate] = tank.most_unstable();
  auto model = tank.model();
  model.perturb(1e-9, 21);
  std::vector<double> times, energy;
  const double efold = 1 / rate;
  while (model.time() < 45 * efold) {
    model.run(20);
    times.push_back(model.time()), energy.push_back(model.energy());
  }
  j.num("most_unstable", m).num("rate", rate).num("dominant", fluids::RotatingTank::dominant_wavenumber(model, 1));
  j.arr("time", times).arr("energy", energy).arr("psi1", model.streamfunction(1)).arr("psi2", model.streamfunction(2));
  j.num("nx", q.nx).num("ny", q.ny).num("gap", q.gap).num("radius", q.mean_radius);
  save("rotating", j);
}

static void faraday() {
  Json j;
  fluids::FaradayLayer layer({});
  const double kr = layer.resonant_wavenumber();
  std::vector<double> k, A, small;
  for (int i = 0; i <= 40; i++) {
    const double kk = kr * (0.7 + 0.015 * i);
    k.push_back(kk), A.push_back(layer.threshold(kk)), small.push_back(layer.small_damping_threshold(kk));
  }
  j.arr("k", k).arr("threshold", A).arr("small_damping", small).num("resonant_k", kr);
  std::vector<double> hz, onset_a, onset_l;
  for (double f : {15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0}) {
    fluids::FaradayOptions o;
    o.drive_frequency = 2 * std::numbers::pi * f;
    const auto on = fluids::FaradayLayer(o).onset();
    hz.push_back(f), onset_a.push_back(on.acceleration), onset_l.push_back(on.wavelength);
  }
  j.arr("drive_hz", hz).arr("onset_acceleration", onset_a).arr("onset_wavelength", onset_l);
  save("faraday", j);
}

static void walkers() {
  Json j;
  std::vector<double> beta, speed;
  for (double b : {1.0, 1.5, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 5.0, 6.0}) {
    fluids::WalkerOptions o;
    o.beta = b;
    fluids::Walker w(o);
    w.run(80);
    beta.push_back(b), speed.push_back(w.speed());
  }
  std::vector<double> bt, ut;
  for (int i = 0; i <= 120; i++) bt.push_back(1.0 + 0.05 * i), ut.push_back(fluids::Walker::steady_speed(1.0 + 0.05 * i));
  j.arr("beta", beta).arr("speed", speed).arr("beta_theory", bt).arr("speed_theory", ut);
  // In a harmonic well: F = −k x.
  fluids::WalkerOptions o;
  o.beta = 4.5;
  const double kw = 0.02;
  o.external = [kw](const fluids::Vec2& x, const fluids::Vec2&) { return fluids::Vec2{-kw * x[0], -kw * x[1]}; };
  fluids::Walker w(o, {4, 0}, {0, 1});
  std::vector<double> X, Y;
  for (int s = 0; s < 20000; s++) {
    w.step();
    if (s % 10 == 0) X.push_back(w.position()[0]), Y.push_back(w.position()[1]);
  }
  j.arr("orbit_x", X).arr("orbit_y", Y).num("well", kw).num("orbit_beta", 4.5);
  save("walkers", j);
}

static void tank() {
  Json j;
  fluids::WaveTankOptions o;
  o.motors = {{0.08, 0.15, 0.012}, {0.32, 0.15, 0.012}};
  fluids::WaveTank t(o);
  t.reset(1);
  const double h = t.basin().step_size();
  for (std::size_t k = 0; k < static_cast<std::size_t>(14.0 / h); k++) {
    const std::vector<double> u{std::sin(2 * std::numbers::pi * 2.0 * k * h), std::sin(2 * std::numbers::pi * 2.6 * k * h)};
    t.apply(u, h);
  }
  j.arr("surface", t.basin().surface()).arr("camera", t.read()).arr("binary", t.binary_image());
  j.num("nx", o.nx).num("ny", o.ny).num("cx", o.camera_x).num("cy", o.camera_y).num("length", o.length).num("width", o.width);
  // The echo state property: two tanks, one from rest and one from a random surface of the same volume.
  fluids::WaveTank a(o), b(o);
  a.reset(1), b.reset(2);
  calx::random::Stream rng(3, 0, 0, 0);
  std::vector<double> start(o.nx * o.ny);
  double mean = 0;
  for (double& v : start) v = 1e-3 * rng.normal(), mean += v;
  for (double& v : start) v -= mean / static_cast<double>(start.size());
  b.basin().set(start);
  std::vector<double> time, gap;
  for (std::size_t k = 0; k < static_cast<std::size_t>(25.0 / h); k++) {
    const std::vector<double> u{std::sin(1.3 * k * h), std::cos(0.7 * k * h)};
    a.apply(u, h), b.apply(u, h);
    if (k % 10 == 0) {
      double d = 0;
      for (std::size_t c = 0; c < start.size(); c++) d += std::pow(a.basin().surface()[c] - b.basin().surface()[c], 2);
      time.push_back(a.time()), gap.push_back(std::sqrt(d));
    }
  }
  j.arr("esp_time", time).arr("esp_gap", gap).num("damping", o.damping);
  save("tank", j);
}

static void capillary() {
  Json j;
  const std::size_t W = 20;
  const auto maze = fluids::CapillaryNetwork::maze(W, W, 11);
  const auto sol = maze.flow(0, W * W - 1, 100);
  std::vector<double> a, b, q;
  for (std::size_t k = 0; k < maze.tubes().size(); k++) a.push_back(maze.tubes()[k].a), b.push_back(maze.tubes()[k].b), q.push_back(sol.flow[k]);
  j.num("width", W).arr("a", a).arr("b", b).arr("flow", q);
  save("capillary", j);
}

static void foam() {
  Json j;
  fluids::FoamOptions o;
  o.size = 128, o.cells = 500, o.seed = 6;
  fluids::Foam f(o);
  auto states = [&]() {
    std::vector<double> s;
    for (auto v : f.lattice().states_of_sites()) s.push_back(static_cast<double>(v));
    return s;
  };
  j.arr("s0", states()).num("cells0", f.cells().size());
  f.coarsen(40);
  const auto before = f.cells();
  j.arr("s1", states()).num("cells1", before.size()).num("sweeps1", f.sweeps());
  f.coarsen(200);
  const auto after = f.cells();
  j.arr("s2", states()).num("cells2", after.size()).num("sweeps2", f.sweeps());
  std::vector<double> n, rate, count;
  for (const auto& [sides, r] : fluids::Foam::area_rates(before, after, 200)) n.push_back(sides), rate.push_back(r.first), count.push_back(r.second);
  j.arr("sides", n).arr("rate", rate).arr("count", count).num("size", o.size);
  save("foam", j);
}

static void bubbles() {
  Json j;
  std::vector<double> table;
  for (auto [x, y] : {std::pair{false, false}, {true, false}, {false, true}, {true, true}}) {
    const auto g = fluids::BubbleCircuit::evaluate(x, y);
    table.push_back(x), table.push_back(y), table.push_back(g.or_out), table.push_back(g.and_out);
  }
  j.arr("truth", table);
  save("bubbles", j);
}

static void phillips() {
  Json j;
  fluids::PhillipsOptions o;
  fluids::PhillipsMachine m(o);
  m.reset(0);
  const double tau = m.time_constant();
  std::vector<double> t, Y, I;
  const double h = tau / 50;
  for (int k = 0; k < 800; k++) {
    const double inj = k < 300 ? 1e-4 : (k < 550 ? 2e-4 : 0.5e-4);
    m.apply(std::vector<double>{inj}, h);
    t.push_back(m.time()), Y.push_back(m.income()), I.push_back(inj);
  }
  j.arr("time", t).arr("income", Y).arr("injections", I).num("propensity", o.propensity).num("tau", tau);
  save("phillips", j);
}

// ---- chemical ---------------------------------------------------------------------------------------------------------
static void oregonator() {
  Json j;
  chemical::Oregonator o;
  const auto hopf = chemical::hopf_points_in_f(o);
  j.arr("hopf", hopf).num("epsilon", o.epsilon).num("q", o.q);
  std::vector<double> f, us, stable, umin, umax;
  for (int k = 0; k <= 90; k++) {
    chemical::Oregonator p = o;
    p.f = 0.2 + 0.035 * k;
    const double u0 = p.steady_u();
    const auto [tr, det] = p.trace_determinant();
    pd::State x{u0 * 1.5 + 0.01, u0};
    pd::integrate(p, 0.0, 150.0, x, 2e-4);
    double lo = 1e300, hi = -1e300;
    pd::integrate(p, 0.0, 30.0, x, 2e-4, [&](double, const pd::State& y) { lo = std::min(lo, y[0]), hi = std::max(hi, y[0]); });
    f.push_back(p.f), us.push_back(u0), stable.push_back(tr < 0 && det > 0 ? 1 : 0), umin.push_back(lo), umax.push_back(hi);
  }
  j.arr("f", f).arr("steady", us).arr("stable", stable).arr("u_min", umin).arr("u_max", umax);
  // One oscillation in time, between the Hopf points.
  chemical::Oregonator p = o;
  p.f = 1.0;
  pd::State x{0.1, 0.1};
  pd::integrate(p, 0.0, 50.0, x, 2e-4);
  std::vector<double> T, U, V;
  std::size_t n = 0;
  pd::integrate(p, 0.0, 15.0, x, 2e-4, [&](double t, const pd::State& y) { if (n++ % 25 == 0) T.push_back(t), U.push_back(y[0]), V.push_back(y[1]); });
  j.arr("t", T).arr("u", U).arr("v", V);
  save("oregonator", j);
}

static void spiral() {
  Json j;
  // A spiral needs a medium wider than its wavelength, about speed × recovery time, which spans some T/ε front
  // widths: at ε = 0.04 that is a few hundred cells, more than a picture is worth; at ε = 0.12 a dish of 160 × 160
  // holds one, at a spacing that resolves the fronts. The dish has no-flux walls, where a plane wave leaves, so that
  // activity still there at t = 40 is the rotating spiral (on a torus a plane wave can circulate for ever).
  chemical::Oregonator o;
  o.epsilon = 0.12;
  const auto hopf = chemical::hopf_points_in_f(o);
  o.f = hopf[1] + 0.05;
  const std::size_t n = 160;
  const double spacing = 0.15;
  auto m = chemical::BZMedium::dish(o, n, n, spacing);
  m.break_wave(0.03);
  const double h = m.stable_step();
  std::vector<std::vector<double>> frames;
  std::vector<double> times, excited;
  for (double t_end : {4.0, 14.0, 27.0, 40.0}) {
    m.run(t_end - m.time(), h);
    frames.push_back(m.activator()), times.push_back(m.time());
    double on = 0;
    for (double u : frames.back()) on += u > 0.1;
    excited.push_back(on);
  }
  j.arr2("frames", frames).arr("times", times).arr("excited", excited).num("n", n).num("spacing", spacing).num("f", o.f).num("epsilon", o.epsilon).arr("catalyst", m.catalyst());
  save("spiral", j);
}

static void array() {
  Json j;
  chemical::Oregonator base;
  const auto hopf = chemical::hopf_points_in_f(base);
  chemical::Oregonator a = base, b = base;
  a.f = 0.5 * (hopf[0] + hopf[1]), b.f = a.f * 1.05;
  for (auto [name, coupling] : {std::pair{"uncoupled", 0.0}, std::pair{"coupled", 2.0}}) {
    chemical::BZArray arr({a, b}, {{0, 1, coupling}});
    arr.run(40, 2e-3);
    std::vector<double> T, va, vb;
    for (int k = 0; k < 1500; k++) {
      arr.run(0.02, 2e-3);
      T.push_back(arr.time()), va.push_back(arr.state()[1]), vb.push_back(arr.state()[3]);
    }
    j.arr(std::string(name) + "_t", T).arr(std::string(name) + "_a", va).arr(std::string(name) + "_b", vb);
  }
  save("array", j);
}

static void programs() {
  using morphe::arithmetic::Integer;
  using morphe::arithmetic::Rational;
  Json j;
  const auto P = physis::chemistry::lorenz_polynomials(10, 28, Rational(Integer(8), Integer(3)));
  const auto D = chemical::lower_flow(P);
  pd::State c{1, 0, 1, 0, 20, 0}, L{1, 1, 20};
  auto crn = [&](double, const pd::State& y, pd::State& dy) { physis::chemistry::mass_action(D.network, y, dy); };
  std::vector<double> t, xp, xm, xl, err;
  for (int k = 0; k < 400; k++) {
    pd::integrate(crn, 0.0, 0.01, c, 1e-4);
    pd::integrate(pd::Lorenz{}, 0.0, 0.01, L, 1e-4);
    t.push_back((k + 1) * 0.01), xp.push_back(c[0]), xm.push_back(c[1]), xl.push_back(L[0]), err.push_back(std::fabs(c[0] - c[1] - L[0]));
  }
  j.arr("t", t).arr("x_plus", xp).arr("x_minus", xm).arr("x_lorenz", xl).arr("error", err).num("reactions", D.network.reactions.size());
  // The max program's law of absorption, exactly.
  std::vector<double> table;
  for (long long x1 = 0; x1 <= 4; x1++)
    for (long long x2 = 0; x2 <= 4; x2++) {
      const auto law = chemical::run_exactly(chemical::max_program(), {x1, x2});
      table.push_back(static_cast<double>(x1)), table.push_back(static_cast<double>(x2));
      table.push_back(static_cast<double>(law.begin()->first)), table.push_back(law.begin()->second.to_double());
    }
  j.arr("max_table", table);
  save("programs", j);
}

static void electrodes() {
  Json j;
  std::vector<double> K, r, theory;
  const double gamma = 0.05;
  for (int k = 0; k <= 14; k++) {
    const double ratio = 0.2 + 0.2 * k;
    chemical::ElectrodeArray a(600, 2 * std::numbers::pi * 0.4, gamma);
    const double Kv = ratio * a.critical_coupling();
    K.push_back(ratio), r.push_back(a.run(Kv, 400, 100)), theory.push_back(a.predicted_order(Kv));
  }
  j.arr("K_over_Kc", K).arr("r", r).arr("theory", theory);
  save("electrodes", j);
}

// ---- molecular --------------------------------------------------------------------------------------------------------
static void strands() {
  using morphe::arithmetic::Rational;
  Json j;
  physis::chemistry::Network formal;
  formal.species = {"X", "Y", "Z"};
  formal.add("X → Y", {{"X", 1}}, {{"Y", 1}}, 1);
  formal.add("X + Y → Z", {{"X", 1}, {"Y", 1}}, {{"Z", 1}}, 1);
  auto formal_at = [&](double T) {
    pd::State c{1, 0, 0};
    auto f = [&](double, const pd::State& y, pd::State& dy) { physis::chemistry::mass_action(formal, y, dy); };
    pd::integrate(f, 0.0, T, c, 1e-4);
    return c;
  };
  auto dsd_at = [&](long long cmax, long long qmax, double T) {
    const auto D = molecular::compile_strand_displacement(formal, Rational(cmax), Rational(qmax));
    std::vector<double> c(D.initial.begin(), D.initial.end());
    c[0] = 1;
    return physis::chemistry::integrate_stiff(c.size(), [&](const double* y, double* dy) {
      std::vector<double> yy(y, y + D.network.species.size()), dd(yy.size());
      physis::chemistry::mass_action(D.network, yy, dd);
      for (std::size_t i = 0; i < dd.size(); i++) dy[i] = dd[i];
    }, c, T, 1e-3);
  };
  std::vector<double> cm, err;
  for (long long C : {30, 100, 300, 1000, 3000}) {
    const auto d = dsd_at(C, 10 * C, 2.0);
    const auto f = formal_at(2.0);
    double e = 0;
    for (int i = 0; i < 3; i++) e = std::max(e, std::fabs(d[i] - f[i]));
    cm.push_back(static_cast<double>(C)), err.push_back(e);
  }
  j.arr("c_max", cm).arr("error", err);
  std::vector<double> t, fx, fy, fz, dx, dy, dz;
  for (int k = 1; k <= 40; k++) {
    const double T = 0.1 * k;
    const auto f = formal_at(T), d = dsd_at(30, 300, T);
    t.push_back(T), fx.push_back(f[0]), fy.push_back(f[1]), fz.push_back(f[2]), dx.push_back(d[0]), dy.push_back(d[1]), dz.push_back(d[2]);
  }
  j.arr("t", t).arr("formal_x", fx).arr("formal_y", fy).arr("formal_z", fz).arr("dsd_x", dx).arr("dsd_y", dy).arr("dsd_z", dz);
  save("strands", j);
}

static void seesaw() {
  Json j;
  for (auto [name, theta] : {std::pair{"or", 0.6}, std::pair{"and", 1.2}}) {
    const auto g = molecular::two_input_gate(theta);
    std::vector<std::vector<double>> curves;
    for (auto [a, b] : {std::pair{false, false}, {true, false}, {false, true}, {true, true}}) {
      std::vector<double> y;
      for (int k = 0; k <= 60; k++) y.push_back(molecular::run_gate(g, a, b, 2.0 * k + 1e-9, 1e-3));
      curves.push_back(y);
    }
    j.arr2(std::string(name), curves);
  }
  std::vector<double> t;
  for (int k = 0; k <= 60; k++) t.push_back(2.0 * k);
  j.arr("t", t);
  save("seesaw", j);
}

static void pen() {
  Json j;
  molecular::PenPredatorPrey p;
  const double KH = p.hopf_template();
  j.num("hopf_template", KH);
  for (auto [name, ratio] : {std::pair{"below", 0.8}, std::pair{"above", 1.4}}) {
    auto L = p.law;
    L.K = ratio * KH;
    pd::State x{L.prey_at_coexistence() * 1.6, L.predator_at_coexistence() * 0.6};
    std::vector<double> N, P;
    std::size_t n = 0;
    pd::integrate(L, 0.0, 400.0, x, 0.01, [&](double, const pd::State& y) { if (n++ % 20 == 0) N.push_back(y[0]), P.push_back(y[1]); });
    j.arr(std::string(name) + "_prey", N).arr(std::string(name) + "_predator", P).num(std::string(name) + "_K", L.K);
  }
  auto F = p.prey_front(800, 0.1, 1.0);
  for (std::size_t i = 0; i < 20; i++) F.state()[i] = p.law.K;
  std::vector<std::vector<double>> profiles;
  const double dt = 0.2 * F.diffusive_step();
  for (double T : {5.0, 10.0, 15.0, 20.0, 25.0}) {
    F.run(T - F.time(), dt);
    profiles.push_back(F.state());
  }
  j.arr2("profiles", profiles).num("spacing", 0.1).num("capacity", p.law.K);
  save("pen", j);
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: substrates <results directory> [section ...]\n");
    return 2;
  }
  out_dir = argv[1];
  std::filesystem::create_directories(out_dir);
  std::set<std::string> want;
  for (int i = 2; i < argc; i++) want.insert(argv[i]);
  const std::vector<std::pair<std::string, void (*)()>> sections{
      {"convection", convection}, {"rotating", rotating}, {"faraday", faraday}, {"walkers", walkers}, {"tank", tank}, {"capillary", capillary},
      {"foam", foam}, {"bubbles", bubbles}, {"phillips", phillips}, {"oregonator", oregonator}, {"spiral", spiral}, {"array", array},
      {"programs", programs}, {"electrodes", electrodes}, {"strands", strands}, {"seesaw", seesaw}, {"pen", pen}};
  for (const auto& [name, f] : sections) {
    if (!want.empty() && !want.count(name)) continue;
    const auto t0 = std::chrono::steady_clock::now();
    f();
    std::printf("  %s: %.1f s\n", name.c_str(), std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    std::fflush(stdout);
  }
}
