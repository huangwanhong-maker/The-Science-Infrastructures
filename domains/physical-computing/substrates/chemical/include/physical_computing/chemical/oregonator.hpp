// The Belousov–Zhabotinsky medium by the Oregonator (substrates/chemical/README.md, simulators 1 and 2).
//
// The kinetics. The two-variable Oregonator of Tyson and Fife (1980), reduced from Field, Körös and Noyes's
// mechanism (1972) by Field and Noyes (1974), with u the activator (HBrO₂) and v the oxidised catalyst, scaled:
//
//   ε u' = u(1 − u) − (f v + φ)(u − q)/(u + q),   v' = u − v,
//
// with ε ≪ 1 the ratio of time scales, q ≪ 1 a ratio of rate constants, f the stoichiometric factor, and φ the
// bromide produced by light in the photosensitive, ruthenium-catalysed medium (Krug, Pohlmann and Kuhnert 1990;
// Kádár, Amemiya and Showalter 1997), the control by which the medium is programmed. The steady state has
// v* = u*, with u* the positive root of u(1 − u) = (f u + φ)(u − q)/(u + q), found by Newton's method; for φ = 0
// it is closed:
//
//   u* = ½ [(1 − q − f) + √((1 − q − f)² + 4q(1 + f))].
//
// The steady state loses stability at Hopf points where the Jacobian's trace vanishes (its determinant being
// positive); between them the stirred medium oscillates, beyond them it is excitable: at rest, but a large enough
// disturbance runs a pulse. In space (physis::chemistry::FieldKinetics), with D_u = 1 and an immobile catalyst
// (D_v = 0) as in a gel, an excitable medium carries trigger waves, and a broken wave curls into a spiral.
//
// Limits (declared): two variables; the scaled parameters are not fitted to a recipe here (Tyson 1985 gives
// their relation to concentrations); light enters as a bromide source only. Checks: tests/check_oregonator.cpp.
#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <calx/solvers/roots.hpp>
#include <physis/chemistry/field_kinetics.hpp>
#include <physis/dynamics/continuation.hpp>
#include <physis/dynamics/flows.hpp>

#include <physical_computing/runtime/backend.hpp>

namespace physical_computing::chemical {

struct Oregonator {
  double epsilon = 0.04, q = 0.002, f = 1.0, light = 0;  // ε, q, f, φ

  void kinetics(const double* c, double* dc) const {
    const double u = c[0], v = c[1];
    dc[0] = (u * (1 - u) - (f * v + light) * (u - q) / (u + q)) / epsilon;
    dc[1] = u - v;
  }
  void operator()(double, const physis::dynamics::State& x, physis::dynamics::State& dx) const { kinetics(x.data(), dx.data()); }

  // The steady state, by Newton from the closed form at φ = 0 (exact there).
  double steady_u() const {
    const double a = 1 - q - f;
    double u = 0.5 * (a + std::sqrt(a * a + 4 * q * (1 + f)));
    for (int it = 0; it < 60 && light != 0; it++) {
      const double g = u * (1 - u) - (f * u + light) * (u - q) / (u + q);
      const double dg = 1 - 2 * u - f * (u - q) / (u + q) - (f * u + light) * 2 * q / ((u + q) * (u + q));
      const double du = g / dg;
      u -= du;
      if (std::fabs(du) < 1e-15) break;
    }
    return u;
  }
  // The trace and determinant of the Jacobian at the steady state.
  std::pair<double, double> trace_determinant() const {
    const double u = steady_u(), v = u;
    const double j11 = (1 - 2 * u - f * v * 2 * q / ((u + q) * (u + q)) - light * 2 * q / ((u + q) * (u + q))) / epsilon;
    const double j12 = -f * (u - q) / (u + q) / epsilon;
    return {j11 - 1, -j11 - j12};
  }
};

// The values of f at which the steady state's trace vanishes with a positive determinant: the Hopf points, by a
// scan and Brent's method over f ∈ [f_lo, f_hi].
inline std::vector<double> hopf_points_in_f(Oregonator o, double f_lo = 0.05, double f_hi = 4.0, int scan = 4000) {
  std::vector<double> out;
  auto trace = [&](double f) { o.f = f; return o.trace_determinant().first; };
  double a = f_lo, ta = trace(a);
  for (int k = 1; k <= scan; k++) {
    const double b = f_lo + (f_hi - f_lo) * k / scan, tb = trace(b);
    if ((ta < 0) != (tb < 0)) {
      const double root = calx::solvers::brent(trace, a, b, 1e-13).x;
      o.f = root;
      if (o.trace_determinant().second > 0) out.push_back(root);
    }
    a = b, ta = tb;
  }
  return out;
}

// The BZ medium on a periodic grid: physis's FieldKinetics with the Oregonator's kinetics.
class BZMedium {
 public:
  BZMedium(const Oregonator& o, std::size_t nx, std::size_t ny, double spacing, double diffusion_u = 1.0, double diffusion_v = 0.0)
      : BZMedium(o, nx, ny, physis::chemistry::SparseLaplacian::periodic_grid(nx, ny, spacing), diffusion_u, diffusion_v) {}
  // A dish: an nx × ny grid with no-flux walls, where a plane wave leaves and only a rotating spiral sustains itself.
  static BZMedium dish(const Oregonator& o, std::size_t nx, std::size_t ny, double spacing, double diffusion_u = 1.0, double diffusion_v = 0.0) {
    return BZMedium(o, nx, ny, physis::chemistry::SparseLaplacian::grid(nx, ny, spacing), diffusion_u, diffusion_v);
  }
  // A chain of n sites of the given spacing, with no-flux ends: one dimension, for waves and their speeds.
  static BZMedium chain(const Oregonator& o, std::size_t n, double spacing, double diffusion_u = 1.0, double diffusion_v = 0.0) {
    std::vector<std::tuple<std::size_t, std::size_t, double>> edges;
    for (std::size_t i = 0; i + 1 < n; i++) edges.push_back({i, i + 1, 1 / (spacing * spacing)});
    return BZMedium(o, n, 1, physis::chemistry::SparseLaplacian::graph(n, edges), diffusion_u, diffusion_v);
  }
  BZMedium(const Oregonator& o, std::size_t nx, std::size_t ny, physis::chemistry::SparseLaplacian laplacian, double diffusion_u, double diffusion_v)
      : o_(o), nx_(nx), ny_(ny), field_(2, [o](const double* c, double* dc) { o.kinetics(c, dc); }, std::move(laplacian), {diffusion_u, diffusion_v}) {
    const double u = o.steady_u();
    for (std::size_t v = 0; v < nx * ny; v++) field_.state()[2 * v] = u, field_.state()[2 * v + 1] = u;
    field_.set_non_negative();  // u and v are concentrations
  }
  // A broken wave: a band excited from the lower wall to the middle, refractory behind it, so that its free end, at the
  // middle, curls about the refractory region into a spiral (in a dish; on a torus the end can also close the front
  // around the periodic direction).
  void break_wave(double band = 0.05) {
    for (std::size_t j = 0; j < ny_; j++)
      for (std::size_t i = 0; i < nx_; i++) {
        const std::size_t v = j * nx_ + i;
        const double x = static_cast<double>(i) / nx_, y = static_cast<double>(j) / ny_;
        if (y < 0.5 && std::fabs(x - 0.5) < band) field_.state()[2 * v] = 0.8;
        if (y < 0.5 && x < 0.5 && x > 0.5 - 3 * band) field_.state()[2 * v + 1] = 0.6;
      }
  }
  void excite(std::size_t i, std::size_t j, double u = 0.8) { field_.state()[2 * (j * nx_ + i)] = u; }
  // The kinetics are stiff in a wave's refractory tail (∂u'/∂u ≈ 2fv/(qε)), so the medium is stepped by physis's
  // split step, implicit in the kinetics; the step is bounded by the diffusion and, for accuracy, by ε/4.
  void run(double duration, double step) { field_.run_split(duration, step); }
  double stable_step() const { return std::min(0.5 * field_.diffusive_step(), 0.25 * o_.epsilon); }
  std::vector<double> activator() const {
    std::vector<double> a(nx_ * ny_);
    for (std::size_t v = 0; v < a.size(); v++) a[v] = field_.state()[2 * v];
    return a;
  }
  std::vector<double> catalyst() const {
    std::vector<double> a(nx_ * ny_);
    for (std::size_t v = 0; v < a.size(); v++) a[v] = field_.state()[2 * v + 1];
    return a;
  }
  double time() const { return field_.time(); }
  std::size_t nx() const { return nx_; }
  std::size_t ny() const { return ny_; }

  runtime::TargetDescription describe() const {
    runtime::TargetDescription d;
    d.family = "chemical", d.substrate = "BZ medium";
    d.accepts = {tasks::Kind::phenomenon, tasks::Kind::filter};
    d.inputs = {{"light φ", "1"}};
    d.outputs = {{"activator u", "1"}, {"oxidised catalyst v", "1"}};
    d.materials.materials = {{"malonic acid", "the organic substrate"}, {"sodium bromate", "the oxidant"}, {"sulfuric acid", "the acid"}, {"ferroin or a ruthenium complex", "the catalyst"}};
    d.materials.hazards = {"bromate: a strong oxidant, toxic", "sulfuric acid: corrosive"};
    d.materials.reuse = computer::Reuse::consumed_per_run;
    d.materials.end_of_life = computer::EndOfLife::neither;
    d.model = "the two-variable Oregonator (Tyson and Fife 1980) on a grid (physis::chemistry::FieldKinetics)";
    d.limits = {"two variables", "scaled parameters, not a recipe", "light as a bromide source only"};
    return d;
  }

 private:
  Oregonator o_;
  std::size_t nx_, ny_;
  physis::chemistry::FieldKinetics field_;
};

}  // namespace physical_computing::chemical
