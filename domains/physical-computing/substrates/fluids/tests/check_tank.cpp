// Checks of the wave tank (tank.hpp).
//
// - Its wave speed is √(g h): 0.2215 m/s at 5 mm.
// - Resonance: a motor driven at the tank's mode (1, 0) excites it more than five times as much as at 0.6 of
//   its frequency (lightly damped, γ = 0.1/s; the mode's amplitude by projection).
// - The echo state property: two runs from different surfaces of equal volume under the same inputs approach each
//   other as e^{−γt} (the gap after 30 s at most 3 e^{−γt} of the first); the mean level, the volume, is kept.
// - The runtime drives it as any backend: the camera's pixels, thresholded images of 0 and 1, the provenance
//   of a simulator; its materials name copper (critical) for the motors.

#include <physical_computing/fluids/tank.hpp>
#include <physical_computing/runtime/run.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numbers>
#include <vector>

using namespace physical_computing;

int main() {
  calx::verification::Suite s;
  fluids::WaveTankOptions o;
  s.check("tank: the wave speed √(g h) at 5 mm (m/s)", fluids::WaveTank(o).wave_speed() - std::sqrt(9.80665 * 0.005), 1e-15, "m/s");

  {
    fluids::WaveTankOptions q = o;
    q.damping = 0.1;
    q.motors = {{0.08, 0.15, 0.012}};
    auto amplitude = [&](double omega) {
      fluids::WaveTank t(q);
      t.reset(1);
      const auto mode = t.basin().mode(1, 0);
      double nn = 0;
      for (double v : mode) nn += v * v;
      const double h = t.basin().step_size();
      double peak = 0;
      const std::size_t steps = static_cast<std::size_t>(80.0 / h);
      for (std::size_t k = 0; k < steps; k++) {
        const double u = std::sin(omega * k * h);
        t.apply(std::vector<double>{u}, h);
        if (k * h > 50) {
          double a = 0;
          for (std::size_t c = 0; c < mode.size(); c++) a += t.basin().surface()[c] * mode[c];
          peak = std::max(peak, std::fabs(a / nn));
        }
      }
      return peak;
    };
    fluids::WaveTank t(q);
    const double w10 = t.basin().grid_frequency(1, 0);
    const double ratio = amplitude(w10) / amplitude(0.6 * w10);
    s.require("tank: a motor at the mode's frequency excites it over five times more than at 0.6 of it", ratio > 5, "ratio " + std::to_string(ratio));
  }

  {
    fluids::WaveTank a(o), b(o);
    a.reset(1), b.reset(2);
    calx::random::Stream rng(3, 0, 0, 0);
    std::vector<double> start(o.nx * o.ny);
    double mean = 0;
    for (double& v : start) v = 1e-3 * rng.normal(), mean += v;
    mean /= static_cast<double>(start.size());
    for (double& v : start) v -= mean;  // of equal volume: the mean level is never forgotten
    b.basin().set(start);
    auto gap = [&]() {
      double d = 0;
      for (std::size_t c = 0; c < start.size(); c++) d += std::pow(a.basin().surface()[c] - b.basin().surface()[c], 2);
      return std::sqrt(d);
    };
    const double g0 = gap();
    const double h = a.basin().step_size();
    double worst = 0;
    for (std::size_t k = 0; k < static_cast<std::size_t>(30.0 / h); k++) {
      const std::vector<double> u{std::sin(1.3 * k * h), std::cos(0.7 * k * h)};
      a.apply(u, h), b.apply(u, h);
      if (k * h > 20) worst = std::max(worst, gap() / (g0 * std::exp(-o.damping * a.time())));
    }
    s.check("tank: the echo state property, the gap at most 3 e^{−γt} of the first (ratio)", worst, 3);
  }

  {
    fluids::WaveTank t(o);
    schedule::Schedule sch;
    sch.hold = 0.05, sch.washout = 10;
    std::vector<std::vector<double>> u;
    for (int k = 0; k < 60; k++) u.push_back({std::sin(0.3 * k), 0.5});
    const auto d = runtime::drive(t, u, sch);
    s.require("tank: the runtime drives it, 16 × 12 pixels a reading", d.runs[0].readings.size() == 50 && d.runs[0].readings[0].size() == 192);
    bool binary = true;
    for (double v : t.binary_image()) binary = binary && (v == 0 || v == 1);
    s.require("tank: a thresholded image is of 0 and 1", binary);
    s.require("tank: its readings are a simulator's", !d.runs[0].provenance.observed() && d.runs[0].provenance.family == "fluids");
    const auto crit = t.describe().materials.critical();
    s.require("tank: its critical materials are the motors' copper", crit.size() == 1 && crit[0] == "copper");
  }
  return s.done();
}
