// Checks of the modalities' physics and forward models (modalities.hpp; ADR 0070).
//
// - Hounsfield units: water 0, air −1000, and the inverse.
// - The spoiled gradient echo: M_z(0) = M₀ after no pulse, the steady state after many, and the Ernst angle
//   arccos E₁ maximising the steady signal; TOF's inflow: faster blood brighter, saturated blood darker than fresh.
// - The CT forward model: a uniform region keeps its value (within the noise's standard error); an edge's rise
//   from 10% to 90% is 2.563 σ of the PSF (a Gaussian's), within 10%.

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>
#include <medical_imaging/modalities/modalities.hpp>

using calx::verification::Suite;
namespace mm = medical_imaging::modalities;
namespace eg = eikon::grid;

int main() {
  Suite s;
  calx::random::Stream rng(20261104, 21, 0, 0);
  {
    const double mw = 0.19, ma = 0.0002;  // per cm, illustrative
    s.check("Hounsfield units: water 0", mm::hounsfield(mw, mw, ma), 1e-12);
    s.check("air −1000", mm::hounsfield(ma, mw, ma) + 1000, 1e-9);
    s.check("and back to attenuation", mm::attenuation(mm::hounsfield(0.21, mw, ma), mw, ma) - 0.21, 1e-15);
    const auto w = mm::window_named("soft tissue");
    s.require("the soft-tissue window maps 40 HU to the middle and clips beyond", w(40) == 0.5 && w(-1000) == 0 && w(1000) == 1);
  }
  {
    mm::SpoiledGradientEcho seq;
    const double t1 = 1.2;
    s.check("spoiled GRE: no pulse leaves M₀", mm::longitudinal(0, t1, seq) - 1, 1e-15);
    const double E1 = std::exp(-seq.tr / t1);
    s.check("many pulses reach M₀(1 − E₁)/(1 − E₁ cos α)", mm::longitudinal(1e6, t1, seq) - (1 - E1) / (1 - E1 * std::cos(seq.alpha)), 1e-15);
    double best = 0, best_alpha = 0;
    for (double a = 0.001; a < 1.0; a += 0.0005) {
      mm::SpoiledGradientEcho q{seq.tr, a};
      if (mm::steady_signal(t1, q) > best) best = mm::steady_signal(t1, q), best_alpha = a;
    }
    s.check("the Ernst angle arccos E₁ maximises the steady signal", best_alpha - std::acos(E1), 0.001, "rad");
    const double slow = mm::inflow_signal(30, 50, 1.6, seq), fast = mm::inflow_signal(30, 400, 1.6, seq), sat = mm::inflow_signal(30, 400, 1.6, seq, 0);
    s.require("TOF: faster inflow brighter, saturated inflow darker than fresh", fast > slow && sat < 0.5 * fast && slow > mm::steady_signal(1.6, seq),
              std::to_string(slow) + " < " + std::to_string(fast) + ", saturated " + std::to_string(sat));
  }
  {
    eg::Geometry g;
    g.size = {80, 40, 30};
    g.spacing = {0.25, 0.25, 0.5};
    g.origin = {-10, -5, -7.5};
    eg::Image<float> occ(g);
    for (std::size_t v = 0; v < g.count(); v++) occ.data[v] = g.world(double(g.unlinear(v)[0]), 0, 0)[0] >= 0 ? 1.0f : 0.0f;
    mm::CtProtocol p;
    p.grid = g;
    p.fwhm_inplane = 1.5;
    p.noise = 0;
    p.background = 0;
    const auto img = mm::ct({{occ, 100}}, p, rng);
    // the edge's 10–90% rise along x on the middle row
    double x10 = 0, x90 = 0;
    for (std::size_t i = 1; i < 80; i++) {
      const double a = img.at(i - 1, 20, 15), b = img.at(i, 20, 15);
      const double xa = g.world(double(i - 1), 0, 0)[0];
      if (a < 10 && b >= 10) x10 = xa + (10 - a) / (b - a) * 0.25;
      if (a < 90 && b >= 90) x90 = xa + (90 - a) / (b - a) * 0.25;
    }
    const double sigma = 1.5 / (2 * std::sqrt(2 * std::log(2.0)));
    s.check("CT: an edge rises from 10% to 90% over 2.563 σ of the PSF, relatively", (x90 - x10) / (2.563 * sigma) - 1, 0.1);
    p.noise = 20;
    const auto noisy = mm::ct({{occ, 100}}, p, rng);
    double sum = 0;
    std::size_t n = 0;
    for (std::size_t k = 0; k < 30; k++)
      for (std::size_t j = 0; j < 40; j++)
        for (std::size_t i = 60; i < 80; i++) sum += noisy.at(i, j, k), n++;
    s.check("a uniform region keeps its value, within 4 standard errors", sum / double(n) - 100, 4 * 20 / std::sqrt(double(n)), "HU");
  }
  return s.done();
}
