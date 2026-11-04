// The modalities of medical imaging as physics: what each measures, how a vessel appears in it, and forward
// models that make an image of a known body, so that every reconstruction can be checked against its truth
// (ADR 0070).
//
// CT measures the linear attenuation coefficient μ of X-rays, reported in Hounsfield units, HU = 1000 (μ − μ_w) /
// (μ_w − μ_a), so that water is 0 and air −1000 (Hounsfield 1973, Br. J. Radiol. 46, 1016). Blood is bright only
// where an iodinated contrast agent fills it; in an arterial phase the arteries, in a portal-venous phase the
// veins as well. A reconstructed CT is the true map blurred by the system's point-spread function and
// noisy; the forward model here takes a Gaussian PSF of stated widths in-plane and along z (slice profile), and
// white Gaussian noise. *Caveat.* Real CT noise is correlated (its power spectrum follows the reconstruction
// kernel); white noise makes the checks' images harder to segment at high frequency and easier at low.
//
// Time-of-flight MR angiography is a spoiled gradient echo whose stationary tissue is saturated by repeated
// pulses, while blood flowing into the slab brings fresh magnetisation. With E₁ = e^{−TR/T₁}, the longitudinal
// magnetisation after n pulses of flip angle α, from M_z(0), is M_n = M_ss + (M_z(0) − M_ss)(E₁ cos α)^n, with the
// steady state M_ss = M₀(1 − E₁)/(1 − E₁ cos α); the signal is ∝ M_n sin α (the Bloch equations with perfect
// spoiling; Haacke et al. 1999, Magnetic Resonance Imaging: Physical Principles and Sequence Design, Wiley; the
// Ernst angle cos α_E = E₁ maximises the steady state's signal). Blood that has travelled a depth d into the slab
// at speed v has received n = d/(v TR) pulses. A saturation band upstream sets M_z(0) = 0 for blood entering from
// that side, which is how TOF suppresses veins: arteries bright, veins dark, the contrast the integration of
// modalities uses. Magnitude images carry Rician noise (eikon's phantoms).
//
// Windows display a range of values: [level − width/2, level + width/2] mapped to [0, 1]. The named windows are
// the radiological convention for CT (soft tissue 40/400, lung −600/1500, bone 300/2000 HU), a display choice
// that changes no value. Checks: tests/check_modalities.cpp.

#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <calx/random/stream.hpp>
#include <eikon/filters/filters.hpp>
#include <eikon/grid/grid.hpp>
#include <eikon/phantoms/phantoms.hpp>

namespace medical_imaging::modalities {

using eikon::grid::Geometry;
using eikon::grid::Image;
using eikon::grid::Vec3;

// ---- CT ----------------------------------------------------------------------------------------------------

inline double hounsfield(double mu, double mu_water, double mu_air = 0) { return 1000 * (mu - mu_water) / (mu_water - mu_air); }
inline double attenuation(double hu, double mu_water, double mu_air = 0) { return mu_water + hu * (mu_water - mu_air) / 1000; }

struct Window {
  double level = 40, width = 400;
  double operator()(double v) const { return std::clamp((v - (level - width / 2)) / width, 0.0, 1.0); }
};
inline Window window_named(const std::string& name) {
  if (name == "lung") return {-600, 1500};
  if (name == "bone") return {300, 2000};
  return {40, 400};  // soft tissue
}

// A tissue of the body for the forward models: its occupancy on a fine grid and its value.
struct Tissue {
  Image<float> occupancy;  // fraction of each voxel, on the truth's grid
  double value = 0;        // HU for CT; for MR, unused (see tof)
};

struct CtProtocol {
  Geometry grid;             // the scanner's reconstruction grid
  double fwhm_inplane = 0.8; // mm
  double fwhm_z = 1.2;       // mm, the slice profile
  double noise = 15;         // HU, standard deviation
  double background = -1000; // HU outside every tissue (air)
};

// The image a CT gives of tissues (later tissues drawn over earlier ones within their occupancy).
inline Image<float> ct(const std::vector<Tissue>& tissues, const CtProtocol& p, calx::random::Stream& rng) {
  const auto& g = tissues.at(0).occupancy.geometry;
  Image<float> truth(g, float(p.background));
  for (const auto& t : tissues)
    for (std::size_t v = 0; v < truth.data.size(); v++) truth.data[v] += float((t.value - truth.data[v]) * t.occupancy.data[v]);
  // the PSF, separable: σ = FWHM/(2√(2 ln 2)) per axis of the truth's grid (assumed aligned with the scanner's)
  const double k = 2 * std::sqrt(2 * std::log(2.0));
  const double si = p.fwhm_inplane / k, sz = p.fwhm_z / k;
  Image<float> blurred = truth;
  const auto kern = [&](double sigma, double h) { return eikon::filters::gaussian_kernel(sigma / h, 0); };
  blurred = eikon::filters::convolve_axis(blurred, 0, kern(si, g.spacing[0]));
  blurred = eikon::filters::convolve_axis(blurred, 1, kern(si, g.spacing[1]));
  blurred = eikon::filters::convolve_axis(blurred, 2, kern(sz, g.spacing[2]));
  auto out = eikon::grid::resample(blurred, p.grid, nullptr, false, p.background);
  for (auto& x : out.data) x = float(x + p.noise * rng.normal());
  return out;
}

// ---- MR: time of flight ------------------------------------------------------------------------------------

struct SpoiledGradientEcho {
  double tr = 0.025;   // s
  double alpha = 0.35; // rad
};

// M_z after n pulses from M_z(0) = m0 (in units of M₀).
inline double longitudinal(double n, double t1, const SpoiledGradientEcho& s, double m0 = 1) {
  const double E1 = std::exp(-s.tr / t1), ca = std::cos(s.alpha);
  const double mss = (1 - E1) / (1 - E1 * ca);
  return mss + (m0 - mss) * std::pow(E1 * ca, n);
}
inline double steady_signal(double t1, const SpoiledGradientEcho& s) { return std::sin(s.alpha) * longitudinal(INFINITY, t1, s); }

// The signal of blood that has travelled depth d (mm) into the slab at speed v (mm/s), entering fresh (m0 = 1) or
// saturated (m0 = 0).
inline double inflow_signal(double depth, double speed, double t1, const SpoiledGradientEcho& s, double m0 = 1) {
  const double n = speed > 0 ? std::max(0.0, depth) / (speed * s.tr) : INFINITY;
  return std::sin(s.alpha) * longitudinal(n, t1, s, m0);
}

// A vessel for TOF: its occupancy, the blood's T1, its speed, and whether it enters through a saturation band.
struct FlowingBlood {
  Image<float> occupancy;
  double t1 = 1.6, speed = 300;
  bool saturated = false;
};

struct TofProtocol {
  Geometry grid;               // the acquisition grid
  SpoiledGradientEcho sequence;
  Vec3 slab_entry{0, 0, 0};    // a point on the slab's entry face (mm)
  Vec3 slab_axis{0, 0, 1};     // the direction blood travels into the slab (unit)
  double tissue_t1 = 1.0;      // s, the stationary tissue
  double noise = 0.01;         // in units of M₀, per channel (Rician)
  double fwhm = 0.8;           // mm, isotropic blur
};

// The TOF image: stationary tissue at its steady state where `tissue` occupies, blood by its inflow.
inline Image<float> tof(const Image<float>& tissue, const std::vector<FlowingBlood>& blood, const TofProtocol& p, calx::random::Stream& rng) {
  const auto& g = tissue.geometry;
  Image<float> truth(g, 0.0f);
  const double st = steady_signal(p.tissue_t1, p.sequence);
  for (std::size_t v = 0; v < truth.data.size(); v++) truth.data[v] = float(st * tissue.data[v]);
  for (const auto& b : blood)
    for (std::size_t v = 0; v < truth.data.size(); v++) {
      if (b.occupancy.data[v] <= 0) continue;
      const auto ijk = g.unlinear(v);
      const auto x = g.world(double(ijk[0]), double(ijk[1]), double(ijk[2]));
      const double depth = (x[0] - p.slab_entry[0]) * p.slab_axis[0] + (x[1] - p.slab_entry[1]) * p.slab_axis[1] + (x[2] - p.slab_entry[2]) * p.slab_axis[2];
      const double sig = inflow_signal(depth, b.speed, b.t1, p.sequence, b.saturated ? 0 : 1);
      truth.data[v] += float((sig - truth.data[v]) * b.occupancy.data[v]);
    }
  const double sg = p.fwhm / (2 * std::sqrt(2 * std::log(2.0)));
  auto out = eikon::grid::resample(eikon::filters::smooth(truth, sg), p.grid, nullptr, false, 0);
  eikon::phantoms::rician_noise(out, p.noise, rng);
  return out;
}

}  // namespace medical_imaging::modalities
