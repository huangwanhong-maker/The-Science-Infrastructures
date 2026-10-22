// Relative vorticity and divergence of the wind, in C++: the algorithm of
// vorticity.js and vorticity.py, as a program the server compiles and runs
// (plugins/cpp/plugin.hpp). Centred differences on the sphere,
//
//   ζ = (1 / (a cos φ)) [ ∂v/∂λ − ∂(u cos φ)/∂φ ],    δ = (1 / (a cos φ)) [ ∂u/∂λ + ∂(v cos φ)/∂φ ],
//
// rows running north to south; the poles and, on a grid that does not wrap,
// the edge columns are left without a value.

#include "cpp/plugin.hpp"

struct Vorticity : plugin::Diagnostic {
  static constexpr double RADIUS = 6371008.8, DEG = 3.14159265358979323846 / 180;

  Vorticity() {
    id = "vorticity-cpp";
    title = "Vorticity and divergence (C++)";
    description = "Relative vorticity and divergence of the wind, by centred differences on the sphere, computed by a C++ program in the server.";
    needs = {{"wind", "Wind", {"eastward_wind", "northward_wind"}, "", 0}};
    outputs = {{"zeta", "Relative vorticity", "atmosphere_relative_vorticity", "s-1", "", false},
               {"delta", "Divergence of the wind", "divergence_of_wind", "s-1", "", false}};
    show = {{"surface", "zeta"}};
  }

  std::map<std::string, std::vector<float>> compute(const plugin::Inputs& inputs, const plugin::Params&, int index) override {
    const plugin::Input& wind = inputs.at("wind");
    const plugin::Grid& g = wind.grid;
    const float* u = wind.slice(0, index);
    const float* v = wind.slice(1, index);
    const float nan = std::nanf("");
    std::vector<float> zeta(g.size(), nan), delta(g.size(), nan);
    const double dx = 2 * g.dlon * DEG * RADIUS, dy = 2 * g.dlat * DEG * RADIUS;
    for (int r = 1; r + 1 < g.nlat; r++) {
      const double lat = g.lat0 - r * g.dlat, cosl = std::cos(lat * DEG);
      if (cosl < 1e-6) continue;
      const double cosN = std::cos((lat + g.dlat) * DEG), cosS = std::cos((lat - g.dlat) * DEG);
      for (int c = 0; c < g.nlon; c++) {
        const int e = c + 1 < g.nlon ? c + 1 : (g.periodic ? 0 : -1), w = c > 0 ? c - 1 : (g.periodic ? g.nlon - 1 : -1);
        if (e < 0 || w < 0) continue;
        const int i = r * g.nlon + c, iN = (r - 1) * g.nlon + c, iS = (r + 1) * g.nlon + c;
        const double dvdl = (v[r * g.nlon + e] - v[r * g.nlon + w]) / dx, dudl = (u[r * g.nlon + e] - u[r * g.nlon + w]) / dx;
        const double ducos = (u[iN] * cosN - u[iS] * cosS) / dy, dvcos = (v[iN] * cosN - v[iS] * cosS) / dy;
        zeta[i] = float((dvdl - ducos) / cosl);
        delta[i] = float((dudl + dvcos) / cosl);
      }
    }
    return {{"zeta", zeta}, {"delta", delta}};
  }
};

PLUGIN_MAIN(Vorticity)
