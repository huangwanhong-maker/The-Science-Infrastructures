// The geoid as a grid of undulations N above the ellipsoid, and heights above it: the orthometric
// (mean-sea-level) height H = h − N of a point at ellipsoidal height h.
//
// The grid is EGM96's at 15′ as GeographicLib distributes it (the data library's product
// geoid/egm96-15.pgm, from NGA's EGM96 to degree and order 360; Lemoine et al. 1998, NASA
// TP-1998-206861): a binary PGM of 1440 × 721 16-bit values from 90° N and 0° E, point-registered,
// with N = Offset + Scale × value in metres (Offset −108, Scale 0.003: a quantum of 3 mm). It is
// interpolated bilinearly, as PROJ does with NGA's grid (the checks hold the two together), or by
// cubic convolution over the 4 × 4 nodes about the point (Keys 1981, IEEE Trans. ASSP 29, 1153,
// with a = −1/2), smooth across the nodes. Neither reproduces the spherical-harmonic model itself
// between the nodes, where the geoid bends sharply (GeographicLib's documentation of its Geoid class
// tabulates the errors of interpolating each of its grids); the 3 mm quantum is the grid's own.

#pragma once

#include <cmath>
#include <cstdint>
#include <fstream>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace geodesy::geoid {

enum class Interpolation { bilinear, cubic };

class GeoidGrid {
 public:
  // A grid read from a binary PGM with GeographicLib's Offset and Scale comments.
  static GeoidGrid read_pgm(const std::string& path, std::string name = "EGM96, 15′ grid") {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::invalid_argument("geodesy: cannot read the geoid grid " + path);
    std::string data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    GeoidGrid g;
    g.name_ = std::move(name);
    std::size_t pos = 0;
    std::vector<std::string> tokens;
    double offset = NAN, scale = NAN;
    while (tokens.size() < 4) {
      const std::size_t end = data.find('\n', pos);
      if (end == std::string::npos) throw std::invalid_argument("geodesy: the geoid grid's header is cut short");
      std::string line = data.substr(pos, end - pos);
      pos = end + 1;
      if (!line.empty() && line[0] == '#') {
        std::istringstream c(line.substr(1));
        std::string key;
        double v;
        if (c >> key >> v) {
          if (key == "Offset") offset = v;
          if (key == "Scale") scale = v;
        }
        continue;
      }
      std::istringstream t(line);
      std::string w;
      while (t >> w) tokens.push_back(w);
    }
    if (tokens[0] != "P5") throw std::invalid_argument("geodesy: the geoid grid is not a binary PGM");
    g.width_ = std::stoi(tokens[1]);
    g.height_ = std::stoi(tokens[2]);
    if (std::isnan(offset) || std::isnan(scale)) throw std::invalid_argument("geodesy: the geoid grid names no Offset and Scale");
    if (data.size() < pos + 2ull * g.width_ * g.height_) throw std::invalid_argument("geodesy: the geoid grid is cut short");
    g.n_.resize(std::size_t(g.width_) * g.height_);
    const auto* b = reinterpret_cast<const unsigned char*>(data.data() + pos);
    for (std::size_t i = 0; i < g.n_.size(); i++) g.n_[i] = offset + scale * double((b[2 * i] << 8) | b[2 * i + 1]);
    g.step_ = 360.0 / g.width_;
    if (std::fabs(g.step_ * (g.height_ - 1) - 180.0) > 1e-9) throw std::invalid_argument("geodesy: the geoid grid is not a whole sphere");
    return g;
  }

  const std::string& name() const { return name_; }
  int width() const { return width_; }
  int height() const { return height_; }
  double step_degrees() const { return step_; }

  // The node (row from 90° N, column from 0° E).
  double node(int row, int col) const {
    row = std::max(0, std::min(height_ - 1, row));
    col = ((col % width_) + width_) % width_;
    return n_[std::size_t(row) * width_ + col];
  }

  // The undulation N (m) at a latitude and longitude (rad).
  double undulation(double lat, double lon, Interpolation how = Interpolation::cubic) const {
    const double deg = 180 / std::numbers::pi;
    double x = std::fmod(lon * deg, 360.0);
    if (x < 0) x += 360;
    const double fx = x / step_, fy = (90 - lat * deg) / step_;
    const int c = int(std::floor(fx)), r = std::min(height_ - 2, std::max(0, int(std::floor(fy))));
    const double u = fx - c, v = std::min(1.0, std::max(0.0, fy - r));
    if (how == Interpolation::bilinear)
      return (1 - v) * ((1 - u) * node(r, c) + u * node(r, c + 1)) + v * ((1 - u) * node(r + 1, c) + u * node(r + 1, c + 1));
    // Cubic convolution (Keys 1981, a = −1/2); beyond a pole the rows are reflected across it, at the
    // longitude half a turn away.
    auto at = [&](int rr, int cc) {
      if (rr < 0) return node(-rr, cc + width_ / 2);
      if (rr > height_ - 1) return node(2 * (height_ - 1) - rr, cc + width_ / 2);
      return node(rr, cc);
    };
    auto w = [](double t) {
      t = std::fabs(t);
      if (t <= 1) return (1.5 * t - 2.5) * t * t + 1;
      if (t < 2) return ((-0.5 * t + 2.5) * t - 4) * t + 2;
      return 0.0;
    };
    double sum = 0;
    for (int j = -1; j <= 2; j++) {
      double row = 0;
      for (int i = -1; i <= 2; i++) row += w(u - i) * at(r + j, c + i);
      sum += w(v - j) * row;
    }
    return sum;
  }

  // The orthometric height H = h − N (m) of a point at ellipsoidal height h.
  double orthometric_height(double lat, double lon, double h, Interpolation how = Interpolation::cubic) const {
    return h - undulation(lat, lon, how);
  }

 private:
  std::string name_;
  int width_ = 0, height_ = 0;
  double step_ = 0;
  std::vector<double> n_;
};

}  // namespace geodesy::geoid
