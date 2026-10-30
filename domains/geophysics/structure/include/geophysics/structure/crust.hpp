// The Earth's crust beneath a place, from CRUST1.0 (Laske, Masters, Ma and Pasyanos 2013): its layers
// from the surface to the Moho, with their P- and S-wave speeds and densities, and the top of the
// mantle beneath. The C++ of the domain's `structure` beneath its JavaScript (js/crust.js), held to the
// same checks.
//
// CRUST1.0 describes each cell of 1° by nine layers: water, ice, upper, middle and lower sediments,
// upper, middle and lower crystalline crust, and the mantle below the Moho. A cell's values are its
// average, so a place takes the values of the cell that contains it, not an interpolation between
// cells (the model's readme: "please use the midpoint of this cell"). A layer's base is the top of the
// next layer present; a layer of zero thickness is absent. Heights are above sea level, negative
// below; the model does not say whether its latitudes are geodetic or geocentric, a difference of at
// most 0.19°, within its cells.
//
// The data library's product crust/1 (data/products/crust_model.py) holds the model as 16-bit
// integers: for each of the quantities top (units of 10 m), vp, vs (m/s) and rho (kg/m³), in that
// order, nine blocks, one for each layer, of 180 × 360 cells, rows from 89.5° N and columns from
// 179.5° W, little-endian. Crust1::read_product reads it from a path that is passed in.

#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace geophysics::structure {

inline constexpr std::array<std::string_view, 9> CRUST1_LAYERS{
    "water", "ice", "upper sediments", "middle sediments", "lower sediments",
    "upper crystalline crust", "middle crystalline crust", "lower crystalline crust", "mantle"};

enum class CrustQuantity { top, vp, vs, rho };

// One layer present in a column: heights of its top and bottom above sea level (m), the speeds of P
// and S waves (m/s), its density (kg/m³), and its index among CRUST1_LAYERS.
struct CrustLayer {
  int index = 0;
  std::string_view name;
  double top_m = 0, bottom_m = 0;
  double vp_m_s = 0, vs_m_s = 0, density_kg_m3 = 0;
  double thickness_m() const { return top_m - bottom_m; }
};

struct MantleTop {
  double vp_m_s = 0, vs_m_s = 0, density_kg_m3 = 0;
};

// The column of a cell.
struct CrustColumn {
  std::vector<CrustLayer> layers;  // the layers present, from the top
  double surface_m = 0;            // the column's top: of the water, the ice or the solid ground
  double solid_m = 0;              // the solid surface: the top of the upper sediments
  double moho_depth_m = 0;         // the Moho's depth below sea level
  double crust_m = 0;              // the thickness of the crust, the solid surface to the Moho
  double sediments_m = 0, water_m = 0, ice_m = 0;
  MantleTop mantle;
  // Oceanic crust: under more than 500 m of water, and thinner than 15 km of rock.
  bool oceanic = false;
  // The layer that holds a height above sea level (m): an index into `layers`, or −1 above the
  // column's top, or layers.size() below the Moho.
  int layer_at(double height_m) const {
    if (height_m > surface_m) return -1;
    for (std::size_t k = 0; k < layers.size(); k++)
      if (height_m >= layers[k].bottom_m) return int(k);
    return int(layers.size());
  }
};

// The column from a getter of the model's values at a cell, get(quantity, layer) in the model's own
// units (km above sea level, km/s, g/cm³), as js/crust.js's column(get).
inline CrustColumn crust_column(const std::function<double(CrustQuantity, int)>& get) {
  std::array<double, 9> tops{};
  for (int k = 0; k < 9; k++) tops[k] = get(CrustQuantity::top, k);
  CrustColumn c;
  for (int k = 0; k < 8; k++) {
    const double thickness = tops[k] - tops[k + 1];
    if (thickness <= 1e-6) continue;
    c.layers.push_back({k, CRUST1_LAYERS[k], tops[k] * 1e3, tops[k + 1] * 1e3, get(CrustQuantity::vp, k) * 1e3, get(CrustQuantity::vs, k) * 1e3,
                        get(CrustQuantity::rho, k) * 1e3});
  }
  c.surface_m = tops[0] * 1e3;
  c.solid_m = tops[2] * 1e3;
  c.moho_depth_m = -tops[8] * 1e3;
  c.crust_m = (tops[2] - tops[8]) * 1e3;
  c.sediments_m = (tops[2] - tops[5]) * 1e3;
  c.water_m = (tops[0] - tops[1]) * 1e3;
  c.ice_m = (tops[1] - tops[2]) * 1e3;
  c.mantle = {get(CrustQuantity::vp, 8) * 1e3, get(CrustQuantity::vs, 8) * 1e3, get(CrustQuantity::rho, 8) * 1e3};
  c.oceanic = c.water_m > 500 && c.crust_m < 15e3;
  return c;
}

// The grid of CRUST1.0: 1° cells, rows from 89.5° N southward, columns from 179.5° W eastward.
struct CrustGrid {
  double lon0_deg = -179.5, dlon_deg = 1, lat0_deg = 89.5, dlat_deg = 1;
  int nlon = 360, nlat = 180;
  // The index of the cell containing (lat, lon), degrees, as js/crust.js's cellIndex: the poles and
  // the antimeridian fall in the cells beside them.
  int cell(double lat_deg, double lon_deg) const {
    const int r = std::min(nlat - 1, std::max(0, int(std::floor((lat0_deg + dlat_deg / 2 - lat_deg) / dlat_deg))));
    const int c0 = int(std::floor((lon_deg - (lon0_deg - dlon_deg / 2)) / dlon_deg));
    const int c = ((c0 % nlon) + nlon) % nlon;
    return r * nlon + c;
  }
  // The centre of a cell (degrees).
  std::array<double, 2> centre(int cell) const {
    return {lat0_deg - (cell / nlon) * dlat_deg, lon0_deg + (cell % nlon) * dlon_deg};
  }
};

class Crust1 {
 public:
  static constexpr int cells = 180 * 360;

  // The product crust/1 at `path` (the .bin beside crust-crust1.json).
  static Crust1 read_product(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("geophysics::structure: cannot read " + path);
    std::vector<std::int16_t> raw(std::size_t(4) * 9 * cells);
    in.read(reinterpret_cast<char*>(raw.data()), std::streamsize(raw.size() * sizeof(std::int16_t)));
    if (in.gcount() != std::streamsize(raw.size() * sizeof(std::int16_t)) || in.peek() != std::char_traits<char>::eof())
      throw std::runtime_error("geophysics::structure: " + path + " is not CRUST1.0's product crust/1 (4 × 9 × 64,800 16-bit integers)");
    Crust1 m;
    constexpr std::array<double, 4> scale{0.01, 0.001, 0.001, 0.001};  // to km, km/s, km/s, g/cm³
    m.values_.resize(raw.size());
    for (std::size_t i = 0; i < raw.size(); i++) {
      std::int16_t v = raw[i];
      if constexpr (std::endian::native == std::endian::big) v = std::int16_t((std::uint16_t(v) >> 8) | (std::uint16_t(v) << 8));
      m.values_[i] = v * scale[i / (std::size_t(9) * cells)];
    }
    return m;
  }

  const CrustGrid& grid() const { return grid_; }

  // A value at a cell, in the model's own units (km above sea level, km/s, g/cm³).
  double value(CrustQuantity q, int layer, int cell) const {
    return values_[(std::size_t(static_cast<int>(q)) * 9 + std::size_t(layer)) * cells + std::size_t(cell)];
  }

  CrustColumn column_of_cell(int cell) const {
    return crust_column([&](CrustQuantity q, int k) { return value(q, k, cell); });
  }
  // The column beneath a place (degrees).
  CrustColumn column(double lat_deg, double lon_deg) const { return column_of_cell(grid_.cell(lat_deg, lon_deg)); }

  std::string_view source() const {
    return "CRUST1.0 (Laske, Masters, Ma and Pasyanos 2013), 1° cells averaged";
  }

 private:
  CrustGrid grid_;
  std::vector<double> values_;
};

}  // namespace geophysics::structure
