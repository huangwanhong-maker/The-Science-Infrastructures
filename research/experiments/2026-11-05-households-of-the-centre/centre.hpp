// The households of the centre: the City of Melbourne's dwellings, rents and census read from the data library's
// product (data/products/housing_records.py), and the city, the households and the simulation built on them with
// urban-science's housing.
//
// Conventions. Time in quarters: tick 0 is 2002Q3 (the first quarter after CLUE's first census year's middle), tick 92
// is 2025Q3 (the last quarter of Homes Victoria's report). A CLUE census year's stock is placed at the year's middle
// and interpolated between years. Rents are A$ a week, nominal. Places in the product's local frame (metres east and
// north of its origin).

#pragma once

#include <urban/housing/calibration.hpp>
#include <urban/housing/households.hpp>
#include <urban/housing/market.hpp>

#include <calx/random/categorical.hpp>
#include <calx/random/stream.hpp>
#include <science/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace centre {

namespace uh = urban::housing;
using science::Json;

constexpr int first_year = 2002, last_year = 2024, years = last_year - first_year + 1;
constexpr int ticks = 93;     // 2002Q3 … 2025Q3
constexpr int q0 = 10;        // the report's quarter index of tick 0 (2000Q1 is 0)
constexpr int census_tick = 76;  // 2021Q3: the Census of 10 August 2021

inline int quarter_of_tick(int t) { return q0 + t; }
inline std::string quarter_label(int q) { return std::to_string(2000 + q / 4) + "Q" + std::to_string(q % 4 + 1); }
inline double quarter_time(int q) { return 2000 + q / 4 + (q % 4 + 0.5) / 4; }  // its middle, in years

inline Json read_json(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return Json::parse(ss.str());
}

inline void write_text(const std::string& path, const std::string& text) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << text;
}

// ---- geometry ---------------------------------------------------------------------------------------------------

using Ring = std::vector<std::array<double, 2>>;

inline bool inside(const Ring& r, double x, double y) {
  bool in = false;
  for (std::size_t i = 0, j = r.size() - 1; i < r.size(); j = i++) {
    const auto& a = r[i];
    const auto& b = r[j];
    if ((a[1] > y) != (b[1] > y) && x < a[0] + (y - a[1]) * (b[0] - a[0]) / (b[1] - a[1])) in = !in;
  }
  return in;
}

inline std::vector<Ring> rings_of(const Json& j) {
  std::vector<Ring> out;
  for (const auto& ring : j.as_array()) {
    Ring r;
    for (const auto& p : ring.as_array()) r.push_back({p[0].as_number(), p[1].as_number()});
    out.push_back(std::move(r));
  }
  return out;
}

// ---- the data ---------------------------------------------------------------------------------------------------

struct Block {
  int id = 0;
  std::size_t area = 0;
  double x = 0, y = 0;
  std::vector<Ring> rings;
};

struct BuildingRecord {
  std::string property;
  std::size_t zone = 0, area = 0;
  double x = 0, y = 0;
  std::array<std::array<int, years>, uh::types> stock{};  // by type and census year
  bool development = false;   // a building of the Development Activity Monitor's pipeline, not in CLUE
};

struct Series {
  std::vector<double> count, median;  // by report quarter; NaN where none
};

struct Sa1 {
  std::string code, sa2;
  std::size_t zone = 0, area = 0;
  double x = 0, y = 0;
  std::vector<Ring> rings;
  std::vector<double> g02, g33, g35, g36, g37, g40, g44;
};

struct Development {
  std::string status;
  int year = 0;  // completed, or 0
  std::size_t building = 0;
  double dwellings = 0, student = 0, x = 0, y = 0;
};

struct Data {
  std::vector<std::string> areas;            // CLUE's small areas, sorted
  std::vector<Block> blocks;                 // the zones, by block id
  std::vector<BuildingRecord> buildings;
  std::vector<std::string> quarters;
  std::map<std::string, std::map<std::string, Series>> rents;  // group → sheet → series
  std::vector<std::string> group_of_area;    // the report's suburb group of each area
  std::vector<std::string> census_columns_g33, census_columns_g37, census_columns_g36, census_columns_g44, census_columns_g35, census_columns_g02;
  std::vector<std::vector<double>> g33_bands;
  std::vector<Sa1> sa1;
  std::vector<Development> developments;
  Json product_sources;

  std::size_t area_index(const std::string& name) const {
    const auto it = std::find(areas.begin(), areas.end(), name);
    if (it == areas.end()) throw std::runtime_error("an unknown small area: " + name);
    return static_cast<std::size_t>(it - areas.begin());
  }
  std::size_t zone_of_point(double x, double y) const {
    for (std::size_t z = 0; z < blocks.size(); z++)
      for (const auto& r : blocks[z].rings)
        if (inside(r, x, y)) return z;
    std::size_t best = 0;
    double d = 1e300;
    for (std::size_t z = 0; z < blocks.size(); z++) {
      const double e = std::hypot(blocks[z].x - x, blocks[z].y - y);
      if (e < d) d = e, best = z;
    }
    return best;
  }
};

inline double num(const Json& j) { return j.is_null() ? std::nan("") : j.as_number(); }

inline Data load(const std::string& path) {
  const Json p = read_json(path);
  Data D;
  std::set<std::string> names;
  for (const auto& b : p["blocks"]["rows"].as_array()) names.insert(b[1].as_string());
  D.areas.assign(names.begin(), names.end());
  std::map<int, std::size_t> zone_of_block;
  for (const auto& b : p["blocks"]["rows"].as_array()) {
    Block k;
    k.id = static_cast<int>(b[0].as_number()), k.area = D.area_index(b[1].as_string()), k.x = b[2].as_number(), k.y = b[3].as_number();
    k.rings = rings_of(b[5]);
    zone_of_block[k.id] = D.blocks.size();
    D.blocks.push_back(std::move(k));
  }
  std::map<std::string, std::size_t> building_of_property;
  for (const auto& b : p["buildings"]["rows"].as_array()) {
    BuildingRecord r;
    r.property = b[0].is_string() ? b[0].as_string() : std::to_string(static_cast<long long>(b[0].as_number()));
    r.x = b[4].as_number(), r.y = b[5].as_number();
    const auto it = b[2].is_null() ? zone_of_block.end() : zone_of_block.find(static_cast<int>(b[2].as_number()));
    r.zone = it != zone_of_block.end() ? it->second : D.zone_of_point(r.x, r.y);
    r.area = D.blocks[r.zone].area;  // a building's area is its block's
    for (std::size_t t = 0; t < uh::types; t++)
      for (int y = 0; y < years; y++) r.stock[t][static_cast<std::size_t>(y)] = static_cast<int>(b[6][t][static_cast<std::size_t>(y)].as_number());
    building_of_property[r.property] = D.buildings.size();
    D.buildings.push_back(std::move(r));
  }
  for (const auto& q : p["rents"]["quarters"].as_array()) D.quarters.push_back(q.as_string());
  for (const auto& [g, sheets] : p["rents"]["groups"].as_object())
    for (const auto& [sheet, s] : sheets.as_object()) {
      Series x;
      for (const auto& v : s["count"].as_array()) x.count.push_back(num(v));
      for (const auto& v : s["median"].as_array()) x.median.push_back(num(v));
      D.rents[g][sheet] = std::move(x);
    }
  D.group_of_area.resize(D.areas.size());
  for (std::size_t a = 0; a < D.areas.size(); a++) D.group_of_area[a] = p["rents"]["group_of_small_area"][D.areas[a]].as_string();
  const auto& C = p["census"];
  for (const auto& c : C["G33"].as_array()) D.census_columns_g33.push_back(c.as_string());
  for (const auto& c : C["G37"].as_array()) D.census_columns_g37.push_back(c.as_string());
  for (const auto& c : C["G36"].as_array()) D.census_columns_g36.push_back(c.as_string());
  for (const auto& c : C["G44"].as_array()) D.census_columns_g44.push_back(c.as_string());
  for (const auto& c : C["G35"].as_array()) D.census_columns_g35.push_back(c.as_string());
  for (const auto& c : C["G02"].as_array()) D.census_columns_g02.push_back(c.as_string());
  for (const auto& band : C["G33_bands"].as_array()) D.g33_bands.push_back({num(band[0]), num(band[1])});
  for (const auto& r : C["rows"].as_array()) {
    Sa1 s;
    s.code = r[0].as_string(), s.sa2 = r[1].as_string();
    const auto it = zone_of_block.find(static_cast<int>(r[2].as_number()));
    s.zone = it->second, s.area = D.blocks[s.zone].area, s.x = r[4].as_number(), s.y = r[5].as_number();
    s.rings = rings_of(r[6]);
    const auto vec = [&](const Json& j) {
      std::vector<double> v;
      for (const auto& x : j.as_array()) v.push_back(num(x));
      return v;
    };
    s.g02 = vec(r[7]), s.g33 = vec(r[8]), s.g35 = vec(r[9]), s.g36 = vec(r[10]), s.g37 = vec(r[11]), s.g40 = vec(r[12]), s.g44 = vec(r[13]);
    D.sa1.push_back(std::move(s));
  }
  // The pipeline's buildings: a development's first property id that CLUE holds, else a building of its own.
  for (const auto& r : p["developments"]["rows"].as_array()) {
    Development d;
    d.status = r[0].as_string();
    d.year = r[1].is_null() ? 0 : static_cast<int>(r[1].as_number());
    d.x = r[5].as_number(), d.y = r[6].as_number();
    d.dwellings = std::isfinite(num(r[8])) ? num(r[8]) : 0, d.student = std::isfinite(num(r[13])) ? num(r[13]) : 0;
    if (!(d.dwellings > 0 || d.student > 0)) continue;
    const bool pipeline = d.status == "UNDER CONSTRUCTION" || d.status == "APPROVED" || (d.status == "COMPLETED" && d.year >= 2025);
    if (!pipeline) continue;
    std::size_t b = D.buildings.size();
    for (const auto& id : r[4].as_array()) {
      const auto it = building_of_property.find(id.as_string());
      if (it != building_of_property.end()) { b = it->second; break; }
    }
    if (b == D.buildings.size()) {
      BuildingRecord n;
      n.property = r[4].size() ? "dam:" + r[4][0].as_string() : "dam:" + std::to_string(D.developments.size());
      n.x = d.x, n.y = d.y, n.development = true;
      const auto it = r[3].is_null() ? zone_of_block.end() : zone_of_block.find(static_cast<int>(r[3].as_number()));
      n.zone = it != zone_of_block.end() ? it->second : D.zone_of_point(n.x, n.y);
      n.area = D.blocks[n.zone].area;
      D.buildings.push_back(std::move(n));
    }
    d.building = b;
    D.developments.push_back(d);
  }
  D.product_sources = p["sources"];
  return D;
}

// The interpolated stock of a building's type at a time (years): CLUE's census year y at y + ½.
inline double stock_at(const BuildingRecord& b, std::size_t t, double time) {
  const double u = std::clamp(time - 0.5 - first_year, 0.0, static_cast<double>(years - 1));
  const std::size_t i = static_cast<std::size_t>(std::floor(u));
  const double f = u - static_cast<double>(i);
  if (i + 1 >= static_cast<std::size_t>(years)) return b.stock[t][static_cast<std::size_t>(years - 1)];
  return (1 - f) * b.stock[t][i] + f * b.stock[t][i + 1];
}

// ---- the report's indices ---------------------------------------------------------------------------------------

// A fixed-weight geometric index of a group's medians over the sheets given: each sheet weighted by its lettings over
// the whole report, renormalised over the sheets that publish a median in the quarter; NaN where none does.
inline std::vector<double> rent_index(const Data& D, const std::string& group, const std::vector<std::string>& sheets, double* lettings = nullptr) {
  const std::size_t Q = D.quarters.size();
  std::vector<double> w;
  const auto& g = D.rents.at(group);
  for (const auto& s : sheets) {
    double total = 0;
    for (std::size_t q = 0; q < Q; q++)
      if (std::isfinite(g.at(s).count[q]) && std::isfinite(g.at(s).median[q])) total += g.at(s).count[q];
    w.push_back(total);
  }
  std::vector<double> out(Q, std::nan(""));
  for (std::size_t q = 0; q < Q; q++) {
    double num = 0, den = 0;
    for (std::size_t k = 0; k < sheets.size(); k++) {
      const double m = g.at(sheets[k]).median[q];
      if (std::isfinite(m) && m > 0 && w[k] > 0) num += w[k] * std::log(m), den += w[k];
    }
    if (den > 0) out[q] = std::exp(num / den);
  }
  if (lettings) {
    double total = 0;
    for (const double x : w) total += x;
    *lettings = total / static_cast<double>(Q) * 4;  // a year
  }
  return out;
}

inline const std::vector<std::string> flats = {"1 bedroom flat", "2 bedroom flat", "3 bedroom flat"};
inline const std::vector<std::string> houses = {"2 bedroom house", "3 bedroom house", "4 bedroom house"};

// ---- the census -------------------------------------------------------------------------------------------------

inline std::size_t column(const std::vector<std::string>& cols, const std::string& name) {
  const auto it = std::find(cols.begin(), cols.end(), name);
  if (it == cols.end()) throw std::runtime_error("no census column " + name);
  return static_cast<std::size_t>(it - cols.begin());
}

// Income classes by weekly household income: low below A$650, middle A$650-1,999, high A$2,000 or more.
constexpr std::size_t classes = 3;
inline std::size_t class_of_band(double lo) { return lo < 650 ? 0 : (lo < 2000 ? 1 : 2); }

struct CensusTotals {
  std::array<double, classes> households{};   // by income class (stated incomes)
  std::array<double, classes> median_income{};  // A$ a week, interpolated within the class's bands
  double owners = 0, renters = 0, other = 0;
  double owners_house = 0, renters_house = 0, owners_flat = 0, renters_flat = 0;
  double one_person = 0, all_sizes = 0;
  double moved = 0, stayed = 0;               // persons at a different, the same address a year before
  double occupied = 0, unoccupied = 0;
};

// The households of SA1 s by income class (stated incomes only).
inline std::array<double, classes> income_classes(const Data& D, const Sa1& s) {
  std::array<double, classes> c{};
  for (std::size_t k = 0; k < D.g33_bands.size(); k++) {
    const double v = s.g33[k];
    if (std::isfinite(v)) c[class_of_band(D.g33_bands[k][0])] += v;
  }
  return c;
}

inline CensusTotals census_totals(const Data& D) {
  CensusTotals T;
  std::vector<double> bands(D.g33_bands.size(), 0.0);
  const auto get = [](const std::vector<double>& v, std::size_t i) { return std::isfinite(v[i]) ? v[i] : 0.0; };
  for (const auto& s : D.sa1) {
    for (std::size_t k = 0; k < bands.size(); k++) bands[k] += get(s.g33, k);
    const auto& c37 = D.census_columns_g37;
    T.owners += get(s.g37, column(c37, "O_OR_Total")) + get(s.g37, column(c37, "O_MTG_Total"));
    T.renters += get(s.g37, column(c37, "R_Tot_Total"));
    T.other += get(s.g37, column(c37, "Oth_ten_type_Total"));
    T.owners_house += get(s.g37, column(c37, "O_OR_DS_Sep_house")) + get(s.g37, column(c37, "O_OR_DS_SemiD_ro_or_tce_h_th")) + get(s.g37, column(c37, "O_MTG_DS_Sep_house")) +
                      get(s.g37, column(c37, "O_MTG_DS_SemiD_ro_or_tce_h_th"));
    T.renters_house += get(s.g37, column(c37, "R_Tot_DS_Sep_house")) + get(s.g37, column(c37, "R_Tot_DS_SemiD_ro_or_tce_h_th"));
    T.owners_flat += get(s.g37, column(c37, "O_OR_DS_Flat_apart")) + get(s.g37, column(c37, "O_MTG_DS_Flat_apart"));
    T.renters_flat += get(s.g37, column(c37, "R_Tot_DS_Flat_apart"));
    T.one_person += get(s.g35, column(D.census_columns_g35, "Num_Psns_UR_1_Total"));
    T.all_sizes += get(s.g35, column(D.census_columns_g35, "Total_Total"));
    T.moved += get(s.g44, column(D.census_columns_g44, "Difnt_Usl_add_1_yr_ago_Tot_P"));
    T.stayed += get(s.g44, column(D.census_columns_g44, "Sme_Usl_ad_1_yr_ago_as_2021_P"));
    T.occupied += get(s.g36, column(D.census_columns_g36, "OPDs_Tot_OPDs_Dwellings"));
    T.unoccupied += get(s.g36, column(D.census_columns_g36, "Unoccupied_PDs_Dwgs"));
  }
  for (std::size_t k = 0; k < bands.size(); k++) T.households[class_of_band(D.g33_bands[k][0])] += bands[k];
  // Each class's median income, linear within its bands (nil and negative incomes at 0; the open band A$4,000 or more
  // taken as A$4,000-6,000, an assumption).
  for (std::size_t c = 0; c < classes; c++) {
    double half = T.households[c] / 2, run = 0;
    for (std::size_t k = 0; k < bands.size(); k++) {
      if (class_of_band(D.g33_bands[k][0]) != c) continue;
      const double lo = std::max(0.0, D.g33_bands[k][0]), hi = std::isfinite(D.g33_bands[k][1]) ? D.g33_bands[k][1] + 1 : 6000;
      if (run + bands[k] >= half && bands[k] > 0) {
        T.median_income[c] = lo + (hi - lo) * (half - run) / bands[k];
        break;
      }
      run += bands[k];
    }
  }
  return T;
}

// ---- the city ---------------------------------------------------------------------------------------------------

// The centre's anchor, Melbourne Town Hall at the corner of Swanston and Collins Streets (37.81502 S, 144.96644 E;
// an assumption of the model, as the place to which distance is reckoned), in the product's frame by the
// equirectangular approximation (under a metre's error within the municipality).
inline std::array<double, 2> anchor() {
  const double lat0 = -37.8085, lon0 = 144.9525, R = 6378137.0, k = 3.14159265358979323846 / 180;
  return {(144.96644 - lon0) * k * R * std::cos(lat0 * k), (-37.81502 - lat0) * k * R};
}

inline uh::City make_city(const Data& D) {
  uh::City c;
  c.areas = D.areas.size(), c.zones = D.blocks.size();
  const auto a = anchor();
  for (const auto& b : D.buildings) c.buildings.push_back({b.area, b.zone, b.x, b.y, std::hypot(b.x - a[0], b.y - a[1]) / 1000});
  return c;
}

// The capacities at tick t: CLUE's stock interpolated, with the pipeline's dwellings added by `extra` (by building,
// type).
inline std::vector<std::array<int, uh::types>> capacities(const Data& D, int t, const std::vector<std::array<double, uh::types>>* extra = nullptr) {
  const double time = quarter_time(quarter_of_tick(t));
  std::vector<std::array<int, uh::types>> cap(D.buildings.size());
  for (std::size_t b = 0; b < D.buildings.size(); b++)
    for (std::size_t k = 0; k < uh::types; k++) {
      double n = D.buildings[b].development ? 0.0 : stock_at(D.buildings[b], k, time);
      if (extra) n += (*extra)[b][k];
      cap[b][k] = static_cast<int>(std::lround(n));
    }
  return cap;
}

}  // namespace centre
