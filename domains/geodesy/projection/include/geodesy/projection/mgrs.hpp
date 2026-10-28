// The Military Grid Reference System (MGRS): a position on the universal grids written as its grid
// zone, its 100 km square and its digits (NGA.SIG.0012_2.0.0_UTMUPS, 2014, §3; DMA TM 8358.1,
// 1990, §3).
//
// In UTM's range the grid zone is the UTM zone and the latitude band (32V); the 100 km square is a
// column letter and a row letter. Columns cycle through three sets of eight letters (A–H, J–R, S–Z;
// I and O omitted) with the zone number modulo three; the column of easting E is the set's letter
// ⌊E / 100 km⌋ (1 to 8). Rows cycle through twenty letters (A–V without I and O) every 2000 km of
// northing, offset by five letters in the even zones (the "AA" lettering of WGS 84). Beyond UTM's
// range, in UPS, the grid zone is A or B (south, west and east of the meridian 0°/180°) or Y or Z
// (north), with columns and rows of their own lettering from fixed false origins (the tables below,
// held by the checks to NGA's GEOTRANS through the mgrs package). The digits are the easting and
// northing within the square, truncated (a reference names the square that holds the point) to
// 1 to 5 digits each: precision 10 km to 1 m.

#pragma once

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

#include <geodesy/projection/utm.hpp>

namespace geodesy::projection {

namespace mgrs_tables {
inline constexpr const char* COLUMN_SETS[3] = {"ABCDEFGH", "JKLMNPQR", "STUVWXYZ"};
inline constexpr const char* UTM_ROWS = "ABCDEFGHJKLMNPQRSTUV";
// UPS: columns of the four polar zones from their first easting, rows of the two caps from their first northing.
inline constexpr const char* UPS_COLUMNS_A = "JKLPQRSTUXYZ";  // south, west: from 800 km
inline constexpr const char* UPS_COLUMNS_B = "ABCFGHJKLPQR";  // south, east: from 2000 km
inline constexpr const char* UPS_COLUMNS_Y = "RSTUXYZ";       // north, west: from 1300 km
inline constexpr const char* UPS_COLUMNS_Z = "ABCFGHJ";       // north, east: from 2000 km
inline constexpr const char* UPS_ROWS_SOUTH = "ABCDEFGHJKLMNPQRSTUVWXYZ";  // from 800 km
inline constexpr const char* UPS_ROWS_NORTH = "ABCDEFGHJKLMNP";            // from 1300 km
inline constexpr double TILE = 100000.0;
}  // namespace mgrs_tables

class Mgrs {
 public:
  explicit Mgrs(const UniversalGrids& grids = UniversalGrids()) : g_(grids) {}

  // The MGRS reference of a point (rad), with `digits` digits for each of easting and northing (0–5).
  std::string forward(double lat, double lon, int digits = 5) const {
    if (digits < 0 || digits > 5) throw std::invalid_argument("geodesy: MGRS precision is 0 to 5 digits");
    const UtmUps p = g_.forward(lat, lon);
    return format(p, digits);
  }

  std::string format(const UtmUps& p, int digits = 5) const {
    using namespace mgrs_tables;
    std::string out;
    char col, row;
    double e0, n0;
    if (p.zone == 0) {
      const char z = p.band;
      const char* cols = z == 'A' ? UPS_COLUMNS_A : z == 'B' ? UPS_COLUMNS_B : z == 'Y' ? UPS_COLUMNS_Y : UPS_COLUMNS_Z;
      const double cstart = z == 'A' ? 800000 : z == 'Y' ? 1300000 : 2000000;
      const char* rows = p.north ? UPS_ROWS_NORTH : UPS_ROWS_SOUTH;
      const double rstart = p.north ? 1300000 : 800000;
      const int ci = int(std::floor((p.easting - cstart) / TILE)), ri = int(std::floor((p.northing - rstart) / TILE));
      if (ci < 0 || ci >= int(std::char_traits<char>::length(cols)) || ri < 0 || ri >= int(std::char_traits<char>::length(rows)))
        throw std::invalid_argument("geodesy: the point is outside UPS's lettered squares");
      col = cols[ci];
      row = rows[ri];
      e0 = cstart + ci * TILE;
      n0 = rstart + ri * TILE;
      out += z;
    } else {
      const int set = (p.zone - 1) % 3;
      const int ci = int(std::floor(p.easting / TILE));  // 1..8
      if (ci < 1 || ci > 8) throw std::invalid_argument("geodesy: an easting outside its UTM zone's squares");
      col = COLUMN_SETS[set][ci - 1];
      const long rn = long(std::floor(p.northing / TILE));
      row = UTM_ROWS[((rn + (p.zone % 2 == 0 ? 5 : 0)) % 20 + 20) % 20];
      e0 = ci * TILE;
      n0 = rn * TILE;
      char zone[4];
      std::snprintf(zone, sizeof zone, "%02d", p.zone);
      out += zone;
      out += p.band;
    }
    out += col;
    out += row;
    if (digits > 0) {
      const double unit = std::pow(10.0, 5 - digits);
      const long de = long(std::floor((p.easting - e0) / unit + 1e-9)), dn = long(std::floor((p.northing - n0) / unit + 1e-9));
      char buf[16];
      std::snprintf(buf, sizeof buf, "%0*ld%0*ld", digits, de, digits, dn);
      out += buf;
    }
    return out;
  }

  // The point an MGRS reference names: the centre of its square at its precision unless `corner`, in
  // which case the square's south-west corner. Spaces are ignored; letters may be in either case.
  struct Parsed {
    UtmUps grid;          // zone, hemisphere, easting and northing (the centre or the corner)
    int digits = 0;
    geodesy::figure::Geodetic geodetic;
  };
  Parsed reverse(std::string text, bool corner = false) const {
    using namespace mgrs_tables;
    std::string s;
    for (char c : text) if (!std::isspace(static_cast<unsigned char>(c))) s += char(std::toupper(static_cast<unsigned char>(c)));
    std::size_t i = 0;
    int zone = 0;
    while (i < s.size() && i < 2 && std::isdigit(static_cast<unsigned char>(s[i]))) zone = zone * 10 + (s[i++] - '0');
    if (i + 3 > s.size()) throw std::invalid_argument("geodesy: not an MGRS reference: " + text);
    const char band = s[i++], col = s[i++], row = s[i++];
    const std::string digits = s.substr(i);
    if (digits.size() % 2 || digits.size() > 10 || digits.find_first_not_of("0123456789") != std::string::npos)
      throw std::invalid_argument("geodesy: an MGRS reference has an even number of digits, at most ten: " + text);
    const int n = int(digits.size() / 2);
    const double unit = std::pow(10.0, 5 - n);
    const double de = n ? std::stod(digits.substr(0, n)) * unit : 0, dn = n ? std::stod(digits.substr(n)) * unit : 0;
    const double half = corner ? 0 : unit / 2;
    Parsed out;
    out.digits = n;
    UtmUps& p = out.grid;
    auto index = [&](const char* letters, char c) {
      const char* at = std::strchr(letters, c);
      if (!at || !c) throw std::invalid_argument(std::string("geodesy: the MGRS letter ") + c + " is not used here: " + text);
      return int(at - letters);
    };
    if (zone == 0) {
      if (band != 'A' && band != 'B' && band != 'Y' && band != 'Z') throw std::invalid_argument("geodesy: a UPS zone is A, B, Y or Z: " + text);
      p.zone = 0;
      p.band = band;
      p.north = band == 'Y' || band == 'Z';
      const char* cols = band == 'A' ? UPS_COLUMNS_A : band == 'B' ? UPS_COLUMNS_B : band == 'Y' ? UPS_COLUMNS_Y : UPS_COLUMNS_Z;
      const double cstart = band == 'A' ? 800000 : band == 'Y' ? 1300000 : 2000000;
      const char* rows = p.north ? UPS_ROWS_NORTH : UPS_ROWS_SOUTH;
      const double rstart = p.north ? 1300000 : 800000;
      p.easting = cstart + index(cols, col) * TILE + de + half;
      p.northing = rstart + index(rows, row) * TILE + dn + half;
      double conv = 0, scale = 1;
      out.geodetic = g_.reverse(0, p.north, p.easting, p.northing, &conv, &scale);
      p.convergence = conv;
      p.scale = scale;
      return out;
    }
    if (zone < 1 || zone > 60) throw std::invalid_argument("geodesy: a UTM zone is 1 to 60: " + text);
    const char* bands = "CDEFGHJKLMNPQRSTUVWX";
    const int b = index(bands, band);
    p.zone = zone;
    p.band = band;
    p.north = band >= 'N';
    p.easting = (index(COLUMN_SETS[(zone - 1) % 3], col) + 1) * TILE + de + half;
    const int r = ((index(UTM_ROWS, row) - (zone % 2 == 0 ? 5 : 0)) % 20 + 20) % 20;
    // The row letters repeat every 2000 km: take the northing whose latitude falls in the band.
    const double lat_lo = (-80.0 + 8 * b) * DEG, lat_hi = (b == 19 ? 84.0 : -80.0 + 8 * (b + 1)) * DEG;
    const double base = r * TILE + dn + half;
    double best = 0, best_miss = 1e9;
    for (int k = 0; k < 10; k++) {
      const double north = base + k * 2000000.0;
      const geodesy::figure::Geodetic g = g_.reverse(zone, p.north, p.easting, north);
      const double miss = g.latitude < lat_lo ? lat_lo - g.latitude : g.latitude > lat_hi ? g.latitude - lat_hi : 0;
      if (miss < best_miss) { best_miss = miss; best = north; }
    }
    if (best_miss > 1.0 * DEG) throw std::invalid_argument("geodesy: the MGRS square does not lie in its latitude band: " + text);
    p.northing = best;
    double conv = 0, scale = 1;
    out.geodetic = g_.reverse(zone, p.north, p.easting, p.northing, &conv, &scale);
    p.convergence = conv;
    p.scale = scale;
    return out;
  }

 private:
  UniversalGrids g_;
};

}  // namespace geodesy::projection
