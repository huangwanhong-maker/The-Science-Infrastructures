// Earth orientation parameters: the observed irregularities of the Earth's rotation.
//
// The rotation of the Earth is not predictable to the precision of its measurement. The IERS
// measures and predicts it daily, in five parameters (IERS Conventions 2010, ch. 5; the
// Earth Orientation Centre's explanatory notes to finals2000A):
//   - x_p, y_p: the coordinates of the Celestial Intermediate Pole in the terrestrial frame,
//     polar motion (arcseconds; a few tenths, about 10 m on the ground);
//   - UT1 − UTC: the Earth's rotation angle as a time, the departure of UT1 from UTC
//     (seconds; kept within ±0.9 s by leap seconds, about 400 m at the equator);
//   - LOD: the excess of the length of the day over 86400 s (seconds);
//   - dX, dY: the observed offsets of the celestial pole from the IAU 2006/2000A model
//     (arcseconds; a few tenths of a milliarcsecond).
//
// The file finals2000A (Bulletin A, daily since 1973, with a year of predictions) is read by
// `parse_finals2000a`; the data library fetches it, and the domain is handed its text. Values
// are interpolated linearly in UTC between the daily rows; UT1 − UTC is interpolated as
// UT1 − TAI across a leap second, so that the jump of one second is not smeared over the day
// before it (IERS Gazette 13; the error of linear interpolation is below 30 µs in UT1 and
// 0.1 mas in the pole, against 4-point Lagrange interpolation with tidal terms).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace metrology::earth_orientation {

inline constexpr double ARCSEC = 3.14159265358979323846 / 648000.0;

struct Values {
  double mjd = 0;          // UTC, the day the row is for (0h)
  double xp = 0, yp = 0;   // polar motion, arcseconds
  double dut1 = 0;         // UT1 − UTC, seconds
  double lod = 0;          // excess length of day, seconds
  double dx = 0, dy = 0;   // celestial pole offsets from IAU 2006/2000A, arcseconds
  bool predicted = false;  // a prediction rather than a measurement (flag "P")
  bool has_nutation = false;
};

class Series {
 public:
  // The IERS finals2000A file (all, data or daily): fixed columns, one day per line. The
  // Bulletin A values are used; lines without polar motion are skipped.
  static Series parse_finals2000a(std::string_view text, std::string source = "IERS finals2000A") {
    Series s;
    s.source_ = std::move(source);
    std::istringstream in{std::string(text)};
    std::string line;
    auto field = [](const std::string& l, std::size_t a, std::size_t b, bool& ok) {
      // Columns a..b, 1-based and inclusive, as the IERS's description gives them.
      if (l.size() < b) { ok = false; return 0.0; }
      const std::string f = l.substr(a - 1, b - a + 1);
      if (f.find_first_not_of(' ') == std::string::npos) { ok = false; return 0.0; }
      char* end = nullptr;
      const double v = std::strtod(f.c_str(), &end);
      if (end == f.c_str()) ok = false;
      return v;
    };
    while (std::getline(in, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      bool ok = true;
      Values v;
      v.mjd = field(line, 8, 15, ok);
      v.xp = field(line, 19, 27, ok);
      v.yp = field(line, 38, 46, ok);
      v.dut1 = field(line, 59, 68, ok);
      if (!ok) continue;  // the empty rows at the end of the predictions
      v.predicted = line.size() > 57 && line[57] == 'P';
      bool lod_ok = true;
      v.lod = field(line, 80, 86, lod_ok) * 1e-3;
      if (!lod_ok) v.lod = 0;
      bool nut_ok = true;
      v.dx = field(line, 98, 106, nut_ok) * 1e-3;
      v.dy = field(line, 117, 125, nut_ok) * 1e-3;
      v.has_nutation = nut_ok;
      if (!nut_ok) v.dx = v.dy = 0;
      s.rows_.push_back(v);
    }
    if (s.rows_.size() < 2) throw std::invalid_argument("finals2000A: fewer than two rows of Earth orientation parameters");
    for (std::size_t i = 1; i < s.rows_.size(); i++)
      if (!(s.rows_[i].mjd > s.rows_[i - 1].mjd)) throw std::invalid_argument("finals2000A: rows out of order");
    return s;
  }
  // A series given row by row (for checks and for other sources).
  static Series from_rows(std::vector<Values> rows, std::string source) {
    Series s;
    s.rows_ = std::move(rows);
    s.source_ = std::move(source);
    if (s.rows_.size() < 2) throw std::invalid_argument("earth_orientation: fewer than two rows");
    return s;
  }

  bool covers(double mjd_utc) const { return !rows_.empty() && mjd_utc >= rows_.front().mjd && mjd_utc <= rows_.back().mjd; }
  double first_mjd() const { return rows_.front().mjd; }
  double last_mjd() const { return rows_.back().mjd; }
  // The first day whose values are predictions; last_mjd() + 1 when none are.
  double first_predicted_mjd() const {
    for (const Values& v : rows_)
      if (v.predicted) return v.mjd;
    return rows_.back().mjd + 1;
  }

  // The parameters at a UTC MJD, interpolated. Throws outside the series.
  Values at(double mjd_utc) const {
    if (!covers(mjd_utc))
      throw std::domain_error("metrology::earth_orientation: MJD " + std::to_string(mjd_utc) + " is outside the Earth orientation series (" +
                              std::to_string(first_mjd()) + " to " + std::to_string(last_mjd()) + ")");
    auto it = std::upper_bound(rows_.begin(), rows_.end(), mjd_utc, [](double m, const Values& v) { return m < v.mjd; });
    if (it == rows_.end()) return rows_.back();
    const Values& b = *it;
    const Values& a = *(it - 1);
    const double w = (mjd_utc - a.mjd) / (b.mjd - a.mjd);
    auto lerp = [w](double p, double q) { return p + w * (q - p); };
    Values v;
    v.mjd = mjd_utc;
    v.xp = lerp(a.xp, b.xp);
    v.yp = lerp(a.yp, b.yp);
    // A leap second between the rows: interpolate as UT1 − TAI (remove the jump), then
    // express the result in the earlier row's UTC.
    const double jump = std::round(b.dut1 - a.dut1);
    v.dut1 = lerp(a.dut1, b.dut1 - jump);
    v.lod = lerp(a.lod, b.lod);
    v.dx = lerp(a.dx, b.dx);
    v.dy = lerp(a.dy, b.dy);
    v.predicted = a.predicted || b.predicted;
    v.has_nutation = a.has_nutation && b.has_nutation;
    return v;
  }
  const std::vector<Values>& rows() const { return rows_; }
  const std::string& source() const { return source_; }

 private:
  std::vector<Values> rows_;
  std::string source_;
};

}  // namespace metrology::earth_orientation
