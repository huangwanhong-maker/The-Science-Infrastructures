// How a place is used relative to the rest of the city, and when that departs from its habit: the Urban Space Usage Index of the
// World Bank's eca-resilience (I(h, d) = U(h, d) / U(d), the day's share at each place, and Z(h, d) = (I - mu(h)) / sigma(h), the
// departure from a baseline period), written here from that statement. Dividing by the day's total removes what moves the whole
// record: there it is the data's volume, here it is that some days are busier everywhere.
//
// Places are any units with a record of counts by day (hexagons of GPS pings there; pedestrian counters here). A day on which a place has
// no record is NaN, and the day's total is over the places in the panel that do report: a panel that changes shifts the shares, so that
// `balanced` keeps only the places reporting on every day of a period.

#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace urban::usage {

using Matrix = std::vector<std::vector<double>>;   // [place][day]; NaN where there is no record

// The day's share at each place; NaN where the place has no record, or where the day holds nothing at all.
inline Matrix shares(const Matrix& counts) {
  if (counts.empty()) return {};
  const std::size_t days = counts.front().size();
  Matrix out(counts.size(), std::vector<double>(days, std::numeric_limits<double>::quiet_NaN()));
  for (std::size_t d = 0; d < days; d++) {
    double total = 0;
    for (const auto& row : counts) {
      if (row.size() != days) throw std::invalid_argument("usage: every place has the same days");
      if (std::isfinite(row[d])) total += row[d];
    }
    if (!(total > 0)) continue;
    for (std::size_t h = 0; h < counts.size(); h++) if (std::isfinite(counts[h][d])) out[h][d] = counts[h][d] / total;
  }
  return out;
}

// The places reporting on every day in [first, last).
inline std::vector<std::size_t> balanced(const Matrix& counts, std::size_t first, std::size_t last) {
  std::vector<std::size_t> out;
  for (std::size_t h = 0; h < counts.size(); h++) {
    bool all = true;
    for (std::size_t d = first; d < last && all; d++) all = std::isfinite(counts[h][d]);
    if (all) out.push_back(h);
  }
  return out;
}

struct Baseline {
  std::vector<double> mean, sd;
};

// A baseline over the days marked true: each place's mean and sample standard deviation of its index.
inline Baseline baseline_of(const Matrix& index, const std::vector<bool>& in) {
  Baseline b;
  for (const auto& row : index) {
    double n = 0, m = 0, q = 0;
    for (std::size_t d = 0; d < row.size(); d++) {
      if (!in[d] || !std::isfinite(row[d])) continue;
      n += 1;
      const double e = row[d] - m;
      m += e / n;
      q += e * (row[d] - m);
    }
    b.mean.push_back(n > 0 ? m : std::numeric_limits<double>::quiet_NaN());
    b.sd.push_back(n > 1 ? std::sqrt(q / (n - 1)) : std::numeric_limits<double>::quiet_NaN());
  }
  return b;
}

// Z(h, d); NaN where the index or the baseline is missing, or the baseline has no spread.
inline Matrix departures(const Matrix& index, const Baseline& b) {
  Matrix out = index;
  for (std::size_t h = 0; h < index.size(); h++) for (auto& v : out[h]) v = b.sd[h] > 0 ? (v - b.mean[h]) / b.sd[h] : std::numeric_limits<double>::quiet_NaN();
  return out;
}

}  // namespace urban::usage
