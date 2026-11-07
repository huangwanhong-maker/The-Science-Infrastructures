// Events counted on hexagons, through time: each event of the record shown placed in its H3 cell (physis's H3) and in its step of time, and the counts given back as cells that the
// workbench raises and lowers on the map. An example of a program in C++: it is compiled here on first use, against the foundations and the domains, so that it may use any of them.
//
// Limits: a count per step depends on the step and the resolution chosen; a cell's area varies by a few percent over the Earth, and these are counts, not densities.

#include "spacetime_plugin.hpp"

#include <physis/computational_geometry/h3.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <stdexcept>

namespace st = science::spacetime;
namespace H = physis::computational_geometry::h3;

struct HexBins : st::Analysis {
  science::Json describe() const override {
    return science::Json::object({{"id", "hex-bins"}, {"title", "Events on hexagons, through time"},
                                  {"description", "Each event in its H3 cell and its step of time: the counts as columns rising and falling."},
                                  {"params", science::Json::array({st::param("resolution", "H3 resolution", 3, "integer", 0, 9), st::param("steps", "Steps of time", 60, "integer", 2, 400)})}});
  }

  std::vector<science::Json> run(const std::vector<st::Record>& records, const st::Params& p, const st::Context&, st::Progress& progress) override {
    const auto found = std::find_if(records.begin(), records.end(), [](const st::Record& r) { return r.kind == "events"; });
    if (found == records.end()) throw std::invalid_argument("This needs events: show a record of events first.");
    const st::Events& e = found->events;
    const int res = p.integer("resolution", 3), steps = p.integer("steps", 60);
    constexpr double DEG = std::numbers::pi / 180;
    double t0 = 1e300, t1 = -1e300;
    for (const double t : e.time) if (std::isfinite(t)) { t0 = std::min(t0, t); t1 = std::max(t1, t); }
    if (t0 > t1) throw std::invalid_argument("These events have no times.");
    const double width = std::max(1.0, (t1 - t0) / steps);
    std::map<std::uint64_t, std::vector<double>> counts;
    for (std::size_t k = 0; k < e.size(); k++) {
      if (!std::isfinite(e.time[k]) || !std::isfinite(e.lat[k]) || !std::isfinite(e.lon[k])) continue;
      auto& row = counts[H::cell_at({e.lat[k] * DEG, e.lon[k] * DEG}, res)];
      if (row.empty()) row.assign(static_cast<std::size_t>(steps), 0.0);
      row[static_cast<std::size_t>(std::min(steps - 1, static_cast<int>((e.time[k] - t0) / width)))] += 1;
      if (k % 20000 == 0) progress(static_cast<double>(k) / static_cast<double>(e.size()), "binning");
    }
    std::vector<st::Cell> cells;
    std::vector<std::vector<double>> frames(static_cast<std::size_t>(steps));
    for (auto& f : frames) f.reserve(counts.size());
    for (const auto& [cell, row] : counts) {
      st::Cell c;
      c.id = H::to_string(cell);
      double previous = 1e9;
      for (const auto& q : H::boundary(cell)) {
        double lon = q.lng / DEG;
        if (previous < 1e8) { while (lon - previous > 180) lon -= 360; while (lon - previous < -180) lon += 360; }
        previous = lon;
        c.ring.push_back({lon, q.lat / DEG});
      }
      cells.push_back(std::move(c));
      for (std::size_t s = 0; s < row.size(); s++) frames[s].push_back(row[s]);
    }
    std::vector<double> time;
    for (int s = 0; s < steps; s++) time.push_back(t0 + (s + 0.5) * width);
    return {st::cells_result("On hexagons: " + found->title, time, cells, {{"events", frames}}, science::Json::object({{"map", "plum"}, {"scale", "sqrt"}}),
                             "Counts on H3 cells of resolution " + std::to_string(res) + ", each step " + std::to_string(static_cast<long>(width / st::day_ms)) + " days."),
            st::note(std::to_string(cells.size()) + " cells of resolution " + std::to_string(res) + " hold these events.")};
  }
};

SCIENCE_SPACETIME_PLUGIN(HexBins)
