// The toolbox of quantitative economics (ADR 0056, moved here from the Economy Atlas by ADR 0080, so that any application's
// engine may answer it): each tool takes the columns of a table, by the roles the page gave them, and parameters, and
// answers blocks the page draws:
//
//   {"kind": "figures", "title", "items": [[label, value, note?]]}       key numbers (or words)
//   {"kind": "table", "title", "columns": [...], "rows": [[...]], "note"?}
//   {"kind": "lines", "title", "x": [...], "xlabel", "ylabel", "series": [{"name", "y", "lower"?, "upper"?, "dashed"?}]}
//   {"kind": "bars", "title", "labels": [...], "series": [{"name", "values"}], "band"?: [lo, hi], "stacked"?}
//   {"kind": "scatter", "title", "x", "y", "xlabel", "ylabel", "labels"?, "line"?: {"x", "y"}}
//   {"kind": "matrix", "title", "rows": [...], "columns": [...], "values": [[...]]}
//   {"kind": "map", "title", "values": {id: number}, "scale": "diverging"|"sequential", "unit"}   to the map, as a layer
//   {"kind": "rows", "title", "rows": [i], "values": [...], "scale", "unit"}   a value for each row kept, drawn where the rows are places (asked by "by_row", or a tool's own)
//   {"kind": "each", "title", "columns", "keys", "measures", "left"}           the tool's numbers for each part ("each": [column, ...]; toolbox/each.hpp)
//   {"kind": "note", "text"}
//
// The methods are the foundations' and the domain's: tyche (probability, statistics, regression,
// time series), nomos (econometrics, measurement, index numbers, games, matching, choice), and the
// economics domain (growth, macro, trade); the engine only reads, calls and writes. The tools are
// in toolbox/, a header for each group.

#pragma once

#include <functional>
#include <map>
#include <string>

#include "toolbox/answer.hpp"
#include "toolbox/each.hpp"
#include "toolbox/econometrics.hpp"
#include "toolbox/games.hpp"
#include "toolbox/growth.hpp"
#include "toolbox/measurement.hpp"
#include "toolbox/series.hpp"
#include "toolbox/spatial.hpp"
#include "toolbox/trade.hpp"

namespace science::toolbox {

inline const std::map<std::string, std::function<Json(const Request&)>>& registry() {
  static const std::map<std::string, std::function<Json(const Request&)>> R = {
      {"summary", summary},
      {"ols", ols},
      {"glm", glm},
      {"panel", panel},
      {"iv", iv},
      {"did", did},
      {"convergence", convergence},
      {"series", autocorrelation},
      {"unit-root", unit_root},
      {"arima", arima_tool},
      {"hp-filter", hp_filter},
      {"var", var_tool},
      {"cointegration", cointegration},
      {"solow", solow},
      {"islm", islm},
      {"growth-accounting", growth_accounting},
      {"inequality", inequality},
      {"poverty", poverty},
      {"concentration", concentration},
      {"price-index", price_index},
      {"specialisation", specialisation},
      {"moran", moran},
      {"game", game},
      {"zero-sum", zero_sum},
      {"matching", matching},
      {"voting", voting},
      {"trade-concentration", trade_concentration},
      {"trade-intensity", trade_intensity},
      {"trade-gravity", trade_gravity},
  };
  return R;
}

inline Json run(const Json& request) {
  const Json* t = request.find("tool");
  if (!t || !t->is_string()) throw std::invalid_argument("name the tool");
  const auto& R = registry();
  const auto it = R.find(t->as_string());
  if (it == R.end()) throw std::invalid_argument("the engine has no tool \"" + t->as_string() + "\" yet");
  // For each place (toolbox/each.hpp), when a table's columns name what to split by.
  if (const Json* e = request.find("each"); e && e->is_array() && !e->as_array().empty()) return each(request, it->second);
  return it->second(Request(request));
}

}  // namespace science::toolbox
