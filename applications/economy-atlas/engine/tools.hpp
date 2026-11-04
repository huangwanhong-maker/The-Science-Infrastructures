// The toolbox of the Economy Atlas's engine (ADR 0048): each tool takes the columns of a table, by
// the roles the page gave them, and parameters, and answers blocks the page draws:
//
//   {"kind": "figures", "title", "items": [[label, value, note?]]}       key numbers (or words)
//   {"kind": "table", "title", "columns": [...], "rows": [[...]], "note"?}
//   {"kind": "lines", "title", "x": [...], "xlabel", "ylabel", "series": [{"name", "y", "lower"?, "upper"?, "dashed"?}]}
//   {"kind": "bars", "title", "labels": [...], "series": [{"name", "values"}], "band"?: [lo, hi], "stacked"?}
//   {"kind": "scatter", "title", "x", "y", "xlabel", "ylabel", "labels"?, "line"?: {"x", "y"}}
//   {"kind": "matrix", "title", "rows": [...], "columns": [...], "values": [[...]]}
//   {"kind": "map", "title", "values": {id: number}, "scale": "diverging"|"sequential", "unit"}   to the map, as a layer
//   {"kind": "note", "text"}
//
// The methods are the foundations' and the domain's: tyche (probability, statistics, regression,
// time series), nomos (econometrics, measurement, index numbers, games, matching, choice), and the
// economics domain (growth, macro, trade); the engine only reads, calls and writes. The tools are
// in tools/, a header for each group.

#pragma once

#include <functional>
#include <map>
#include <string>

#include "tools/answer.hpp"
#include "tools/econometrics.hpp"
#include "tools/games.hpp"
#include "tools/growth.hpp"
#include "tools/measurement.hpp"
#include "tools/series.hpp"
#include "tools/trade.hpp"

namespace tools {

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
  return it->second(Request(request));
}

}  // namespace tools
