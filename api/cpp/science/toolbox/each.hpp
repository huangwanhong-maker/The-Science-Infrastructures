// A tool run for each place: the rows of the table split by the columns "each" names (a place, or a place and a time), the
// tool run on each part, and the numbers of its answers read back as measures by key, which the page can draw. Since the
// tools answer in blocks that say what they are, every tool on a table gains this at once, and none need know of it.
//
//   {"kind": "each", "title", "columns": [...], "keys": [[...]], "measures": [{"name", "values": [...]}], "left": [{"reason", "count"}]}
//
// A part a tool cannot answer (too few rows, a value that does not vary) is left out and counted by its reason, and is not
// given a zero. A part is run on the same tool with the same parameters; its standard errors are its own, and nothing
// here corrects for testing many parts at once.

#pragma once

#include <functional>
#include <unordered_map>

#include "answer.hpp"

namespace science::toolbox {

// The numbers in an answer, by name: figures by label, and from short tables each row's estimate and p (or every cell).
inline void harvest(const Json& answer_json, std::vector<std::pair<std::string, double>>& out) {
  const Json* blocks = answer_json.find("blocks");
  if (!blocks || !blocks->is_array()) return;
  for (const auto& b : blocks->as_array()) {
    const Json* kind = b.find("kind");
    if (!kind || !kind->is_string()) continue;
    if (kind->as_string() == "figures") {
      for (const auto& item : b["items"].as_array()) {
        const auto& pair = item.as_array();
        if (pair.size() >= 2 && pair[0].is_string() && pair[1].is_number()) out.emplace_back(pair[0].as_string(), pair[1].as_number());
      }
    } else if (kind->as_string() == "table") {
      const auto& columns = b["columns"].as_array();
      const auto& rows = b["rows"].as_array();
      if (rows.size() > 12 || columns.size() < 2) continue;
      std::vector<std::size_t> take;
      for (std::size_t c = 1; c < columns.size(); c++)
        if (columns[c].is_string() && (columns[c].as_string() == "estimate" || columns[c].as_string() == "p")) take.push_back(c);
      if (take.empty())
        for (std::size_t c = 1; c < columns.size(); c++) take.push_back(c);
      for (const auto& row : rows) {
        const auto& cells = row.as_array();
        if (cells.empty() || !cells[0].is_string()) continue;
        for (std::size_t c : take)
          if (c < cells.size() && cells[c].is_number() && columns[c].is_string()) out.emplace_back(cells[0].as_string() + " · " + columns[c].as_string(), cells[c].as_number());
      }
    }
  }
}

inline Json each(const Json& request, const std::function<Json(const Request&)>& tool) {
  const Request R(request);
  std::vector<std::string> names;
  for (const auto& e : request["each"].as_array())
    if (e.is_string()) names.push_back(e.as_string());
  if (names.empty()) throw std::invalid_argument("say what to run the tool for each of");
  std::vector<std::vector<std::string>> key_columns;
  for (const auto& n : names) key_columns.push_back(R.texts(n));
  const std::size_t n = key_columns[0].size();

  // The rows of each part, in the order the parts first appear.
  std::vector<std::vector<std::string>> keys;
  std::vector<std::vector<std::size_t>> parts;
  std::unordered_map<std::string, std::size_t> index;
  for (std::size_t i = 0; i < n; i++) {
    std::string joined;
    for (const auto& c : key_columns) joined += c[i] + '\x1f';
    const auto [it, fresh] = index.try_emplace(joined, keys.size());
    if (fresh) {
      std::vector<std::string> key;
      for (const auto& c : key_columns) key.push_back(c[i]);
      keys.push_back(std::move(key));
      parts.emplace_back();
    }
    parts[it->second].push_back(i);
  }
  if (keys.size() > 6000) throw std::invalid_argument("too many parts (" + std::to_string(keys.size()) + ", the limit is 6000): choose fewer rows");

  std::vector<std::string> order;
  std::map<std::string, std::vector<double>> measures;
  std::map<std::string, double> left;
  const Json& all = request["columns"];
  for (std::size_t p = 0; p < parts.size(); p++) {
    Obj columns;
    for (const auto& [name, column] : all.as_object()) {
      Arr cut;
      const auto& cells = column.as_array();
      cut.reserve(parts[p].size());
      for (std::size_t i : parts[p]) cut.push_back(cells[i]);
      columns.push_back({name, Json(std::move(cut))});
    }
    Obj sub;
    for (const auto& [name, value] : request.as_object())
      if (name != "each" && name != "columns" && name != "by_row" && name != "places") sub.push_back({name, value});
    sub.push_back({"columns", Json(std::move(columns))});
    const Json part(std::move(sub));
    std::vector<std::pair<std::string, double>> numbers;
    try {
      harvest(tool(Request(part)), numbers);
    } catch (const std::exception& error) {
      left[error.what()] += 1;
      continue;
    }
    for (const auto& [name, value] : numbers) {
      const auto [it, fresh] = measures.try_emplace(name, std::vector<double>(keys.size(), NaN));
      if (fresh) order.push_back(name);
      it->second[p] = value;
    }
  }
  if (measures.empty()) {
    std::string why = left.empty() ? "the tool gave no numbers" : left.begin()->first;
    throw std::invalid_argument("no part could be answered: " + why);
  }

  Arr jkeys;
  for (const auto& k : keys) jkeys.push_back(strs(k));
  Arr jmeasures;
  for (const auto& name : order) jmeasures.push_back(Json(Obj{{"name", Json(name)}, {"values", nums(measures[name])}}));
  Arr jleft;
  for (const auto& [reason, count] : left) jleft.push_back(Json(Obj{{"reason", Json(reason)}, {"count", Json(count)}}));
  std::string joined;
  for (std::size_t k = 0; k < names.size(); k++) joined += (k ? " and " : "") + names[k];
  Arr blocks;
  blocks.push_back(Json(Obj{{"kind", Json("each")}, {"title", Json("For each " + joined)}, {"columns", strs(names)}, {"keys", Json(std::move(jkeys))}, {"measures", Json(std::move(jmeasures))}, {"left", Json(std::move(jleft))}}));
  return answer("For each " + joined, std::move(blocks));
}

}  // namespace science::toolbox
