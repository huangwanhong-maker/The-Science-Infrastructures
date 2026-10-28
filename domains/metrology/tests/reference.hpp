// The reference values of the metrology domain's C++ checks: lines of tests/fixtures/
// erfa-reference.txt (ERFA through pyerfa, make_fixtures.py), read by name.

#pragma once

#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace metrology::checks {

struct Row {
  std::vector<std::string> fields;  // after the name
  double number(std::size_t i) const { return std::strtod(fields.at(i).c_str(), nullptr); }
  int integer(std::size_t i) const { return std::atoi(fields.at(i).c_str()); }
};

// Rows by name; empty when the file is missing.
inline std::map<std::string, std::vector<Row>> read_reference(const std::string& path) {
  std::map<std::string, std::vector<Row>> out;
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream f(line);
    std::string name, x;
    f >> name;
    Row r;
    while (f >> x) r.fields.push_back(x);
    out[name].push_back(r);
  }
  return out;
}

}  // namespace metrology::checks
