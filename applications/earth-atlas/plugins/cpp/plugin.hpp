// Plug-ins of the Earth Atlas in C++ (ADR 0029): the same contract as the
// JavaScript and Python classes, for a program the server compiles and runs.
// A plug-in derives from Diagnostic, fills in its identity, needs, parameters
// and outputs, writes compute(), and ends with PLUGIN_MAIN(Type):
//
//     struct Mine : plugin::Diagnostic { Mine() { id = "mine"; ... } std::map<std::string, std::vector<float>> compute(...) override; };
//     PLUGIN_MAIN(Mine)
//
// The host (applications/earth-atlas/py/host.py) runs the program twice: with
// --describe, which prints the plug-in's description as JSON, and with
// --run <dir>, where dir holds input.txt and one float32 file per input
// variable, and receives output.txt and one float32 file per output. Fields are
// row-major on a grid whose rows run north to south, NaN where there is no
// value, one block per time; the outputs lie on the first need's grid and times.
// Only the standard library is used, so any C++20 compiler builds a plug-in.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace plugin {

struct Grid { double lon0 = 0, dlon = 1; int nlon = 0; double lat0 = 90, dlat = 1; int nlat = 0; bool periodic = true;
              int size() const { return nlon * nlat; } };
struct Need { std::string name, label; std::vector<std::string> standard_names; std::string level; int times = 0; };
struct Output { std::string name, long_name, standard_name, units, level; bool difference = false; };
struct Parameter { std::string name, label; double value = 0, min = 0, max = 0, step = 0; std::string unit; };

// One resolved need: the product chosen, its grid and times, and its variables as float32 blocks.
struct Input {
  std::string product, model, level;
  Grid grid;
  std::vector<std::string> times, names;
  std::vector<std::vector<float>> data;      // per variable: times × nlat × nlon
  const float* slice(int variable, int index) const {
    int last = std::max(0, int(times.size()) - 1);
    index = std::max(0, std::min(last, index));
    return data[variable].data() + size_t(index) * grid.size();
  }
};
using Inputs = std::map<std::string, Input>;
using Params = std::map<std::string, double>;

struct Plugin {
  std::string id = "plugin", title = "A plug-in", description, domain = "meteorology";
  std::vector<Need> needs;
  std::vector<Parameter> parameters;
  virtual ~Plugin() = default;
};

// A field computed from fields, at every time of the first need, on its grid.
struct Diagnostic : Plugin {
  std::vector<Output> outputs;
  std::map<std::string, std::string> show;   // what to select: {"surface", name}, {"flow", vector}, {"lines", name}
  virtual std::map<std::string, std::vector<float>> compute(const Inputs& inputs, const Params& params, int index) = 0;
};

// ---- JSON out, for --describe ---------------------------------------------------------

inline std::string quoted(const std::string& s) {
  std::string out = "\"";
  for (unsigned char c : s) {
    if (c == '"') out += "\\\""; else if (c == '\\') out += "\\\\"; else if (c == '\n') out += "\\n";
    else if (c < 0x20) { char buf[8]; std::snprintf(buf, sizeof buf, "\\u%04x", c); out += buf; } else out += char(c);
  }
  return out + "\"";
}

inline std::string describe(const Diagnostic& p) {
  std::ostringstream o;
  o << "{\"id\":" << quoted(p.id) << ",\"title\":" << quoted(p.title) << ",\"description\":" << quoted(p.description)
    << ",\"domain\":" << quoted(p.domain) << ",\"parameters\":{";
  for (size_t i = 0; i < p.parameters.size(); i++) {
    const auto& q = p.parameters[i];
    o << (i ? "," : "") << quoted(q.name) << ":{\"label\":" << quoted(q.label) << ",\"default\":" << q.value
      << ",\"min\":" << q.min << ",\"max\":" << q.max << ",\"step\":" << q.step << ",\"unit\":" << quoted(q.unit) << "}";
  }
  o << "},\"needs\":[";
  for (size_t i = 0; i < p.needs.size(); i++) {
    const auto& n = p.needs[i];
    o << (i ? "," : "") << "{\"name\":" << quoted(n.name) << ",\"label\":" << quoted(n.label) << ",\"standard_names\":[";
    for (size_t k = 0; k < n.standard_names.size(); k++) o << (k ? "," : "") << quoted(n.standard_names[k]);
    o << "]";
    if (!n.level.empty()) o << ",\"level\":" << quoted(n.level);
    if (n.times) o << ",\"times\":" << n.times;
    o << "}";
  }
  o << "],\"outputs\":[";
  for (size_t i = 0; i < p.outputs.size(); i++) {
    const auto& v = p.outputs[i];
    o << (i ? "," : "") << "{\"name\":" << quoted(v.name) << ",\"long_name\":" << quoted(v.long_name)
      << ",\"standard_name\":" << quoted(v.standard_name) << ",\"units\":" << quoted(v.units) << "}";
  }
  o << "],\"show\":{";
  size_t k = 0;
  for (const auto& [key, value] : p.show) o << (k++ ? "," : "") << quoted(key) << ":" << quoted(value);
  o << "}}\n";
  return o.str();
}

// ---- the files of a run -----------------------------------------------------------------

inline std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> out;
  std::string part;
  std::istringstream in(line);
  while (std::getline(in, part, '\t')) out.push_back(part);
  return out;
}

inline std::vector<float> read_f32(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  in.seekg(0, std::ios::end);
  std::vector<float> out(size_t(in.tellg()) / sizeof(float));
  in.seekg(0);
  in.read(reinterpret_cast<char*>(out.data()), std::streamsize(out.size() * sizeof(float)));
  return out;
}

inline void write_f32(const std::string& path, const std::vector<float>& v) {
  std::ofstream out(path, std::ios::binary);
  out.write(reinterpret_cast<const char*>(v.data()), std::streamsize(v.size() * sizeof(float)));
}

struct Run { Inputs inputs; Params params; long long time_ms = 0; bool pinned = false; double pin_lon = 0, pin_lat = 0; std::string first; };

inline Run load(const std::string& dir) {
  Run run;
  std::ifstream in(dir + "/input.txt");
  if (!in) throw std::runtime_error("no input.txt in " + dir);
  std::string line;
  while (std::getline(in, line)) {
    auto f = split(line);
    if (f.empty()) continue;
    if (f[0] == "time" && f.size() > 1) run.time_ms = std::stoll(f[1]);
    else if (f[0] == "pinned" && f.size() > 2) { run.pinned = true; run.pin_lon = std::stod(f[1]); run.pin_lat = std::stod(f[2]); }
    else if (f[0] == "param" && f.size() > 2) { try { run.params[f[1]] = std::stod(f[2]); } catch (...) {} }
    else if (f[0] == "need" && f.size() > 12) {
      Input& inp = run.inputs[f[1]];
      if (run.first.empty()) run.first = f[1];
      inp.product = f[2]; inp.model = f[3]; inp.level = f[4];
      inp.grid = { std::stod(f[6]), std::stod(f[7]), std::stoi(f[8]), std::stod(f[9]), std::stod(f[10]), std::stoi(f[11]), f[12] == "1" };
    } else if (f[0] == "times" && f.size() > 1) {
      run.inputs[f[1]].times.assign(f.begin() + 2, f.end());
    } else if (f[0] == "data" && f.size() > 3) {
      Input& inp = run.inputs[f[1]];
      inp.names.push_back(f[2]);
      inp.data.push_back(read_f32(dir + "/" + f[3]));
    }
  }
  return run;
}

template <class T>
int run_main(int argc, char** argv) {
  T p;
  try {
    if (argc >= 2 && std::string(argv[1]) == "--describe") { std::cout << describe(p); return 0; }
    if (argc >= 3 && std::string(argv[1]) == "--run") {
      std::string dir = argv[2];
      Run run = load(dir);
      const Input& first = run.inputs.at(run.first);
      int ntimes = std::max<int>(1, int(first.times.size()));
      std::map<std::string, std::vector<float>> buffers;
      for (const auto& o : p.outputs) buffers[o.name].reserve(size_t(ntimes) * first.grid.size());
      for (int k = 0; k < ntimes; k++) {
        auto values = p.compute(run.inputs, run.params, k);
        for (const auto& o : p.outputs) {
          auto it = values.find(o.name);
          if (it == values.end() || int(it->second.size()) != first.grid.size()) throw std::runtime_error("compute() did not return " + o.name + " on the grid");
          buffers[o.name].insert(buffers[o.name].end(), it->second.begin(), it->second.end());
        }
      }
      std::ofstream out(dir + "/output.txt");
      for (const auto& o : p.outputs) {
        std::string file = "out." + o.name + ".f32";
        write_f32(dir + "/" + file, buffers[o.name]);
        out << "field\t" << o.name << '\t' << o.long_name << '\t' << o.standard_name << '\t' << o.units << '\t'
            << (o.level.empty() ? first.level : o.level) << '\t' << (o.difference ? 1 : 0) << '\t' << file << '\n';
      }
      for (const auto& [key, value] : p.show) out << "show\t" << key << '\t' << value << '\n';
      return 0;
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
  std::cerr << "usage: " << argv[0] << " --describe | --run <dir>\n";
  return 2;
}

}  // namespace plugin

#define PLUGIN_MAIN(Type) int main(int argc, char** argv) { return plugin::run_main<Type>(argc, argv); }
