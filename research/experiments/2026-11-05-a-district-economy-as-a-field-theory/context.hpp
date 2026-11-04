// The inputs of the study and its calibration of the people's field (README §3): the shared problem, the counters on
// and around the district, and universal kriging line by line with its parameters chosen by leave-one-site-out.

#pragma once

#include <fstream>
#include <memory>
#include <sstream>

#include "model.hpp"

namespace field {

using economy::Json;

inline std::string read_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream b;
  b << in.rdbuf();
  return b.str();
}
inline void write_file(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  f << text;
  if (!f) throw std::runtime_error("cannot write " + path);
}

struct Context {
  Json inputs, problem_json, sensors, counters;
  std::unique_ptr<Problem> problem;
  std::vector<Counter> inside, outside;  // the counters on the district's network, and those around it
};

inline Context load(const std::string& inputs, const std::string& problem, const std::string& sensors, const std::string& counters) {
  Context C;
  C.inputs = Json::parse(read_file(inputs));
  C.problem_json = Json::parse(read_file(problem));
  C.sensors = Json::parse(read_file(sensors));
  if (!counters.empty()) C.counters = Json::parse(read_file(counters));
  C.problem = std::make_unique<Problem>(C.inputs, &C.problem_json);
  const auto& box = C.inputs["district"]["district_box"];
  std::size_t row = 0;
  for (const auto& s : C.sensors["sensors"].as_array()) {
    Counter c;
    c.row = row++;
    c.id = static_cast<long long>(s["id"].as_number());
    c.name = s["name"].text();
    c.at = {s["x"].as_number(), s["y"].as_number()};
    for (const auto& v : s["mean_week"].as_array()) c.week.push_back(v.as_number());
    const bool in = c.at[0] >= box["x0"].as_number() && c.at[0] <= box["x1"].as_number() && c.at[1] >= box["y0"].as_number() && c.at[1] <= box["y1"].as_number();
    (in ? C.inside : C.outside).push_back(std::move(c));
  }
  return C;
}

// A counter's foot on the base network: its nearest segment's centreline; a counter of a crossing stands on the way it
// crosses to, the nearest segment of another name.
inline Point counter_foot(const Graph& G, const Counter& c) {
  auto o = observe(G, c.at);
  if (c.name.find("Crossing") != std::string::npos) {
    const std::string first = G.name[o.edge];
    o = observe(G, c.at, [&](std::size_t e) { return !G.link[e] && G.name[e] != first; });
  }
  return o.on;
}

// ---- calibration --------------------------------------------------------------------------------------------------------

struct Skill {
  double nsse = 0;          // the mean over counters of Σ_h (ŷ − y)² / Σ_h y²
  double median_r2 = 0;     // the median over counters of 1 − Σ(ŷ − y)²/Σ(y − ȳ)²
  double log_total_rms = 0; // the RMS of ln(predicted / observed weekly total)
  double total_corr = 0;    // the correlation of ln weekly totals, predicted and observed
  std::vector<double> r2, ratio;
  std::vector<std::vector<double>> predicted;  // per counter, its 168 hours
};

inline Skill skill_of(const std::vector<Counter>& cs, const std::vector<std::vector<double>>& pred) {
  Skill S;
  S.predicted = pred;
  double lt = 0, lp = 0, n = 0;
  std::vector<double> a, b;
  for (std::size_t c = 0; c < cs.size(); c++) {
    double se = 0, yy = 0, mean = 0, sp = 0, sy = 0;
    for (double y : cs[c].week) mean += y / week;
    for (std::size_t h = 0; h < week; h++) {
      const double y = cs[c].week[h], p = pred[c][h];
      se += (p - y) * (p - y), yy += y * y, sp += p, sy += y;
    }
    double vv = 0;
    for (double y : cs[c].week) vv += (y - mean) * (y - mean);
    S.nsse += se / yy / static_cast<double>(cs.size());
    S.r2.push_back(1 - se / vv);
    S.ratio.push_back(sp / sy);
    const double r = std::log(std::max(1e-9, sp) / sy);
    S.log_total_rms += r * r;
    a.push_back(std::log(std::max(1e-9, sp))), b.push_back(std::log(sy));
    lt += b.back(), lp += a.back(), n++;
  }
  S.log_total_rms = std::sqrt(S.log_total_rms / n);
  auto r2 = S.r2;
  std::sort(r2.begin(), r2.end());
  S.median_r2 = r2.size() % 2 ? r2[r2.size() / 2] : 0.5 * (r2[r2.size() / 2 - 1] + r2[r2.size() / 2]);
  double sab = 0, saa = 0, sbb = 0;
  for (std::size_t i = 0; i < a.size(); i++) sab += (a[i] - lp / n) * (b[i] - lt / n), saa += (a[i] - lp / n) * (a[i] - lp / n), sbb += (b[i] - lt / n) * (b[i] - lt / n);
  S.total_corr = sab / std::sqrt(saa * sbb);
  return S;
}

// Sites: counters within 30 m of one another stand at one site, and are left out together.
inline std::vector<std::vector<std::size_t>> sites_of(const std::vector<Counter>& cs) {
  std::vector<int> site(cs.size(), -1);
  std::vector<std::vector<std::size_t>> out;
  for (std::size_t i = 0; i < cs.size(); i++) {
    if (site[i] >= 0) continue;
    site[i] = static_cast<int>(out.size());
    out.push_back({i});
    for (std::size_t j = i + 1; j < cs.size(); j++)
      if (site[j] < 0 && std::hypot(cs[i].at[0] - cs[j].at[0], cs[i].at[1] - cs[j].at[1]) < 30) site[j] = site[i], out.back().push_back(j);
  }
  return out;
}

struct Fitted {
  Physics phys;
  bool edge_from_counters = true, trams = false;
  Skill loo, regression_loo, in_sample;
  std::vector<LineFit> fits;  // all counters, line by line
};

// Leave-one-site-out predictions (kriged, and the groups' regression alone) of a kriging on a set of counters.
inline std::pair<Skill, Skill> leave_out(const Kriging& K, const std::vector<Counter>& cs) {
  const auto sites = sites_of(cs);
  std::vector<std::vector<Complex>> pk(cs.size(), std::vector<Complex>(lines)), pr(cs.size(), std::vector<Complex>(lines));
  for (const auto& site : sites) {
    std::vector<std::size_t> set;
    for (std::size_t i = 0; i < cs.size(); i++)
      if (std::find(site.begin(), site.end(), i) == site.end()) set.push_back(i);
    for (std::size_t m = 0; m < lines; m++) {
      const auto F = K.fit(m, set);
      for (auto q : site) {
        pk[q][m] = K.predict(m, F, q);
        Complex s{0, 0};
        for (std::size_t g = 0; g < K.groups; g++) s += K.X[m][q][g] * F.beta[g];
        pr[q][m] = s;
      }
    }
  }
  const auto& W = the_week();
  std::vector<std::vector<double>> a, b;
  for (std::size_t q = 0; q < cs.size(); q++) a.push_back(W.inverse(pk[q])), b.push_back(W.inverse(pr[q]));
  return {skill_of(cs, a), skill_of(cs, b)};
}

}  // namespace field
