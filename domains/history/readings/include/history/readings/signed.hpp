// The signed ties of a history read as a ℤ₂ gauge field in temporal gauge (theory note §3).
//
// Each worldline is a path, so the worldlines form a forest and every configuration is gauge
// equivalent to one whose temporal links are the identity; the ties' signs at their meetings
// are then the spatial links. In this gauge:
//
//   the holonomy of a lens is the change of its tie between two meetings, s(q)·s(p): a flip
//   where it is −1;
//   a triangle of ties at an occasion has the product of their signs: frustrated where −1;
//   on a closed shell (three lenses between the same two meetings of three mutually tied
//   beings) the Bianchi identity holds: product(q) = product(p)·(−1)^flips.
//
// Each reading is a class function, so it is unchanged when every sign of one being's ties
// is reversed (switching, Zaslavsky 1982): the gauge invariance the checks verify.
//
// Persistence. For ℤ₂ with the weight e^{β s s'} on a lens, the correlation of a tie's signs
// across one lens is tanh β (the character expansion: e^{βs} = cosh β (1 + s tanh β)). Its
// estimate is the memory m = 1 − 2p, p the fraction of lenses with both ends attested that
// flip; a sign carried k lenses from its last attestation is expected to be s·m^k. Where
// flips outnumber persistence (m ≤ 0) β is not defined, and the memory used is 0.

#pragma once

#include <history/construction/complex.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace history {

struct Persistence {
  std::string colour;             // or "all"
  std::size_t lenses = 0, flips = 0;
  double p = 0, memory = 1;       // p = flips / lenses; memory = 1 − 2p = tanh β
  std::optional<double> beta;     // artanh(memory), when 0 < memory < 1
};

inline Persistence persistence_of(std::string colour, std::size_t lenses, std::size_t flips) {
  Persistence P;
  P.colour = std::move(colour);
  P.lenses = lenses;
  P.flips = flips;
  if (lenses) {
    P.p = static_cast<double>(flips) / static_cast<double>(lenses);
    P.memory = 1 - 2 * P.p;
    if (P.memory > 0 && P.memory < 1) P.beta = std::atanh(P.memory);
  }
  return P;
}

// The persistence of each colour that has lenses, in the codebook's order, then of all.
inline std::vector<Persistence> persistence(const HistoryComplex& K) {
  std::map<std::string, std::array<std::size_t, 2>> count;
  std::size_t n = 0, f = 0;
  for (std::size_t i = 0; i < K.lenses.size(); i++) {
    const Lens& L = K.lenses[i];
    const auto& g = K.rungs[L.tie];
    if (!g[L.rung_from].attested || !g[L.rung_to].attested) continue;
    auto& c = count[K.record.ties[L.tie].colour];
    c[0]++;
    c[1] += K.flip(i);
    n++;
    f += K.flip(i);
  }
  std::vector<Persistence> out;
  for (const auto& colour : colours())
    if (count.count(colour)) out.push_back(persistence_of(colour, count[colour][0], count[colour][1]));
  out.push_back(persistence_of("all", n, f));
  return out;
}

// The memory of each colour: its own where it has attested lenses, else that of all.
struct Memory {
  std::map<std::string, double> by_colour;
  double all = 1;
  double of(const std::string& colour) const {
    const auto it = by_colour.find(colour);
    return std::clamp(it == by_colour.end() ? all : it->second, 0.0, 1.0);
  }
};

inline Memory memory(const std::vector<Persistence>& P) {
  Memory m;
  for (const auto& p : P) {
    if (p.colour == "all") m.all = p.lenses ? p.memory : 1;
    else if (p.lenses) m.by_colour[p.colour] = p.memory;
  }
  return m;
}

inline double expected(const Rung& g, double m) { return g.sign * std::pow(m, static_cast<double>(g.lenses_since)); }

// The faces whose tie flips across them.
inline std::vector<std::size_t> flips(const HistoryComplex& K) {
  std::vector<std::size_t> out;
  for (std::size_t f = 0; f < K.lenses.size(); f++)
    if (K.flip(f)) out.push_back(f);
  return out;
}

namespace detail {

// One tie between two beings, chosen among several: the one attested most recently.
inline std::map<std::array<std::size_t, 2>, std::pair<std::size_t, const Rung*>> pairs_at(const HistoryComplex& K, std::size_t occasion) {
  std::map<std::array<std::size_t, 2>, std::pair<std::size_t, const Rung*>> out;
  for (std::size_t t = 0; t < K.rungs.size(); t++)
    for (const auto& g : K.rungs[t]) {
      if (g.occasion != occasion) continue;
      const Tie& tie = K.record.ties[t];
      const std::array<std::size_t, 2> key{std::min(tie.a, tie.b), std::max(tie.a, tie.b)};
      auto it = out.find(key);
      if (it == out.end() || (g.attested && !it->second.second->attested)) out[key] = {t, &g};
    }
  return out;
}

}  // namespace detail

// A triangle of ties at an occasion: three beings present, pairwise tied, all three ties live.
struct Triangle {
  std::size_t occasion = 0;
  std::array<std::size_t, 3> beings{}, ties{};
  int product = 1;      // of the signs there, attested or carried
  double expected = 1;  // of the product, carried signs weighed by their memory
  bool attested = false;
};

inline std::vector<Triangle> triangles(const HistoryComplex& K, const Memory& m) {
  std::vector<Triangle> out;
  for (std::size_t o : K.order) {
    const auto pairs = detail::pairs_at(K, o);
    std::vector<std::size_t> beings;
    for (const auto& p : K.record.occasions[o].present) beings.push_back(p.being);
    std::sort(beings.begin(), beings.end());
    for (std::size_t i = 0; i < beings.size(); i++)
      for (std::size_t j = i + 1; j < beings.size(); j++)
        for (std::size_t k = j + 1; k < beings.size(); k++) {
          const std::array<std::array<std::size_t, 2>, 3> keys{{{beings[i], beings[j]}, {beings[j], beings[k]}, {beings[i], beings[k]}}};
          if (!pairs.count(keys[0]) || !pairs.count(keys[1]) || !pairs.count(keys[2])) continue;
          Triangle T;
          T.occasion = o;
          T.beings = {beings[i], beings[j], beings[k]};
          T.attested = true;
          for (std::size_t q = 0; q < 3; q++) {
            const auto& [t, g] = pairs.at(keys[q]);
            T.ties[q] = t;
            T.product *= g->sign;
            T.expected *= expected(*g, m.of(K.record.ties[t].colour));
            T.attested = T.attested && g->attested;
          }
          out.push_back(T);
        }
  }
  return out;
}

// A closed shell: the lenses of three mutually tied beings between the same two meetings.
struct Shell {
  std::array<std::size_t, 3> beings{}, ties{}, faces{};
  std::size_t from = 0, to = 0;
  int product_from = 1, product_to = 1;
  std::size_t flips = 0;
};

inline std::vector<Shell> shells(const HistoryComplex& K) {
  std::map<std::array<std::size_t, 2>, std::vector<std::size_t>> between;  // (from, to) -> faces
  for (std::size_t f = 0; f < K.lenses.size(); f++) between[{K.lenses[f].from, K.lenses[f].to}].push_back(f);
  std::vector<Shell> out;
  for (const auto& [ends, faces] : between) {
    for (std::size_t i = 0; i < faces.size(); i++)
      for (std::size_t j = i + 1; j < faces.size(); j++)
        for (std::size_t k = j + 1; k < faces.size(); k++) {
          const std::array<std::size_t, 3> fs{faces[i], faces[j], faces[k]};
          std::map<std::size_t, int> degree;
          for (std::size_t f : fs) {
            degree[K.record.ties[K.lenses[f].tie].a]++;
            degree[K.record.ties[K.lenses[f].tie].b]++;
          }
          if (degree.size() != 3 || std::any_of(degree.begin(), degree.end(), [](const auto& d) { return d.second != 2; })) continue;
          Shell S;
          S.from = ends[0];
          S.to = ends[1];
          std::size_t q = 0;
          for (const auto& d : degree) S.beings[q++] = d.first;
          for (std::size_t n = 0; n < 3; n++) {
            const Lens& L = K.lenses[fs[n]];
            S.faces[n] = fs[n];
            S.ties[n] = L.tie;
            S.product_from *= K.rungs[L.tie][L.rung_from].sign;
            S.product_to *= K.rungs[L.tie][L.rung_to].sign;
            S.flips += K.flip(fs[n]);
          }
          out.push_back(S);
        }
  }
  return out;
}

// The network of the moment at a day: the beings and live ties, each at its last value.
struct Link {
  std::size_t tie = 0, a = 0, b = 0;
  int sign = 1;
  double expected = 1;
  std::size_t since = 0, lenses_since = 0;
};

struct MomentTriangle {
  std::array<std::size_t, 3> beings{}, ties{};
  int product = 1;
  double expected = 1;
};

struct Moment {
  long jdn = 0;
  std::vector<std::size_t> nodes;
  std::vector<Link> links;
  std::vector<MomentTriangle> triangles;
  std::size_t frustrated = 0;
  double expected = 0;  // the expected number of frustrated triangles
};

inline Moment moment(const HistoryComplex& K, long jdn, const Memory& m) {
  Moment M;
  M.jdn = jdn;
  const auto& occ = K.record.occasions;
  std::map<std::array<std::size_t, 2>, std::size_t> by_pair;  // -> index into links
  for (std::size_t t = 0; t < K.rungs.size(); t++) {
    if (K.rungs[t].empty()) continue;
    if (occ[K.life[t][0]].date.first > jdn || occ[K.life[t][1]].date.last < jdn) continue;
    const Rung* last = nullptr;
    for (const auto& g : K.rungs[t])
      if (occ[g.occasion].date.first <= jdn) last = &g;
    if (!last) continue;
    const Tie& tie = K.record.ties[t];
    Link l{t, tie.a, tie.b, last->sign, expected(*last, m.of(tie.colour)), last->since, last->lenses_since};
    const std::array<std::size_t, 2> key{std::min(tie.a, tie.b), std::max(tie.a, tie.b)};
    const auto it = by_pair.find(key);
    if (it == by_pair.end()) {
      by_pair[key] = M.links.size();
      M.links.push_back(l);
    } else if (occ[l.since].date.first > occ[M.links[it->second].since].date.first) {
      M.links[it->second] = l;  // of several ties between two beings, the most recently attested
    }
  }
  std::map<std::size_t, std::vector<std::size_t>> neighbours;
  for (const auto& l : M.links) {
    neighbours[l.a].push_back(l.b);
    neighbours[l.b].push_back(l.a);
  }
  for (const auto& [b, _] : neighbours) M.nodes.push_back(b);
  for (const auto& [key, i] : by_pair) {
    const std::size_t a = key[0], b = key[1];
    for (std::size_t c : neighbours[b]) {
      if (c <= b) continue;
      const auto ac = by_pair.find({a, c});
      if (ac == by_pair.end()) continue;
      const Link &x = M.links[i], &y = M.links[by_pair.at({b, c})], &z = M.links[ac->second];
      MomentTriangle T{{a, b, c}, {x.tie, y.tie, z.tie}, x.sign * y.sign * z.sign, x.expected * y.expected * z.expected};
      M.frustrated += T.product < 0;
      M.expected += (1 - T.expected) / 2;
      M.triangles.push_back(T);
    }
  }
  return M;
}

// The network of the moment at the end of each year of the record, the years read in the
// historians' convention (chronology.hpp): Julian before the reform of 1582, Gregorian after.
struct YearMoment {
  long year = 0, jdn = 0;
  std::size_t ties = 0, triangles = 0, frustrated = 0;
  double expected = 0;
};

inline std::vector<YearMoment> moments_by_year(const HistoryComplex& K, const Memory& m) {
  std::vector<YearMoment> out;
  if (K.order.empty()) return out;
  long first = historical_of(K.record.occasions[K.order.front()].date.first).year;
  long last = first;
  for (const auto& o : K.record.occasions) last = std::max(last, historical_of(o.date.last).year);
  for (long y = first; y <= last; y++) {
    const long jdn = jdn_historical(y, 12, 31);
    const Moment M = moment(K, jdn, m);
    out.push_back({y, jdn, M.links.size(), M.triangles.size(), M.frustrated, M.expected});
  }
  return out;
}

}  // namespace history
