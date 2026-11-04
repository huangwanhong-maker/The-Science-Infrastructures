// What a city remembers of its uses (ADR 0067; the design, §3.2 and §4; stage E4): the uses of a district's buildings
// and premises observed at every census, read as chains on a classified set (physis's nested_chain.hpp), with the
// record's probability on its event complex (physis's event_chains.hpp) and the rates of its kinds of events.
//
// The classifications.
//   ANZSIC 2006 (the Australian and New Zealand Standard Industrial Classification): the whole; its divisions; its
//   subdivisions (two digits); its groups (three); its classes (four). CLUE's vacant space (0000) is a division of its
//   own: a premises falling vacant has left every division.
//   CLUE's space uses: the whole; a family, the text before " - " where CLUE writes one ("Retail - Shop" and "Retail
//   - Showroom" are Retail; "Unoccupied - Unused" and "Unoccupied - Under Renovation" are Unoccupied), else the use
//   alone; the use. No family is invented beyond what CLUE writes.
// The law π of each is the uses' frequency over the panel's observations.
//
// The panels. A unit (a base property for its building's use; a premises with one establishment for its class) seen
// at consecutive censuses gives one transition each time, x then y a number of years later. The fit is physis's
// fit_rates; its memories μ_m say how fast each level is forgotten, ln 2 / μ_m how long it is half kept. The flat chain
// (one level, a redraw from anywhere) is the null: twice the difference of log-likelihoods against χ² on L − 1
// degrees of freedom says whether the classification's levels are remembered apart.
//
// On the complex. Each building's use, observed at its censuses, is a chain on the segments of its worldline in K: at
// its first vertex its law π, at each of its censuses the observation (diag of the use seen), at its other events the
// identity, at its last vertex the sum. Its probability by physis's contraction is the product of its chain's
// transitions exactly (the design's Theorem 4.2, with Theorem 4.5 for the events that only cut the worldline), which
// the checks hold on the real complex.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <physis/statistical_mechanics/event_chains.hpp>
#include <physis/statistical_mechanics/nested_chain.hpp>
#include <urban/events/complex.hpp>
#include <urban/events/record.hpp>

namespace urban::events {

namespace sm = physis::statistical_mechanics;

struct UseClassification {
  std::vector<std::string> uses;    // the states
  std::vector<std::string> levels;  // the levels' names, from the whole to the uses
  std::vector<std::vector<std::string>> block_names;  // block_names[ℓ][x]
  sm::Classification C;
  std::map<std::string, std::size_t> index;
  std::size_t of(const std::string& use) const {
    const auto it = index.find(use);
    if (it == index.end()) throw std::invalid_argument("urban::events::UseClassification: no use \"" + use + "\"");
    return it->second;
  }
};

namespace detail {

// A classification from each use's names at the levels between the whole and the use, and the uses' counts.
inline UseClassification classify(std::vector<std::string> level_names, const std::map<std::string, std::vector<std::string>>& path,
                                  const std::map<std::string, double>& counts) {
  UseClassification U;
  U.levels = std::move(level_names);
  for (const auto& [use, n] : counts)
    if (n > 0 && path.count(use)) U.index[use] = U.uses.size(), U.uses.push_back(use);
  const std::size_t N = U.uses.size(), inner = U.levels.size() - 2;
  if (N == 0) throw std::invalid_argument("urban::events::classify: no use observed");
  U.block_names.assign(U.levels.size(), std::vector<std::string>(N));
  U.C.block.assign(U.levels.size(), std::vector<std::size_t>(N, 0));
  for (std::size_t x = 0; x < N; x++) {
    const auto& p = path.at(U.uses[x]);
    if (p.size() != inner) throw std::invalid_argument("urban::events::classify: a use's path is not one name a level");
    U.block_names[0][x] = "all";
    for (std::size_t l = 0; l < inner; l++) {
      // A block is named by its whole path, so that two subdivisions of one name in two divisions stay two.
      std::string key;
      for (std::size_t k = 0; k <= l; k++) key += (k ? "/" : "") + p[k];
      U.block_names[l + 1][x] = key;
    }
    U.block_names.back()[x] = U.uses[x];
  }
  for (std::size_t l = 0; l < U.levels.size(); l++) {
    std::map<std::string, std::size_t> ids;
    for (std::size_t x = 0; x < N; x++) {
      const auto& k = U.block_names[l][x];
      if (!ids.count(k)) ids[k] = ids.size();
      U.C.block[l][x] = ids[k];
    }
  }
  double total = 0;
  for (const auto& u : U.uses) total += counts.at(u);
  U.C.law.resize(N);
  for (std::size_t x = 0; x < N; x++) U.C.law[x] = counts.at(U.uses[x]) / total;
  sm::validate(U.C);
  return U;
}

}  // namespace detail

// ANZSIC: each code with its division, subdivision and group (vacant space, division none, its own division).
struct AnzsicCode {
  std::string code, division, subdivision, group;
};

inline UseClassification anzsic_classification(const std::vector<AnzsicCode>& codes, const std::map<std::string, double>& counts) {
  std::map<std::string, std::vector<std::string>> path;
  for (const auto& c : codes) {
    const std::string division = c.division.empty() ? "vacant" : c.division;
    path[c.code] = {division, c.subdivision, c.group};
  }
  return detail::classify({"all", "division", "subdivision", "group", "class"}, path, counts);
}

inline std::string family_of(const std::string& use) {
  const auto k = use.find(" - ");
  return k == std::string::npos ? use : use.substr(0, k);
}

inline UseClassification prefix_classification(const std::map<std::string, double>& counts) {
  std::map<std::string, std::vector<std::string>> path;
  for (const auto& [use, n] : counts) path[use] = {family_of(use)};
  return detail::classify({"all", "family", "use"}, path, counts);
}

// An observation: a unit seen with a use at a census.
struct Observation {
  std::string unit, use;
  int year = 0;
};

// The uses' counts over the observations (the law's weights).
inline std::map<std::string, double> use_counts(const std::vector<Observation>& O) {
  std::map<std::string, double> n;
  for (const auto& o : O) n[o.use] += 1;
  return n;
}

// Each unit's consecutive observations as transitions, aggregated by (x, y, years between).
inline std::vector<sm::Transition> transitions_of(const UseClassification& U, std::vector<Observation> O) {
  std::sort(O.begin(), O.end(), [](const Observation& a, const Observation& b) { return a.unit != b.unit ? a.unit < b.unit : a.year < b.year; });
  std::map<std::tuple<std::size_t, std::size_t, int>, double> agg;
  for (std::size_t i = 0; i + 1 < O.size(); i++) {
    if (O[i].unit != O[i + 1].unit || O[i + 1].year <= O[i].year) continue;
    if (!U.index.count(O[i].use) || !U.index.count(O[i + 1].use)) continue;
    agg[{U.of(O[i].use), U.of(O[i + 1].use), O[i + 1].year - O[i].year}] += 1;
  }
  std::vector<sm::Transition> T;
  for (const auto& [k, n] : agg) T.push_back({std::get<0>(k), std::get<1>(k), static_cast<double>(std::get<2>(k)), n});
  return T;
}

// The deepest level at which two uses share a block (the last level: the same use).
inline std::size_t agreement(const UseClassification& U, std::size_t x, std::size_t y) {
  std::size_t l = 0;
  while (l + 1 < U.C.block.size() && U.C.block[l + 1][x] == U.C.block[l + 1][y]) l++;
  return l;
}

struct UseMemory {
  std::vector<std::string> levels;
  sm::NestedFit fit;
  std::vector<double> half_life, half_life_error;  // years, for m = 1 … L (index m − 1)
  std::vector<double> agree;                        // transitions whose two uses share level ℓ and no deeper, ℓ = 0 … L
  double flat_log_likelihood = 0, ratio = 0;        // the one-level chain; twice the gain over it
  std::size_t degrees = 0;                          // of freedom of the ratio
  double transitions = 0;
  std::size_t units = 0, uses = 0;
};

// The memory of a classification fitted to its transitions.
inline UseMemory remember(const UseClassification& U, const std::vector<sm::Transition>& T, std::size_t units = 0) {
  UseMemory M;
  M.levels = U.levels;
  M.uses = U.uses.size(), M.units = units;
  const std::size_t L = U.C.levels();
  M.agree.assign(L + 1, 0.0);
  for (const auto& t : T) M.agree[agreement(U, t.x, t.y)] += t.count, M.transitions += t.count;
  // A start from the counts: a change that keeps level ℓ and no deeper, per transition-year, as λ_ℓ.
  double years = 0;
  for (const auto& t : T) years += t.count * t.tau;
  std::vector<double> start(L, 0.0);
  for (std::size_t l = 0; l < L; l++) start[l] = std::max(1e-4, M.agree[l] / std::max(1.0, years));
  M.fit = sm::fit_rates(U.C, T, start);
  for (std::size_t m = 1; m <= L; m++) {
    const double mu = M.fit.mu[m];
    M.half_life.push_back(mu > 0 ? std::log(2.0) / mu : std::numeric_limits<double>::infinity());
    M.half_life_error.push_back(mu > 0 ? M.half_life.back() * M.fit.mu_error[m] / mu : std::nan(""));
  }
  // The flat chain: one redraw from anywhere at the rate that best fits.
  sm::Classification flat;
  flat.law = U.C.law;
  flat.block = {U.C.block.front(), U.C.block.back()};
  const auto F = sm::fit_rates(flat, T, {std::max(1e-4, (M.transitions - M.agree[L]) / std::max(1.0, years))});
  M.flat_log_likelihood = F.log_likelihood;
  M.ratio = 2 * (M.fit.log_likelihood - F.log_likelihood);
  M.degrees = L - 1;
  return M;
}

// A building's use between two of its censuses under the fitted chain: the expected changes nobody recorded, and the
// weights with which each level's mode carries the transition (the design's foam).
struct UseSpan {
  int from = 0, to = 0;
  std::string first, last;
  double expected_changes = 0;
  std::vector<double> weights;  // m = 0 … L
};

inline std::vector<UseSpan> spans(const UseClassification& U, const sm::NestedChain& X, std::vector<Observation> O) {
  std::sort(O.begin(), O.end(), [](const Observation& a, const Observation& b) { return a.year < b.year; });
  std::vector<UseSpan> out;
  for (std::size_t i = 0; i + 1 < O.size(); i++) {
    if (!U.index.count(O[i].use) || !U.index.count(O[i + 1].use) || O[i + 1].year <= O[i].year) continue;
    const std::size_t x = U.of(O[i].use), y = U.of(O[i + 1].use);
    const double tau = O[i + 1].year - O[i].year;
    out.push_back({O[i].year, O[i + 1].year, O[i].use, O[i + 1].use, X.expected_changes(x, y, tau), X.mode_weights(x, y, tau)});
  }
  return out;
}

// The chains of the buildings' uses on the event complex: for each coupling with observations (its use at a census
// vertex, by the record's event), the segments of its worldline as cells, its law at its first vertex, each census an
// observation of the use seen, its other events the identity, its last vertex the sum. Time in years.
inline sm::EventChains chains_on_complex(const EventComplex& X, const sm::NestedChain& chain,
                                         const std::map<std::size_t, std::map<std::size_t, std::size_t>>& seen) {
  // seen: coupling ↦ (vertex ↦ use observed there)
  sm::EventChains M;
  const std::size_t N = chain.C.states();
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> cell;  // (coupling, segment) ↦ cell
  for (const auto& [c, obs] : seen) {
    (void)obs;
    for (auto e : X.K.segments_of[c]) {
      const auto& s = X.K.segments[e];
      cell[{c, e}] = M.cells.size();
      M.cells.push_back({N, chain.kernel_matrix((X.time[s.to] - X.time[s.from]) / 365.2425)});
    }
  }
  std::vector<double> identity(N * N, 0.0), sum(N, 1.0);
  for (std::size_t i = 0; i < N; i++) identity[i * N + i] = 1;
  for (std::size_t v = 0; v < X.K.cw.vertices; v++)
    for (auto c : X.encounters.present[v]) {
      const auto it = seen.find(c);
      if (it == seen.end()) continue;
      const auto& w = X.K.worldline[c];
      const auto& segs = X.K.segments_of[c];
      const std::size_t k = static_cast<std::size_t>(std::find(w.begin(), w.end(), v) - w.begin());
      std::vector<std::size_t> in, out;
      if (k > 0) in.push_back(cell.at({c, segs[k - 1]}));
      if (k + 1 < w.size()) out.push_back(cell.at({c, segs[k]}));
      if (in.empty() && out.empty()) continue;
      const auto o = it->second.find(v);
      std::vector<double> K;
      if (in.empty()) K = chain.C.law;                      // the first vertex: the law
      else if (out.empty()) K = sum;                         // the last: the sum
      else K = identity;                                     // an event that only cuts the worldline
      if (o != it->second.end()) {                           // a census: the use seen
        if (in.empty()) for (std::size_t x = 0; x < N; x++) K[x] *= x == o->second ? 1.0 : 0.0;
        else if (out.empty()) for (std::size_t x = 0; x < N; x++) K[x] *= x == o->second ? 1.0 : 0.0;
        else for (std::size_t x = 0; x < N; x++) K[x * N + x] = x == o->second ? 1.0 : 0.0;
      }
      M.pieces.push_back({in, out, K});
    }
  return M;
}

// The rates of the record's kinds of events, per year of exposure: the events of a kind touching couplings of a type,
// over the years those couplings stood in the window; Poisson's errors.
struct EventRate {
  std::string type, kind, of;  // the event's type and kind; the type of coupling it touches
  double events = 0, exposure = 0, rate = 0, error = 0;
};

inline std::vector<EventRate> event_rates(const Record& R, const EventComplex& X) {
  std::map<std::string, double> exposure;  // coupling type ↦ coupling-years
  for (std::size_t c = 0; c < R.couplings.size(); c++) {
    const auto& t = R.couplings[c].type;
    const std::string kind = t == "building" ? "building" : t == "junction" ? "junction" : "street";
    exposure[kind] += (X.time[X.life_c[c][1]] - X.time[X.life_c[c][0]]) / 365.2425;
  }
  std::map<std::tuple<std::string, std::string, std::string>, double> n;
  for (const auto& e : R.events) {
    if (e.type == "observation" || e.type == "tick") continue;
    std::set<std::string> of;
    for (auto c : e.couplings) of.insert(R.couplings[c].type == "building" ? "building" : R.couplings[c].type == "junction" ? "junction" : "street");
    for (auto c : e.begins_couplings) of.insert(R.couplings[c].type == "building" ? "building" : "street");
    if (of.empty() && !e.primitives.empty()) of.insert("street");
    for (const auto& o : of) n[{e.type, e.kind, o}] += 1;
  }
  std::vector<EventRate> out;
  for (const auto& [k, m] : n) {
    const double E = exposure[std::get<2>(k)];
    out.push_back({std::get<0>(k), std::get<1>(k), std::get<2>(k), m, E, E > 0 ? m / E : 0, E > 0 ? std::sqrt(m) / E : 0});
  }
  std::sort(out.begin(), out.end(), [](const EventRate& a, const EventRate& b) { return a.events > b.events; });
  return out;
}

}  // namespace urban::events
