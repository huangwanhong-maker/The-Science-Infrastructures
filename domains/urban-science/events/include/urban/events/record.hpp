// The record of a city's events (the infrastructure's ADR 0067; the design, development/ideas/2026-11-04-event-
// based-gis-on-a-2-complex.md, §1–2): couplings, the persistent things that hold primitives together (volumes:
// buildings, street spaces, junction spaces); primitives, the relations or interfaces that persist between exactly
// two couplings (party walls, frontages, the cross-sections of streets); and events, each with its label in time,
// its type and what it touches, begins and ends.
//
// Time is a Modified Julian Date in UTC (days; MJD 0 = 1858-11-17T00:00Z), and a label is a span: a year of
// construction is the year, a census the year of the census, a permit its first and last day. Events are ordered by
// the start of their labels, then by their order in the record, which is the source's (the design's axiom A4: the
// order is structure, clock time a label).
//
// What a record must satisfy (the design's §1.5) is checked by `problems`:
//   - every index in range; a primitive between two distinct couplings (a loop is subdivided before it comes here);
//   - a cell begun at most once and ended at most once, and not ended before it is begun;
//   - a primitive's life within both its couplings' lives;
//   - an event touching only cells alive at it (begun at or before it, not yet ended), within the record's window;
//   - an event touching something: a tick of the clock names the couplings whose worldlines it cuts.
// A cell no event begins is alive from the window's start; one no event ends, to its end.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace urban::events {

struct Span {
  double from = 0, to = 0;  // MJD, UTC
};

// The MJD of a civil date (proleptic Gregorian), 0 h UTC.
inline double mjd(int year, int month = 1, int day = 1) {
  // Fliegel and Van Flandern (1968), the Julian Day Number of a Gregorian date, less 2400001 for the MJD.
  const long a = (14 - month) / 12, y = year + 4800 - a, m = month + 12 * a - 3;
  const long jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;
  return static_cast<double>(jdn - 2400001);
}
inline Span year_span(int year) { return {mjd(year), mjd(year + 1)}; }

struct Coupling {
  std::string id, type, name;   // type: building, street, junction, …
  std::array<double, 2> at{};   // its anchor, local east and north (m)
  double base = 0, top = 0;     // the heights it stands between (m), for drawing
};

struct Primitive {
  std::string id, type;         // type: party-wall, frontage, cross-section, …
  std::size_t a = 0, b = 0;     // its two couplings
  double area = 0;              // m², an interface's
  std::array<double, 3> normal{};  // unit, from a toward b (zero where it has none)
};

struct Event {
  std::string id, type, kind;   // type: birth, death, alteration, regime, institution, occurrence, observation, tick
  Span when;
  std::vector<std::size_t> couplings;   // the couplings it touches
  std::vector<std::size_t> primitives;  // the primitives it meets
  std::vector<std::size_t> begins_couplings, ends_couplings, begins_primitives, ends_primitives;
  std::string source, provenance = "recorded";  // provenance: recorded, inferred, modelled, synthetic
};

struct Record {
  Span window;
  std::vector<Coupling> couplings;
  std::vector<Primitive> primitives;
  std::vector<Event> events;
};

// The events in the order of the record: by the start of their labels, then as they were given.
inline std::vector<std::size_t> event_order(const Record& R) {
  std::vector<std::size_t> o(R.events.size());
  std::iota(o.begin(), o.end(), 0);
  std::stable_sort(o.begin(), o.end(), [&](std::size_t i, std::size_t j) { return R.events[i].when.from < R.events[j].when.from; });
  return o;
}

// Each cell's life as positions in the order of the events: (begun, ended), where −1 is the window's start and
// the number of events its end.
struct Lives {
  std::vector<std::array<long, 2>> couplings, primitives;
};

inline Lives lives(const Record& R, const std::vector<std::size_t>& order) {
  const long n = static_cast<long>(order.size());
  Lives L;
  L.couplings.assign(R.couplings.size(), {-1, n});
  L.primitives.assign(R.primitives.size(), {-1, n});
  for (long k = 0; k < n; k++) {
    const auto& e = R.events[order[static_cast<std::size_t>(k)]];
    for (auto c : e.begins_couplings) if (c < L.couplings.size()) L.couplings[c][0] = k;
    for (auto c : e.ends_couplings) if (c < L.couplings.size()) L.couplings[c][1] = k;
    for (auto p : e.begins_primitives) if (p < L.primitives.size()) L.primitives[p][0] = k;
    for (auto p : e.ends_primitives) if (p < L.primitives.size()) L.primitives[p][1] = k;
  }
  return L;
}

// What is wrong with a record, each problem named with its event or cell; empty when the record is sound.
inline std::vector<std::string> problems(const Record& R) {
  std::vector<std::string> out;
  const std::size_t C = R.couplings.size(), P = R.primitives.size();
  for (std::size_t p = 0; p < P; p++) {
    const auto& q = R.primitives[p];
    if (q.a >= C || q.b >= C) out.push_back("primitive " + q.id + ": a coupling out of range");
    else if (q.a == q.b) out.push_back("primitive " + q.id + ": its two couplings are one (a loop)");
  }
  std::vector<int> begun_c(C, 0), ended_c(C, 0), begun_p(P, 0), ended_p(P, 0);
  for (const auto& e : R.events) {
    const auto range = [&](const std::vector<std::size_t>& v, std::size_t n, const char* what) {
      for (auto i : v)
        if (i >= n) out.push_back("event " + e.id + ": " + what + " out of range");
    };
    range(e.couplings, C, "a coupling"), range(e.begins_couplings, C, "a coupling begun"), range(e.ends_couplings, C, "a coupling ended");
    range(e.primitives, P, "a primitive"), range(e.begins_primitives, P, "a primitive begun"), range(e.ends_primitives, P, "a primitive ended");
    for (auto c : e.begins_couplings) if (c < C) begun_c[c]++;
    for (auto c : e.ends_couplings) if (c < C) ended_c[c]++;
    for (auto p : e.begins_primitives) if (p < P) begun_p[p]++;
    for (auto p : e.ends_primitives) if (p < P) ended_p[p]++;
    if (e.when.from < R.window.from || e.when.from > R.window.to) out.push_back("event " + e.id + ": outside the record's window");
    if (e.when.to < e.when.from) out.push_back("event " + e.id + ": its label ends before it begins");
    if (e.couplings.empty() && e.primitives.empty() && e.begins_couplings.empty() && e.ends_couplings.empty() && e.begins_primitives.empty() && e.ends_primitives.empty())
      out.push_back("event " + e.id + ": touches nothing (a tick names the cells it cuts)");
  }
  for (std::size_t c = 0; c < C; c++) {
    if (begun_c[c] > 1) out.push_back("coupling " + R.couplings[c].id + ": begun twice");
    if (ended_c[c] > 1) out.push_back("coupling " + R.couplings[c].id + ": ended twice");
  }
  for (std::size_t p = 0; p < P; p++) {
    if (begun_p[p] > 1) out.push_back("primitive " + R.primitives[p].id + ": begun twice");
    if (ended_p[p] > 1) out.push_back("primitive " + R.primitives[p].id + ": ended twice");
  }
  if (!out.empty()) return out;
  const auto order = event_order(R);
  const auto L = lives(R, order);
  for (std::size_t c = 0; c < C; c++)
    if (L.couplings[c][1] < L.couplings[c][0]) out.push_back("coupling " + R.couplings[c].id + ": ended before it is begun");
  for (std::size_t p = 0; p < P; p++) {
    const auto& q = R.primitives[p];
    const auto& l = L.primitives[p];
    if (l[1] < l[0]) out.push_back("primitive " + q.id + ": ended before it is begun");
    for (auto c : {q.a, q.b})
      if (l[0] < L.couplings[c][0] || l[1] > L.couplings[c][1]) out.push_back("primitive " + q.id + ": lives outside its coupling " + R.couplings[c].id);
  }
  for (long k = 0; k < static_cast<long>(order.size()); k++) {
    const auto& e = R.events[order[static_cast<std::size_t>(k)]];
    const auto alive = [&](const std::array<long, 2>& l) { return l[0] <= k && k <= l[1]; };
    for (auto c : e.couplings) if (!alive(L.couplings[c])) out.push_back("event " + e.id + ": touches coupling " + R.couplings[c].id + " outside its life");
    for (auto p : e.primitives) if (!alive(L.primitives[p])) out.push_back("event " + e.id + ": meets primitive " + R.primitives[p].id + " outside its life");
  }
  return out;
}

}  // namespace urban::events
