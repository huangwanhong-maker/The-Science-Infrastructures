// The engine's protocol: one request, one answer, each a JSON object on one line.
//
//   {"id": 7, "op": "genome.tracks", ...}            a request
//   {"id": 7, "ok": true, "seconds": 0.012, "result": {...}}
//   {"id": 7, "ok": false, "error": "why"}
//
// With --serve the engine reads requests line by line until its input ends, so that the server
// keeps workers alive between requests (and a cluster keeps them on other machines, over ssh,
// with the same pipes); without it, it answers the one request on its input and ends. The
// operations are registered by name; each takes the request and the worker's state (its
// thread pool and its cache of loaded genomes) and returns the result.

#pragma once

#include <science/json.hpp>

#include <calx/scheduler/thread_pool.hpp>

#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace workbench {

using science::Json;

// A request's fault, reported to the page as it is written.
struct Fault : std::runtime_error {
  explicit Fault(const std::string& what) : std::runtime_error(what) {}
};

inline double num(const Json& j, const char* key, double fallback) {
  const Json& v = j[key];
  return v.is_number() ? v.as_number() : fallback;
}
inline long long integer(const Json& j, const char* key, long long fallback) {
  const Json& v = j[key];
  return v.is_number() ? static_cast<long long>(std::llround(v.as_number())) : fallback;
}
inline std::string str(const Json& j, const char* key, const std::string& fallback = {}) { return j[key].text(fallback); }
inline bool flag(const Json& j, const char* key, bool fallback) {
  const Json& v = j[key];
  return v.is_bool() ? v.as_bool() : fallback;
}
// An integer within bounds, or a fault naming the key.
inline long long bounded(const Json& j, const char* key, long long fallback, long long lo, long long hi) {
  const long long v = integer(j, key, fallback);
  if (v < lo || v > hi) throw Fault(std::string(key) + " must lie between " + std::to_string(lo) + " and " + std::to_string(hi));
  return v;
}

// Numbers rounded for the page: a JSON array of values to `digits` significant decimals.
template <class Range>
Json numbers(const Range& xs, int decimals = 4) {
  Json a = Json::array();
  const double scale = std::pow(10.0, decimals);
  for (const auto& x : xs) {
    const double v = static_cast<double>(x);
    a.push(std::isfinite(v) ? std::round(v * scale) / scale : Json());
  }
  return a;
}

struct Worker;
using Operation = std::function<Json(const Json& request, Worker& worker)>;

inline std::map<std::string, Operation>& operations() {
  static std::map<std::string, Operation> ops;
  return ops;
}
// Registers an operation at static initialisation: BENCH_OP("genome.tracks") { ... }.
struct Register {
  Register(const char* name, Operation op) { operations()[name] = std::move(op); }
};
#define BENCH_CONCAT2(a, b) a##b
#define BENCH_CONCAT(a, b) BENCH_CONCAT2(a, b)
// __COUNTER__ is expanded once per operation, so that operations in different headers never
// share a name (as __LINE__ would when two of them stand on the same line).
#define BENCH_OP(name) BENCH_OP_WITH(name, __COUNTER__)
#define BENCH_OP_WITH(name, n) BENCH_OP_NAMED(name, n)
#define BENCH_OP_NAMED(name, n)                                                                \
  static science::Json BENCH_CONCAT(bench_op_, n)(const science::Json&, workbench::Worker&);         \
  static workbench::Register BENCH_CONCAT(bench_reg_, n)(name, BENCH_CONCAT(bench_op_, n));          \
  static science::Json BENCH_CONCAT(bench_op_, n)([[maybe_unused]] const science::Json& request, [[maybe_unused]] workbench::Worker& worker)

}  // namespace workbench
