// The client of the API in C++ (api/README.md §3.5, api/openapi.json):
// header-only, C++17 (it builds as C++20 too), on the standard library and the
// platform's sockets, so that a model, a sensor's program on Linux or a research
// code reads fields and observations, posts readings, follows the stream in real
// time and calls the domains' functions, as the Python client does
// (api/python/science/client.py), in C++'s own spelling.
//
//     #include "science.hpp"                                      // with -I api/cpp; link -pthread only for threads of your own
//     science::Client api("http://127.0.0.1:8800");               // or $SCIENCE_API_URL; the token from $SCIENCE_API_TOKEN
//     science::Field gfs = api.field("gfs-forecast");
//     double t = gfs.sample("tmp2m", 139.69, 35.69, "now+1d");    // K; NaN where there is no value
//     science::Series s = gfs.series("tmp2m", 139.69, 35.69);     // s.times (ISO 8601, UTC), s.values
//     science::Grid g = gfs.grid("tmp2m", "0", science::Box{120, 20, 150, 50});   // float32, rows north to south
//     science::Json sun = api.call("astronomy", "sunPosition", 2460000.5);       // a domain's function, run by the service
//     api.stream("garden", [](const science::Event& e) { std::puts(e.readings.dump().c_str()); return true; });
//
// The service may be the API's own (api/service/server.py, port 8800) or an
// application's server that mounts it (the Earth Atlas's, port 8790): give
// either's address, with or without /api/v1. Times are UTC, as ISO 8601 text
// (an instant without an offset is UTC), "now" or "now+6h"; an empty time is the
// service's default, now. Positions are longitude and latitude in degrees on
// WGS 84; values are in the units the product states, SI unless it says
// otherwise. Errors are thrown as science::ApiError, with the HTTP status (0
// when no answer came) and the service's message; JSON that does not parse, as
// science::JsonError.
//
// Each call opens its own connection (HTTP/1.1 with Connection: close, plain
// http:// only: science/http.hpp says why and what to do for HTTPS). A Client
// holds only its address, its token and its time limit, so one Client may serve
// several threads at once, and copies are cheap.

#pragma once

#include "science/http.hpp"
#include "science/json.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace science {

inline constexpr const char* VERSION = "0.1.0";
inline constexpr const char* DEFAULT_URL = "http://127.0.0.1:8800";
inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

using Query = std::vector<std::pair<std::string, std::string>>;
using Headers = std::vector<std::pair<std::string, std::string>>;

// A box in degrees on WGS 84; a west greater than its east crosses the antimeridian.
struct Box {
  double west, south, east, north;
};

// The value at a place at every time of a product.
struct Series {
  std::string variable, units;
  std::vector<std::string> times;  // ISO 8601, UTC
  std::vector<double> values;      // NaN where there is no value
};

// A variable on its grid, or a box of it: row-major, rows north to south, NaN where there is no value.
struct Grid {
  double lon0 = 0, dlon = 0;
  int nlon = 0;
  double lat0 = 0, dlat = 0;
  int nlat = 0;
  bool periodic = false;
  std::string variable, units, time;  // time: the instant, ISO 8601 UTC
  std::vector<float> values;
  float at(int row, int col) const { return values[static_cast<std::size_t>(row) * static_cast<std::size_t>(nlon) + static_cast<std::size_t>(col)]; }
  double lon(int col) const { return lon0 + col * dlon; }
  double lat(int row) const { return lat0 - row * dlat; }
};

// A file the service made (GeoTIFF, netCDF, CSV, JSON): its bytes, its name and its media type.
struct Download {
  std::string data, filename, content_type;
};

// A reading of a sensor, as a device posts it.
struct Reading {
  std::string station;
  double lon = NaN, lat = NaN;  // omitted when NaN: a station already described needs no place
  std::string time;             // ISO 8601 UTC; empty: now, by this computer's clock
  std::vector<std::pair<std::string, double>> values;  // {variable: value}, in the units the source declared
  Json json() const;
};

// A message of the stream: new readings of a source, as {source, readings} (data), with its event id.
struct Event {
  std::string id, source;
  Json readings, data;
};

struct StreamOptions {
  std::string last_event_id;   // start after this event (it and those before it are not sent again)
  bool reconnect = true;       // after a dropped connection or a silence, connect again, asking for what was missed
  double retry = 3.0;          // seconds before reconnecting (the service may ask for another delay)
  double silence = 45.0;       // seconds without a byte (heartbeats come every 15) after which the connection is dead
  std::function<bool()> idle;  // called at each heartbeat; false stops a quiet stream
};

// An answer as it came: status, headers and body.
struct Response {
  int status = 0;
  Headers headers;
  std::string body;
  std::string header(std::string_view name) const {
    std::string want = http::lower(std::string(name));
    for (const auto& h : headers)
      if (http::lower(h.first) == want) return h.second;
    return {};
  }
};

// ---- text: numbers, lists, boxes and instants ---------------------------------------------------------

// A number as the API reads it: the fewest digits that read back to the same double.
inline std::string number_text(double x) { return Json::number_text(x); }

inline std::string list_text(const std::vector<double>& xs) {
  std::string out;
  for (std::size_t i = 0; i < xs.size(); i++) out += (i ? "," : "") + number_text(xs[i]);
  return out;
}

inline std::string box_text(const std::optional<Box>& b) {
  if (!b) return {};
  return list_text({b->west, b->south, b->east, b->north});
}

namespace detail {
// Days from 1970-01-01 to a date of the proleptic Gregorian calendar, and back
// (H. Hinnant, "chrono-Compatible Low-Level Date Algorithms", 2013: days_from_civil, civil_from_days).
inline long long days_from_civil(long long y, unsigned m, unsigned d) {
  y -= m <= 2;
  const long long era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;
  const unsigned doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
  return era * 146097 + static_cast<long long>(doe) - 719468;
}
inline void civil_from_days(long long z, long long& y, unsigned& m, unsigned& d) {
  z += 719468;
  const long long era = (z >= 0 ? z : z - 146096) / 146097;
  const unsigned doe = static_cast<unsigned>(z - era * 146097);
  const unsigned yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
  const unsigned doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
  const unsigned mp = (5 * doy + 2) / 153;
  d = doy - (153 * mp + 2) / 5 + 1;
  m = mp < 10 ? mp + 3 : mp - 9;
  y = static_cast<long long>(yoe) + era * 400 + (m <= 2);
}
inline std::string env(const char* name) {
#ifdef _MSC_VER
  char* value = nullptr;
  std::size_t n = 0;
  if (_dupenv_s(&value, &n, name) != 0 || !value) return {};
  std::string out(value);
  std::free(value);
  return out;
#else
  const char* value = std::getenv(name);
  return value ? value : "";
#endif
}
}  // namespace detail

// An instant as ISO 8601 UTC to the second: 2026-10-23T06:00:00Z.
inline std::string utc_text(std::chrono::system_clock::time_point t) {
  long long s = std::chrono::duration_cast<std::chrono::seconds>(t.time_since_epoch()).count();
  if (t.time_since_epoch() < std::chrono::seconds(s)) s--;  // to the second below, before 1970 too
  long long days = s >= 0 ? s / 86400 : -((-s + 86399) / 86400);
  long long rest = s - days * 86400;
  long long y;
  unsigned m, d;
  detail::civil_from_days(days, y, m, d);
  char buf[40];
  std::snprintf(buf, sizeof buf, "%04lld-%02u-%02uT%02lld:%02lld:%02lldZ", y, m, d, rest / 3600, rest / 60 % 60, rest % 60);
  return buf;
}
inline std::string utc_now() { return utc_text(std::chrono::system_clock::now()); }

// Seconds since 1970-01-01T00:00:00Z of an ISO 8601 instant (YYYY-MM-DD, then
// THH:MM[:SS[.fraction]], then Z or ±HH:MM; no offset is UTC); NaN when the text is not one.
inline double unix_seconds(std::string_view s) {
  std::size_t k = 0;
  auto num = [&](int width, long long& out) {
    out = 0;
    for (int i = 0; i < width; i++, k++) {
      if (k >= s.size() || s[k] < '0' || s[k] > '9') return false;
      out = out * 10 + (s[k] - '0');
    }
    return true;
  };
  auto sep = [&](char c) { return k < s.size() && s[k] == c ? (k++, true) : false; };
  long long y, mo, d, h = 0, mi = 0, sec = 0;
  if (!num(4, y) || !sep('-') || !num(2, mo) || !sep('-') || !num(2, d)) return NaN;
  if (mo < 1 || mo > 12 || d < 1 || d > 31) return NaN;
  double fraction = 0;
  if (k < s.size() && (s[k] == 'T' || s[k] == 't' || s[k] == ' ')) {
    k++;
    if (!num(2, h) || !sep(':') || !num(2, mi)) return NaN;
    if (sep(':')) {
      if (!num(2, sec)) return NaN;
      if (k < s.size() && (s[k] == '.' || s[k] == ',')) {
        k++;
        double scale = 0.1;
        std::size_t from = k;
        for (; k < s.size() && s[k] >= '0' && s[k] <= '9'; k++, scale /= 10) fraction += (s[k] - '0') * scale;
        if (k == from) return NaN;
      }
    }
    if (h > 24 || mi > 59 || sec > 60) return NaN;
  }
  long long offset = 0;
  if (k < s.size()) {
    if (s[k] == 'Z' || s[k] == 'z') {
      k++;
    } else if (s[k] == '+' || s[k] == '-') {
      int sign = s[k++] == '-' ? -1 : 1;
      long long oh, om = 0;
      if (!num(2, oh)) return NaN;
      sep(':');
      if (k < s.size() && !num(2, om)) return NaN;
      offset = sign * (oh * 3600 + om * 60);
    }
  }
  if (k != s.size()) return NaN;
  long long days = detail::days_from_civil(y, static_cast<unsigned>(mo), static_cast<unsigned>(d));
  return static_cast<double>(days * 86400 + h * 3600 + mi * 60 + sec - offset) + fraction;
}

inline Json Reading::json() const {
  Json r = Json::object({{"station", station}});
  if (!std::isnan(lon) && !std::isnan(lat)) {
    r.set("lon", lon);
    r.set("lat", lat);
  }
  r.set("time", time.empty() ? utc_now() : time);
  Json v = Json::object();
  for (const auto& [name, value] : values) v.set(name, value);
  r.set("values", std::move(v));
  return r;
}

class Field;

// A service of the API: fields, observations, events, tracks, regions,
// verification, the compute bridge and plug-ins.
class Client {
 public:
  // A service at `url` (with or without /api/v1; empty: $SCIENCE_API_URL, else
  // http://127.0.0.1:8800), with the token its writes need (empty:
  // $SCIENCE_API_TOKEN), and the time limit in seconds on connecting and on each
  // wait for an answer's bytes. No request is made here.
  explicit Client(std::string url = {}, std::string token = {}, double timeout = 60.0) : timeout_(timeout) {
    if (url.empty()) url = detail::env("SCIENCE_API_URL");
    if (url.empty()) url = DEFAULT_URL;
    while (!url.empty() && url.back() == '/') url.pop_back();
    if (url.size() >= 4 && url.compare(url.size() - 4, 4, "/api") == 0) url += "/v1";
    else if (url.size() < 7 || url.compare(url.size() - 7, 7, "/api/v1") != 0) url += "/api/v1";
    base_ = http::parse_url(url);
    url_ = url;
    token_ = token.empty() ? detail::env("SCIENCE_API_TOKEN") : std::move(token);
  }

  const std::string& url() const { return url_; }  // the base, ending in /api/v1
  double timeout() const { return timeout_; }

  // ---- requests -----------------------------------------------------------------------------------

  // Any request of the API, by its path under /api/v1, with its query and a JSON
  // body; the answer as it came. An answer of status 300 or more is thrown as ApiError.
  Response request(const std::string& method, const std::string& path, const Query& query = {}, const Json* body = nullptr,
                   const Headers& headers = {}) const {
    http::Socket socket;
    http::Message message(socket);
    std::string text = body ? body->dump() : std::string();
    open(socket, message, method, path, query, body ? &text : nullptr, headers, timeout_);
    Response r;
    r.status = message.status;
    r.headers = message.headers;
    r.body = message.read_all(timeout_);
    return r;
  }

  // Any GET of the API, as JSON.
  Json get(const std::string& path, const Query& query = {}) const { return Json::parse(request("GET", path, query).body); }
  // Any POST of the API with a JSON body, as JSON.
  Json post(const std::string& path, const Json& body) const { return Json::parse(request("POST", path, {}, &body).body); }

  // ---- the service and its catalogue ----------------------------------------------------------------

  // The version, the endpoints and the write policy.
  Json info() const { return get("/"); }
  // The registry's datasets (id, kind, title, licence, parameters) and the index of products.
  Json catalog() const { return get("/catalog"); }

  // ---- fields ------------------------------------------------------------------------------------------

  // The field products: [{name, title, variables, times: {count, first, last}, grid, ...}].
  Json fields() const { return get("/fields")["fields"]; }
  // One field product, described.
  Field field(const std::string& name) const;

  // ---- observations and their stream -----------------------------------------------------------------

  // The observation sets and sensor sources.
  Json observations() const { return get("/observations"); }
  // A source's readings ({readings: [...]}) or a set's rows ({fields, rows: [[...]]}), within a span, a box and a limit (0: none).
  Json observations(const std::string& source, const std::string& since = {}, const std::string& until = {},
                    std::optional<Box> bbox = std::nullopt, int limit = 0) const {
    return get("/observations/" + http::encode(source),
               {{"since", since}, {"until", until}, {"bbox", box_text(bbox)}, {"limit", limit > 0 ? std::to_string(limit) : ""}});
  }

  // Sends readings [{station, lon, lat, time, values: {variable: value}}] to a
  // sensor source, defining it first when `define` is {title, variables,
  // stations?}; the answer {accepted, rejected: [{index, reason}]}. A write: a
  // service with a token needs it.
  Json post_observations(const std::string& source, const Json& readings, const Json& define = {}) const {
    Json body = Json::object({{"readings", readings}});
    if (!define.is_null()) body.set("define", define);
    return post("/observations/" + http::encode(source), body);
  }
  Json post_observations(const std::string& source, const std::vector<Reading>& readings, const Json& define = {}) const {
    Json list = Json::array();
    for (const auto& r : readings) list.push(r.json());
    return post_observations(source, list, define);
  }

  // Follows new readings as they arrive (Server-Sent Events): of one source, or
  // of every source when `source` is empty. Calls on_event for each, until it
  // returns false (or options.idle does, at a heartbeat). After a dropped
  // connection or a silence it reconnects, asking by the last event's id for
  // what it missed (the service keeps the last thousand); an answer of an error
  // (404, 503) is thrown.
  void stream(const std::string& source, const std::function<bool(const Event&)>& on_event, const StreamOptions& options = {}) const;

  // ---- events, tracks, regions, verification -------------------------------------------------------------

  // The event sets (no name), or one: {fields, events: [[row]], window, source, ...}.
  Json events(const std::string& name = {}, const std::string& since = {}, const std::string& until = {},
              std::optional<Box> bbox = std::nullopt, int limit = 0) const {
    if (name.empty()) return get("/events")["events"];
    return get("/events/" + http::encode(name),
               {{"since", since}, {"until", until}, {"bbox", box_text(bbox)}, {"limit", limit > 0 ? std::to_string(limit) : ""}});
  }
  Json tracks(const std::string& name = {}) const { return name.empty() ? get("/tracks")["tracks"] : get("/tracks/" + http::encode(name)); }
  Json regions(const std::string& name = {}) const { return name.empty() ? get("/regions")["regions"] : get("/regions/" + http::encode(name)); }
  Json verification(const std::string& name = {}) const {
    return name.empty() ? get("/verification")["verification"] : get("/verification/" + http::encode(name));
  }

  // ---- the compute bridge ---------------------------------------------------------------------------------

  // The domains and foundations that can be called, with their functions; or one unit.
  Json units(const std::string& unit = {}) const { return unit.empty() ? get("/compute")["units"] : get("/compute/" + http::encode(unit)); }

  // A pure exported function of a domain or foundation, called by the service
  // with `args`, the JSON array of its arguments; its result.
  Json compute(const std::string& unit, const std::string& function, const Json& args) const {
    if (!args.is_array()) throw std::invalid_argument("the arguments of a call are a JSON array");
    return post("/compute/" + http::encode(unit) + "/" + http::encode(function), Json::object({{"args", args}}))["result"];
  }
  // The same, with the arguments as they are: call("astronomy", "sunPosition", 2460000.5).
  template <class... A>
  Json call(const std::string& unit, const std::string& function, A&&... args) const {
    return compute(unit, function, Json::array({Json(std::forward<A>(args))...}));
  }

  // ---- plug-ins and jobs --------------------------------------------------------------------------------

  // The plug-ins the service can run, by application.
  Json plugins() const { return get("/plugins")["applications"]; }
  Json describe_plugin(const std::string& app, const std::string& file) const { return get("/plugins/describe", {{"app", app}, {"file", file}}); }
  Json job(const std::string& id) const { return get("/jobs/" + http::encode(id)); }

  // Runs a server plug-in on its sources `chosen` ({need: "product|level"}; the
  // described defaults otherwise) with its parameters, and returns its job:
  // when `wait`, once it is done, asking every `poll` seconds, for at most
  // `limit` seconds (0: no limit). A failed run throws ApiError (status 0) with the job.
  Json run_plugin(const std::string& app, const std::string& file, const Json& chosen = Json::object(),
                  const Json& params = Json::object(), bool wait = true, double poll = 1.0, double limit = 0,
                  const Json& context = Json::object()) const {
    Json job = post("/plugins/run", Json::object({{"app", app}, {"file", file}, {"chosen", chosen}, {"params", params}, {"context", context}}));
    if (!wait) return job;
    auto start = std::chrono::steady_clock::now();
    while (job["state"].text() == "starting" || job["state"].text() == "running") {
      if (limit > 0 && std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() > limit)
        throw ApiError(0, "job " + job["id"].text() + " is still running after " + number_text(limit) + " s", {}, job);
      std::this_thread::sleep_for(std::chrono::duration<double>(poll));
      job = this->job(job["id"].text());
    }
    if (job["state"].text() == "failed") throw ApiError(0, job["message"].text("the run failed"), {}, job);
    return job;
  }

  // The address of a path and query under /api/v1.
  std::string address(const std::string& path, const Query& query = {}) const {
    return "http://" + authority() + target(path, query);
  }

 private:
  http::Url base_;
  std::string url_, token_;
  double timeout_;

  std::string authority() const {
    std::string host = base_.host.find(':') != std::string::npos ? "[" + base_.host + "]" : base_.host;
    return base_.port == "80" ? host : host + ":" + base_.port;
  }

  std::string target(const std::string& path, const Query& query) const {
    std::string out = base_.path + path;
    char joint = '?';
    for (const auto& [key, value] : query) {
      if (value.empty()) continue;
      out += joint;
      out += http::encode(key) + "=" + http::encode(value);
      joint = '&';
    }
    return out;
  }

  // Sends a request and reads the answer's head, leaving its body to be read;
  // an answer of status 300 or more is thrown as ApiError, with the service's message.
  void open(http::Socket& socket, http::Message& message, const std::string& method, const std::string& path, const Query& query,
            const std::string* body, const Headers& headers, double seconds) const {
    std::string where = address(path, query);
    std::string text = method + " " + target(path, query) + " HTTP/1.1\r\nHost: " + authority() +
                       "\r\nUser-Agent: science-cpp/" + VERSION + "\r\nConnection: close\r\n";
    bool accept = false;
    for (const auto& [key, value] : headers) {
      text += key + ": " + value + "\r\n";
      accept = accept || http::lower(key) == "accept";
    }
    if (!accept) text += "Accept: application/json\r\n";
    if (!token_.empty()) text += "Authorization: Bearer " + token_ + "\r\n";
    if (body) text += "Content-Type: application/json\r\nContent-Length: " + std::to_string(body->size()) + "\r\n";
    text += "\r\n";
    if (body) text += *body;
    try {
      socket.connect(base_.host, base_.port, seconds);
      socket.write(text, seconds);
      message.read_head(seconds);
    } catch (const http::Timeout& e) {
      throw http::Timeout(0, e.message, where);
    } catch (const ApiError& e) {
      throw ApiError(0, e.message, where);
    }
    if (message.status < 300) return;
    std::string raw;
    try {
      raw = message.read_all(seconds);
    } catch (const ApiError&) {  // the status says enough
    }
    std::string why = message.reason;
    try {
      Json j = Json::parse(raw);
      if (j["error"].is_string()) why = j["error"].as_string();
    } catch (const JsonError&) {
      if (!raw.empty()) why = raw.substr(0, 500);
    }
    if (message.status < 400) why = "redirected (" + std::to_string(message.status) + " " + message.reason + " to " + message.header("Location") + "), which this client does not follow";
    throw ApiError(message.status, why, where);
  }
};

// A field product of a service: its description, and its values at places,
// over time and on its grid, computed by the service as the Earth Atlas
// computes them (bilinear in space, linear in time; api/README.md §2.1).
class Field {
 public:
  Field(Client client, Json meta) : client_(std::move(client)), meta_(std::move(meta)) {
    name_ = meta_["name"].text();
    const Json& times = meta_["times"];
    if (times.is_array())
      for (const auto& t : times.as_array()) times_.push_back(t.text());
  }

  const std::string& name() const { return name_; }
  const Json& meta() const { return meta_; }                    // grid, times, variables with their units and quanta, source
  const std::vector<std::string>& times() const { return times_; }  // ISO 8601, UTC
  // A variable's description, or null.
  const Json& variable(const std::string& var) const {
    for (const auto& v : meta_["variables"].is_array() ? meta_["variables"].as_array() : empty())
      if (v["name"].text() == var) return v;
    return meta_["no such key"];
  }

  // Values at places and an instant; NaN where there is none. More than 500
  // places are asked by POST, which takes up to 100 000 at once.
  std::vector<double> sample(const std::string& var, const std::vector<double>& lon, const std::vector<double>& lat,
                             const std::string& time = {}) const {
    Json r;
    if (std::max(lon.size(), lat.size()) > 500) {
      Json body = Json::object({{"var", var}, {"lon", Json(lon)}, {"lat", Json(lat)}});
      if (!time.empty()) body.set("time", time);
      r = client_.post(path() + "/sample", body);
    } else {
      r = client_.get(path() + "/sample", {{"var", var}, {"lon", list_text(lon)}, {"lat", list_text(lat)}, {"time", time}});
    }
    return numbers(r["values"]);
  }
  // The value at one place.
  double sample(const std::string& var, double lon, double lat, const std::string& time = {}) const {
    std::vector<double> v = sample(var, std::vector<double>{lon}, std::vector<double>{lat}, time);
    return v.empty() ? NaN : v[0];
  }

  // The value at a place at every time of the product.
  Series series(const std::string& var, double lon, double lat) const {
    Json r = client_.get(path() + "/series", {{"var", var}, {"lon", number_text(lon)}, {"lat", number_text(lat)}});
    Series s{var, r["units"].text(), {}, numbers(r["values"])};
    for (const auto& t : r["times"].as_array()) s.times.push_back(t.text());
    return s;
  }

  // The variable on its grid, or a box of it, every `stride`-th point, at a
  // stored time (an index, as its digits: "0") or an instant; as float32.
  Grid grid(const std::string& var, const std::string& time = {}, std::optional<Box> bbox = std::nullopt, int stride = 1) const {
    static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559, "float must be IEEE 754 binary32");
    Response r = client_.request("GET", path() + "/grid",
                                 {{"var", var}, {"time", time}, {"bbox", box_text(bbox)},
                                  {"stride", stride != 1 ? std::to_string(stride) : ""}, {"format", "f32"}},
                                 nullptr, {{"Accept", "application/octet-stream"}});
    Json described = Json::parse(r.header("X-Grid"));
    const Json& g = described["grid"];
    Grid out;
    out.lon0 = g["lon0"].as_number();
    out.dlon = g["dlon"].as_number();
    out.nlon = static_cast<int>(g["nlon"].as_number());
    out.lat0 = g["lat0"].as_number();
    out.dlat = g["dlat"].as_number();
    out.nlat = static_cast<int>(g["nlat"].as_number());
    out.periodic = g["periodic"].is_bool() && g["periodic"].as_bool();
    out.variable = var;
    out.units = described["units"].text();
    out.time = described["time"].text();
    std::size_t cells = static_cast<std::size_t>(out.nlon) * static_cast<std::size_t>(out.nlat);
    if (r.body.size() != 4 * cells)
      throw ApiError(0, "a grid of " + std::to_string(cells) + " cells came with " + std::to_string(r.body.size()) + " bytes",
                     client_.address(path() + "/grid"));
    out.values.resize(cells);
    const auto* b = reinterpret_cast<const unsigned char*>(r.body.data());
    for (std::size_t k = 0; k < cells; k++, b += 4) {  // little-endian, whatever this computer's order
      std::uint32_t u = std::uint32_t(b[0]) | std::uint32_t(b[1]) << 8 | std::uint32_t(b[2]) << 16 | std::uint32_t(b[3]) << 24;
      std::memcpy(&out.values[k], &u, 4);
    }
    return out;
  }

  // The variable at a stored time as a file: GeoTIFF (tif), netCDF (nc), CSV or
  // JSON; time "all" gives every time, in netCDF or JSON.
  Download download(const std::string& var, const std::string& time = {}, const std::string& format = "tif",
                    std::optional<Box> bbox = std::nullopt, int stride = 1) const {
    Response r = client_.request("GET", path() + "/export",
                                 {{"var", var}, {"time", time}, {"format", format}, {"bbox", box_text(bbox)},
                                  {"stride", stride != 1 ? std::to_string(stride) : ""}},
                                 nullptr, {{"Accept", "*/*"}});
    Download d{std::move(r.body), name_ + "-" + var + "." + format, r.header("Content-Type")};
    std::string disposition = r.header("Content-Disposition");
    std::size_t at = disposition.find("filename=");
    if (at != std::string::npos) {
      std::string f = disposition.substr(at + 9);
      if (!f.empty() && f[0] == '"') f = f.substr(1, f.find('"', 1) - 1);
      else f = f.substr(0, f.find(';'));
      if (!f.empty()) d.filename = f;
    }
    return d;
  }
  // The same, written to a file at `path`; the number of bytes written.
  std::size_t export_to(const std::string& file, const std::string& var, const std::string& time = {},
                        const std::string& format = "tif", std::optional<Box> bbox = std::nullopt, int stride = 1) const {
    Download d = download(var, time, format, bbox, stride);
    std::ofstream out(file, std::ios::binary);
    out.write(d.data.data(), static_cast<std::streamsize>(d.data.size()));
    if (!out) throw std::runtime_error("cannot write " + file);
    return d.data.size();
  }

 private:
  Client client_;
  Json meta_;
  std::string name_;
  std::vector<std::string> times_;

  std::string path() const { return "/fields/" + http::encode(name_); }
  static const Json::Array& empty() {
    static const Json::Array none;
    return none;
  }
  static std::vector<double> numbers(const Json& list) {
    std::vector<double> out;
    if (!list.is_array()) return out;
    out.reserve(list.size());
    for (const auto& v : list.as_array()) out.push_back(v.as_number());  // null: NaN
    return out;
  }
};

inline Field Client::field(const std::string& name) const { return Field(*this, get("/fields/" + http::encode(name))); }

inline void Client::stream(const std::string& source, const std::function<bool(const Event&)>& on_event,
                           const StreamOptions& options) const {
  const std::string path = source.empty() ? "/stream" : "/observations/" + http::encode(source) + "/stream";
  std::string last = options.last_event_id;
  double retry = options.retry;
  for (;;) {
    try {
      http::Socket socket;
      http::Message message(socket);
      Headers headers{{"Accept", "text/event-stream"}};
      if (!last.empty()) headers.emplace_back("Last-Event-ID", last);
      open(socket, message, "GET", path, {}, nullptr, headers, timeout_);
      http::EventStream events;
      events.last_id = last;
      std::string chunk;
      bool go = true;
      while (go) {
        chunk.clear();
        if (!message.read_body(chunk, options.silence)) break;  // the service ended the stream
        go = events.feed(
            chunk,
            [&](const http::ServerEvent& e) {
              last = e.id;
              if (e.type != "readings") return true;
              Event out;
              out.id = e.id;
              try {
                out.data = Json::parse(e.data);
              } catch (const JsonError&) {
                return true;  // an event that is not JSON is passed over
              }
              out.source = out.data["source"].text();
              out.readings = out.data["readings"];
              return on_event(out);
            },
            [&](std::string_view) { return options.idle ? options.idle() : true; });
        last = events.last_id;
        if (events.retry_ms >= 0) retry = events.retry_ms / 1000.0;
      }
      if (!go) return;
    } catch (const http::Timeout&) {
      if (!options.reconnect) throw;
    } catch (const ApiError& e) {
      if (e.status != 0 || !options.reconnect) throw;  // an answer of an error is not worth retrying
    }
    if (!options.reconnect) return;
    std::this_thread::sleep_for(std::chrono::duration<double>(retry));
  }
}

}  // namespace science
