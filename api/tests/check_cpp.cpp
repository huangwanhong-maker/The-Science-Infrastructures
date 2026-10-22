// Checks of the C++ client (api/cpp/science.hpp), built by api/Makefile and run
// by api/tests/run-tests.mjs. Offline, its JSON (numbers with exponents and at
// the limits of doubles, escapes, surrogate pairs, nesting, malformed texts,
// the C locale) and its HTTP (chunked bodies, framing, interim answers, the
// events of a stream, instants, addresses). Against a running service, its
// answers held to the Python client's for the same requests, read on standard
// input as check_cpp.py writes them: samples and series exactly, grids exactly
// as float32, calls of the domains exactly; then readings posted and received on
// a stream within a second, a reconnection that receives what it missed, a
// plug-in run to its end, and the errors. One JSON object a line, as
// check_client.py reports.
//
//     check_cpp --offline
//     check_cpp --probe <url>                          exit 0 when the service answers
//     check_cpp <url> <sensor source> [token] < reference.json

#include "science.hpp"

#include <atomic>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iterator>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

using science::Json;

namespace {

const char* PREFIX = "c++ client: ";

void emit(const Json& record) {
  std::string line = record.dump();
  std::fwrite(line.data(), 1, line.size(), stdout);
  std::fputc('\n', stdout);
  std::fflush(stdout);
}
void check(const std::string& name, double value, double limit, const std::string& unit = "") {
  emit(Json::object({{"kind", "check"}, {"name", PREFIX + name}, {"value", value}, {"limit", limit}, {"unit", unit}}));
}
void affirm(const std::string& name, bool ok, const std::string& detail = "") {
  emit(Json::object({{"kind", "assert"}, {"name", PREFIX + name}, {"ok", ok}, {"detail", detail}}));
}
void skip(const std::string& what) { emit(Json::object({{"kind", "skip"}, {"name", PREFIX + what}})); }

// The largest difference of two lists, NaN matching NaN (else infinite, which reports as a failure).
double difference(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) return INFINITY;
  double worst = 0;
  for (std::size_t i = 0; i < a.size(); i++) {
    if (std::isnan(a[i]) || std::isnan(b[i])) {
      if (std::isnan(a[i]) != std::isnan(b[i])) return INFINITY;
      continue;
    }
    worst = std::max(worst, std::fabs(a[i] - b[i]));
  }
  return worst;
}

std::vector<double> numbers(const Json& list) {
  std::vector<double> out;
  for (const auto& v : list.as_array()) out.push_back(v.as_number());
  return out;
}

bool refused(const std::string& text) {
  try {
    Json::parse(text);
    return false;
  } catch (const science::JsonError&) {
    return true;
  }
}

// ---- offline: JSON ------------------------------------------------------------------------------

void json_checks() {
  // Numbers: each text reads as strtod reads it, and each double writes in digits that read back to it.
  const char* texts[] = {"0", "-0", "1", "-1", "0.1", "1e5", "1E5", "1e+5", "1e-5", "-1.5e-7", "6.02214076e23", "2.5E-3",
                         "1.7976931348623157e308", "2.2250738585072014e-308", "5e-324", "4.9406564584124654e-324",
                         "9007199254740992", "9007199254740993", "123456789012345678901234567890", "0.30000000000000004",
                         "1761234567890", "-0.0e0", "1e400", "-1e400", "1e-400"};
  int read_right = 0, round_trip = 0, total = 0;
  std::string wrong;
  for (const char* t : texts) {
    total++;
    double x = Json::parse(t).as_number(), want = std::strtod(t, nullptr);
    if (x == want && std::signbit(x) == std::signbit(want)) read_right++;
    else wrong += std::string(" ") + t;
    if (!std::isfinite(x)) {
      round_trip += Json(x).dump() == "null";
      continue;
    }
    double back = Json::parse(Json(x).dump()).as_number();
    if (back == x && std::signbit(back) == std::signbit(x)) round_trip++;
    else wrong += std::string(" ") + t + "->" + Json(x).dump();
  }
  affirm("json: numbers with exponents and at the limits of doubles read as strtod reads them", read_right == total, wrong);
  affirm("json: every double written reads back to itself; NaN and infinities write as null", round_trip == total, wrong);
  affirm("json: the fewest digits", Json(0.1).dump() == "0.1" && Json(1.0 / 3).dump() == "0.3333333333333333" &&
                                        Json(1761234567890.0).dump() == "1761234567890" && Json(1e22).dump() == "1e+22",
         Json(0.1).dump() + " " + Json(1.0 / 3).dump() + " " + Json(1761234567890.0).dump());
  affirm("json: NaN is written as null, and null reads as NaN",
         Json(std::nan("")).dump() == "null" && std::isnan(Json::parse("null").as_number()));

  // Strings: escapes, \u to UTF-8, surrogate pairs joined, lone surrogates as U+FFFD.
  struct Case {
    const char* json;
    std::string want;
  } strings[] = {
      {R"("\"\\\/\b\f\n\r\t")", "\"\\/\b\f\n\r\t"},
      {R"("\u00e9\u00E9")", "\xC3\xA9\xC3\xA9"},
      {R"("\u65e5\u672c")", "\xE6\x97\xA5\xE6\x9C\xAC"},
      {R"("\ud83d\ude00")", "\xF0\x9F\x98\x80"},
      {R"("\uD834\uDD1E clef")", "\xF0\x9D\x84\x9E clef"},
      {R"("\ud83d")", "\xEF\xBF\xBD"},
      {R"("\ude00x")", "\xEF\xBF\xBDx"},
      {R"("\ud83d\u0041")", "\xEF\xBF\xBD" "A"},
      {R"("\ud83d\ud83d\ude00")", "\xEF\xBF\xBD\xF0\x9F\x98\x80"},
      {"\"\xE6\x97\xA5\xE6\x9C\xAC raw\"", "\xE6\x97\xA5\xE6\x9C\xAC raw"},
      {R"("\u0000")", std::string(1, '\0')},
  };
  int ok = 0, n = 0;
  std::string bad;
  for (const auto& c : strings) {
    n++;
    try {
      if (Json::parse(c.json).as_string() == c.want) ok++;
      else bad += std::string(" ") + c.json;
    } catch (const science::JsonError& e) {
      bad += std::string(" ") + c.json + " (" + e.what() + ")";
    }
  }
  affirm("json: escapes, \\u to UTF-8, surrogate pairs joined, a lone surrogate as U+FFFD", ok == n, bad);
  std::string written = Json(std::string("a\"b\\c\n\x01\x1f\xC3\xA9")).dump();
  affirm("json: strings are written with their escapes, and UTF-8 as it is", written == "\"a\\\"b\\\\c\\n\\u0001\\u001f\xC3\xA9\"", written);
  std::string again = Json::parse(written).as_string();
  affirm("json: a string written reads back to itself", again == "a\"b\\c\n\x01\x1f\xC3\xA9");

  // Nesting, order and lookup.
  Json j = Json::parse(R"( {"z": 1, "a": {"b": [1, {"c": [true, false, null]}, "x"]}, "e": [], "o": {}} )");
  affirm("json: nested objects and arrays", j["a"]["b"][1]["c"][0].as_bool() && j["a"]["b"][1]["c"][2].is_null() && j["a"]["b"][2].text() == "x" &&
                                               j["e"].is_array() && j["e"].size() == 0 && j["o"].is_object() && j["o"].size() == 0);
  affirm("json: objects keep the order of their keys", j.dump() == R"({"z":1,"a":{"b":[1,{"c":[true,false,null]},"x"]},"e":[],"o":{}})", j.dump());
  affirm("json: an absent key and an index beyond the end read as null", j["missing"].is_null() && j["a"]["b"][9].is_null() && j["z"]["deeper"].is_null());
  affirm("json: objects are equal whatever the order of their keys", Json::parse(R"({"a":1,"b":[2,{"c":3,"d":4}]})") == Json::parse(R"({"b":[2,{"d":4,"c":3}],"a":1})"));
  std::string deep(512, '['), deeper(513, '[');
  deep += "1" + std::string(512, ']');
  deeper += "1" + std::string(513, ']');
  affirm("json: 512 levels deep are read, 513 refused", !refused(deep) && refused(deeper));

  const char* malformed[] = {"", " ", "[1,]", "{\"a\":1,}", "01", "-01", "1.", ".5", "+1", "1e", "1e+", "-", "[1 2]", "{a:1}", "'x'",
                             "\"abc", "\"a\x01\"", "nul", "tru", "[1]x", "{\"a\" 1}", "\"\\x\"", "\"\\u12\"", "\"\\u12g4\"", "NaN",
                             "Infinity", "[", "{", "{\"a\":}", "[,1]", "//", "0x10"};
  int refusals = 0, count = 0;
  std::string accepted;
  for (const char* m : malformed) {
    count++;
    if (refused(m)) refusals++;
    else accepted += std::string(" |") + m + "|";
  }
  affirm("json: every malformed text is refused", refusals == count, accepted);
  try {
    Json::parse("{\"a\": [1, 2,, 3]}");
  } catch (const science::JsonError& e) {
    affirm("json: a refusal says where", std::string(e.what()).find("byte 12") != std::string::npos, e.what());
  }
  bool threw = false;
  try {
    Json::parse("\"x\"").as_number();
  } catch (const science::JsonError&) {
    threw = true;
  }
  affirm("json: a value read as what it is not throws", threw);

  Json built = Json::object({{"var", "tmp2m"}, {"lon", Json(std::vector<double>{139.69, -0.13})}, {"n", 3}, {"ok", true}, {"none", nullptr}});
  built.set("lon", Json::array({1.5})).set("extra", "é");
  affirm("json: values built and changed", built.dump() == R"({"var":"tmp2m","lon":[1.5],"n":3,"ok":true,"none":null,"extra":"é"})", built.dump());

  // A program that sets a locale with a decimal comma still reads and writes JSON's point.
  const char* locales[] = {"de_DE.UTF-8", "de_DE.utf8", "fr_FR.UTF-8", "German_Germany.1252"};
  const char* set = nullptr;
  for (const char* l : locales)
    if ((set = std::setlocale(LC_NUMERIC, l))) break;
  if (set && *std::localeconv()->decimal_point == ',') {
    bool fine = Json::parse("1.5").as_number() == 1.5 && Json(2.25).dump() == "2.25";
    std::setlocale(LC_NUMERIC, "C");
    affirm("json: under a locale with a decimal comma, numbers still read and write with a point", fine);
  } else {
    std::setlocale(LC_NUMERIC, "C");
    skip("json: no locale with a decimal comma is installed, so its check is not made");
  }
}

// ---- offline: HTTP, the stream's events, instants -------------------------------------------------------

// A text handed out a few bytes at a time, as a slow network would.
class Text : public science::http::Input {
 public:
  explicit Text(std::string s, std::size_t step = 3) : s_(std::move(s)), step_(step) {}
  std::size_t read(char* out, std::size_t n, double) override {
    std::size_t k = std::min({n, step_, s_.size() - at_});
    std::memcpy(out, s_.data() + at_, k);
    at_ += k;
    return k;
  }

 private:
  std::string s_;
  std::size_t step_, at_ = 0;
};

// An answer's body, read from its text in pieces of `step` bytes; its status and headers in head.
struct Head {
  int status = 0;
  std::string type, grid;
};
std::string body_of(const std::string& raw, std::size_t step, Head* head = nullptr) {
  Text in(raw, step);
  science::http::Message m(in);
  m.read_head(1);
  std::string body = m.read_all(1);
  if (head) *head = Head{m.status, m.header("CONTENT-TYPE"), m.header("x-grid")};
  return body;
}

void http_checks() {
  using namespace science;
  const std::string chunked =
      "HTTP/1.1 200 OK\r\nContent-Type: text/plain\r\ntransfer-encoding: Chunked\r\nX-Grid: {\"a\": 1}\r\n\r\n"
      "5\r\nhello\r\n7;name=value\r\n, world\r\n1A\r\n abcdefghijklmnopqrstuvwxy\r\n0\r\nTrailer: x\r\n\r\n";
  bool all = true;
  for (std::size_t step : {1u, 2u, 3u, 7u, 4096u}) all = all && body_of(chunked, step) == "hello, world abcdefghijklmnopqrstuvwxy";
  affirm("http: a chunked body is decoded, whatever the pieces it arrives in", all);
  Head head;
  body_of(chunked, 5, &head);
  affirm("http: headers are found in any case", head.grid == "{\"a\": 1}" && head.type == "text/plain");
  affirm("http: a body of a Content-Length",
         body_of("HTTP/1.0 200 OK\r\nContent-Length: 4\r\n\r\nabcdEXTRA", 2) == "abcd");
  affirm("http: a body to the end of the connection", body_of("HTTP/1.0 200 OK\nServer: x\n\nto the end", 4) == "to the end");
  affirm("http: an interim 100 Continue before the answer",
         body_of("HTTP/1.1 100 Continue\r\n\r\nHTTP/1.1 201 Created\r\nContent-Length: 2\r\n\r\nok", 3, &head) == "ok" && head.status == 201);
  bool cut = false;
  try {
    body_of("HTTP/1.1 200 OK\r\nContent-Length: 10\r\n\r\nshort", 3);
  } catch (const ApiError& e) {
    cut = e.status == 0;
  }
  affirm("http: a body cut short is an error", cut);
  bool notHttp = false;
  try {
    body_of("SSH-2.0-OpenSSH\r\n\r\n", 3);
  } catch (const ApiError&) {
    notHttp = true;
  }
  affirm("http: an answer that is not HTTP is an error", notHttp);

  // The events of a stream, fed in pieces, with comments, CR LF, multi-line data, ids and retry.
  const std::string stream =
      ": readings of every source, a heartbeat every 15 s\nretry: 3000\n\n"
      "id: 41\r\nevent: readings\r\ndata: {\"a\"\r\ndata: :1}\r\n\r\n"
      ": heartbeat 2026-10-23T06:00:00Z\n\n"
      "event: other\ndata\n\n"
      "data: no id given\n\n"
      "id: 42\nevent: readings\n\n"
      "event: readings\ndata:{\"b\":2}\n\n";
  for (std::size_t step : {1u, 5u, 1000u}) {
    http::EventStream parser;
    std::vector<http::ServerEvent> got;
    int comments = 0;
    for (std::size_t k = 0; k < stream.size(); k += step)
      parser.feed(std::string_view(stream).substr(k, step), [&](const http::ServerEvent& e) { got.push_back(e); return true; },
                  [&](std::string_view) { comments++; return true; });
    bool fine = got.size() == 4 && got[0].type == "readings" && got[0].data == "{\"a\"\n:1}" && got[0].id == "41" &&
                got[1].type == "other" && got[1].data.empty() && got[2].type == "message" && got[2].id == "41" &&
                got[3].data == "{\"b\":2}" && got[3].id == "42" && comments == 2 && parser.retry_ms == 3000 && parser.last_id == "42";
    if (step == 1) affirm("sse: events, comments, CR LF, data over lines, ids that persist and retry, fed byte by byte", fine);
    else if (step == 5) affirm("sse: the same, fed five bytes at a time", fine);
    else affirm("sse: the same, fed at once", fine);
  }
  http::EventStream stopper;
  int seen = 0;
  bool went = stopper.feed("data: 1\n\ndata: 2\n\n", [&](const http::ServerEvent&) { return ++seen < 1; }, [](std::string_view) { return true; });
  affirm("sse: a callback that returns false stops the stream", !went && seen == 1);

  // Instants: Unix seconds of ISO 8601 and back.
  affirm("time: Unix seconds of ISO 8601 instants",
         unix_seconds("1970-01-01T00:00:00Z") == 0 && unix_seconds("2000-01-01T12:00:00Z") == 946728000 &&
             unix_seconds("2026-10-23T06:00:00Z") == 1792735200 && unix_seconds("2026-10-23T15:00:00+09:00") == 1792735200 &&
             unix_seconds("2026-10-23T06:00Z") == 1792735200 && unix_seconds("2026-10-23") == 1792713600 &&
             unix_seconds("1969-12-31T23:59:59Z") == -1 && unix_seconds("2024-02-29T00:00:00.25Z") == 1709164800.25 &&
             std::isnan(unix_seconds("yesterday")) && std::isnan(unix_seconds("2026-13-01")));
  bool back = true;
  for (long long s : {0LL, -1LL, 951782400LL, 1709164800LL, 4107542400LL, 4107456000LL + 86399, 1792735200LL, -2208988800LL}) {
    std::string t = utc_text(std::chrono::system_clock::time_point(std::chrono::seconds(s)));
    back = back && unix_seconds(t) == static_cast<double>(s);
  }
  affirm("time: ISO 8601 of instants reads back, across leap days and centuries", back && utc_text(std::chrono::system_clock::time_point{}) == "1970-01-01T00:00:00Z");

  // Addresses and their encoding.
  http::Url u = http::parse_url("http://[::1]:8800/x/api/v1");
  affirm("http: addresses, IPv6 literals among them", u.host == "::1" && u.port == "8800" && u.path == "/x/api/v1" &&
                                                         http::parse_url("http://localhost").port == "80");
  affirm("http: query values are percent-encoded, '+' among them", http::encode("now+1d") == "now%2B1d" && http::encode("a b,c/é") == "a%20b%2Cc%2F%C3%A9");
  bool https = false;
  try {
    Client("https://example.org");
  } catch (const std::invalid_argument& e) {
    https = std::string(e.what()).find("TLS proxy") != std::string::npos;
  }
  affirm("http: an https address is refused, saying what to do", https);
  Client c("http://example.org:8790/");
  affirm("client: the base gains /api/v1", c.url() == "http://example.org:8790/api/v1" && Client("http://h/api").url() == "http://h/api/v1" &&
                                               c.address("/fields/x/sample", {{"time", "now+1d"}, {"bbox", ""}}) ==
                                                   "http://example.org:8790/api/v1/fields/x/sample?time=now%2B1d");
}

// ---- against a service --------------------------------------------------------------------------------

void fields_checks(const science::Client& api, const Json& ref) {
  if (!ref.has("sample")) {
    skip("fields: gfs-forecast is not in the cache");
    return;
  }
  science::Field gfs = api.field("gfs-forecast");
  affirm("fields: a description with its times", gfs.times().size() == ref["series"]["times"].size() && gfs.variable("tmp2m")["units"].text() == "K");
  const Json& s = ref["sample"];
  std::vector<double> got = gfs.sample(s["var"].as_string(), numbers(s["lon"]), numbers(s["lat"]), s["time"].as_string());
  check("fields: samples equal the Python client's (ten places, the poles and the antimeridian)", difference(got, numbers(s["values"])), 0, "K");
  double one = gfs.sample("tmp2m", s["lon"][0].as_number(), s["lat"][0].as_number(), s["time"].as_string());
  check("fields: one place gives one value, the same", std::fabs(one - s["values"][0].as_number()), 0, "K");
  const Json& m = ref["many"];
  std::vector<double> many = gfs.sample(m["var"].as_string(), numbers(m["lon"]), numbers(m["lat"]), m["time"].as_string());
  check("fields: 600 places, asked by POST, equal the Python client's by GET", difference(many, numbers(m["values"])), 0, "%");

  const Json& r = ref["series"];
  science::Series series = gfs.series(r["var"].as_string(), r["lon"].as_number(), r["lat"].as_number());
  check("fields: a series equals the Python client's", difference(series.values, numbers(r["values"])), 0, "K");
  bool times = series.times.size() == r["times"].size() && series.units == r["units"].text();
  for (std::size_t k = 0; times && k < series.times.size(); k++)
    times = series.times[k] == r["times"][k].text() && science::unix_seconds(series.times[k]) == r["times_unix"][k].as_number();
  affirm("fields: ... at the same times, read as the same Unix seconds", times);

  for (const char* which : {"grid", "blend"}) {
    const Json& g = ref[which];
    const Json& b = g["bbox"];
    science::Grid grid = gfs.grid(g["var"].as_string(), g["time"].as_string(),
                                  science::Box{b[0].as_number(), b[1].as_number(), b[2].as_number(), b[3].as_number()},
                                  static_cast<int>(g["stride"].as_number()));
    const Json& d = g["grid"];
    bool shape = grid.nlon == static_cast<int>(d["nlon"].as_number()) && grid.nlat == static_cast<int>(d["nlat"].as_number()) &&
                 grid.lon0 == d["lon0"].as_number() && grid.lat0 == d["lat0"].as_number() && grid.dlon == d["dlon"].as_number() &&
                 grid.dlat == d["dlat"].as_number() && grid.time == g["at"].text() && grid.units == "K";
    std::vector<double> values(grid.values.begin(), grid.values.end());
    if (std::string(which) == "grid") {
      affirm("fields: a grid across the antimeridian, every second point, has the Python client's shape", shape,
             std::to_string(grid.nlon) + " x " + std::to_string(grid.nlat));
      check("fields: ... and its float32 values", difference(values, numbers(g["values"])), 0, "K");
    } else {
      affirm("fields: a grid between stored times has the Python client's shape", shape);
      check("fields: ... and its blended float32 values", difference(values, numbers(g["values"])), 0, "K");
    }
  }

  try {
    science::Download tif = gfs.download("tmp2m", "0", "tif", science::Box{122, 24, 150, 46});
    affirm("fields: an export as GeoTIFF, named", tif.data.size() > 8 && tif.data.compare(0, 4, std::string("II*\0", 4)) == 0 &&
                                                       tif.filename.rfind("gfs-forecast-tmp2m", 0) == 0 && tif.filename.find(".tif") != std::string::npos,
           tif.filename);
  } catch (const science::ApiError& e) {
    if (e.status == 501) skip("fields: export answers 501 until data/formats/export.py exists");
    else affirm("fields: an export as GeoTIFF", false, e.what());
  }

  int status = -1;
  try {
    api.field("no-such-field");
  } catch (const science::ApiError& e) {
    status = e.status;
  }
  affirm("errors: a missing field is an ApiError of status 404", status == 404, std::to_string(status));
  std::string message;
  try {
    gfs.sample("no-such-variable", 0, 0);
  } catch (const science::ApiError& e) {
    status = e.status;
    message = e.message;
  }
  affirm("errors: ... and a missing variable, with the service's message", status == 404 && message.find("no-such-variable") != std::string::npos, message);
}

void compute_checks(const science::Client& api, const Json& ref) {
  if (!ref.has("compute")) {
    skip("compute: the Python client could not call the domains: " + ref["compute_error"].text());
    return;
  }
  int equal = 0;
  std::string differ;
  for (const auto& c : ref["compute"].as_array()) {
    Json got = api.compute(c["unit"].as_string(), c["function"].as_string(), c["args"]);
    if (got == c["result"]) equal++;
    else differ += " " + c["function"].as_string() + ": " + got.dump() + " != " + c["result"].dump();
  }
  affirm("compute: " + std::to_string(ref["compute"].size()) + " calls (astronomy, meteorology, geodesy, physis) equal the Python client's",
         equal == static_cast<int>(ref["compute"].size()), differ);
  Json sun = api.call("astronomy", "sunPosition", 2460000.5);
  affirm("compute: call() with plain arguments equals compute() with a JSON array", sun == ref["compute"][0]["result"], sun.dump());
  int status = 0;
  try {
    api.call("astronomy", "noSuchFunction");
  } catch (const science::ApiError& e) {
    status = e.status;
  }
  affirm("compute: a function that is not there answers 404", status == 404, std::to_string(status));
}

void product_checks(const science::Client& api, const Json& ref) {
  Json info = api.info();
  affirm("the service's index", info["api"].text() == "v1" && info["version"] == ref["version"], info["version"].dump());
  Json names = Json::array();
  for (const auto& f : api.fields().as_array()) names.push(f["name"]);
  affirm("the fields are the Python client's", names == ref["fields"]);
  if (ref.has("events")) {
    Json e = api.events(ref["events"]["name"].as_string(), {}, {}, std::nullopt, 5);
    affirm("the event set " + ref["events"]["name"].as_string() + ", its latest five, equals the Python client's", e["events"] == ref["events"]["events"]);
  } else {
    skip("no event set in the cache");
  }
  Json v = Json::array();
  for (const auto& x : api.verification().as_array()) v.push(x["name"]);
  affirm("the verifications are the Python client's", v == ref["verification"]);
  Json catalog = api.catalog();
  affirm("the catalogue lists datasets and products", catalog.is_object() && catalog.size() > 0);
}

// The plug-ins: listed, described, and one run through to its field, as the JavaScript client's check runs it.
void plugin_checks(const science::Client& api) {
  Json listed;
  try {
    listed = api.plugins();
  } catch (const science::ApiError& e) {
    skip("plug-ins: " + e.message);
    return;
  }
  bool found = false;
  for (const auto& p : listed["earth-atlas"].is_array() ? listed["earth-atlas"].as_array() : Json::Array{})
    found = found || p["file"].text() == "plugins/vorticity.py";
  affirm("plug-ins: the Earth Atlas's vorticity plug-in is listed", found);
  if (!found) return;
  Json d = api.describe_plugin("earth-atlas", "plugins/vorticity.py");
  affirm("plug-ins: its description, with its needs", d["id"].text() == "vorticity-py" && d["needs"][0]["name"].text() == "wind", d["id"].dump());
  if (d["needs"][0]["chosen"].is_null()) {
    skip("plug-ins: no wind in the cache to run the plug-in on");
    return;
  }
  auto t0 = std::chrono::steady_clock::now();
  Json job = api.run_plugin("earth-atlas", "plugins/vorticity.py", Json::object(), Json::object(), true, 0.3, 300);
  double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::string field;
  for (const auto& r : job["results"].as_array())
    if (r["kind"].text() == "field") field = r["field"].text();
  affirm("plug-ins: a run waited for to its end (" + science::number_text(std::round(seconds * 10) / 10) + " s), with a field among its results",
         job["state"].text() == "done" && !field.empty(), job["state"].text());
  if (!field.empty()) affirm("plug-ins: the run's field is served among the fields", api.field(field).variable("zeta").is_object(), field);
}

// Readings posted, received on a stream within a second, and again after a reconnection.
void stream_checks(const science::Client& api, const std::string& source) {
  using clock = std::chrono::steady_clock;
  std::atomic<bool> connected{false};
  std::mutex lock;
  std::vector<std::pair<clock::time_point, science::Event>> got;
  std::string failure;
  auto deadline = clock::now() + std::chrono::seconds(8);
  std::thread listener([&] {
    science::StreamOptions options;
    options.idle = [&] {
      connected = true;
      return clock::now() < deadline;
    };
    try {
      api.stream(source, [&](const science::Event& e) {
        std::lock_guard<std::mutex> hold(lock);
        got.emplace_back(clock::now(), e);
        return false;
      }, options);
    } catch (const std::exception& e) {
      std::lock_guard<std::mutex> hold(lock);
      failure = e.what();
    }
  });
  while (!connected && clock::now() < deadline) std::this_thread::sleep_for(std::chrono::milliseconds(20));
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  Json define = Json::object({{"title", "API check (C++)"}, {"variables", Json::array({"air_temperature", "relative_humidity"})}});
  science::Reading reading{"cpp1", 139.76, 35.68, science::utc_text(std::chrono::system_clock::now() - std::chrono::seconds(1)),
                           {{"air_temperature", 291.4}, {"relative_humidity", 61}}};
  auto t0 = clock::now();
  Json posted = api.post_observations(source, std::vector<science::Reading>{reading}, define);
  listener.join();
  affirm("observations: a source defined and a reading posted", posted["accepted"].as_number() == 1 && posted["defined"].as_bool(), posted.dump());
  std::string first_id;
  if (!got.empty()) {
    const science::Event& e = got[0].second;
    check("observations: a posted reading reaches the stream", std::chrono::duration<double>(got[0].first - t0).count(), 1.0, "s");
    affirm("observations: the event carries its id, its source and the reading",
           !e.id.empty() && e.source == source && e.readings[0]["station"].text() == "cpp1" &&
               e.readings[0]["values"]["air_temperature"].as_number() == 291.4,
           e.data.dump());
    first_id = e.id;
  } else {
    affirm("observations: a posted reading reaches the stream", false, failure.empty() ? "nothing within 8 s" : failure);
  }

  // Missed while no one listened: received after reconnecting with the last id seen.
  science::Reading second{"cpp2", 139.7, 35.6, science::utc_now(), {{"air_temperature", 290.0}}};
  Json again = api.post_observations(source, std::vector<science::Reading>{second});
  affirm("observations: a second reading, of a new station, is accepted", again["accepted"].as_number() == 1, again.dump());
  if (!first_id.empty()) {
    science::Event missed;
    science::StreamOptions options;
    options.last_event_id = first_id;
    auto until = clock::now() + std::chrono::seconds(5);
    options.idle = [&] { return clock::now() < until; };
    api.stream({}, [&](const science::Event& e) {
      missed = e;
      return false;
    }, options);
    affirm("observations: a stream opened with the last id seen receives what it missed", missed.readings[0]["station"].text() == "cpp2", missed.data.dump());
  }
  Json back = api.observations(source, {}, {}, std::nullopt, 1);
  affirm("observations: the latest reading reads back", back["readings"].size() == 1 && back["readings"][0]["station"].text() == "cpp2", back.dump().substr(0, 300));
  science::Reading hot{"cpp1", 139.76, 35.68, science::utc_now(), {{"air_temperature", 1000}}};
  Json answer = api.post_observations(source, std::vector<science::Reading>{hot});
  affirm("observations: an implausible reading is refused, with its reason",
         answer["accepted"].as_number() == 0 && !answer["rejected"][0]["reason"].text().empty(), answer.dump());
}

void unreachable_checks() {
  int status = -1;
  std::string why;
  auto t0 = std::chrono::steady_clock::now();
  try {
    science::Client("http://127.0.0.1:1", {}, 2).info();
  } catch (const science::ApiError& e) {
    status = e.status;
    why = e.message;
  }
  double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  // Refused at once, or (behind some firewalls, and in WSL's mirrored network) unanswered until the limit.
  affirm("errors: a service that is not there is an ApiError of status 0, within the time limit", status == 0 && seconds < 3, why);
  status = -1;
  try {
    science::Client("http://no-such-host.invalid", {}, 5).info();
  } catch (const science::ApiError& e) {
    status = e.status;
    why = e.message;
  }
  affirm("errors: a name that does not resolve is an ApiError of status 0", status == 0, why);
}

int online(const std::string& url, const std::string& source, const std::string& token) {
  std::string text((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
  Json ref;
  try {
    ref = Json::parse(text);
  } catch (const science::JsonError& e) {
    affirm("the Python client's reference values are read", false, e.what());
    return 1;
  }
  science::Client api(url, token);
  struct Part {
    const char* name;
    std::function<void()> run;
  } parts[] = {
      {"the service and its products", [&] { product_checks(api, ref); }},
      {"fields", [&] { fields_checks(api, ref); }},
      {"compute", [&] { compute_checks(api, ref); }},
      {"observations and the stream", [&] { stream_checks(api, source); }},
      {"plug-ins", [&] { plugin_checks(api); }},
      {"errors", [&] { unreachable_checks(); }},
  };
  for (const auto& p : parts) {
    try {
      p.run();
    } catch (const std::exception& e) {
      affirm(std::string(p.name) + " ran to the end", false, e.what());
    }
  }
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  std::string first = argc > 1 ? argv[1] : "";
  if (first == "--probe") {
    try {
      science::Client(argc > 2 ? argv[2] : "", {}, 3).info();
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "%s\n", e.what());
      return 3;
    }
  }
  json_checks();
  http_checks();
  if (first == "--offline" || first.empty()) return 0;
  return online(first, argc > 2 ? argv[2] : "api-check-cpp", argc > 3 ? argv[3] : "");
}
