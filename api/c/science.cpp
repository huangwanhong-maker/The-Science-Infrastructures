// The C interface of the API's client (science.h), implemented over the C++
// header (api/cpp/science.hpp). Every function catches what the C++ throws and
// returns it as a code, with its message and status kept per thread; nothing
// thrown crosses into C. Memory handed to the caller is malloc's, so that
// science_free is free.

#include "science.h"

#include "../cpp/science.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <new>
#include <string>

struct science_client {
  science::Client api;
};

namespace {

thread_local std::string last_error;
thread_local int last_status = 0;

// A failure of the caller's own input, told apart from the service's.
struct Argument : std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

int fail(int code, const std::string& why, int status = 0) {
  last_error = why;
  last_status = status;
  return code;
}

// Runs f, turning what it throws into a code; SCIENCE_OK, with no error, otherwise.
template <class F>
int guard(F&& f) {
  try {
    int code = f();
    if (code == SCIENCE_OK) {
      last_error.clear();
      last_status = 0;
    }
    return code;
  } catch (const science::ApiError& e) {
    if (e.status) return fail(SCIENCE_ERROR_HTTP, e.what(), e.status);
    return fail(e.job.is_null() ? SCIENCE_ERROR_UNREACHABLE : SCIENCE_ERROR_FAILED, e.what());
  } catch (const Argument& e) {
    return fail(SCIENCE_ERROR_ARGUMENT, e.what());
  } catch (const std::invalid_argument& e) {
    return fail(SCIENCE_ERROR_ARGUMENT, e.what());
  } catch (const science::JsonError& e) {
    return fail(SCIENCE_ERROR_ANSWER, std::string("the answer is not what the API sends: ") + e.what());
  } catch (const std::bad_alloc&) {
    return fail(SCIENCE_ERROR_FAILED, "out of memory");
  } catch (const std::exception& e) {
    return fail(SCIENCE_ERROR_FAILED, e.what());
  } catch (...) {
    return fail(SCIENCE_ERROR_FAILED, "an unknown failure");
  }
}

char* copy(const std::string& s) {
  char* out = static_cast<char*>(std::malloc(s.size() + 1));
  if (!out) throw std::bad_alloc();
  std::memcpy(out, s.data(), s.size());
  out[s.size()] = '\0';
  return out;
}

template <class T>
T* allocate(std::size_t n) {
  T* out = static_cast<T*>(std::malloc((n ? n : 1) * sizeof(T)));
  if (!out) throw std::bad_alloc();
  return out;
}

std::string text(const char* s) { return s ? s : ""; }

science::Json json_argument(const char* text, const char* what) {
  if (!text) throw Argument(std::string(what) + " is NULL");
  try {
    return science::Json::parse(text);
  } catch (const science::JsonError& e) {
    throw Argument(std::string(what) + " is " + e.what());
  }
}

void need(const void* p, const char* what) {
  if (!p) throw Argument(std::string(what) + " is NULL");
}

void bounded_copy(char* out, std::size_t size, const std::string& s) {
  std::size_t n = std::min(size - 1, s.size());
  std::memcpy(out, s.data(), n);
  out[n] = '\0';
}

// The value at a path of keys and indices joined by '.'; ARGUMENT when there is none.
const science::Json& at_path(const science::Json& root, const char* path) {
  const science::Json* v = &root;
  std::string_view rest = path ? path : "";
  std::string_view whole = rest;
  while (!rest.empty()) {
    std::size_t dot = rest.find('.');
    std::string_view part = rest.substr(0, dot);
    rest = dot == std::string_view::npos ? std::string_view() : rest.substr(dot + 1);
    if (v->is_array()) {
      if (part.empty() || part.size() > 9 || !std::all_of(part.begin(), part.end(), [](unsigned char c) { return c >= '0' && c <= '9'; }))
        throw Argument("at " + std::string(whole) + ": an array is indexed by a number, not '" + std::string(part) + "'");
      std::size_t i = static_cast<std::size_t>(std::stoul(std::string(part)));
      if (i >= v->size()) throw Argument("at " + std::string(whole) + ": index " + std::string(part) + " is beyond the array's " + std::to_string(v->size()));
      v = &(*v)[i];
    } else if (v->is_object()) {
      v = v->find(part);
      if (!v) throw Argument("at " + std::string(whole) + ": no key '" + std::string(part) + "'");
    } else {
      throw Argument("at " + std::string(whole) + ": '" + std::string(part) + "' is asked of a value that has no parts");
    }
  }
  return *v;
}

}  // namespace

extern "C" {

science_client* science_connect(const char* url, const char* token, double timeout_seconds) {
  try {
    auto* client = new science_client{science::Client(text(url), text(token), timeout_seconds > 0 ? timeout_seconds : 60.0)};
    last_error.clear();
    last_status = 0;
    return client;
  } catch (const science::ApiError& e) {  // an address this client cannot use
    fail(SCIENCE_ERROR_ARGUMENT, e.what());
  } catch (const std::exception& e) {
    fail(SCIENCE_ERROR_FAILED, e.what());
  }
  return nullptr;
}

void science_close(science_client* client) { delete client; }

const char* science_last_error(void) { return last_error.c_str(); }
int science_last_status(void) { return last_status; }
void science_free(void* memory) { std::free(memory); }

int science_get(science_client* client, const char* path, char** json) {
  return guard([&] {
    need(client, "the client");
    need(path, "the path");
    need(json, "the place for the answer");
    *json = nullptr;
    // The query is the caller's, already encoded: it is sent as given.
    science::Response r = client->api.request("GET", path);
    *json = copy(r.body);
    return SCIENCE_OK;
  });
}

int science_info(science_client* client, char** json) { return science_get(client, "/", json); }

int science_sample(science_client* client, const char* field, const char* variable, const char* time, const double* lon,
                   const double* lat, size_t count, double* values) {
  return guard([&] {
    need(client, "the client");
    need(field, "the field");
    need(variable, "the variable");
    need(lon, "lon");
    need(lat, "lat");
    need(values, "the place for the values");
    if (count == 0) return SCIENCE_OK;
    science::Field f(client->api, science::Json::object({{"name", field}}));
    std::vector<double> got = f.sample(variable, std::vector<double>(lon, lon + count), std::vector<double>(lat, lat + count), text(time));
    if (got.size() != count) throw science::JsonError("the service gave " + std::to_string(got.size()) + " values for " + std::to_string(count) + " places");
    std::copy(got.begin(), got.end(), values);
    return SCIENCE_OK;
  });
}

int science_series(science_client* client, const char* field, const char* variable, double lon, double lat, size_t* count,
                   double** times, double** values) {
  return guard([&] {
    need(client, "the client");
    need(field, "the field");
    need(variable, "the variable");
    need(count, "the place for the count");
    need(times, "the place for the times");
    need(values, "the place for the values");
    *count = 0;
    *times = *values = nullptr;
    science::Field f(client->api, science::Json::object({{"name", field}}));
    science::Series s = f.series(variable, lon, lat);
    double* t = allocate<double>(s.values.size());
    double* v = allocate<double>(s.values.size());
    for (std::size_t k = 0; k < s.values.size(); k++) {
      t[k] = k < s.times.size() ? science::unix_seconds(s.times[k]) : science::NaN;
      v[k] = s.values[k];
    }
    *count = s.values.size();
    *times = t;
    *values = v;
    return SCIENCE_OK;
  });
}

int science_grid_read(science_client* client, const char* field, const char* variable, const char* time, const double* bbox,
                      int stride, science_grid* grid) {
  if (grid) std::memset(grid, 0, sizeof *grid);
  return guard([&] {
    need(client, "the client");
    need(field, "the field");
    need(variable, "the variable");
    need(grid, "the grid");
    std::optional<science::Box> box;
    if (bbox) box = science::Box{bbox[0], bbox[1], bbox[2], bbox[3]};
    science::Field f(client->api, science::Json::object({{"name", field}}));
    science::Grid g = f.grid(variable, text(time), box, stride > 1 ? stride : 1);
    float* v = allocate<float>(g.values.size());
    std::copy(g.values.begin(), g.values.end(), v);
    grid->lon0 = g.lon0;
    grid->dlon = g.dlon;
    grid->lat0 = g.lat0;
    grid->dlat = g.dlat;
    grid->nlon = static_cast<size_t>(g.nlon);
    grid->nlat = static_cast<size_t>(g.nlat);
    grid->values = v;
    bounded_copy(grid->time, sizeof grid->time, g.time);
    bounded_copy(grid->units, sizeof grid->units, g.units);
    return SCIENCE_OK;
  });
}

void science_grid_free(science_grid* grid) {
  if (!grid) return;
  std::free(grid->values);
  grid->values = nullptr;
}

int science_define(science_client* client, const char* source, const char* title, const char* const* variables, size_t count) {
  return guard([&] {
    need(client, "the client");
    need(source, "the source");
    if (count) need(variables, "the variables");
    science::Json names = science::Json::array();
    for (std::size_t k = 0; k < count; k++) {
      need(variables[k], "a variable");
      names.push(variables[k]);
    }
    science::Json define = science::Json::object({{"title", title ? title : source}, {"variables", names}});
    client->api.post_observations(source, science::Json::array(), define);
    return SCIENCE_OK;
  });
}

int science_post_reading(science_client* client, const char* source, const char* station, double lon, double lat, const char* time,
                         const char* const* names, const double* values, size_t count) {
  return guard([&] {
    need(client, "the client");
    need(source, "the source");
    need(station, "the station");
    if (count) {
      need(names, "the names");
      need(values, "the values");
    }
    science::Reading r;
    r.station = station;
    r.lon = lon;
    r.lat = lat;
    r.time = text(time);
    for (std::size_t k = 0; k < count; k++) {
      need(names[k], "a name");
      r.values.emplace_back(names[k], values[k]);
    }
    science::Json answer = client->api.post_observations(source, std::vector<science::Reading>{r});
    if (answer["accepted"].as_number() >= 1) return SCIENCE_OK;
    std::string why = answer["rejected"][0]["reason"].text("the store accepted none");
    return fail(SCIENCE_ERROR_REFUSED, why);
  });
}

int science_post(science_client* client, const char* source, const char* body, char** answer) {
  return guard([&] {
    need(client, "the client");
    need(source, "the source");
    if (answer) *answer = nullptr;
    science::Json b = json_argument(body, "the body");
    science::Json a = client->api.post("/observations/" + science::http::encode(source), b);
    if (answer) *answer = copy(a.dump());
    return SCIENCE_OK;
  });
}

int science_stream(science_client* client, const char* source, const char* last_event_id, science_event_callback callback, void* user) {
  return guard([&] {
    need(client, "the client");
    if (!callback) throw Argument("the callback is NULL");
    science::StreamOptions options;
    options.last_event_id = text(last_event_id);
    options.idle = [&] { return callback(nullptr, nullptr, nullptr, user) != 0; };
    client->api.stream(text(source), [&](const science::Event& e) {
      std::string readings = e.readings.dump();
      return callback(e.id.c_str(), e.source.c_str(), readings.c_str(), user) != 0;
    }, options);
    return SCIENCE_OK;
  });
}

int science_compute(science_client* client, const char* unit, const char* function, const char* args, char** result) {
  return guard([&] {
    need(client, "the client");
    need(unit, "the unit");
    need(function, "the function");
    need(result, "the place for the result");
    *result = nullptr;
    science::Json a = json_argument(args ? args : "[]", "the arguments");
    if (!a.is_array()) throw Argument("the arguments are a JSON array, as [2460000.5]");
    *result = copy(client->api.compute(unit, function, a).dump());
    return SCIENCE_OK;
  });
}

int science_json_number(const char* json, const char* path, double* value) {
  return guard([&] {
    need(value, "the place for the number");
    science::Json j = json_argument(json, "the JSON text");
    const science::Json& v = at_path(j, path);
    if (!v.is_number() && !v.is_null()) throw Argument(std::string("at ") + text(path) + ": not a number");
    *value = v.as_number();
    return SCIENCE_OK;
  });
}

int science_json_numbers(const char* json, const char* path, double* values, size_t capacity, size_t* count) {
  return guard([&] {
    need(count, "the place for the count");
    if (capacity) need(values, "the place for the numbers");
    science::Json j = json_argument(json, "the JSON text");
    const science::Json& v = at_path(j, path);
    if (!v.is_array()) throw Argument(std::string("at ") + text(path) + ": not an array");
    for (std::size_t k = 0; k < v.size() && k < capacity; k++) {
      if (!v[k].is_number() && !v[k].is_null()) throw Argument(std::string("at ") + text(path) + "." + std::to_string(k) + ": not a number");
      values[k] = v[k].as_number();
    }
    *count = v.size();
    return SCIENCE_OK;
  });
}

int science_json_text(const char* json, const char* path, char** out) {
  return guard([&] {
    need(out, "the place for the text");
    *out = nullptr;
    science::Json j = json_argument(json, "the JSON text");
    const science::Json& v = at_path(j, path);
    *out = copy(v.is_string() ? v.as_string() : v.dump());
    return SCIENCE_OK;
  });
}

int science_json_length(const char* json, const char* path, size_t* length) {
  return guard([&] {
    need(length, "the place for the length");
    science::Json j = json_argument(json, "the JSON text");
    *length = at_path(j, path).size();
    return SCIENCE_OK;
  });
}

}  // extern "C"
