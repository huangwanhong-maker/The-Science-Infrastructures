// A sensor that posts its readings in a loop, from C++: it defines its source
// once (its variables, with the units it measures in, which the store converts
// to SI on arrival, and its station), then posts a reading every few seconds,
// and keeps going while the network or the service is away. Run cpp-stream in a
// second terminal to watch the readings arrive.
//
// The readings are made up (a daily cycle and some noise), so that the example
// runs anywhere; a real program posts what its device measures, the same way.
// The C example sensor.c is the template for a device running Linux.
//
//   api/build/cpp-sensor [source] [seconds between readings] [count, 0 for ever] [url]
//   (a service with a token needs it: $SCIENCE_API_TOKEN)

#include "science.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <thread>

int main(int argc, char** argv) {
  const std::string source = argc > 1 ? argv[1] : "example-garden";
  const double every = argc > 2 ? std::atof(argv[2]) : 5.0;
  const int count = argc > 3 ? std::atoi(argv[3]) : 0;
  science::Client api(argc > 4 ? argv[4] : "");

  using science::Json;
  const Json define = Json::object({
      {"title", "Example garden station (C++)"},
      {"variables", Json::array({"air_temperature:degree_C", "relative_humidity:%"})},
      {"stations", Json::array({Json::object({{"id", "garden-1"}, {"name", "The garden"}, {"lon", 139.69}, {"lat", 35.69}, {"elevation", 40}})})},
  });
  std::mt19937 random(std::random_device{}());
  std::normal_distribution<double> noise(0.0, 1.0);
  bool defined = false;

  for (int sent = 0; count == 0 || sent < count; sent++) {
    // A daily cycle peaking in the afternoon, in local solar time at 139.69° E.
    auto now = std::chrono::system_clock::now();
    double hours = std::fmod(std::chrono::duration<double>(now.time_since_epoch()).count() / 3600.0 + 139.69 / 15.0, 24.0);
    double phase = (hours - 9.0) / 24.0 * 2.0 * std::acos(-1.0);
    double celsius = 17.0 + 5.0 * std::sin(phase) + 0.2 * noise(random);
    double humidity = 65.0 - 15.0 * std::sin(phase) + noise(random);

    science::Reading r;
    r.station = "garden-1";  // described in define: no place is needed
    r.time = science::utc_text(now);
    r.values = {{"air_temperature", std::round(celsius * 100) / 100}, {"relative_humidity", std::round(humidity * 10) / 10}};
    try {
      Json answer = api.post_observations(source, std::vector<science::Reading>{r}, defined ? Json() : define);
      defined = true;
      std::printf("%s  %5.2f °C  %4.1f %%   accepted %g", r.time.c_str(), celsius, humidity, answer["accepted"].as_number());
      if (answer["rejected"].size()) std::printf(", refused: %s", answer["rejected"][0]["reason"].text().c_str());
      std::printf("\n");
    } catch (const science::ApiError& e) {
      if (e.status == 401 || e.status == 403) {  // no token, or the wrong one: trying again will not help
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
      }
      std::fprintf(stderr, "%s  not sent (%s); trying again\n", r.time.c_str(), e.what());
    }
    std::fflush(stdout);
    std::this_thread::sleep_for(std::chrono::duration<double>(every));
  }
  return 0;
}
