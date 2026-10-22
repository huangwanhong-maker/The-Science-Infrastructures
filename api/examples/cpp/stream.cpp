// Following the stream of readings from C++: each reading posted to the
// service (by cpp-sensor, c-sensor, a device, or the Python and JavaScript
// examples) is printed as it arrives, within milliseconds on one computer.
// The stream is Server-Sent Events: after a dropped connection the client
// reconnects by itself and receives, by the last event's id, what it missed.
//
//   api/build/cpp-stream [source, or "" for every source] [count, 0 for ever] [url]

#include "science.hpp"

#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
  const std::string source = argc > 1 ? argv[1] : "";
  const int count = argc > 2 ? std::atoi(argv[2]) : 0;
  science::Client api(argc > 3 ? argv[3] : "");
  std::printf("Listening to %s at %s\n", source.empty() ? "every source" : source.c_str(), api.url().c_str());
  int seen = 0;
  try {
    api.stream(source, [&](const science::Event& e) {
      for (const auto& r : e.readings.as_array()) {
        std::printf("[%s] %s %s at %s:", e.id.c_str(), e.source.c_str(), r["station"].text().c_str(), r["time"].text().c_str());
        for (const auto& [name, value] : r["values"].as_object()) std::printf("  %s %g", name.c_str(), value.as_number());
        std::printf("\n");  // values in SI units, as the store keeps them
      }
      std::fflush(stdout);
      return count == 0 || ++seen < count;  // false ends the stream
    });
  } catch (const std::exception& e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
  return 0;
}
