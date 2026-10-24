// The Genome Workbench's engine: the genomics and evolution domains' C++, behind the JSON-lines
// protocol of protocol.hpp. Built by engine.py (g++ 11 or later, -O3 -march=native) into
// data/cache/engines/genome/.
//
//   genome --serve        answer requests line by line until the input ends (a worker)
//   genome                answer the one request on standard input
//
// The operations are in ops_*.hpp, each compiled in a unit of its own (ops_*.cpp) and registered
// at static initialisation in the table of protocol.hpp; this unit is the loop that answers.
//
// Every answer carries the seconds it took; an operation's fault is an answer with "ok": false,
// never a crash, so that a worker survives a bad request.

#include "protocol.hpp"
#include "genomes.hpp"

#include <chrono>
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>

namespace {

workbench::Json answer(const std::string& line, workbench::Worker& worker) {
  using workbench::Json;
  Json id;
  const auto t0 = std::chrono::steady_clock::now();
  try {
    const Json request = Json::parse(line);
    id = request["id"];
    const std::string op = request["op"].text();
    auto& ops = workbench::operations();
    auto it = ops.find(op);
    if (it == ops.end()) throw workbench::Fault("no operation \"" + op + "\"");
    Json result = it->second(request, worker);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return Json::object({{"id", id}, {"ok", true}, {"seconds", std::round(seconds * 1e6) / 1e6}, {"result", std::move(result)}});
  } catch (const std::exception& e) {
    return Json::object({{"id", id}, {"ok", false}, {"error", e.what()}});
  }
}

}  // namespace

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);
  bool serve = false;
  for (int i = 1; i < argc; i++) serve = serve || std::strcmp(argv[i], "--serve") == 0;
  workbench::Worker worker;
  if (serve) {
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line.empty() || line == "\r") continue;
      std::cout << answer(line, worker).dump() << '\n' << std::flush;
    }
    return 0;
  }
  std::ostringstream all;
  all << std::cin.rdbuf();
  std::cout << answer(all.str(), worker).dump() << '\n';
  return 0;
}
