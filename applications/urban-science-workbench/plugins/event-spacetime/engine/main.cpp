// The engine of the plug-in "Event Spacetime": the urban-science domain's events module on the district the page
// sends and the record of its events the data library compiles (ADR 0067).
//
// It answers requests, one JSON object a line ({"op": …, …}), with one JSON line each
// ({"ok": true, "result": …} or {"ok": false, "error": …}): with --serve until its input ends, as the
// workbench's server keeps it; otherwise the one request given as an argument. --data names the directory of the
// data library's products, of which it reads city-events-<area>.json.
//
//   record  {area, junctions, segments, buildings[, record]}
//                                                    the district's volumes and interfaces, its record of events (the
//                                                    one given, else the area's product) and
//                                                    its event complex: couplings, interfaces and events for drawing,
//                                                    the topology of K and K̂, what was left out and why
//   cut     {t[, …the record's request]}             the spin network of the moment at t (MJD): of the record the request
//                                                    names (built again only if it differs), else of the one last built
//   status  {}                                       the engine answers

#include <cstring>
#include <iostream>
#include <string>

#include "ops.hpp"

namespace {

using science::Json;

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "record") return es::op_record(r);
  if (op == "cut") return es::op_cut(r);
  if (op == "status") return Json::object({{"ok", true}});
  throw std::invalid_argument("no operation \"" + op + "\"");
}

std::string respond(const std::string& line) {
  Json reply = Json::object();
  try {
    const Json r = Json::parse(line);
    if (r.has("id")) reply.set("id", r["id"]);
    reply.set("ok", true);
    reply.set("result", answer(r));
  } catch (const std::exception& e) {
    reply = Json::object({{"ok", false}, {"error", e.what()}});
  }
  return reply.dump();
}

}  // namespace

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);
  bool serve = false;
  std::string request;
  for (int i = 1; i < argc; i++) {
    if (std::strcmp(argv[i], "--serve") == 0) serve = true;
    else if (std::strcmp(argv[i], "--data") == 0 && i + 1 < argc) es::data_dir() = argv[++i];
    else request = argv[i];
  }
  if (serve) {
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      std::cout << respond(line) << "\n" << std::flush;
    }
    return 0;
  }
  if (request.empty()) std::getline(std::cin, request);
  std::cout << respond(request) << "\n";
  return 0;
}
