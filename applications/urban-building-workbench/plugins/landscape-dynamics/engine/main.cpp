// The engine of the plug-in "Landscape Dynamics": the urban-science domain's landscape module, on the district
// and the building the page sends.
//
// It answers requests, one JSON object a line ({"op": …, …}), with one JSON line each
// ({"ok": true, "result": …} or {"ok": false, "error": …}): with --serve until its input ends, as the
// workbench's server keeps it; otherwise the one request given as an argument. It reads no files.
//
//   district  {junctions, segments, stores}            the district as cells: geometry, homology, the stores, the lane's plan
//   pairs     {…the district, within?, least?, limit?} of the stores sent, the pairs a lane could join, nearest first
//   lane      {…the district}                          the lane built between the first two stores and closed again
//   run       {…the district, params, entries, altitude} the landscape running for a day: events, frames, strata, order
//   skeleton  {alphaA, …, grant, …, gain, seed, hours} the landscape's events composed exactly: law, irreversibility,
//                                                       what coarser descriptions forget, interventions, evidence
//   building  {x, y, heading, width, depth, storey, nx, ny, floors, altitude, azimuth, airByHour, events}
//                                                       a building's rooms, their heat through a day, its events
//   scales    {…the district}                          the district's Hodge spectra and harmonic flows; memory; operators

#include <cstring>
#include <iostream>
#include <string>

#include "ops.hpp"

namespace {

using science::Json;

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "district") return ld::op_district(r);
  if (op == "pairs") return ld::op_pairs(r);
  if (op == "lane") return ld::op_lane(r);
  if (op == "run") return ld::op_run(r);
  if (op == "skeleton") return ld::op_skeleton(r);
  if (op == "building") return ld::op_building(r);
  if (op == "scales") return ld::op_scales(r);
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
    else if (std::strcmp(argv[i], "--data") == 0 && i + 1 < argc) ++i;  // the workbench names the data; this engine reads none
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
