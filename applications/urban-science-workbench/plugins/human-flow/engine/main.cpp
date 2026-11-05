// The engine of the plug-in "Human flow": the urban-science domain's walkers as agents (urban::movement, on nomos's agents
// and tyche's decisions), on the district the page sends, calibrated to the City of Melbourne's counters each way.
//
// It answers requests, one JSON object a line ({"op": …, …}), with one JSON line each ({"ok": true, "result": …} or
// {"ok": false, "error": …}): with --serve until its input ends, as the workbench's server keeps it; otherwise the one
// request given as an argument. --data names the directory of the data library's products, the only files it reads.
// A run and the advice are jobs on the engine's worker thread: asked again, they answer how far they have come, then
// what they found.
//
// Every operation but status takes the district: {box, junctions, segments, buildings, sensors: [{id, at, way}]}
// (src/district.js's network and buildings; each sensor's place and nearest way in the area's local frame). A run is
// named by {weekday, from, to, seed?, regime?: {street | edges, kind: open | one-way | closed, entry?, from, to},
// event?: {street | edges, crowd, from, to, dwell}}.
//
//   status                               which products are present
//   prepare   {…district}                the walkway, counters, entrances, streets, the streams and their calibration
//   run       {…district, …run}          a run started (or found): how far it has come
//   frame     {…district, …run, hour}    the quarter hour's walkers every ten seconds for a minute, densities, queues
//   result    {…district, …run}          a run's totals, hours, counters (observed, expected, simulated), segments
//   compare   {…district, a, b, over, until}   two runs of one window compared, segment by segment
//   permits   {…district}                the City's event permits on the district's streets, and the counts on their dates
//   advice    {…district, …run, epochs, epochMinutes, simulations, truth}   POMCP's regime for the event's street
//   record    {…district, …run, mjd0}    a run written as a record of events, every one simulated

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "ops.hpp"

namespace {

using science::Json;

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "status") {
    Json p = Json::object();
    for (const char* f : {"street-sensors.json", "street-holidays-victoria.json", "street-terrain-melbourne.json", "city-events-melbourne.json"}) p.set(f, hf::product_exists(f));
    return p;
  }
  if (op == "prepare") return hf::op_prepare(r);
  if (op == "run") return hf::op_run(r);
  if (op == "frame") return hf::op_frame(r);
  if (op == "result") return hf::op_result(r);
  if (op == "compare") return hf::op_compare(r);
  if (op == "permits") return hf::op_permits(r);
  if (op == "advice") return hf::op_advice(r);
  if (op == "record") return hf::op_record(r);
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
    else if (std::strcmp(argv[i], "--data") == 0 && i + 1 < argc) hf::data_dir() = argv[++i];
    else request = argv[i];
  }
  if (serve) {
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      std::cout << respond(line) << "\n" << std::flush;
    }
    std::_Exit(0);   // the worker thread may be running a job
  }
  if (request.empty()) std::getline(std::cin, request);
  std::cout << respond(request) << "\n" << std::flush;
  std::_Exit(0);
}
