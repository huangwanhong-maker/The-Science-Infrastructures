// The engine of the plug-in "Street Oceanography": the urban-science domain on the data library's records.
//
// It answers requests, one JSON object a line ({"op": …, …}), with one JSON line each
// ({"ok": true, "result": …} or {"ok": false, "error": …}): with --serve until its input ends, as the
// workbench's server keeps it; otherwise the one request given as an argument. --data names the directory of
// the data library's products, the only files it reads.
//
//   status                                  which products are present
//   sensors | streets                       Melbourne's sensors; OpenStreetMap's named streets
//   street      {street, walk?, place?}     a named street: holdings, mixing, terrain length, hours, ascent and descent
//   visits      {street, types}             the visits of declared kinds of people (a representation)
//   terrain     {sensor, year?, walk?}      a sensor's surroundings in CLUE: uses, mixing, terrain length, their record
//   dominance   {year?}                     every sensor's terrain: its dominance order and the rankings of weightings
//   establishments {year?}                  every establishment about the sensors, as offsets from its sensor
//   profiles    {era?}                      every outdoor sensor's mean count at each civil hour of the week
//   affinity    {where, cell?, radius?}     co-location of uses (Melbourne, or a named street) against two nulls
//   spacing     {j11…, du, dv, street?, use?}   predicted (Turing) and observed (Lomb–Scargle) spacing
//   rhythm      {sensor, era?, nw?, order?} a record's spectrum, lines, harmonic constants, selectivity, profile
//   clock       {sensor}                    the civil share of its clock at each change of daylight saving
//   bicoherence {sensor, era?}              the bicoherence of what the calendars leave
//   constructed                             the paper's Figure 1, and the case it misreads
//   response    {sensor, era?}              susceptibility to rain, heat and holidays; holidays paired with their neighbours
//   summary                                 every outdoor sensor: selectivity, compression, rain, holidays
//   plane                                   the classification of the sensors in the plane
//   pathway     {…}                         the bridal store and the ring store (a representation)
//   streetday   {…}                         the crowd on a street through a day, with and without a rule opening its carriageway (a representation)
//   district    {…}                         movement in a district: circulation and transit (a representation)
//   persontime  {stays, window}             Little's identity on given stays
//   tides                                   every sensor's tide: high and low water, range, form, its mean week
//   sea         {date}                      a civil date: each sensor's predicted level, observed counts and counts each way
//   currents                                every sensor that counts each way: its current meter
//   streetwave  {hour, blockLengths, …}     the crowd on a street second by second about an hour (a representation)
//   plume       {length, bridal, ring, …}   a visitor's law along its street after leaving a store (a representation)

#include <cstring>
#include <iostream>
#include <string>

#include "ops_models.hpp"
#include "ops_rhythm.hpp"
#include "ops_sea.hpp"
#include "ops_street.hpp"

namespace {

using science::Json;

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "status") return so::op_status();
  if (op == "sensors") return so::op_sensors();
  if (op == "streets") return so::op_streets();
  if (op == "street") return so::op_street(r);
  if (op == "visits") return so::op_visits(r);
  if (op == "terrain") return so::op_terrain(r);
  if (op == "establishments") return so::op_establishments(r);
  if (op == "profiles") return so::op_profiles(r);
  if (op == "dominance") return so::op_dominance(r);
  if (op == "affinity") return so::op_affinity(r);
  if (op == "spacing") return so::op_spacing(r);
  if (op == "rhythm") return so::op_rhythm(r);
  if (op == "clock") return so::op_clock(r);
  if (op == "bicoherence") return so::op_bicoherence(r);
  if (op == "constructed") return so::op_constructed();
  if (op == "response") return so::op_response(r);
  if (op == "summary") return so::op_sensors_summary();
  if (op == "plane") return so::op_plane();
  if (op == "pathway") return so::op_pathway(r);
  if (op == "streetday") return so::op_streetday(r);
  if (op == "district") return so::op_district(r);
  if (op == "persontime") return so::op_person_time(r);
  if (op == "tides") return so::op_tides();
  if (op == "sea") return so::op_sea(r);
  if (op == "currents") return so::op_currents();
  if (op == "streetwave") return so::op_streetwave(r);
  if (op == "plume") return so::op_plume(r);
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
    else if (std::strcmp(argv[i], "--data") == 0 && i + 1 < argc) so::data_dir() = argv[++i];
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
