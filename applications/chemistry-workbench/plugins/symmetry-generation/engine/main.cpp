// The engine of the plug-in "symmetry and generation": the manuscript's two formalisms, executed.
//
// It answers requests, one JSON object a line ({"op": …, …}), with one JSON line each
// ({"ok": true, "result": …} or {"ok": false, "error": …}): with --serve until its input ends, as
// the workbench's server keeps it; otherwise the one request given as an argument or on stdin.
//
//   catalogue                       the example molecules, their point groups and character tables
//   grammar   {example}             ethanol | glycine | stereocentre | proton: a grammar's derivation space
//   complexes                       the six-cycle, its filled complex, the two-vertex complex
//   field     {group, psi?}         the tetrahedral interaction over 1, Z2, Z3 or S3
//   wick      {n, group}            its Feynman graphs at order n (1–3) with their 2-complexes
//   frank     {k, ki, k0, …}        Frank's autocatalysis: mass action and the Markov chain
//   pitchfork {b}                   the cubic normal form's equilibria
//   landau    {mu, lambda}          the action along the field's modes
//   checks                          every checked claim of the manuscript

#include <iostream>
#include <string>

#include "checks.hpp"

namespace {

using science::Json;

double number(const Json& r, const char* key, double fallback) { return r[key].is_number() ? r[key].as_number() : fallback; }

Json answer(const Json& r) {
  const std::string op = r["op"].text();
  if (op == "catalogue") {
    Json list = Json::array();
    for (const auto& m : chem::catalogue()) {
      Json j = chem::molecule_json(m);
      j.set("pointGroup", chem::point_group_json(m));
      list.push(j);
    }
    return list;
  }
  if (op == "grammar") {
    const std::string ex = r["example"].text("stereocentre");
    if (ex == "ethanol") return chem::ethanol();
    if (ex == "glycine") return chem::glycine(static_cast<std::size_t>(number(r, "units", 3)));
    if (ex == "stereocentre") return chem::stereocentre();
    if (ex == "proton") return chem::proton_transfer();
    throw std::invalid_argument("no grammar " + ex);
  }
  if (op == "complexes") return chem::complexes();
  if (op == "field") return chem::field(r["group"].text("Z2"), r["psi"]);
  if (op == "wick") return chem::wick(static_cast<std::size_t>(number(r, "n", 2)), r["group"].text("Z2"));
  if (op == "frank") {
    chem::FrankParameters P;
    P.k = number(r, "k", P.k), P.ki = number(r, "ki", P.ki), P.k0 = number(r, "k0", P.k0);
    P.a0 = number(r, "a0", P.a0), P.l0 = number(r, "l0", P.l0), P.r0 = number(r, "r0", P.r0), P.t_end = number(r, "t", P.t_end);
    P.precursor = static_cast<long long>(number(r, "precursor", static_cast<double>(P.precursor)));
    P.seed_l = static_cast<long long>(number(r, "seedL", 1)), P.seed_r = static_cast<long long>(number(r, "seedR", 1));
    if (P.precursor < 1 || P.precursor > 80) throw std::invalid_argument("the precursor count is 1 to 80");
    return chem::frank(P);
  }
  if (op == "pitchfork") return chem::pitchfork(number(r, "b", 1));
  if (op == "landau") return chem::landau(number(r, "mu", -1), number(r, "lambda", 0.01));
  if (op == "checks") return chem::checks();
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
  const std::string mode = argc > 1 ? argv[1] : "";
  if (mode == "--serve") {
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      std::cout << respond(line) << "\n";
      std::cout.flush();
    }
    return 0;
  }
  std::string request = mode;
  if (request.empty() || request == "-") {
    std::string line;
    request.clear();
    while (std::getline(std::cin, line)) request += line;
  }
  std::cout << respond(request) << "\n";
  return 0;
}
