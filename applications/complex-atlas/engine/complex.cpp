// The engine of the Complex Atlas: one request in JSON on standard input, one answer in JSON
// on standard output (the application's README, §4). It composes the foundations: morphe's
// topology for complexes and their invariants, physis's gauge theory and quantum gravity for
// the fields and state sums on them. Errors are answered as {"error": "..."} with exit
// status 1.
//
//   {"op": "catalogue"}                                   the families of complexes
//   {"op": "build", "family": id, "params": {...}}        a complex of a family
//   {"op": "analyse", "complex": {...}}                   its invariants
//   {"op": "gauge", "complex": {...}, "group": {...}, ...} a lattice gauge field on it, by Monte Carlo
//   {"op": "turaev_viro", "complex": {...}, "r": 4}       the Turaev–Viro invariant of a triangulated 3-manifold

#include <exception>
#include <iostream>
#include <iterator>
#include <string>

#include "complexes.hpp"
#include "fields.hpp"
#include "generation.hpp"
#include "theories.hpp"

int main() {
  using science::Json;
  std::ios::sync_with_stdio(false);
  const std::string text((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
  try {
    const Json request = Json::parse(text);
    const std::string op = request["op"].text();
    Json answer;
    if (op == "catalogue") answer = atlas::catalogue();
    else if (op == "build") answer = Json::object({{"complex", atlas::to_json(atlas::build(request["family"].text(), request["params"]))}});
    else if (op == "analyse") answer = Json::object({{"analysis", atlas::analyse(atlas::from_json(request["complex"]))}});
    else if (op == "gauge") answer = atlas::gauge(atlas::from_json(request["complex"]), request);
    else if (op == "turaev_viro") answer = atlas::turaev_viro(atlas::from_json(request["complex"]), request);
    else if (op == "theories") answer = atlas::theories();
    else if (op == "expand") answer = atlas::expand(request);
    else if (op == "large_n") answer = atlas::large_n(request);
    else if (op == "graph_complex") answer = atlas::graph_complex(request);
    else if (op == "condensate") answer = atlas::condensate(request);
    else if (op == "dimensions") answer = atlas::dimensions(atlas::from_json(request["complex"]), request);
    else if (op == "random_geometry") answer = atlas::random_geometry(request);
    else if (op == "spin_foam") answer = atlas::spin_foam(atlas::from_json(request["complex"]), request);
    else throw std::invalid_argument("no operation \"" + op + "\"");
    std::cout << answer.dump() << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cout << Json::object({{"error", e.what()}}).dump() << "\n";
    return 1;
  }
}
