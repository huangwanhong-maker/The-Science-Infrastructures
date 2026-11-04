// Fabrica's engine: BodyParts3D's body read through the anatomy domain, measured once, and asked for
// sections, compounds and profiles (ADR 0066).
//
// It reads the data library's product (data/products/body_anatomy.py) from the directory passed with
// --data, welds each element's mesh, and measures it with the domain (anatomy::measures): topology by
// morphe, the solid by physis. It speaks JSON lines with --serve, a request a line and an answer a line:
//
//   {"op": "catalogue"}
//       the systems; for every element its concept, system, side, copy, measures and topology; for every
//       concept its parents in both relations and the elements below it; the pairs of parts that touch
//       (within 1.5 mm); the vertebral levels. Seconds to compute, so the server keeps it on disk.
//   {"op": "section", "plane": "transverse", "level": 1137, "convention": "radiological", "hidden": [..]}
//       every visible part cut by the plane at that level (mm along its axis in LPS), as loops in the
//       display's coordinates (mm, to 0.1), with each piece's area; for a transverse plane, its vertebral
//       level.
//   {"op": "part", "concept": "FMA7197", "relation": "partof"}
//       a compound's measures from its elements' (copies left out), its extents among them; when not all its
//       elements bound solids, its volume is null and the volume of those that do is given apart.
//   {"op": "profile", "element": 812, "axis": "long", "samples": 120}
//   {"op": "cases"} | {"op": "reconstruct", "case": …} | {"op": "slice", "case": …, "k": …}   (reconstruct.hpp)
// Every request may name its body: "body": "body-bp3d" (BodyParts3D's man, the default) or "body-hra-f" (the
// Human Reference Atlas's woman); each is read and measured when first asked for.
//       the part's cross-sectional area along its long principal axis (or "superior"), and the diameter of
//       the circle of that area: the calibre of a vessel along its course.
//
// Usage: fabrica --serve --data data/cache/products   |   fabrica --catalogue --data data/cache/products

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <map>
#include <numbers>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <anatomy/body/body.hpp>
#include <anatomy/frames/frames.hpp>
#include <anatomy/measures/measures.hpp>
#include <anatomy/parts/parts.hpp>
#include <anatomy/sections/sections.hpp>
#include <science/json.hpp>

#include "reconstruct.hpp"

namespace ab = anatomy::body;
namespace af = anatomy::frames;
namespace am = anatomy::measures;
namespace ap = anatomy::parts;
namespace as = anatomy::sections;
namespace cg = physis::computational_geometry;
using science::Json;

namespace {

constexpr double touch_distance = 1.5;  // mm

struct State {
  std::filesystem::path data;
  std::string name = "body-bp3d";  // the body's product: BodyParts3D's man, or the Atlas's woman (body-hra-f)
  bool loaded = false;
  ab::Body body;
  std::vector<cg::TriangleSurface> welded;  // per element
  std::vector<am::Measures> measures;
  std::vector<std::size_t> system;          // per element
  std::vector<cg::Box3> bounds;
  std::vector<as::Level> levels;
  fabrica::Reconstructions recon;  // the vessels reconstructed from images, by case (ADR 0070)
};

double r1(double x) { return std::round(x * 10) / 10; }
double r3(double x) { return std::round(x * 1000) / 1000; }
Json vec(const std::array<double, 3>& v, double (*round)(double) = r1) { return Json::array({round(v[0]), round(v[1]), round(v[2])}); }

void load(State& S) {
  if (S.loaded) return;
  S.body = ab::read_body(S.data, S.name);
  const std::size_t n = S.body.elements.size();
  S.welded.resize(n);
  S.measures.resize(n);
  S.system.resize(n);
  S.bounds.resize(n);
  std::atomic<std::size_t> next{0};
  const unsigned threads = std::max(1u, std::min(16u, std::thread::hardware_concurrency()));
  std::vector<std::thread> pool;
  for (unsigned t = 0; t < threads; t++)
    pool.emplace_back([&] {
      for (std::size_t e = next++; e < n; e = next++) {
        S.welded[e] = cg::welded(S.body.surface(e), am::weld_tolerance);
        S.measures[e] = am::measure(S.welded[e]);
        S.bounds[e] = S.measures[e].bounds;
      }
    });
  for (auto& t : pool) t.join();
  for (std::size_t e = 0; e < n; e++) S.system[e] = ap::system_of_element(S.body.atlas, e);
  for (const auto& v : as::vertebrae_for(S.body.atlas)) {
    const std::size_t c = S.body.atlas.find(v.fma);
    for (std::size_t e = 0; e < n; e++)
      if (S.body.elements[e].concept_ == c && S.body.elements[e].duplicate_of == ap::none) {
        S.levels.push_back({std::string(v.label), S.measures[e].centre[2]});
        break;
      }
  }
  S.loaded = true;
}

Json measures_json(const am::Measures& M) {
  Json j = Json::object();
  j.set("volume", M.solid() ? Json(r3(M.volume)) : Json());
  j.set("area", r1(M.area));
  j.set("centre", vec(M.centre));
  j.set("axes", Json::array({vec(M.axes[0], r3), vec(M.axes[1], r3), vec(M.axes[2], r3)}));
  j.set("extents", vec(M.extents));
  j.set("bounds", Json::array({r1(M.bounds.lower[0]), r1(M.bounds.lower[1]), r1(M.bounds.lower[2]), r1(M.bounds.upper[0]), r1(M.bounds.upper[1]), r1(M.bounds.upper[2])}));
  j.set("inward", M.inward);
  j.set("topology", Json::object({{"components", M.topology.components}, {"surfaces", M.topology.surfaces}, {"closed", M.topology.closed},
                                  {"consistent", M.topology.consistent}, {"euler", static_cast<double>(M.topology.euler)}, {"genus", M.topology.genus},
                                  {"boundaryCircles", M.topology.boundary_circles}, {"welded", M.topology.welded}, {"reason", M.topology.reason}}));
  return j;
}

Json catalogue(State& S) {
  load(S);
  const auto& B = S.body;
  const auto& A = B.atlas;
  const std::size_t n = B.elements.size();
  Json out = Json::object();
  out.set("frame", Json::object({{"name", "LPS"}, {"unit", "mm"}}));
  Json systems = Json::array();
  std::vector<std::size_t> counts(ap::SYSTEMS.size(), 0);
  for (std::size_t e = 0; e < n; e++)
    if (S.system[e] != ap::none && B.elements[e].duplicate_of == ap::none) counts[S.system[e]]++;
  for (std::size_t k = 0; k < ap::SYSTEMS.size(); k++)
    systems.push(Json::object({{"key", ap::SYSTEMS[k].key}, {"name", ap::SYSTEMS[k].name}, {"count", counts[k]}}));
  out.set("systems", systems);

  Json elements = Json::array();
  for (std::size_t e = 0; e < n; e++) {
    const auto& E = B.elements[e];
    const auto side = ap::side_by_name(E.name);
    Json j = Json::object();
    j.set("i", e);
    j.set("id", E.id);
    j.set("concept", A.concepts[E.concept_].id);
    j.set("name", E.name);
    j.set("conceptFromHeader", E.concept_from_header);
    j.set("system", S.system[e] == ap::none ? Json() : Json(ap::SYSTEMS[S.system[e]].key));
    j.set("side", side == ap::Side::left ? "left" : side == ap::Side::right ? "right" : "");
    j.set("copyOf", E.duplicate_of == ap::none ? Json() : Json(E.duplicate_of));
    j.set("published", std::isfinite(E.published_volume) ? Json(E.published_volume) : Json());
    j.set("level", as::level_at(S.levels, S.measures[e].centre[2]));
    j.set("measures", measures_json(S.measures[e]));
    elements.push(std::move(j));
  }
  out.set("elements", elements);

  // Concepts: those with elements below them in either relation, with their parents and elements.
  std::vector<std::vector<std::size_t>> below_isa(A.concepts.size()), below_partof(A.concepts.size());
  for (std::size_t e = 0; e < n; e++) {
    for (auto [rel, below] : {std::pair{ap::Relation::isa, &below_isa}, std::pair{ap::Relation::partof, &below_partof}}) {
      // Ancestors within the relation's tree only: a concept outside the part-of tree has no place in it.
      auto in_tree = [&](std::size_t x) { return rel == ap::Relation::isa ? A.concepts[x].in_isa : A.concepts[x].in_partof; };
      std::set<std::size_t> up;
      if (in_tree(B.elements[e].concept_))
        for (std::size_t a : A.forest(rel).ancestors(B.elements[e].concept_)) up.insert(a);
      for (std::size_t k : A.listed(rel)[e])
        for (std::size_t a : A.forest(rel).ancestors(k)) up.insert(a);
      for (std::size_t a : up) (*below)[a].push_back(e);
    }
  }
  Json concepts = Json::array();
  for (std::size_t c = 0; c < A.concepts.size(); c++) {
    if (below_isa[c].empty() && below_partof[c].empty()) continue;
    const auto sys = ap::system_of_concept(A, c);
    Json j = Json::object();
    j.set("id", A.concepts[c].id);
    j.set("name", A.concepts[c].name);
    j.set("isa", A.isa.parent[c] == ap::none ? Json() : Json(A.concepts[A.isa.parent[c]].id));
    j.set("partof", A.partof.parent[c] == ap::none ? Json() : Json(A.concepts[A.partof.parent[c]].id));
    j.set("system", sys == ap::none ? Json() : Json(ap::SYSTEMS[sys].key));
    j.set("inIsa", A.concepts[c].in_isa);
    j.set("inPartof", A.concepts[c].in_partof);
    j.set("elementsIsa", Json(below_isa[c]));
    j.set("elementsPartof", Json(below_partof[c]));
    concepts.push(std::move(j));
  }
  out.set("concepts", concepts);

  std::vector<const cg::TriangleSurface*> live(n, nullptr);
  for (std::size_t e = 0; e < n; e++)
    if (B.elements[e].duplicate_of == ap::none) live[e] = &S.welded[e];
  Json touching = Json::array();
  for (const auto& p : am::touching(live, touch_distance)) touching.push(Json::array({p.a, p.b, r3(p.distance)}));
  out.set("touching", touching);
  out.set("touchDistance", touch_distance);

  Json levels = Json::array();
  for (const auto& l : S.levels) levels.push(Json::object({{"label", l.label}, {"z", r1(l.height)}}));
  out.set("levels", levels);
  double top = -1e9, bottom = 1e9;
  for (const auto& M : S.measures) top = std::max(top, M.bounds.upper[2]), bottom = std::min(bottom, M.bounds.lower[2]);
  out.set("height", Json::object({{"top", r1(top)}, {"bottom", r1(bottom)}}));
  return out;
}

Json section(State& S, const Json& q) {
  load(S);
  const auto plane = af::plane_named(q["plane"].text("transverse"));
  const auto convention = q["convention"].text("radiological") == "neurological" ? af::Convention::neurological : af::Convention::radiological;
  const auto V = af::view(plane, convention);
  const double level = q["level"].as_number();
  const std::size_t n = S.body.elements.size();
  std::vector<char> cut(n, 1);
  for (std::size_t e = 0; e < n; e++)
    if (S.body.elements[e].duplicate_of != ap::none) cut[e] = 0;
  if (q["hidden"].is_array())
    for (std::size_t k = 0; k < q["hidden"].size(); k++) {
      const auto e = static_cast<std::size_t>(q["hidden"][k].as_number());
      if (e < n) cut[e] = 0;
    }
  std::vector<const cg::TriangleSurface*> surfaces(n);
  for (std::size_t e = 0; e < n; e++) surfaces[e] = &S.welded[e];
  const auto pieces = as::section(surfaces, S.bounds, V, level, cut);
  Json out = Json::object();
  out.set("plane", af::name(plane));
  out.set("convention", convention == af::Convention::radiological ? "radiological" : "neurological");
  out.set("level", level);
  out.set("right", af::name(V.right));
  out.set("up", af::name(V.up));
  out.set("seenFrom", af::name(af::seen_from(V)));
  out.set("u", vec(V.u));
  out.set("v", vec(V.v));
  out.set("w", vec(V.w));
  out.set("vertebra", plane == af::Plane::transverse ? as::level_at(S.levels, level) : "");
  Json list = Json::array();
  for (const auto& p : pieces) {
    auto flat = [](const std::vector<cg::Point2>& l) {
      Json a = Json::array();
      for (const auto& q2 : l) a.push(r1(q2[0])), a.push(r1(q2[1]));
      return a;
    };
    Json loops = Json::array(), chains = Json::array();
    for (const auto& l : p.loops) loops.push(flat(l));
    for (const auto& c : p.chains) chains.push(flat(c));
    list.push(Json::object({{"e", p.part}, {"area", r1(p.area)}, {"loops", loops}, {"chains", chains}}));
  }
  out.set("pieces", list);
  return out;
}

Json part(State& S, const Json& q) {
  load(S);
  const auto& A = S.body.atlas;
  const std::size_t c = A.find(q["concept"].text());
  if (c == ap::none) throw std::invalid_argument("no concept " + q["concept"].text());
  const auto rel = q["relation"].text("partof") == "isa" ? ap::Relation::isa : ap::Relation::partof;
  std::vector<const am::Measures*> M;
  std::vector<const cg::TriangleSurface*> W;
  Json elements = Json::array();
  double published = 0, solid_volume = 0;
  std::size_t solids = 0;
  for (std::size_t e : A.elements_of(c, rel)) {
    if (S.body.elements[e].duplicate_of != ap::none) continue;
    if (S.measures[e].solid()) solids++, solid_volume += S.measures[e].volume;
    M.push_back(&S.measures[e]);
    W.push_back(&S.welded[e]);
    elements.push(e);
    if (std::isfinite(S.body.elements[e].published_volume)) published += S.body.elements[e].published_volume;
  }
  Json out = Json::object();
  out.set("concept", A.concepts[c].id);
  out.set("name", A.concepts[c].name);
  out.set("elements", elements);
  out.set("published", published);
  // The volume of the parts that bound solids, when not all do (open vessels among them): a lower bound.
  out.set("solids", solids);
  out.set("solidVolume", r3(solid_volume));
  if (!M.empty()) {
    const auto C = am::combine(M, W);
    out.set("measures", measures_json(C));
    out.set("level", as::level_at(S.levels, C.centre[2]));
  }
  return out;
}

Json profile(State& S, const Json& q) {
  load(S);
  const auto e = static_cast<std::size_t>(q["element"].as_number());
  if (e >= S.body.elements.size()) throw std::invalid_argument("no such element");
  const auto& M = S.measures[e];
  const auto& W = S.welded[e];
  const bool along_body = q["axis"].text("long") == "superior";
  const cg::Point3 N = along_body ? cg::Point3{0, 0, 1} : M.axes[0];
  double lo = 1e300, hi = -1e300;
  for (const auto& p : W.points) {
    const double h = N[0] * p[0] + N[1] * p[1] + N[2] * p[2];
    lo = std::min(lo, h), hi = std::max(hi, h);
  }
  const int samples = std::clamp(static_cast<int>(q["samples"].is_number() ? q["samples"].as_number() : 120), 8, 1000);
  const double dh = (hi - lo) / samples;
  Json points = Json::array();
  double integral = 0;
  for (int i = 0; i < samples; i++) {
    const double h = lo + (i + 0.5) * dh;
    const cg::Plane3 P{{h * N[0], h * N[1], h * N[2]}, N};
    const auto chains = cg::plane_section(W, P);
    std::size_t loops = 0;
    for (const auto& ch : chains) loops += ch.closed;
    const double a = cg::section_area(chains, P);
    integral += a * dh;
    points.push(Json::object({{"s", r1(h - lo)}, {"area", r1(a)}, {"diameter", r3(a > 0 ? 2 * std::sqrt(a / std::numbers::pi) : 0)}, {"loops", loops}}));
  }
  Json out = Json::object();
  out.set("element", e);
  out.set("axis", along_body ? "superior" : "long");
  out.set("direction", vec(N, r3));
  out.set("length", r1(hi - lo));
  out.set("samples", points);
  out.set("integral", r1(integral));  // Cavalieri: ≈ the volume when the mesh bounds a solid
  out.set("volume", M.solid() ? Json(r1(M.volume)) : Json());
  return out;
}

// The bodies the engine knows, by their products' names.
const std::set<std::string> BODIES{"body-bp3d", "body-hra-f"};

Json answer(State& S, const Json& q) {
  const std::string op = q["op"].text();
  if (op == "catalogue") return catalogue(S);
  if (op == "section") return section(S, q);
  if (op == "part") return part(S, q);
  if (op == "profile") return profile(S, q);
  if (op == "ping") return Json::object({{"loaded", S.loaded}});
  if (op == "cases") return fabrica::list_cases(S.recon);
  if (op == "reconstruct") return fabrica::ensure(S.recon, q["case"].text()).result;
  if (op == "slice") return fabrica::slice(S.recon, q);
  throw std::invalid_argument("no operation " + op);
}

}  // namespace

int main(int argc, char** argv) {
  std::map<std::string, State> bodies;  // each body loaded when first asked for
  State& S = bodies["body-bp3d"];
  S.data = "data/cache/products";
  bool serve = false, once = false;
  for (int i = 1; i < argc; i++) {
    const std::string a = argv[i];
    if (a == "--serve") serve = true;
    else if (a == "--catalogue") once = true;
    else if (a == "--data" && i + 1 < argc) S.data = argv[++i];
  }
  S.recon.products = S.data;
  std::ios::sync_with_stdio(false);
  if (once) {
    try {
      std::cout << catalogue(S).dump() << "\n";
    } catch (const std::exception& e) {
      std::cerr << e.what() << "\n";
      return 1;
    }
    return 0;
  }
  if (!serve) {
    std::cerr << "usage: fabrica --serve --data <products>  |  fabrica --catalogue --data <products>\n";
    return 2;
  }
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.empty()) continue;
    Json reply;
    try {
      const Json q = Json::parse(line);
      const std::string body = q["body"].text("body-bp3d");
      if (!BODIES.count(body)) throw std::invalid_argument("no body " + body);
      State& B = bodies[body];
      if (B.data.empty()) B.data = S.data, B.name = body;
      // the reconstructions are the same whatever body is shown
      reply = Json::object({{"ok", true}, {"result", answer(q["op"].text() == "cases" || q["op"].text() == "reconstruct" || q["op"].text() == "slice" ? S : B, q)}});
    } catch (const std::exception& e) {
      reply = Json::object({{"ok", false}, {"error", e.what()}});
    }
    std::cout << reply.dump() << "\n" << std::flush;
  }
  return 0;
}
