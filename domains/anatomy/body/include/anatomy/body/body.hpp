// A body: BodyParts3D's adult man read from the data library's product (data/products/body_anatomy.py), its
// concepts and relations into an Atlas, its elements' meshes as triangle surfaces in the LPS frame, in
// millimetres (ADR 0066).
//
// The product's tables are tab-separated with a header line: -concepts.tsv (concept, name, in_isa,
// in_partof), -isa.tsv and -partof.tsv (parent, child), -compounds-isa.tsv and -compounds-partof.tsv
// (concept, element: BodyParts3D's own lists of the elements composing each compound), and -elements.tsv
// (one line per element: its concept, name, pack and byte offsets, counts, quantisation box, published
// volume, and the element it copies, if any). The packs hold each element's positions as uint16 triples over
// its box, normals as int8 quadruples and triangles as uint16 or uint32 indices, each block 4-aligned.
// Only the directory passed in is read (CLAUDE.md §3.3). Checks: tests/check_body.cpp, on the product when
// it has been built.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <anatomy/parts/parts.hpp>
#include <physis/computational_geometry/surfaces.hpp>

namespace anatomy::body {

using physis::computational_geometry::Point3;
using physis::computational_geometry::TriangleSurface;

struct Element {
  std::string id;          // BodyParts3D's file id, "FJ1252"; "M" marks a part mirrored to the other side
  std::size_t concept_ = parts::none;
  std::string name;
  bool concept_from_header = true;  // false: taken from the most specific compound containing it
  std::size_t pack = 0, position_offset = 0, normal_offset = 0, index_offset = 0;
  std::size_t vertices = 0, triangles = 0, index_bytes = 2;
  Point3 lower{}, upper{};  // the quantisation box, mm
  double published_volume = std::numeric_limits<double>::quiet_NaN();  // cm³, as the database computed it
  std::size_t duplicate_of = parts::none;  // the earlier element this one copies
};

struct Body {
  parts::Atlas atlas;
  std::vector<Element> elements;
  std::vector<std::pair<std::size_t, std::size_t>> compounds_isa, compounds_partof;  // (concept, element), as given
  std::vector<std::vector<unsigned char>> packs;

  // The element's mesh: points in LPS millimetres, triangles in the files' order.
  TriangleSurface surface(std::size_t e) const {
    const Element& E = elements.at(e);
    const auto& buf = packs.at(E.pack);
    const std::size_t need = E.index_offset + 3 * E.triangles * E.index_bytes;
    if (need > buf.size() || E.position_offset + 6 * E.vertices > buf.size()) throw std::out_of_range("anatomy::body::Body::surface: an element beyond its pack");
    TriangleSurface S;
    S.points.resize(E.vertices);
    for (std::size_t v = 0; v < E.vertices; v++)
      for (int k = 0; k < 3; k++) {
        std::uint16_t q;
        std::memcpy(&q, buf.data() + E.position_offset + 6 * v + 2 * k, 2);
        S.points[v][k] = E.lower[k] + (E.upper[k] - E.lower[k] > 0 ? (E.upper[k] - E.lower[k]) : 1e-3) * (q / 65535.0);
      }
    S.triangles.resize(E.triangles);
    for (std::size_t t = 0; t < E.triangles; t++)
      for (int k = 0; k < 3; k++) {
        std::size_t i;
        if (E.index_bytes == 2) {
          std::uint16_t x;
          std::memcpy(&x, buf.data() + E.index_offset + 2 * (3 * t + k), 2);
          i = x;
        } else {
          std::uint32_t x;
          std::memcpy(&x, buf.data() + E.index_offset + 4 * (3 * t + k), 4);
          i = x;
        }
        if (i >= E.vertices) throw std::out_of_range("anatomy::body::Body::surface: a corner beyond the element's points");
        S.triangles[t][k] = i;
      }
    return S;
  }
};

namespace detail {
inline std::vector<std::vector<std::string>> read_tsv(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  if (!in) throw std::runtime_error("anatomy::body: cannot read " + p.string());
  std::vector<std::vector<std::string>> rows;
  std::string line;
  bool header = true;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (header) { header = false; continue; }
    if (line.empty()) continue;
    std::vector<std::string> cells;
    std::size_t at = 0;
    while (true) {
      const std::size_t tab = line.find('\t', at);
      cells.push_back(line.substr(at, tab == std::string::npos ? std::string::npos : tab - at));
      if (tab == std::string::npos) break;
      at = tab + 1;
    }
    rows.push_back(std::move(cells));
  }
  return rows;
}
inline std::vector<unsigned char> read_bytes(const std::filesystem::path& p) {
  std::ifstream in(p, std::ios::binary);
  if (!in) throw std::runtime_error("anatomy::body: cannot read " + p.string());
  return std::vector<unsigned char>(std::istreambuf_iterator<char>(in), {});
}
inline double number(const std::string& s) { return s.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(s); }
}  // namespace detail

// The body of the product `name` in `directory` (data/cache/products), with its packs.
inline Body read_body(const std::filesystem::path& directory, const std::string& name = "body-bp3d") {
  Body B;
  auto file = [&](const std::string& suffix) { return directory / (name + suffix); };
  for (const auto& r : detail::read_tsv(file("-concepts.tsv"))) {
    if (r.size() < 4) throw std::runtime_error("anatomy::body: a short line in the concepts");
    B.atlas.add_concept({r[0], r[1], r[2] == "1", r[3] == "1"});
  }
  auto concept_of = [&](const std::string& id) {
    const std::size_t c = B.atlas.find(id);
    return c != parts::none ? c : B.atlas.add_concept({id, id, false, false});
  };
  for (auto [suffix, relation] : {std::pair{"-isa.tsv", parts::Relation::isa}, std::pair{"-partof.tsv", parts::Relation::partof}})
    for (const auto& r : detail::read_tsv(file(suffix))) B.atlas.forest(relation).relate(concept_of(r[1]), concept_of(r[0]));
  std::unordered_map<std::string, std::size_t> element_index;
  std::size_t packs = 0;
  for (const auto& r : detail::read_tsv(file("-elements.tsv"))) {
    if (r.size() < 20) throw std::runtime_error("anatomy::body: a short line in the elements");
    Element E;
    E.id = r[1];
    E.concept_ = concept_of(r[2]);
    E.name = r[3];
    E.concept_from_header = r[4] == "header";
    E.pack = std::stoul(r[5]);
    E.position_offset = std::stoul(r[6]);
    E.normal_offset = std::stoul(r[7]);
    E.index_offset = std::stoul(r[8]);
    E.vertices = std::stoul(r[9]);
    E.triangles = std::stoul(r[10]);
    E.index_bytes = std::stoul(r[11]);
    for (int k = 0; k < 3; k++) E.lower[k] = std::stod(r[12 + k]), E.upper[k] = std::stod(r[15 + k]);
    E.published_volume = detail::number(r[18]);
    if (!r[19].empty()) {
      const auto it = element_index.find(r[19]);
      if (it == element_index.end()) throw std::runtime_error("anatomy::body: a copy of an element not yet read, " + r[19]);
      E.duplicate_of = it->second;
    }
    packs = std::max(packs, E.pack + 1);
    element_index.emplace(E.id, B.elements.size());
    B.atlas.element_concept.push_back(E.concept_);
    B.elements.push_back(std::move(E));
  }
  B.atlas.listed_isa.assign(B.elements.size(), {});
  B.atlas.listed_partof.assign(B.elements.size(), {});
  for (auto [suffix, out, listed] : {std::tuple{"-compounds-isa.tsv", &B.compounds_isa, &B.atlas.listed_isa},
                                      std::tuple{"-compounds-partof.tsv", &B.compounds_partof, &B.atlas.listed_partof}})
    for (const auto& r : detail::read_tsv(file(suffix))) {
      const auto it = element_index.find(r[1]);
      if (it == element_index.end()) throw std::runtime_error("anatomy::body: a compound of an unknown element, " + r[1]);
      const std::size_t c = concept_of(r[0]);
      out->push_back({c, it->second});
      (*listed)[it->second].push_back(c);
    }
  for (std::size_t k = 0; k < packs; k++) B.packs.push_back(detail::read_bytes(file("-" + std::to_string(k) + ".bin")));
  return B;
}

}  // namespace anatomy::body
