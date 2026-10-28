// A taxonomy built from the lineages of GenBank records, and the lowest common ancestor of its
// taxa.
//
// - A GenBank record's ORGANISM field gives the organism's name and its lineage, the names of the
//   taxa above it from the top down (sequence/genbank.hpp reads both). The taxonomy here is the
//   trie of the lineages under one root: a taxon is a name under a parent, so that records whose
//   lineages begin alike share those taxa, and the organism is a leaf below its lineage's last
//   name (that taxon itself when the lineage ends with the organism's name). No taxonomy is typed
//   in or read from elsewhere: the tree knows exactly what its records say.
//   *Caveat.* A flat file's lineage carries no ranks. A taxon's depth here is its position along
//   the lineages, not a Linnaean rank, and lineages of different lengths put taxa of one rank at
//   different depths (Cercopithecinae stands between Cercopithecidae and Macaca, while Hominidae
//   is followed by Homo). Nor does a lineage go below the genus for most organisms: Homo sapiens
//   and Homo sapiens neanderthalensis are two leaves of Homo, siblings, not species and
//   subspecies.
// - The lowest common ancestor lca(a, b) is the deepest taxon above both, each taxon being above
//   itself: the deeper of the two is lifted to the other's depth, then both together until they
//   meet, O(depth). It is the meet of the tree's order (a ≤ b when a is above b), so it is
//   associative, commutative and idempotent, with the root absorbing: the LCA of a set of taxa is
//   the same in any order, which is what lets a database of k-mers be built in any order of its
//   genomes (kraken.hpp).
//
// Taxa are numbered in the order they are first met; the root is 0.

#pragma once

#include <genomics/sequence/genbank.hpp>

#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genomics::classify {

using TaxonId = std::uint32_t;
inline constexpr TaxonId ROOT = 0;
inline constexpr TaxonId NO_TAXON = ~TaxonId{0};

struct Taxon {
  std::string name;
  TaxonId parent = ROOT;  // the root is its own parent
  std::uint32_t depth = 0;
};

class Taxonomy {
 public:
  Taxonomy() { taxa_.push_back({"root", ROOT, 0}); }

  std::size_t size() const { return taxa_.size(); }
  const Taxon& operator[](TaxonId t) const { return taxa_.at(t); }

  // The taxon named `name` under `parent`, added if absent.
  TaxonId child(TaxonId parent, const std::string& name) {
    if (parent >= taxa_.size()) throw std::invalid_argument("genomics::classify::Taxonomy::child: no such parent");
    const auto key = std::make_pair(parent, name);
    const auto it = index_.find(key);
    if (it != index_.end()) return it->second;
    const auto id = static_cast<TaxonId>(taxa_.size());
    taxa_.push_back({name, parent, taxa_[parent].depth + 1});
    index_.emplace(key, id);
    return id;
  }
  // The leaf of an organism under its lineage, from the top down, added as needed.
  TaxonId add(const std::vector<std::string>& lineage, const std::string& organism) {
    TaxonId t = ROOT;
    for (const auto& name : lineage) t = child(t, name);
    if (!organism.empty() && (lineage.empty() || lineage.back() != organism)) t = child(t, organism);
    return t;
  }
  TaxonId add(const sequence::GenBank& record) { return add(record.lineage, record.organism); }

  // The taxon with this name, the first met; NO_TAXON if none.
  TaxonId find(const std::string& name) const {
    for (std::size_t t = 0; t < taxa_.size(); t++)
      if (taxa_[t].name == name) return static_cast<TaxonId>(t);
    return NO_TAXON;
  }

  TaxonId lca(TaxonId a, TaxonId b) const {
    if (a >= taxa_.size() || b >= taxa_.size()) throw std::invalid_argument("genomics::classify::Taxonomy::lca: no such taxon");
    while (taxa_[a].depth > taxa_[b].depth) a = taxa_[a].parent;
    while (taxa_[b].depth > taxa_[a].depth) b = taxa_[b].parent;
    while (a != b) {
      a = taxa_[a].parent;
      b = taxa_[b].parent;
    }
    return a;
  }
  // Whether a is above b or is b.
  bool above(TaxonId a, TaxonId b) const {
    if (a >= taxa_.size() || b >= taxa_.size()) throw std::invalid_argument("genomics::classify::Taxonomy::above: no such taxon");
    while (taxa_[b].depth > taxa_[a].depth) b = taxa_[b].parent;
    return a == b;
  }
  // The taxa from the root down to t.
  std::vector<TaxonId> path(TaxonId t) const {
    std::vector<TaxonId> p(taxa_.at(t).depth + 1);
    for (std::size_t i = p.size(); i-- > 0; t = taxa_[t].parent) p[i] = t;
    return p;
  }
  // The names from the root's child down to t, joined by "; ".
  std::string lineage(TaxonId t) const {
    std::string s;
    for (TaxonId u : path(t))
      if (u != ROOT) s += (s.empty() ? "" : "; ") + taxa_[u].name;
    return s.empty() ? "root" : s;
  }

 private:
  std::vector<Taxon> taxa_;
  std::map<std::pair<TaxonId, std::string>, TaxonId> index_;
};

}  // namespace genomics::classify
