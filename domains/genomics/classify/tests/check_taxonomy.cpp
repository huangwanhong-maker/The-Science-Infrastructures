// Checks of the taxonomy of lineages: the lowest common ancestor against the deepest common taxon
// of the two root paths, by brute force, on random taxonomies; the LCA associative, commutative
// and idempotent, with the root absorbing; "above" against membership of the root path; a
// lineage wrapped over two lines of a GenBank record read as one name; and, when the data cache
// has them, the fourteen RefSeq genomes' lineages joined into one tree: the human and Neanderthal
// mitochondria two leaves of Homo, chimpanzee and bonobo of Pan, the LCA of human and chimpanzee
// Hominidae, of human and macaque Catarrhini, of E. coli and its phage lambda the root.

#include <genomics/classify/taxonomy.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::classify;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

// A random taxonomy: lineages over a small vocabulary, so that they share prefixes.
Taxonomy random_taxonomy(calx::random::Stream& r, std::vector<TaxonId>& leaves) {
  Taxonomy t;
  for (int i = 0; i < 40; i++) {
    std::vector<std::string> lineage;
    const std::size_t depth = 1 + r.below(7);
    for (std::size_t d = 0; d < depth; d++) lineage.push_back("t" + std::to_string(d) + "_" + std::to_string(r.below(3)));
    leaves.push_back(t.add(lineage, "organism " + std::to_string(i)));
  }
  return t;
}

// The LCA by brute force: the deepest taxon on both root paths.
TaxonId brute_lca(const Taxonomy& t, TaxonId a, TaxonId b) {
  const auto pa = t.path(a), pb = t.path(b);
  TaxonId best = ROOT;
  for (TaxonId x : pa)
    if (std::find(pb.begin(), pb.end(), x) != pb.end() && t[x].depth >= t[best].depth) best = x;
  return best;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(2014, 0x5441584f, 0, 0);

  // ---- random taxonomies against brute force ----------------------------------------------------------
  {
    bool lca = true, laws = true, above = true, paths = true;
    for (int trial = 0; trial < 30; trial++) {
      std::vector<TaxonId> leaves;
      const Taxonomy t = random_taxonomy(r, leaves);
      const auto n = static_cast<TaxonId>(t.size());
      for (int q = 0; q < 300; q++) {
        const auto a = static_cast<TaxonId>(r.below(n)), b = static_cast<TaxonId>(r.below(n)), c = static_cast<TaxonId>(r.below(n));
        lca = lca && t.lca(a, b) == brute_lca(t, a, b);
        laws = laws && t.lca(a, b) == t.lca(b, a) && t.lca(t.lca(a, b), c) == t.lca(a, t.lca(b, c)) && t.lca(a, a) == a && t.lca(ROOT, a) == ROOT;
        const auto pb = t.path(b);
        above = above && t.above(a, b) == (std::find(pb.begin(), pb.end(), a) != pb.end()) && (t.lca(a, b) == a) == t.above(a, b);
      }
      for (TaxonId x = 0; x < n; x++) {
        const auto p = t.path(x);
        paths = paths && p.front() == ROOT && p.back() == x && p.size() == t[x].depth + 1;
        for (std::size_t i = 1; i < p.size(); i++) paths = paths && t[p[i]].parent == p[i - 1];
      }
    }
    s.require("on 30 random taxonomies of 40 lineages, the LCA of 9000 pairs is the deepest taxon on both root paths (brute force)", lca);
    s.require("the LCA is commutative, associative and idempotent, and the root absorbs it", laws);
    s.require("a is above b exactly when a lies on b's root path, exactly when lca(a, b) = a", above);
    s.require("every root path runs from the root to its taxon through parents, one taxon a depth", paths);
  }

  // ---- a lineage wrapped over two lines -------------------------------------------------------------
  {
    std::string gb;
    gb += "LOCUS       TEST0002                  60 bp    DNA     circular PHG 24-OCT-2026\n";
    gb += "DEFINITION  A test record whose lineage wraps a name over two lines.\n";
    gb += "ACCESSION   TEST0002\nVERSION     TEST0002.1\n";
    gb += "SOURCE      test phage\n  ORGANISM  Test phage one\n            Viruses; Testviria; Testvirae; Testviricota; Testviricetes;\n";
    gb += "            Testvirales; Testviridae; Testvirus; Testvirus\n            unus.\n";
    gb += "FEATURES             Location/Qualifiers\n     source          1..60\n";
    gb += "ORIGIN      \n        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n//\n";
    const auto recs = sequence::read_genbank(gb);
    const std::vector<std::string> want{"Viruses", "Testviria", "Testvirae", "Testviricota", "Testviricetes", "Testvirales", "Testviridae", "Testvirus", "Testvirus unus"};
    s.require("GenBank: a lineage name wrapped over two lines (Testvirus / unus) is read as one name", recs.size() == 1 && recs[0].lineage == want,
              recs.empty() ? "" : std::to_string(recs[0].lineage.size()) + " names, the last \"" + (recs[0].lineage.empty() ? "" : recs[0].lineage.back()) + "\"");
  }

  // ---- the fourteen genomes ---------------------------------------------------------------------------
  {
    const char* env = std::getenv("SCIENCE_CACHE");
    const std::filesystem::path dir = (env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE)) / "ncbi-genomes";
    const std::vector<std::string> files{"NC_000913.3", "NC_001416.1", "NC_001422.1", "NC_001643.1", "NC_001644.1", "NC_001646.1", "NC_002082.1",
                                         "NC_002083.1", "NC_005943.1", "NC_011120.1", "NC_011137.1", "NC_012920.1", "NC_020006.2", "NC_045512.2"};
    bool present = true;
    for (const auto& f : files) present = present && std::filesystem::exists(dir / (f + ".gb"));
    if (!present) {
      s.skip("the fourteen RefSeq genomes (data/cache/ncbi-genomes) are not all in the data cache: their taxonomy is not built");
    } else {
      Taxonomy t;
      std::map<std::string, TaxonId> leaf;
      std::map<std::string, std::string> organism;
      for (const auto& f : files) {
        const auto recs = sequence::read_genbank(read_file(dir / (f + ".gb")));
        leaf[f] = t.add(recs.at(0));
        organism[f] = recs.at(0).organism;
      }
      const auto name = [&](TaxonId x) { return t[x].name; };
      s.require("the 14 lineages join into one tree of " + std::to_string(t.size()) + " taxa, each genome a distinct leaf",
                [&] { std::vector<TaxonId> l; for (auto& [f, x] : leaf) l.push_back(x); std::sort(l.begin(), l.end()); return std::unique(l.begin(), l.end()) == l.end(); }());
      s.require("Homo sapiens and Homo sapiens neanderthalensis are two leaves of Homo", t[leaf["NC_012920.1"]].parent == t[leaf["NC_011137.1"]].parent && name(t[leaf["NC_012920.1"]].parent) == "Homo");
      s.require("chimpanzee and bonobo are two leaves of Pan, the two orangutans of Pongo",
                name(t.lca(leaf["NC_001643.1"], leaf["NC_001644.1"])) == "Pan" && name(t.lca(leaf["NC_001646.1"], leaf["NC_002083.1"])) == "Pongo");
      s.require("the LCA of human and chimpanzee is Hominidae, of human and gibbon or macaque Catarrhini, of macaque and baboon Cercopithecinae",
                name(t.lca(leaf["NC_012920.1"], leaf["NC_001643.1"])) == "Hominidae" && name(t.lca(leaf["NC_012920.1"], leaf["NC_002082.1"])) == "Catarrhini" &&
                    name(t.lca(leaf["NC_012920.1"], leaf["NC_005943.1"])) == "Catarrhini" && name(t.lca(leaf["NC_005943.1"], leaf["NC_020006.2"])) == "Cercopithecinae");
      s.require("the LCA of E. coli and lambda is the root, of lambda and phiX174 Viruses",
                t.lca(leaf["NC_000913.3"], leaf["NC_001416.1"]) == ROOT && name(t.lca(leaf["NC_001416.1"], leaf["NC_001422.1"])) == "Viruses");
      s.require("phiX174's lineage keeps its wrapped species name, and ends with its organism", t.lineage(leaf["NC_001422.1"]).find("Sinsheimervirus; Sinsheimervirus phiX174; Escherichia phage phiX174") != std::string::npos,
                t.lineage(leaf["NC_001422.1"]));
    }
  }
  return s.done();
}
