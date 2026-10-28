// Checks of Kraken's classification. The database against brute force on small genomes sharing
// segments on both strands: every canonical k-mer mapped to the LCA of the taxa of the genomes
// that contain it, absent k-mers to none, and the same database whatever the order of the
// genomes. The heaviest root-to-leaf path against the classification tree built explicitly, its
// leaves weighed one by one, ties resolved by their LCA; the confidence threshold against the
// label climbed by hand. And, when the data cache has the fourteen RefSeq genomes (E. coli, its
// phages lambda and phiX174, SARS-CoV-2, ten primate mitochondria), reads simulated from each with
// sequencing errors: classified on their own lineage (to their genome's leaf, or to an ancestor
// on k-mers shared with relatives) at the rates measured on 2026-10-24, the same on one thread
// and on four; and, with each genome left out of the database in turn, its reads unclassified when
// no relative shares its 31-mers (E. coli, phiX174, SARS-CoV-2), lambda's on E. coli K-12 or
// unclassified, a primate's never off the primates and mostly within the clade of its deepest
// ancestor still in the database: on that ancestor, or (Kraken's rule allows no less) on a
// relative's own taxon.

#include <genomics/classify/kraken.hpp>
#include <genomics/classify/taxonomy.hpp>
#include <genomics/simulate/genome.hpp>
#include <genomics/simulate/reads.hpp>

#include <calx/random/stream.hpp>
#include <calx/scheduler/thread_pool.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::classify;
using sequence::Codes;
using sequence::Kmer;

namespace {

std::string fixed(double v, int digits = 3) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*f", digits, v);
  return b;
}

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

Codes random_dna(calx::random::Stream& r, std::size_t n) {
  Codes x(n);
  for (auto& c : x) c = static_cast<sequence::Code>(r.below(4));
  return x;
}

std::set<Kmer> canonical_kmers(const Codes& x, int k) {
  std::set<Kmer> out;
  sequence::for_each_kmer(x, k, [&](std::size_t, Kmer f, Kmer r) { out.insert(std::min(f, r)); });
  return out;
}

// The label by the definition: the classification tree (the hit taxa and their ancestors), its
// leaves, each leaf's path weighed, the LCA of the heaviest leaves.
TaxonId brute_label(const Taxonomy& t, const std::vector<std::pair<TaxonId, std::uint32_t>>& hits) {
  if (hits.empty()) return NO_TAXON;
  std::set<TaxonId> nodes;
  for (const auto& [x, c] : hits)
    for (TaxonId u : t.path(x)) nodes.insert(u);
  std::vector<TaxonId> leaves;
  for (TaxonId u : nodes) {
    bool leaf = true;
    for (TaxonId v : nodes)
      if (v != u && t[v].parent == u) leaf = false;
    if (leaf) leaves.push_back(u);
  }
  std::uint64_t best = 0;
  std::vector<TaxonId> top;
  for (TaxonId l : leaves) {
    std::uint64_t w = 0;
    const auto p = t.path(l);
    for (const auto& [x, c] : hits)
      if (std::find(p.begin(), p.end(), x) != p.end()) w += c;
    if (w > best) { best = w; top = {l}; }
    else if (w == best) top.push_back(l);
  }
  TaxonId label = top[0];
  for (TaxonId l : top) label = t.lca(label, l);
  return label;
}

struct Rates {
  std::size_t reads = 0, leaf = 0, ancestor = 0, relative = 0, outside = 0, unclassified = 0;
};

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(2014, 0x4b52414b, 0, 0);

  // ---- the database against brute force ----------------------------------------------------------------
  {
    bool lca = true, absent = true, order = true, size = true;
    for (int trial = 0; trial < 12; trial++) {
      Taxonomy t;
      std::vector<TaxonId> taxa;
      for (int i = 0; i < 10; i++) {
        std::vector<std::string> lineage;
        for (std::size_t d = 0, depth = 1 + r.below(4); d < depth; d++) lineage.push_back("c" + std::to_string(d) + "_" + std::to_string(r.below(2)));
        taxa.push_back(t.add(lineage, "g" + std::to_string(i)));
      }
      // Genomes of fresh sequence and copies of shared segments, some reverse-complemented; some on
      // an inner taxon.
      std::vector<Codes> segments;
      for (int i = 0; i < 6; i++) segments.push_back(random_dna(r, 20 + r.below(40)));
      std::vector<Codes> genomes;
      std::vector<TaxonId> of;
      for (int i = 0; i < 10; i++) {
        Codes g = random_dna(r, 50 + r.below(100));
        for (int j = 0; j < 4; j++) {
          Codes seg = segments[r.below(segments.size())];
          if (r.below(2)) seg = sequence::reverse_complement(seg);
          g.insert(g.end(), seg.begin(), seg.end());
          const Codes more = random_dna(r, r.below(30));
          g.insert(g.end(), more.begin(), more.end());
        }
        genomes.push_back(g);
        of.push_back(r.below(5) == 0 ? t[taxa[static_cast<std::size_t>(i)]].parent : taxa[static_cast<std::size_t>(i)]);
      }
      const int k = 7 + static_cast<int>(r.below(5));
      Database forward(t, k), backward(t, k);
      for (std::size_t i = 0; i < genomes.size(); i++) forward.add(genomes[i], of[i]);
      for (std::size_t i = genomes.size(); i-- > 0;) backward.add(genomes[i], of[i]);
      std::vector<std::set<Kmer>> sets;
      std::set<Kmer> all;
      for (const auto& g : genomes) {
        sets.push_back(canonical_kmers(g, k));
        all.insert(sets.back().begin(), sets.back().end());
      }
      for (Kmer x : all) {
        TaxonId want = NO_TAXON;
        for (std::size_t i = 0; i < genomes.size(); i++)
          if (sets[i].count(x)) want = want == NO_TAXON ? of[i] : t.lca(want, of[i]);
        lca = lca && forward.lookup(x) == want;
        order = order && backward.lookup(x) == want;
      }
      size = size && forward.size() == all.size() && backward.size() == all.size();
      for (int q = 0; q < 2000; q++) {
        const Kmer x = sequence::canonical(r.next64() & sequence::kmer_mask(k), k);
        if (!all.count(x)) absent = absent && forward.lookup(x) == NO_TAXON;
      }
    }
    s.require("on 12 sets of 10 genomes sharing segments on both strands (k = 7 to 11), every canonical k-mer maps to the LCA of the taxa of the genomes containing it (brute force)", lca);
    s.require("the database is the same built in the opposite order of the genomes", order);
    s.require("it holds each distinct canonical k-mer once, and absent k-mers map to no taxon", size && absent);
  }

  // ---- the heaviest root-to-leaf path against the classification tree -------------------------------
  {
    bool same = true, climb = true;
    for (int trial = 0; trial < 3000; trial++) {
      Taxonomy t;
      for (int i = 0; i < 12; i++) {
        std::vector<std::string> lineage;
        for (std::size_t d = 0, depth = 1 + r.below(5); d < depth; d++) lineage.push_back("c" + std::to_string(d) + "_" + std::to_string(r.below(2)));
        t.add(lineage, r.below(3) ? "g" + std::to_string(i) : "");
      }
      std::map<TaxonId, std::uint32_t> m;
      for (std::size_t j = 0, n = 1 + r.below(6); j < n; j++) m[static_cast<TaxonId>(r.below(t.size()))] += 1 + static_cast<std::uint32_t>(r.below(4));
      const std::vector<std::pair<TaxonId, std::uint32_t>> hits(m.begin(), m.end());
      const TaxonId label = heaviest_path(t, hits);
      same = same && label == brute_label(t, hits);
      // The confidence rule by hand: climb while the clade holds less than τ of `kmers`.
      std::uint32_t total = 0;
      for (const auto& [x, c] : hits) total += c;
      const std::uint32_t kmers = total + static_cast<std::uint32_t>(r.below(10));
      const double tau = 0.3 + 0.1 * static_cast<double>(r.below(7));
      TaxonId u = label;
      const auto inside = [&](TaxonId v) { std::uint32_t n = 0; for (const auto& [x, c] : hits) if (t.above(v, x)) n += c; return n; };
      while (inside(u) < tau * kmers && u != ROOT) u = t[u].parent;
      const bool keep = inside(u) >= tau * kmers;
      climb = climb && clade_hits(t, u, hits) == inside(u) && (keep || inside(ROOT) < tau * kmers);
    }
    s.require("on 3000 random sets of hits, the heaviest path's label equals the LCA of the heaviest leaves of the classification tree built explicitly", same);
    s.require("clade hits count the hits at or below a taxon, and the climb to the confidence stops at the first clade that holds enough", climb);
  }

  // ---- the fourteen genomes ------------------------------------------------------------------------------
  const char* env = std::getenv("SCIENCE_CACHE");
  const std::filesystem::path dir = (env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE)) / "ncbi-genomes";
  const std::vector<std::string> files{"NC_000913.3", "NC_001416.1", "NC_001422.1", "NC_045512.2", "NC_012920.1", "NC_011137.1", "NC_001643.1",
                                       "NC_001644.1", "NC_011120.1", "NC_001646.1", "NC_002083.1", "NC_002082.1", "NC_005943.1", "NC_020006.2"};
  bool present = true;
  for (const auto& f : files) present = present && std::filesystem::exists(dir / (f + ".gb"));
  if (!present) {
    s.skip("the fourteen RefSeq genomes (data/cache/ncbi-genomes) are not all in the data cache: reads of real genomes are not classified");
    return s.done();
  }
  Taxonomy tax;
  std::vector<Codes> genomes;
  std::vector<TaxonId> leaf;
  std::vector<std::string> organism;
  for (const auto& f : files) {
    const auto recs = sequence::read_genbank(read_file(dir / (f + ".gb")));
    genomes.push_back(sequence::encode_dna(recs.at(0).sequence));
    leaf.push_back(tax.add(recs.at(0)));
    organism.push_back(recs.at(0).organism);
  }
  // Reads: 600 of 150 bases from each genome, with the simulator's errors (0.1 % to 1 % along the read).
  const std::size_t per = 600;
  std::vector<std::vector<Codes>> reads(genomes.size());
  for (std::size_t g = 0; g < genomes.size(); g++)
    for (const auto& rd : simulate::short_reads(simulate::Source::haploid({genomes[g]}), per, {}, 4600 + g)) reads[g].push_back(rd.bases);

  calx::scheduler::ThreadPool one(1), four(4);
  const auto t0 = std::chrono::steady_clock::now();
  Database db(tax, 31);
  for (std::size_t g = 0; g < genomes.size(); g++) db.add(genomes[g], leaf[g]);
  const double build = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::size_t at_root = 0, inner = 0;
  db.for_each([&](Kmer, TaxonId x) { at_root += x == ROOT; inner += std::find(leaf.begin(), leaf.end(), x) == leaf.end(); });
  {
    Rates all;
    std::string table;
    bool deterministic = true;
    double seconds = 0;
    std::size_t distinct_leaf = 0, distinct_total = 0;
    for (std::size_t g = 0; g < genomes.size(); g++) {
      const auto t1 = std::chrono::steady_clock::now();
      const auto c = classify_all(db, reads[g], four);
      seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
      const auto c1 = classify_all(db, reads[g], one);
      Rates x;
      for (std::size_t i = 0; i < c.size(); i++) {
        deterministic = deterministic && c[i].taxon == c1[i].taxon && c[i].hits == c1[i].hits;
        x.reads++;
        if (!c[i].classified()) x.unclassified++;
        else if (c[i].taxon == leaf[g]) x.leaf++;
        else if (tax.above(c[i].taxon, leaf[g])) x.ancestor++;
        else x.outside++;
      }
      all.reads += x.reads, all.leaf += x.leaf, all.ancestor += x.ancestor, all.outside += x.outside, all.unclassified += x.unclassified;
      // Genomes with no relative in the database closer than a family (all but Homo, Pan and Pongo).
      const std::string& o = organism[g];
      if (o.rfind("Homo", 0) != 0 && o.rfind("Pan ", 0) != 0 && o.rfind("Pongo", 0) != 0) distinct_leaf += x.leaf, distinct_total += x.reads;
      table += (table.empty() ? "" : "; ") + o + " " + fixed(static_cast<double>(x.leaf) / x.reads, 2) + "/" + fixed(static_cast<double>(x.ancestor) / x.reads, 2);
    }
    const double on_lineage = static_cast<double>(all.leaf + all.ancestor) / all.reads;
    s.require("a database of the 14 genomes at k = 31: " + std::to_string(db.size()) + " canonical k-mers, " + std::to_string(inner) + " of them on an inner taxon, " + std::to_string(at_root) + " at the root; built in " + fixed(build, 2) + " s",
              db.size() > 4'000'000 && db.genomes() == 14);
    s.require("8400 simulated reads, 600 a genome: at least 0.99 classified on their own lineage (their genome's leaf or an ancestor)", on_lineage >= 0.99,
              fixed(on_lineage, 4) + " on the lineage (" + fixed(static_cast<double>(all.leaf) / all.reads, 3) + " at the leaf), " + std::to_string(all.outside) + " off it, " + std::to_string(all.unclassified) +
                  " unclassified; leaf/ancestor by genome: " + table + "; " + fixed(seconds, 3) + " s on 4 threads");
    s.require("reads of genomes without a close relative in the database (E. coli, the phages, SARS-CoV-2, gorilla, gibbon, macaque, baboon) go to their own leaf, at least 0.98 of them",
              static_cast<double>(distinct_leaf) / distinct_total >= 0.98, fixed(static_cast<double>(distinct_leaf) / distinct_total, 4));
    s.require("the classifications are the same on one thread and on four", deterministic);
  }

  // ---- each genome left out ------------------------------------------------------------------------------
  {
    // Measured on 2026-10-24 (600 reads a genome): E. coli, phiX174 and SARS-CoV-2 left out, all
    // unclassified; lambda 0.89 unclassified and the rest on E. coli; a primate left out, never off
    // the primates, its classified reads 0.880 within the clade of its deepest ancestor still in the
    // database (that ancestor, or a relative's taxon below it), the rest right only at a higher
    // ancestor (Hominidae or Catarrhini).
    Rates all, primates;
    std::size_t microbes = 0, microbes_unclassified = 0, lambda_reads = 0, lambda_unclassified = 0, lambda_other = 0, off_primates = 0;
    const TaxonId catarrhini = tax.find("Catarrhini");
    std::string table;
    for (std::size_t out = 0; out < genomes.size(); out++) {
      Database loo(tax, 31);
      for (std::size_t g = 0; g < genomes.size(); g++)
        if (g != out) loo.add(genomes[g], leaf[g]);
      // The deepest ancestor of the genome left out whose clade still holds a genome of the database.
      TaxonId shared = ROOT;
      for (std::size_t g = 0; g < genomes.size(); g++)
        if (g != out) {
          const TaxonId a = tax.lca(leaf[g], leaf[out]);
          if (tax[a].depth > tax[shared].depth) shared = a;
        }
      const auto c = classify_all(loo, reads[out], four);
      const bool primate = catarrhini != NO_TAXON && tax.above(catarrhini, leaf[out]);
      Rates x;
      for (const auto& y : c) {
        x.reads++;
        if (!y.classified()) x.unclassified++;
        else if (tax.above(y.taxon, leaf[out])) x.ancestor++;  // the shared ancestor or above
        else if (tax.above(shared, y.taxon)) x.relative++;     // within the shared clade, a relative's taxon
        else x.outside++;
        if (primate && y.classified() && !tax.above(catarrhini, y.taxon) && !tax.above(y.taxon, catarrhini)) off_primates++;
      }
      if (organism[out] == "Escherichia phage Lambda") {
        lambda_reads += x.reads;
        lambda_unclassified += x.unclassified;
        for (const auto& y : c) lambda_other += y.classified() && y.taxon != leaf[0];
      } else if (!primate) {
        microbes += x.reads;
        microbes_unclassified += x.unclassified;
      } else {
        primates.reads += x.reads, primates.ancestor += x.ancestor, primates.relative += x.relative, primates.outside += x.outside, primates.unclassified += x.unclassified;
      }
      all.reads += x.reads, all.ancestor += x.ancestor, all.relative += x.relative, all.outside += x.outside, all.unclassified += x.unclassified;
      table += (table.empty() ? "" : "; ") + organism[out] + " (" + tax[shared].name + ") " + fixed(static_cast<double>(x.ancestor) / x.reads, 2) + "/" + fixed(static_cast<double>(x.relative) / x.reads, 2) + "/" +
               fixed(static_cast<double>(x.outside) / x.reads, 2) + "/" + fixed(static_cast<double>(x.unclassified) / x.reads, 2);
    }
    const std::string detail = "by genome (deepest shared ancestor) ancestor/relative/beyond/unclassified: " + table;
    s.require("each genome left out in turn: E. coli, phiX174 and SARS-CoV-2, with no relative in the database, leave at least 0.99 of their reads unclassified",
              static_cast<double>(microbes_unclassified) / microbes >= 0.99, fixed(static_cast<double>(microbes_unclassified) / microbes, 4));
    s.require("lambda left out: its reads unclassified (at least 0.85) or classified to E. coli K-12, with which it shares 31-mers, and to nothing else",
              static_cast<double>(lambda_unclassified) / lambda_reads >= 0.85 && lambda_other == 0,
              fixed(static_cast<double>(lambda_unclassified) / lambda_reads) + " unclassified, " + std::to_string(lambda_other) + " elsewhere than E. coli");
    const double classified = static_cast<double>(primates.reads - primates.unclassified);
    s.require("a primate mitochondrion left out: none of its reads classified off the primates (Catarrhini)", off_primates == 0, std::to_string(off_primates) + " off; " + detail);
    s.require("and at least 0.85 of its classified reads land within the clade of its deepest ancestor still in the database (that ancestor or a relative's taxon), the rest at a higher ancestor's",
              static_cast<double>(primates.ancestor + primates.relative) / classified >= 0.85,
              fixed(static_cast<double>(primates.ancestor + primates.relative) / classified) + " (ancestor " + fixed(static_cast<double>(primates.ancestor) / classified) + ", a relative's taxon " +
                  fixed(static_cast<double>(primates.relative) / classified) + "), " + fixed(static_cast<double>(primates.unclassified) / primates.reads) + " of the reads unclassified");
  }
  return s.done();
}
