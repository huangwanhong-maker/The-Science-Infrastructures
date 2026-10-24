// Checks of exact repeats against their definitions (Gusfield 1997, §7.12): on short random
// words, with and without N, on one strand and on both, the maximal and supermaximal repeats
// and the longest repeated substring found from the suffix array are those found by brute
// force over every factor and every pair of its occurrences; on a simulated genome, every
// element planted in exact copies on random strands is one maximal repeat holding all its
// copies; and, when the data cache has Escherichia coli K-12, its repeats of 100 bases and more
// on both strands are found and timed.

#include <genomics/annotation/repeats.hpp>
#include <genomics/sequence/genbank.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
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
using namespace genomics::annotation;

namespace {

std::string read_file(const std::filesystem::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::ostringstream s;
  s << f.rdbuf();
  return s.str();
}

using Key = std::pair<std::size_t, std::vector<Occurrence>>;  // a repeat as its length and occurrences

// Every repeat by brute force: the factors of T (no separator, no N) with a maximal pair, mapped
// to the forward strand as the implementation reports them.
struct Brute {
  std::set<Key> maximal, supermaximal;
  std::size_t longest = 0;
};

Brute brute(const Codes& x, bool both) {
  const std::size_t n = x.size();
  std::vector<int> T;  // letters 10 + code, separators 0 and 1
  for (Code c : x) T.push_back(10 + c);
  if (both) {
    T.push_back(0);
    for (Code c : sequence::reverse_complement(x)) T.push_back(10 + c);
    T.push_back(1);
  }
  const std::size_t N = T.size();
  const auto usable = [&](std::size_t p, std::size_t l) {
    for (std::size_t i = p; i < p + l; i++)
      if (T[i] < 10 || T[i] == 14) return false;
    return true;
  };
  std::map<std::vector<int>, std::vector<std::size_t>> occ;
  for (std::size_t p = 0; p < N; p++)
    for (std::size_t l = 1; p + l <= N && usable(p, l); l++) occ[std::vector<int>(T.begin() + static_cast<std::ptrdiff_t>(p), T.begin() + static_cast<std::ptrdiff_t>(p + l))].push_back(p);
  Brute b;
  std::vector<std::vector<int>> maximal_words;
  for (const auto& [w, ps] : occ) {
    if (ps.size() >= 2) b.longest = std::max(b.longest, w.size());
    const std::size_t l = w.size();
    bool maximal = false;
    for (std::size_t a = 0; a < ps.size() && !maximal; a++)
      for (std::size_t c = a + 1; c < ps.size(); c++) {
        const std::size_t i = ps[a], j = ps[c];
        const bool left = i == 0 || j == 0 || T[i - 1] != T[j - 1];
        const bool right = i + l == N || j + l == N || T[i + l] != T[j + l];
        if (left && right) { maximal = true; break; }
      }
    if (maximal) maximal_words.push_back(w);
  }
  const auto key_of = [&](const std::vector<int>& w) -> std::pair<bool, Key> {
    const std::size_t l = w.size();
    if (both) {
      // Canonical: not larger than its reverse complement.
      std::vector<int> rc(l);
      for (std::size_t i = 0; i < l; i++) rc[i] = 10 + (3 - (w[l - 1 - i] - 10));
      if (rc < w) return {false, {}};
    }
    std::vector<Occurrence> o;
    for (std::size_t p : occ[w]) o.push_back(!both || p < n ? Occurrence{p, false} : Occurrence{n - (p - n - 1) - l, true});
    std::sort(o.begin(), o.end());
    o.erase(std::unique(o.begin(), o.end()), o.end());
    return {true, {l, o}};
  };
  for (const auto& w : maximal_words) {
    const auto [keep, k] = key_of(w);
    if (!keep) continue;
    b.maximal.insert(k);
    bool contained = false;
    for (const auto& v : maximal_words)
      if (v.size() > w.size() && std::search(v.begin(), v.end(), w.begin(), w.end()) != v.end()) { contained = true; break; }
    if (!contained) b.supermaximal.insert(k);
  }
  return b;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(20261024, 0x72657074, 0, 0);

  // ---- against brute force ------------------------------------------------------------------------------------
  {
    std::size_t cases = 0, repeats = 0;
    bool maximal_ok = true, super_ok = true, longest_ok = true;
    std::string first_wrong;
    for (int rep = 0; rep < 400; rep++) {
      const bool both = rep % 2;
      const std::size_t n = 2 + r.below(28), sigma = rep % 4 < 2 ? 2 : 4;
      Codes x(n);
      for (auto& c : x) c = static_cast<Code>(r.below(sigma));
      if (rep % 7 == 0) x[r.below(n)] = sequence::N;
      const Brute b = brute(x, both);
      RepeatOptions o;
      o.both_strands = both;
      std::set<Key> found, found_super;
      for (const auto& rp : maximal_repeats(x, o)) {
        found.insert({rp.length, rp.occurrences});
        if (rp.supermaximal) found_super.insert({rp.length, rp.occurrences});
      }
      std::set<Key> only_super;
      for (const auto& rp : supermaximal_repeats(x, o)) only_super.insert({rp.length, rp.occurrences});
      cases++;
      repeats += found.size();
      if (found != b.maximal && first_wrong.empty()) first_wrong = sequence::decode_dna(x) + (both ? " (both strands)" : "");
      maximal_ok = maximal_ok && found == b.maximal;
      super_ok = super_ok && found_super == b.supermaximal && only_super == b.supermaximal;
      const Repeat lr = longest_repeat(x, both);
      bool real = lr.length == 0 || lr.occurrences.size() == 2;
      longest_ok = longest_ok && lr.length == b.longest && real;
    }
    s.require("maximal repeats from the suffix array are those of brute force (" + std::to_string(cases) + " words of up to 29 bases, " + std::to_string(repeats) + " repeats, one strand and both, with N)",
              maximal_ok, first_wrong);
    s.require("supermaximal repeats are the maximal ones contained in no other, as brute force finds them", super_ok);
    s.require("the longest repeated substring has the length brute force finds", longest_ok);
  }

  // ---- a genome with planted repeats -------------------------------------------------------------------------------
  {
    simulate::GenomeSpec spec;
    spec.length = 60000;
    spec.repeats = {3, 4, 500, 0.0};
    const simulate::Genome g = simulate::make_genome(spec, 23);
    RepeatOptions o;
    o.min_length = 400;
    o.both_strands = true;
    const auto reps = maximal_repeats(g.sequence, o);
    bool all = true;
    for (std::size_t e = 0; e < 3; e++) {
      bool found = false;
      for (const auto& rp : reps) {
        std::size_t held = 0;
        for (const auto& c : g.repeats) {
          if (c.element != e) continue;
          for (const auto& oc : rp.occurrences)
            if (oc.position <= c.begin && oc.position + rp.length >= c.end) { held++; break; }
        }
        if (held == 4 && rp.occurrences.size() == 4) found = true;
      }
      all = all && found;
    }
    const Repeat lr = longest_repeat(g.sequence, true);
    s.require("each element planted in four exact copies on random strands is one maximal repeat holding its four copies (both strands)", all, std::to_string(reps.size()) + " repeats of 400 bases or more");
    s.require("the longest repeat is a planted element, 500 bases and a few by chance", lr.length >= 500 && lr.length < 520, std::to_string(lr.length));
  }

  // ---- Escherichia coli K-12 ------------------------------------------------------------------------------------------
  {
    const char* env = std::getenv("SCIENCE_CACHE");
    const std::filesystem::path cache = env ? std::filesystem::path(env) : std::filesystem::path(GENOMICS_CACHE);
    const std::filesystem::path file = cache / "ncbi-genomes" / "NC_000913.3.gb";
    if (!std::filesystem::exists(file)) {
      s.skip("Escherichia coli K-12 is not in the data cache: its repeats are not measured");
    } else {
      const auto recs = sequence::read_genbank(read_file(file));
      const Codes genome = sequence::encode_dna(recs.at(0).sequence);
      RepeatOptions o;
      o.min_length = 100;
      o.both_strands = true;
      const auto t0 = std::chrono::steady_clock::now();
      const auto reps = maximal_repeats(genome, o);
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      const Repeat lr = longest_repeat(genome, true);
      std::size_t super = 0;
      for (const auto& rp : reps) super += rp.supermaximal;
      char detail[200];
      std::snprintf(detail, sizeof detail, "%zu maximal repeats (%zu supermaximal), the longest %zu bases at %zu and %zu; %.2f s for 9.3 Mb of text", reps.size(), super, lr.length,
                    lr.occurrences.empty() ? 0 : lr.occurrences[0].position, lr.occurrences.size() < 2 ? 0 : lr.occurrences[1].position, seconds);
      s.require("E. coli K-12 MG1655: its maximal repeats of 100 bases or more on both strands, the longest first", !reps.empty() && reps.front().length == lr.length, detail);
    }
  }
  return s.done();
}
