// Checks of inter-sequence vectorisation (one target per lane): the best local score of a query
// against each of a batch of targets of mixed lengths, empty ones among them, equals Gotoh's
// recurrences (gotoh_score in local mode) for DNA and protein; lanes whose scores reach the
// 16-bit limit are rescored and still exact.

#include <genomics/align/gotoh.hpp>
#include <genomics/align/inter.hpp>
#include <genomics/align/scoring.hpp>
#include <genomics/simulate/genome.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cstdint>
#include <string>
#include <vector>

using namespace genomics;
using namespace genomics::align;

int main() {
  calx::verification::Suite s;
  calx::random::Stream r(5, 0x494e5452, 0, 0);
  const auto random_codes = [&](std::size_t n, int alphabet) {
    Codes x(n);
    for (auto& c : x) c = static_cast<Code>(r.below(static_cast<std::uint64_t>(alphabet)));
    return x;
  };
  struct Case {
    std::string name;
    Scoring scoring;
    int alphabet;
  };
  for (const Case& c : {Case{"DNA, minimap2 map-ont", minimap2_map_ont(), 4}, Case{"DNA, blastn", blastn(), 5}, Case{"protein, blastp", blastp(), 20}}) {
    const Codes q = random_codes(300, c.alphabet);
    std::vector<Codes> targets;
    for (std::size_t k = 0; k < 75; k++) {
      Codes t = random_codes(k % 9 == 0 ? 0 : r.below(500), c.alphabet);
      if (k % 3 == 0 && !t.empty()) {  // a copy of a piece of the query, for a high score
        const std::size_t a = r.below(150), len = 50 + r.below(100);  // a + len < 300
        t.insert(t.begin() + static_cast<std::ptrdiff_t>(r.below(t.size())), q.begin() + static_cast<std::ptrdiff_t>(a), q.begin() + static_cast<std::ptrdiff_t>(a + len));
      }
      targets.push_back(std::move(t));
    }
    InterStats st;
    const auto got = inter_scores(q, targets, c.scoring, &st);
    int ok = 0;
    for (std::size_t k = 0; k < targets.size(); k++) ok += got[k] == gotoh_score(q, targets[k], c.scoring, Mode::local);
    s.require(c.name + ": the 75 scores of a batch, in " + std::to_string(st.batches) + " batches of " + std::to_string(calx::cpu::Lanes<std::int16_t>::width) + " lanes, are Gotoh's",
              ok == 75, std::to_string(ok) + " of 75");
  }
  {
    // Scores past 32767: a 20 kb query and its copies.
    const Codes q = simulate::random_sequence(20'000, 0.5, 7);
    std::vector<Codes> targets = {q, Codes(q.begin(), q.begin() + 5000), random_codes(3000, 4)};
    InterStats st;
    const auto got = inter_scores(q, targets, minimap2_map_ont(), &st);
    s.require("scores past the 16-bit limit (40000) are rescored and exact", got[0] == 40000 && got[1] == 10000 && got[2] == gotoh_score(q, targets[2], minimap2_map_ont(), Mode::local) && st.rescored == 1,
              std::to_string(got[0]) + ", " + std::to_string(got[1]));
  }
  return s.done();
}
