// Checks of the history as a Feynman graph and as a sample path (theory note §8).
//
// The stranded graph: every strand end at every occasion lies on one corner; each lens has one
// closing and one opening corner and a passing corner wherever it runs past an occasion; the ranks
// of the segments add up to the lengths of the lenses' boundaries. Proposition 8.1: eliminating
// the δ's of the Feynman rules over ℤ₂ leaves one free strand variable per lens and, as each
// kernel's argument, exactly that lens's boundary in K̂. Bubbles: a symmetric vertex has the
// automorphisms it should. The process: on a history whose answer is known by hand, the rates, the
// binding and the martingale checks. All of it on every canon.

#include <history/field/process.hpp>
#include <history/field/stranded.hpp>
#include <history/foam/refined.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

using namespace history;

namespace {

const std::string HEAD = R"({"kind": "canon", "id": "test", "title": "test", "codebook": 1}
{"kind": "source", "id": "s", "cite": "s"}
)";
std::string being(const std::string& id) { return R"({"kind": "being", "id": ")" + id + R"(", "name": ")" + id + R"(", "type": "person"})" + "\n"; }
std::string occasion(const std::string& id, const std::string& date, const std::string& present) {
  return R"({"kind": "occasion", "id": ")" + id + R"(", "date": ")" + date + R"(", "type": "meeting", "present": )" + present + R"(, "sources": ["s"], "confidence": 1})" + "\n";
}
std::string tie(const std::string& a, const std::string& b, const std::string& attested, const std::string& colour = "letters") {
  return R"({"kind": "tie", "id": ")" + a + "~" + b + R"(", "a": ")" + a + R"(", "b": ")" + b + R"(", "type": "friendship", "colour": ")" + colour + R"(", "attested": )" + attested +
         R"(, "sources": ["s"], "confidence": 1})" + "\n";
}

// The structural checks of §8.1 on one record.
void feynman(calx::verification::Suite& s, const std::string& name, const HistoryComplex& K) {
  const RefinedComplex R = refine(K);
  const StrandedGraph G = stranded(K);
  s.require(name + ": every strand end at every occasion lies on exactly one corner", pairing_complete(K, G));
  std::size_t ranks = 0, lengths = 0;
  for (std::size_t e = 0; e < K.segments.size(); e++) ranks += G.rank(e);
  for (const auto& f : K.cw.faces) lengths += f.boundary.size();
  s.require(name + ": the ranks of the segments add up to the lengths of the lenses", ranks == lengths, std::to_string(ranks) + " strands");
  std::vector<int> closing(K.lenses.size(), 0), opening(K.lenses.size(), 0);
  std::size_t corners = 0;
  for (const auto& B : G.bubbles)
    for (const auto& c : B.corners) {
      corners++;
      closing[c.lens] += c.kind == Corner::closing;
      opening[c.lens] += c.kind == Corner::opening;
    }
  bool once = true;
  for (std::size_t f = 0; f < K.lenses.size(); f++) once = once && closing[f] == 1 && opening[f] == 1;
  s.require(name + ": each lens closes once and opens once, and turns at every vertex of its boundary", once && corners == lengths);
  const FeynmanForms F = feynman_forms(K, R, G);
  const auto boundary = lens_boundaries(R, K.lenses.size());
  bool same = F.consistent;
  for (std::size_t f = 0; f < K.lenses.size(); f++) same = same && F.argument[f] == boundary[f];
  s.require(name + ": Prop. 8.1, each kernel's argument is its lens's holonomy on K̂", same);
  s.require(name + ": Prop. 8.1, one free strand variable per lens", F.free_strands == K.lenses.size(), std::to_string(F.free_strands) + " for " + std::to_string(K.lenses.size()) + " lenses");
  std::size_t typed = 0;
  for (const auto& [type, vs] : G.types) typed += vs.size();
  s.require(name + ": the bubble types partition the occasions", typed == K.order.size(), std::to_string(G.types.size()) + " types");
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // A small history with a tie of one meeting, carried signs and passing strands.
  {
    const HistoryComplex K = construct(read_record_text(
        HEAD + being("a") + being("b") + being("c") + being("d") + occasion("o1", "1800", R"([["a", ""], ["b", ""], ["c", ""]])") +
        occasion("o2", "1801", R"([["a", ""], ["b", ""], ["c", ""], ["d", ""]])") + occasion("o3", "1802", R"([["b", ""], ["c", ""], ["d", ""]])") +
        occasion("o4", "1803", R"([["a", ""], ["b", ""], ["c", ""], ["d", ""]])") + occasion("o5", "1804", R"([["a", ""], ["c", ""], ["d", ""]])") +
        tie("a", "b", R"([["o1", 1], ["o4", -1]])") + tie("b", "c", R"([["o1", 1], ["o3", 1]])") + tie("a", "c", R"([["o2", -1], ["o5", -1]])", "music") +
        tie("c", "d", R"([["o2", 1], ["o5", 1]])", "music") + tie("b", "d", R"([["o3", -1], ["o4", -1]])") + tie("a", "d", R"([["o4", 1]])")));
    feynman(s, "the small history", K);
    const StrandedGraph G = stranded(K);
    // At o3, a is absent: a~b and a~c run past it on b's and c's side.
    std::size_t passing = 0;
    for (const auto& c : G.bubbles[2].corners) passing += c.kind == Corner::passing;
    s.require("the small history: two lenses run past o3, each on the side of the being present", passing == 2, std::to_string(passing));
  }

  // Symmetry: two beings appearing together, untied, then three beings with a ring of three ties.
  {
    const HistoryComplex K = construct(read_record_text(
        HEAD + being("a") + being("b") + being("c") + occasion("o1", "1800", R"([["a", ""], ["b", ""], ["c", ""]])") +
        occasion("o2", "1801", R"([["a", ""], ["b", ""], ["c", ""]])") + occasion("o3", "1802", R"([["a", ""], ["b", ""], ["c", ""]])") +
        tie("a", "b", R"([["o1", 1], ["o3", 1]])") + tie("b", "c", R"([["o1", 1], ["o3", 1]])") + tie("a", "c", R"([["o1", 1], ["o3", 1]])")));
    const StrandedGraph G = stranded(K);
    s.require("a ring of three ties of one colour renewed at o2: the bubble's automorphisms are S₃", G.bubbles[1].automorphisms == 6 && G.bubbles[1].canonical,
              std::to_string(G.bubbles[1].automorphisms));
    s.require("the first and last occasions of the ring differ in type from the middle one", G.bubbles[0].type != G.bubbles[1].type && G.bubbles[2].type != G.bubbles[1].type);
    feynman(s, "the ring", K);
  }

  // The process on a history known by hand: a and b meet every year from 1800 to 1810, tied.
  {
    std::string rec = HEAD + being("a") + being("b");
    for (int y = 1800; y <= 1810; y++) rec += occasion("o" + std::to_string(y), std::to_string(y), R"([["a", ""], ["b", ""]])");
    rec += tie("a", "b", R"([["o1800", 1], ["o1810", 1]])");
    const HistoryComplex K = construct(read_record_text(rec));
    const Process P = process_of(K);
    double l0 = -1, l2 = -1, e2 = 0;
    for (const auto& R : P.rates) {
      if (R.k == 0) l0 = R.rate;
      if (R.k == 2) { l2 = R.rate; e2 = R.exposure; }
    }
    s.check("process: two beings ten years in the record, ten meetings: λ̂₂ = 1 a year", l2 - 1, 1e-12, "exposure " + std::to_string(e2));
    s.check("process: the one occasion creating them both, over ten years: λ̂₀ = 0.1", l0 - 0.1, 1e-12);
    s.check("process: tied the whole time, they meet at μ₁ = 1 a year", P.mu_tied - 1, 1e-12);
    s.require("process: no tie begins between continuing beings, and the tie ends once in ten meetings", P.open.trials == 0 && P.close.trials == 10 && P.close.successes == 1);
    double expected = 0;
    for (const auto& C : P.eras) expected += C.expected;
    s.check("process: the compensators of the eras add up to the occasions", expected - 11, 1e-9);
  }

  // Every canon.
  for (const auto& entry : std::filesystem::directory_iterator(HISTORY_CANONS)) {
    if (entry.path().extension() != ".jsonl") continue;
    const std::string name = "canon " + entry.path().filename().string();
    const HistoryComplex K = construct(read_record_file(entry.path().string()));
    feynman(s, name, K);
    const Process P = process_of(K);
    double expected = 0;
    std::size_t occasions = 0;
    for (const auto& C : P.eras) { expected += C.expected; occasions += C.occasions; }
    s.check(name + ": the eras' compensators add up to the occasions", expected - static_cast<double>(occasions), 1e-6,
            "χ² " + std::to_string(P.chi2) + " on " + std::to_string(P.dof) + "; binding " + std::to_string(P.binding));
    s.require(name + ": every occasion falls in an era", occasions == K.order.size());
  }
  return s.done();
}
