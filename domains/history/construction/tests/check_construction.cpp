// Checks of the construction (theory note §2): on records small enough to know their topology,
// a tie between two meetings is a disc; three encounters no tie fills are a circle, told by
// its beings and occasions; three ties between the same two meetings close a sphere; a lens
// follows a being through the occasions it attends between its meetings; signs are carried
// between attestations. On every canon: the complex is a valid CW complex, its edges and faces
// are the worldlines' segments and the ties' lenses, and χ is the alternating sum of the Betti
// numbers.

#include <history/construction/complex.hpp>
#include <history/construction/loops.hpp>
#include <calx/verification/check.hpp>

#include <filesystem>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using namespace history;
namespace mt = morphe::topology;

namespace {

const std::string HEAD = R"({"kind": "canon", "id": "test", "title": "test", "codebook": 1}
{"kind": "source", "id": "s", "cite": "s"}
)";

std::string being(const std::string& id) { return R"({"kind": "being", "id": ")" + id + R"(", "name": ")" + id + R"(", "type": "person"})" + "\n"; }
std::string occasion(const std::string& id, const std::string& date, const std::string& present) {
  return R"({"kind": "occasion", "id": ")" + id + R"(", "date": ")" + date + R"(", "type": "meeting", "present": )" + present + R"(, "sources": ["s"], "confidence": 1})" + "\n";
}
std::string tie(const std::string& a, const std::string& b, const std::string& attested) {
  return R"({"kind": "tie", "id": ")" + a + "~" + b + R"(", "a": ")" + a + R"(", "b": ")" + b + R"(", "type": "friendship", "colour": "letters", "attested": )" + attested +
         R"(, "sources": ["s"], "confidence": 1})" + "\n";
}

std::vector<std::size_t> betti_of(const HistoryComplex& K) {
  std::vector<std::size_t> b;
  for (const auto& h : mt::homology(K.cw, false)) b.push_back(h.rank);
  b.resize(3, 0);
  return b;
}

std::string text(const std::vector<std::size_t>& v) {
  std::string s;
  for (std::size_t x : v) s += (s.empty() ? "" : " ") + std::to_string(x);
  return s;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // A tie between two meetings: one lens, a disc.
  {
    const Record r = read_record_text(HEAD + being("a") + being("b") + occasion("p", "1800", R"([["a", ""], ["b", ""]])") +
                                      occasion("q", "1801", R"([["a", ""], ["b", ""]])") + tie("a", "b", R"([["p", 1], ["q", 1]])"));
    const HistoryComplex K = construct(r);
    s.require("a tie between two meetings: V, E, F = 2, 2, 1", K.cw.vertices == 2 && K.cw.edges.size() == 2 && K.cw.faces.size() == 1);
    s.require("a tie between two meetings: a disc (Betti 1 0 0)", betti_of(K) == std::vector<std::size_t>{1, 0, 0}, text(betti_of(K)));
    s.require("its lens's word is a[p, q] · b[p, q]⁻¹", K.cw.faces[0].boundary.size() == 2 && K.cw.faces[0].boundary[0].sign == 1 && K.cw.faces[0].boundary[1].sign == -1);
  }

  // Three encounters that no tie fills: a circle.
  {
    const Record r = read_record_text(HEAD + being("a") + being("b") + being("c") + occasion("o1", "1801", R"([["a", ""], ["b", ""]])") +
                                      occasion("o2", "1802", R"([["b", ""], ["c", ""]])") + occasion("o3", "1803", R"([["a", ""], ["c", ""]])"));
    const HistoryComplex K = construct(r);
    s.require("three encounters, no tie: a circle (Betti 1 1 0)", betti_of(K) == std::vector<std::size_t>{1, 1, 0}, text(betti_of(K)));
    const auto H = mt::homology(K.cw, true);
    const auto L = loops(K, H);
    std::set<std::size_t> beings(L.empty() ? std::set<std::size_t>{} : std::set<std::size_t>(L[0].beings.begin(), L[0].beings.end()));
    std::set<std::size_t> occasions(L.empty() ? std::set<std::size_t>{} : std::set<std::size_t>(L[0].occasions.begin(), L[0].occasions.end()));
    s.require("its loop follows a, b and c", L.size() == 1 && beings == std::set<std::size_t>{0, 1, 2});
    s.require("and passes between them at the three occasions", occasions == std::set<std::size_t>{0, 1, 2});
    s.require("its walk is closed: three steps", L.size() == 1 && L[0].steps.size() == 3);
  }

  // Three ties between the same two meetings: a sphere.
  {
    const std::string all = R"([["a", ""], ["b", ""], ["c", ""]])";
    const Record r = read_record_text(HEAD + being("a") + being("b") + being("c") + occasion("p", "1800", all) + occasion("q", "1801", all) +
                                      tie("a", "b", R"([["p", 1], ["q", 1]])") + tie("b", "c", R"([["p", 1], ["q", -1]])") + tie("a", "c", R"([["p", -1], ["q", 1]])"));
    const HistoryComplex K = construct(r);
    s.require("three ties between two meetings: a sphere (Betti 1 0 1)", betti_of(K) == std::vector<std::size_t>{1, 0, 1}, text(betti_of(K)));
    s.require("its Euler characteristic is 2", K.cw.euler_characteristic() == 2);
  }

  // A lens follows a being through what it attends between the meetings; carried signs.
  {
    const std::string ab = R"([["a", ""], ["b", ""]])";
    const Record r = read_record_text(HEAD + being("a") + being("b") + being("d") + occasion("p", "1800", ab) + occasion("x", "1800-06", R"([["a", ""], ["d", ""]])") +
                                      occasion("q", "1801", ab) + occasion("r", "1802", ab) + tie("a", "b", R"([["p", 1], ["r", -1]])"));
    const HistoryComplex K = construct(r);
    const std::size_t t = *K.record.tie("a~b");
    s.require("a tie's rungs are its meetings within its life: p, q, r", K.rungs[t].size() == 3);
    s.require("at q its sign is carried from p, one lens since", !K.rungs[t][1].attested && K.rungs[t][1].sign == 1 && K.rungs[t][1].lenses_since == 1 && K.rungs[t][1].since == *K.record.occasion("p"));
    s.require("two lenses; the second flips", K.lenses.size() == 2 && !K.flip(0) && K.flip(1));
    s.require("the first lens follows a through x: a word of three letters", K.cw.faces[0].boundary.size() == 3);
    s.require("two discs pinched at q (Betti 1 0 0)", betti_of(K) == std::vector<std::size_t>{1, 0, 0}, text(betti_of(K)));
  }

  // Overlapping occasions of one being: their order rests on the file, and the construction says so.
  {
    const Record r = read_record_text(HEAD + being("a") + being("b") + occasion("congress", "1808-09-27/1808-10-14", R"([["a", ""], ["b", ""]])") +
                                      occasion("audience", "1808-10-02", R"([["a", ""]])"));
    const HistoryComplex K = construct(r);
    s.require("overlapping occasions of one being: a warning", !K.problems.empty() && K.problems[0].message.find("rests on the file") != std::string::npos);
  }

  // Every canon.
  for (const auto& entry : std::filesystem::directory_iterator(HISTORY_CANONS)) {
    if (entry.path().extension() != ".jsonl") continue;
    const std::string name = entry.path().filename().string();
    const HistoryComplex K = construct(read_record_file(entry.path().string()));
    bool valid = true;
    std::string why;
    try {
      K.cw.validate();
    } catch (const std::exception& e) {
      valid = false;
      why = e.what();
    }
    s.require("canon " + name + ": a valid CW complex", valid, why);
    std::size_t edges = 0, faces = 0;
    for (const auto& w : K.worldline) edges += w.empty() ? 0 : w.size() - 1;
    for (const auto& g : K.rungs) faces += g.empty() ? 0 : g.size() - 1;
    s.require("canon " + name + ": its edges are the worldlines' segments", K.cw.edges.size() == edges);
    s.require("canon " + name + ": its faces are the ties' lenses", K.cw.faces.size() == faces);
    const auto b = betti_of(K);
    const long long alternating = static_cast<long long>(b[0]) - static_cast<long long>(b[1]) + static_cast<long long>(b[2]);
    s.require("canon " + name + ": χ is the alternating sum of the Betti numbers", K.cw.euler_characteristic() == alternating,
              "χ = " + std::to_string(K.cw.euler_characteristic()) + ", Betti " + text(b));
    bool based = true;
    for (std::size_t f = 0; f < K.lenses.size(); f++) based = based && K.cw.faces[f].base == K.vertex[K.lenses[f].from];
    s.require("canon " + name + ": each lens's word starts at its first meeting", based);
  }
  return s.done();
}
