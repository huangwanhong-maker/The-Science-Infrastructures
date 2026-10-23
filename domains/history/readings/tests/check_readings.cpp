// Checks of the readings (theory note §3): on a closed shell, the Bianchi identity; frustrated
// triangles at occasions and in the network of the moment; carried signs weighed by their
// memory; the persistence of flips recovered from ties drawn with a known rate; and, on every
// canon, the identity on every shell and the invariance of every reading under switching a
// being (reversing the signs of all its ties, a gauge transformation; Zaslavsky 1982).

#include <history/readings/colours.hpp>
#include <history/readings/signed.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <cstdint>
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

// Switching a being: every sign of its ties reversed.
Record switched(Record r, std::size_t being) {
  for (auto& t : r.ties)
    if (t.a == being || t.b == being)
      for (auto& a : t.attested) a.sign = -a.sign;
  return r;
}

// Everything the readings say that a gauge transformation must leave alone.
std::string invariants(const HistoryComplex& K) {
  const Memory m = memory(persistence(K));
  std::string s;
  for (std::size_t f : flips(K)) s += "f" + std::to_string(f) + " ";
  for (const auto& T : triangles(K, m)) s += "t" + std::to_string(T.product) + " ";
  for (const auto& S : shells(K)) s += "s" + std::to_string(S.product_from) + std::to_string(S.product_to) + " ";
  for (const auto& Y : moments_by_year(K, m)) s += "y" + std::to_string(Y.frustrated) + " ";
  for (const auto& P : persistence(K)) s += "p" + std::to_string(P.flips) + "/" + std::to_string(P.lenses) + " ";
  return s;
}

// A small, seeded generator (SplitMix64), so that the check needs nothing beyond the domain.
struct SplitMix {
  std::uint64_t x;
  double uniform() {
    std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return static_cast<double>((z ^ (z >> 31)) >> 11) * 0x1.0p-53;
  }
};

}  // namespace

int main() {
  calx::verification::Suite s;
  const std::string all = R"([["a", ""], ["b", ""], ["c", ""]])";

  // A closed shell: ab (+, +), bc (+, −), ac (−, +).
  {
    const Record r = read_record_text(HEAD + being("a") + being("b") + being("c") + occasion("p", "1800", all) + occasion("q", "1801", all) +
                                      tie("a", "b", R"([["p", 1], ["q", 1]])") + tie("b", "c", R"([["p", 1], ["q", -1]])") + tie("a", "c", R"([["p", -1], ["q", 1]])"));
    const HistoryComplex K = construct(r);
    const auto S = shells(K);
    s.require("one closed shell", S.size() == 1);
    if (!S.empty())
      s.require("the Bianchi identity: product at q = product at p · (−1)^flips", S[0].product_from == -1 && S[0].product_to == -1 && S[0].flips == 2);
    const auto P = persistence(K);
    s.require("persistence: 2 flips in 3 attested lenses, no β (memory < 0)", P.back().lenses == 3 && P.back().flips == 2 && !P.back().beta && std::fabs(P.back().memory + 1.0 / 3) < 1e-15);
    const auto T = triangles(K, memory(P));
    s.require("frustrated triangles at both meetings", T.size() == 2 && T[0].product == -1 && T[1].product == -1 && T[0].attested);
    const Moment M = moment(K, jdn_gregorian(1800, 6, 1), memory(P));
    s.require("the moment between them: three links, one frustrated triangle", M.links.size() == 3 && M.triangles.size() == 1 && M.frustrated == 1 && std::fabs(M.expected - 1) < 1e-15);
    const auto Y = moments_by_year(K, memory(P));
    s.require("the moments by year: 1800 and 1801", Y.size() == 2 && Y[0].year == 1800 && Y[1].year == 1801);
    const Colours C = colours_of(K);
    s.require("colours: three lenses of letters, no meeting of colours", C.lenses.at("letters") == 3 && C.meetings.empty() && C.of_being[0] == "letters");
    s.require("switching b changes no reading", invariants(K) == invariants(construct(switched(r, 1))));
  }

  // A sign carried between attestations is weighed by its memory.
  {
    const std::string ab = R"([["a", ""], ["b", ""]])";
    const Record r = read_record_text(HEAD + being("a") + being("b") + occasion("p", "1800", ab) + occasion("q", "1801", ab) + occasion("r", "1802", ab) +
                                      tie("a", "b", R"([["p", 1], ["r", -1]])"));
    const HistoryComplex K = construct(r);
    Memory m;
    m.by_colour["letters"] = 0.5;
    const Rung& q = K.rungs[0][1];
    s.check("a sign carried one lens with memory 0.5 is expected 0.5", expected(q, m.of("letters")) - 0.5, 1e-15);
    s.require("no lens has both ends attested, so the persistence is not estimated", persistence(K).back().lenses == 0);
  }

  // The persistence recovered: ties drawn with flips at a known rate p0 = 0.2 per lens.
  {
    SplitMix rng{20261024};
    std::string text = HEAD;
    const int pairs = 40, meetings = 60;
    for (int i = 0; i < 2 * pairs; i++) text += being("x" + std::to_string(i));
    std::vector<std::string> attested(pairs);
    std::vector<int> sign(pairs, 1);
    for (int k = 0; k < meetings; k++) {
      std::string present = "[";
      for (int i = 0; i < 2 * pairs; i++) present += std::string(i ? ", " : "") + "[\"x" + std::to_string(i) + "\", \"\"]";
      text += occasion("m" + std::to_string(k), std::to_string(1800 + k), present + "]");
      for (int i = 0; i < pairs; i++) {
        if (k && rng.uniform() < 0.2) sign[i] = -sign[i];
        attested[i] += std::string(k ? ", " : "") + "[\"m" + std::to_string(k) + "\", " + std::to_string(sign[i]) + "]";
      }
    }
    for (int i = 0; i < pairs; i++) text += tie("x" + std::to_string(2 * i), "x" + std::to_string(2 * i + 1), "[" + attested[i] + "]", "art");
    const HistoryComplex K = construct(read_record_text(text));
    const auto P = persistence(K);
    const double n = static_cast<double>(P.back().lenses), sigma = std::sqrt(0.2 * 0.8 / n);
    s.require("persistence: 40 ties over 60 meetings, 2360 lenses", P.back().lenses == 2360, std::to_string(P.back().lenses));
    s.check("persistence: the rate of flips recovered within 3σ (p0 = 0.2)", (P.back().p - 0.2) / sigma, 3.0, "σ");
    s.check("persistence: β = artanh(1 − 2p) near artanh 0.6 = 0.6931", P.back().beta.value_or(0) - std::atanh(0.6), 0.1);
  }

  // Every canon: the Bianchi identity on every shell, and switching each being changes nothing.
  for (const auto& entry : std::filesystem::directory_iterator(HISTORY_CANONS)) {
    if (entry.path().extension() != ".jsonl") continue;
    const std::string name = entry.path().filename().string();
    const Record r = read_record_file(entry.path().string());
    const HistoryComplex K = construct(r);
    bool bianchi = true;
    const auto S = shells(K);
    for (const auto& x : S) bianchi = bianchi && x.product_to == x.product_from * (x.flips % 2 ? -1 : 1);
    s.require("canon " + name + ": the Bianchi identity on every shell", bianchi, std::to_string(S.size()) + " shells");
    const std::string base = invariants(K);
    std::string moved;
    for (std::size_t b = 0; b < r.beings.size() && moved.empty(); b++)
      if (invariants(construct(switched(r, b))) != base) moved = r.beings[b].id;
    s.require("canon " + name + ": switching any one being changes no reading", moved.empty(), moved.empty() ? std::to_string(r.beings.size()) + " beings" : "switching " + moved + " does");
    // The year moments end where the record's years end: an occasion written in the calendar of
    // the historians' convention for its time (Julian before 1582, Gregorian after) enters the
    // moments in its own year, never at the end of the year before.
    const auto Y = moments_by_year(K, memory(persistence(K)));
    std::string early;
    std::size_t counted = 0;
    for (const auto& o : r.occasions) {
      const bool conventional = (o.date.calendar == Calendar::julian && o.date.first < GREGORIAN_REFORM) ||
                                (o.date.calendar == Calendar::gregorian && o.date.first >= GREGORIAN_REFORM);
      if (!conventional) continue;
      counted++;
      const long y = historical_of(o.date.first).year;
      for (const auto& m : Y)
        if (m.year == y - 1 && m.jdn >= o.date.first && early.empty()) early = o.id;
    }
    s.require("canon " + name + ": each occasion enters the year moments in its own year", early.empty(), early.empty() ? std::to_string(counted) + " occasions" : early + " enters a year early");
  }
  return s.done();
}
