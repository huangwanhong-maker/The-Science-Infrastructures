// Checks of the record: a small record read whole, its references resolved and its occasions
// ordered in time across calendars; every fault the codebook names found, with where and why;
// the bereaved derived (codebook §8); and every canon of the domain read without an error.

#include <history/record/record.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <filesystem>
#include <string>

using namespace history;

namespace {

const char* SMALL = R"({"kind": "canon", "id": "small", "title": "A small record", "codebook": 1, "period": "1700/1720"}
{"kind": "source", "id": "s1", "cite": "A book, 1999"}
{"kind": "being", "id": "ann", "name": "Ann", "type": "person", "lived": "1680/1730"}
{"kind": "being", "id": "bob", "name": "Bob", "type": "person"}
{"kind": "being", "id": "cat", "name": "Cat", "type": "institution"}
{"kind": "occasion", "id": "o1710", "date": "1710-05-01", "type": "meeting", "label": "Second", "present": [["ann", "host"], ["bob", "guest"]], "sources": ["s1"], "confidence": 0.9}
{"kind": "occasion", "id": "o1705", "date": "1705-02-10", "calendar": "julian-march", "type": "letter", "label": "First, Old Style", "present": [["bob", "author"], ["ann", "recipient"], ["cat", "witness"]], "sources": ["s1"], "confidence": 1}
{"kind": "occasion", "id": "o1715", "date": "1715-06~", "type": "death", "label": "Ann dies", "present": [["ann", "died"]], "sources": ["s1"], "confidence": 1}
{"kind": "tie", "id": "ann~bob", "a": "ann", "b": "bob", "type": "friendship", "colour": "letters", "attested": [["o1710", -1], ["o1705", 1]], "ends": "o1715", "sources": ["s1"], "confidence": 0.8}
{"kind": "tie", "id": "bob~cat", "a": "bob", "b": "cat", "type": "patronage", "colour": "science", "attested": [["o1705", 1]], "sources": ["s1"], "confidence": 0.7}
)";

bool has(const Record& r, bool error, const std::string& where, const std::string& part) {
  return std::any_of(r.problems.begin(), r.problems.end(), [&](const Problem& p) {
    return p.error == error && p.where.find(where) != std::string::npos && p.message.find(part) != std::string::npos;
  });
}

// The small record with one line replaced.
std::string with(const std::string& from, const std::string& to) {
  std::string s = SMALL;
  const auto i = s.find(from);
  return i == std::string::npos ? s : s.replace(i, from.size(), to);
}

}  // namespace

int main() {
  calx::verification::Suite s;

  {
    const Record r = read_record_text(SMALL);
    s.require("a small record reads with no error", r.ok(), r.problems.empty() ? "" : r.problems.front().where + ": " + r.problems.front().message);
    s.require("its canon, sources, beings, occasions and ties", r.canon.id == "small" && r.sources.size() == 1 && r.beings.size() == 3 && r.occasions.size() == 3 && r.ties.size() == 2);
    // 10 February 1705 written Old Style is 10 February 1706, Julian: after 1 May 1710 is false,
    // before it is true; the order of time is o1705, o1710, o1715 whatever the file's order.
    s.require("occasions ordered in time across calendars", r.occasions[r.order[0]].id == "o1705" && r.occasions[r.order[1]].id == "o1710" && r.occasions[r.order[2]].id == "o1715");
    s.require("an Old Style date in January to March falls in the next year", r.occasions[*r.occasion("o1705")].date.first == jdn_julian(1706, 2, 10));
    const Tie& t = r.ties[*r.tie("ann~bob")];
    s.require("attestations sorted in time", t.attested.size() == 2 && r.occasions[t.attested[0].occasion].id == "o1705" && t.attested[0].sign == 1 && t.attested[1].sign == -1);
    s.require("a tie's end resolved", t.ends && r.occasions[*t.ends].id == "o1715");
    s.require("one attestation and no end: a warning that it makes no lens", has(r, false, "bob~cat", "makes no lens"));

    Record d = r;
    add_bereaved(d);
    const Occasion& death = d.occasions[*d.occasion("o1715")];
    const Presence* bob = presence(death, *d.being("bob"));
    s.require("the bereaved: Bob is present at Ann's death, derived", bob && bob->derived && bob->role == "bereaved");
    add_bereaved(d);
    s.require("the derivation is idempotent", d.occasions[*d.occasion("o1715")].present.size() == 2);
  }

  // Faults, each found where it is.
  auto fault = [&](const char* name, const std::string& text, bool error, const std::string& where, const std::string& part) {
    const Record r = read_record_text(text);
    s.require(name, has(r, error, where, part), r.problems.empty() ? "no problem found" : r.problems.front().where + ": " + r.problems.front().message);
  };
  fault("fault: a line that is not JSON", with(R"({"kind": "source", "id": "s1", "cite": "A book, 1999"})", "{not json"), true, "line 2", "not JSON");
  fault("fault: an identifier used twice", with(R"("id": "bob", "name": "Bob")", R"("id": "ann", "name": "Bob")"), true, "ann", "used twice");
  fault("fault: an unknown being present", with(R"(["cat", "witness"])", R"(["dan", "witness"])"), true, "o1705", "\"dan\" is not in the record");
  fault("fault: a tie attested where one of its beings is absent", with(R"([["o1705", 1]], "sources")", R"([["o1710", 1]], "sources")"), true, "bob~cat", "where cat is not present");
  fault("fault: a colour not in the codebook", with(R"("colour": "science")", R"("colour": "alchemy")"), true, "bob~cat", "not in the codebook");
  fault("fault: a confidence below 0.6", with(R"("confidence": 0.7})", R"("confidence": 0.5})"), true, "bob~cat", "0.6 to 1");
  fault("fault: a date that is not a date", with(R"("date": "1710-05-01")", R"("date": "1710-05-41")"), true, "o1710", "the date");
  fault("fault: a calendar not in the codebook", with(R"("calendar": "julian-march")", R"("calendar": "french-republican")"), true, "o1705", "calendar");
  fault("fault: a source not in the record", with(R"("sources": ["s1"], "confidence": 0.9})", R"("sources": ["s9"], "confidence": 0.9})"), true, "o1710", "\"s9\" is not in the record");
  fault("fault: a sign that is not ±1", with(R"([["o1710", -1], ["o1705", 1]])", R"([["o1710", 0], ["o1705", 1]])"), true, "ann~bob", "+1 or -1");
  fault("fault: an end before the last attestation", with(R"("ends": "o1715")", R"("ends": "o1705")"), true, "ann~bob", "before its last attestation");
  fault("fault: an end where only one is present, who neither dies nor departs", with(R"([["ann", "died"]])", R"([["ann", "host"]])"), true, "ann~bob", "neither dies nor departs");
  fault("fault: no canon line", with(R"({"kind": "canon", "id": "small", "title": "A small record", "codebook": 1, "period": "1700/1720"})", ""), true, "the record", "no canon line");
  fault("warning: an occasion type not in the codebook", with(R"("type": "meeting")", R"("type": "salon")"), false, "o1710", "not in the codebook");
  fault("warning: a source no assertion cites", with(R"({"kind": "source", "id": "s1", "cite": "A book, 1999"})", "{\"kind\": \"source\", \"id\": \"s1\", \"cite\": \"A book, 1999\"}\n{\"kind\": \"source\", \"id\": \"s2\", \"cite\": \"Another\"}"), false, "s2", "no assertion cites");

  // The canons of the domain.
  std::size_t canons = 0;
  for (const auto& entry : std::filesystem::directory_iterator(HISTORY_CANONS)) {
    if (entry.path().extension() != ".jsonl") continue;
    canons++;
    const Record r = read_record_file(entry.path().string());
    const std::string name = entry.path().filename().string();
    std::string first;
    for (const auto& p : r.problems)
      if (p.error) { first = p.where + ": " + p.message; break; }
    s.require("canon " + name + " reads with no error", r.ok(), first.empty() ? std::to_string(r.occasions.size()) + " occasions, " + std::to_string(r.ties.size()) + " ties" : first);
    s.require("canon " + name + " follows this codebook", r.canon.codebook == CODEBOOK);
  }
  s.require("the domain has a canon", canons > 0);

  // The join of the canons (the theory note, Definition 1.3): an id identifies across records
  // only where they agree about what it names.
  {
    const std::string dir = HISTORY_CANONS;
    auto path = [&](const char* canon) { return dir + "/" + canon + ".jsonl"; };
    const Record a = read_record_file(path("athens-470-322bc")), e = read_record_file(path("europe-1785-1828"));
    const Record x = read_record_file(path("europe-1871-1914")), v = read_record_file(path("vienna-1890-1938"));
    const Record J = join_records({path("athens-470-322bc"), path("europe-1785-1828"), path("europe-1871-1914"), path("vienna-1890-1938")});
    std::string first;
    for (const auto& p : J.problems)
      if (p.error) { first = p.where + ": " + p.message; break; }
    s.require("join: four canons join with no error", J.ok(), first);
    auto named = [&](const std::string& id) -> const Being* {
      for (const auto& b : J.beings) if (b.id == id) return &b;
      return nullptr;
    };
    const Being* great = named("alexander");
    const Being* tsar = named("alexander@europe-1785-1828");
    s.require("join: Alexander the Great and Alexander I of Russia, one id with two Wikidata items, stay two beings",
              great && tsar && great->wikidata == "Q8409" && tsar->wikidata == "Q15193");
    s.require("join: the Academy of Athens and the Academy of Fine Arts in Vienna, one id with two names, stay two beings",
              named("academy") && named("academy@vienna-1890-1938") && named("academy")->name != named("academy@vienna-1890-1938")->name);
    s.require("join: Franz Joseph, named alike in the alliances and Vienna, is one being", named("franz_joseph") && !named("franz_joseph@vienna-1890-1938"));
    std::size_t tsar_ties = 0, tsar_in_europe = 0;
    for (const auto& t : J.ties) tsar_ties += J.beings[t.a].id == "alexander@europe-1785-1828" || J.beings[t.b].id == "alexander@europe-1785-1828";
    for (const auto& t : e.ties) tsar_in_europe += e.beings[t.a].id == "alexander" || e.beings[t.b].id == "alexander";
    s.require("join: every tie of Alexander I follows him to his new name", tsar_ties == tsar_in_europe && tsar_ties > 0, std::to_string(tsar_ties) + " of " + std::to_string(tsar_in_europe));
    s.require("join: the beings are the union less the two shared",
              J.beings.size() == a.beings.size() + e.beings.size() + x.beings.size() + v.beings.size() - 2,
              std::to_string(J.beings.size()));
  }
  return s.done();
}
