// The record of a history (codebook §1-3, 5-9): assertions read from a file of JSON lines,
// their references resolved, checked against the codebook, and ordered in time.
//
// A record holds what sources say, never facts: occasions with the beings present, and
// dyadic ties with their signs where they are attested. Reading never throws on content:
// every fault becomes a Problem, an error or a warning, with where it is and why, so that a
// record can be shown and corrected whole. Errors make a record unfit for the construction.

#pragma once

#include <history/record/chronology.hpp>
#include <history/record/json.hpp>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <istream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace history {

// The vocabularies of the codebook, version 1.
inline constexpr int CODEBOOK = 1;
inline const std::vector<std::string>& colours() {
  static const std::vector<std::string> c{"science", "art", "music", "letters", "philosophy", "international-relations"};
  return c;
}
inline const std::vector<std::string>& being_types() {
  static const std::vector<std::string> t{"person", "group", "institution", "state", "work"};
  return t;
}
inline const std::vector<std::string>& occasion_types() {
  static const std::vector<std::string> t{"meeting", "letter", "lesson", "publication", "dedication", "performance", "exhibition", "commission",
                                          "appointment", "treaty", "congress", "battle", "ceremony", "journey", "death", "other"};
  return t;
}
inline const std::vector<std::string>& tie_types() {
  static const std::vector<std::string> t{"alliance", "diplomacy", "service", "teaching", "patronage", "commission", "friendship",
                                          "collaboration", "admiration", "rivalry", "kinship", "marriage"};
  return t;
}

struct Provenance {
  std::string drafted_by, checked_by, note;
};

struct Canon {
  std::string id, title, about, period, status, drafted_by, drafted;
  int codebook = 0;
};

struct Source {
  std::string id, cite;
};

struct Being : Provenance {
  std::string id, name, type, lived, wikidata;
};

struct Presence {
  std::size_t being = 0;
  std::string role;
  bool derived = false;  // added by a derivation of the codebook (§8), not written in the record
};

struct Occasion : Provenance {
  std::string id, label, place, type;
  Date date;
  std::vector<Presence> present;
  std::vector<std::string> sources;
  double confidence = 0;
  std::size_t line = 0;  // where it is written, for the order of ties in time
};

struct Attestation {
  std::size_t occasion = 0;
  int sign = 1;
};

struct Tie : Provenance {
  std::string id, type, colour;
  std::size_t a = 0, b = 0;
  std::vector<Attestation> attested;  // in the order of time, once the record is read
  std::optional<std::size_t> ends;
  std::vector<std::string> sources;
  double confidence = 0;
};

struct Problem {
  bool error = true;
  std::string where, message;
};

struct Record {
  Canon canon;
  std::vector<Source> sources;
  std::vector<Being> beings;
  std::vector<Occasion> occasions;
  std::vector<Tie> ties;
  std::vector<Problem> problems;
  std::vector<std::size_t> order;  // the occasions in time: by first day, last day, then line

  std::size_t errors() const { return static_cast<std::size_t>(std::count_if(problems.begin(), problems.end(), [](const Problem& p) { return p.error; })); }
  bool ok() const { return errors() == 0; }

  template <class T>
  static std::optional<std::size_t> find(const std::vector<T>& v, std::string_view id) {
    for (std::size_t i = 0; i < v.size(); i++)
      if (v[i].id == id) return i;
    return std::nullopt;
  }
  std::optional<std::size_t> being(std::string_view id) const { return find(beings, id); }
  std::optional<std::size_t> occasion(std::string_view id) const { return find(occasions, id); }
  std::optional<std::size_t> tie(std::string_view id) const { return find(ties, id); }

  // The rank of each occasion in time.
  std::vector<std::size_t> ranks() const {
    std::vector<std::size_t> r(occasions.size());
    for (std::size_t k = 0; k < order.size(); k++) r[order[k]] = k;
    return r;
  }
};

inline bool is_present(const Occasion& o, std::size_t being) {
  return std::any_of(o.present.begin(), o.present.end(), [&](const Presence& p) { return p.being == being; });
}

inline const Presence* presence(const Occasion& o, std::size_t being) {
  for (const auto& p : o.present)
    if (p.being == being) return &p;
  return nullptr;
}

// Whether a being dies or departs at an occasion, so that its ties may end there (codebook §8).
inline bool ends_there(const Occasion& o, std::size_t being) {
  const Presence* p = presence(o, being);
  return p && (p->role == "died" || p->role == "departing");
}

namespace detail {

inline std::string text(const Json& line, std::string_view key) {
  const Json* v = line.find(key);
  return v && v->is_string() ? v->as_string() : std::string();
}

inline std::vector<std::string> texts(const Json& line, std::string_view key) {
  std::vector<std::string> out;
  if (const Json* v = line.find(key); v && v->is_array())
    for (const auto& x : v->as_array())
      if (x.is_string()) out.push_back(x.as_string());
  return out;
}

inline bool known(const std::vector<std::string>& vocabulary, const std::string& word) {
  return std::find(vocabulary.begin(), vocabulary.end(), word) != vocabulary.end();
}

template <class T>
void provenance(T& x, const Json& line) {
  x.drafted_by = text(line, "drafted_by");
  x.checked_by = text(line, "checked_by");
  x.note = text(line, "note");
}

}  // namespace detail

// Reads a record from JSON lines. Never throws on content; see Record::problems.
inline Record read_record(std::istream& in) {
  Record r;
  auto error = [&](std::string where, std::string message) { r.problems.push_back({true, std::move(where), std::move(message)}); };
  auto warning = [&](std::string where, std::string message) { r.problems.push_back({false, std::move(where), std::move(message)}); };
  using detail::text;
  using detail::texts;

  // The references are written by identifier, and resolved once every line is read.
  struct RawPresence { std::string being, role; };
  struct RawAttestation { std::string occasion; double sign; };
  std::vector<std::vector<RawPresence>> raw_present;
  struct RawTie { std::string a, b, ends; std::vector<RawAttestation> attested; std::string where; };
  std::vector<RawTie> raw_ties;
  std::set<std::string> ids;
  bool have_canon = false;
  std::size_t n = 0, first_content = 0;

  std::string s;
  while (std::getline(in, s)) {
    n++;
    if (!s.empty() && s.back() == '\r') s.pop_back();
    if (s.find_first_not_of(" \t") == std::string::npos) continue;
    if (!first_content) first_content = n;
    const std::string where = "line " + std::to_string(n);
    Json line;
    try {
      line = Json::parse(s);
    } catch (const std::exception& e) {
      error(where, std::string("not JSON: ") + e.what());
      continue;
    }
    if (!line.is_object()) { error(where, "a line is a JSON object"); continue; }
    const std::string kind = text(line, "kind"), id = text(line, "id");
    if (kind.empty()) { error(where, "no \"kind\""); continue; }
    if (id.empty()) { error(where, "no \"id\""); continue; }
    const std::string at = id + " (" + where + ")";
    if (kind != "canon" && !ids.insert(id).second) { error(at, "the identifier is used twice"); continue; }
    auto confidence = [&](const Json& l) {
      const Json* c = l.find("confidence");
      if (!c) { error(at, "no confidence (codebook §7)"); return 0.0; }
      if (!c->is_number() || c->as_number() < 0.6 || c->as_number() > 1) error(at, "the confidence is a number from 0.6 to 1 (codebook §7)");
      return c->is_number() ? c->as_number() : 0.0;
    };

    if (kind == "canon") {
      if (have_canon) { error(at, "a record has one canon line"); continue; }
      have_canon = true;
      if (n != first_content) warning(at, "the canon line comes first (codebook §2)");
      Canon& c = r.canon;
      c.id = id;
      c.title = text(line, "title");
      c.about = text(line, "about");
      c.period = text(line, "period");
      c.status = text(line, "status");
      c.drafted_by = text(line, "drafted_by");
      c.drafted = text(line, "drafted");
      if (const Json* v = line.find("codebook"); v && v->is_number()) c.codebook = static_cast<int>(v->as_number());
      if (c.codebook != CODEBOOK) warning(at, "coded by codebook " + std::to_string(c.codebook) + "; this is codebook " + std::to_string(CODEBOOK));
    } else if (kind == "source") {
      Source x{id, text(line, "cite")};
      if (x.cite.empty()) error(at, "a source has its \"cite\"");
      r.sources.push_back(std::move(x));
    } else if (kind == "being") {
      Being x;
      x.id = id;
      x.name = text(line, "name");
      x.type = text(line, "type");
      x.lived = text(line, "lived");
      x.wikidata = text(line, "wikidata");
      detail::provenance(x, line);
      if (x.name.empty()) error(at, "a being has its \"name\"");
      if (!detail::known(being_types(), x.type)) warning(at, "the type \"" + x.type + "\" is not in the codebook (§3)");
      if (!x.lived.empty()) {
        try {
          parse_date(x.lived);
        } catch (const std::exception& e) {
          error(at, std::string("\"lived\": ") + e.what());
        }
      }
      r.beings.push_back(std::move(x));
    } else if (kind == "occasion") {
      Occasion x;
      x.id = id;
      x.label = text(line, "label");
      x.place = text(line, "place");
      x.type = text(line, "type");
      x.sources = texts(line, "sources");
      x.line = n;
      x.confidence = confidence(line);
      detail::provenance(x, line);
      if (!detail::known(occasion_types(), x.type)) warning(at, "the type \"" + x.type + "\" is not in the codebook (§5)");
      const auto cal = calendar_of(text(line, "calendar"));
      if (!cal) error(at, "the calendar \"" + text(line, "calendar") + "\" is not gregorian, julian or julian-march");
      try {
        x.date = parse_date(text(line, "date"), cal.value_or(Calendar::gregorian));
      } catch (const std::exception& e) {
        error(at, std::string("the date: ") + e.what());
        x.date.text = text(line, "date");
      }
      std::vector<RawPresence> present;
      if (const Json* p = line.find("present"); p && p->is_array()) {
        for (const auto& e : p->as_array()) {
          if (e.is_array() && e.size() >= 1 && e.as_array()[0].is_string())
            present.push_back({e.as_array()[0].as_string(), e.size() > 1 && e.as_array()[1].is_string() ? e.as_array()[1].as_string() : ""});
          else
            error(at, "each presence is [being, role]");
        }
      }
      if (present.empty()) error(at, "an occasion has beings present (codebook §5)");
      raw_present.push_back(std::move(present));
      r.occasions.push_back(std::move(x));
    } else if (kind == "tie") {
      Tie x;
      x.id = id;
      x.type = text(line, "type");
      x.colour = text(line, "colour");
      x.sources = texts(line, "sources");
      x.confidence = confidence(line);
      detail::provenance(x, line);
      if (!detail::known(colours(), x.colour)) error(at, "the colour \"" + x.colour + "\" is not in the codebook (§6)");
      if (!detail::known(tie_types(), x.type)) warning(at, "the type \"" + x.type + "\" is not in the codebook (§6)");
      RawTie raw{text(line, "a"), text(line, "b"), text(line, "ends"), {}, at};
      if (const Json* p = line.find("attested"); p && p->is_array()) {
        for (const auto& e : p->as_array()) {
          if (e.is_array() && e.size() == 2 && e.as_array()[0].is_string() && e.as_array()[1].is_number())
            raw.attested.push_back({e.as_array()[0].as_string(), e.as_array()[1].as_number()});
          else
            error(at, "each attestation is [occasion, sign]");
        }
      }
      if (raw.attested.empty()) error(at, "a tie is attested at one occasion at least");
      raw_ties.push_back(std::move(raw));
      r.ties.push_back(std::move(x));
    } else {
      error(at, "the kind \"" + kind + "\" is not canon, source, being, occasion or tie");
    }
  }
  if (!have_canon) error("the record", "no canon line (codebook §2)");

  // Sources cited.
  std::set<std::string> source_ids, cited;
  for (const auto& x : r.sources) source_ids.insert(x.id);
  auto cite = [&](const std::string& at, const std::vector<std::string>& ss) {
    for (const auto& s2 : ss) {
      cited.insert(s2);
      if (!source_ids.count(s2)) error(at, "the source \"" + s2 + "\" is not in the record");
    }
  };

  // Presences.
  for (std::size_t i = 0; i < r.occasions.size(); i++) {
    Occasion& o = r.occasions[i];
    cite(o.id, o.sources);
    for (const auto& p : raw_present[i]) {
      const auto b = r.being(p.being);
      if (!b) { error(o.id, "the being \"" + p.being + "\" is not in the record"); continue; }
      if (is_present(o, *b)) { error(o.id, "\"" + p.being + "\" is present twice"); continue; }
      o.present.push_back({*b, p.role, false});
    }
  }

  // The order of time: by first day, last day, then where written (codebook §4).
  r.order.resize(r.occasions.size());
  for (std::size_t i = 0; i < r.order.size(); i++) r.order[i] = i;
  std::stable_sort(r.order.begin(), r.order.end(), [&](std::size_t x, std::size_t y) {
    const auto &a = r.occasions[x].date, &b = r.occasions[y].date;
    if (a.first != b.first) return a.first < b.first;
    if (a.last != b.last) return a.last < b.last;
    return r.occasions[x].line < r.occasions[y].line;
  });
  const auto rank = r.ranks();

  // Ties.
  for (std::size_t i = 0; i < r.ties.size(); i++) {
    Tie& t = r.ties[i];
    const RawTie& raw = raw_ties[i];
    cite(t.id, t.sources);
    const auto a = r.being(raw.a), b = r.being(raw.b);
    if (!a) error(t.id, "the being \"" + raw.a + "\" is not in the record");
    if (!b) error(t.id, "the being \"" + raw.b + "\" is not in the record");
    if (!a || !b) continue;
    if (*a == *b) { error(t.id, "a tie joins two beings"); continue; }
    t.a = *a;
    t.b = *b;
    std::set<std::size_t> seen;
    for (const auto& e : raw.attested) {
      const auto o = r.occasion(e.occasion);
      if (!o) { error(t.id, "the occasion \"" + e.occasion + "\" is not in the record"); continue; }
      if (e.sign != 1 && e.sign != -1) { error(t.id, "the sign at \"" + e.occasion + "\" is +1 or -1"); continue; }
      if (!seen.insert(*o).second) { error(t.id, "attested twice at \"" + e.occasion + "\""); continue; }
      if (!is_present(r.occasions[*o], t.a) || !is_present(r.occasions[*o], t.b)) {
        error(t.id, "attested at \"" + e.occasion + "\", where " + (is_present(r.occasions[*o], t.a) ? r.beings[t.b].id : r.beings[t.a].id) + " is not present");
        continue;
      }
      t.attested.push_back({*o, static_cast<int>(e.sign)});
    }
    std::sort(t.attested.begin(), t.attested.end(), [&](const Attestation& x, const Attestation& y) { return rank[x.occasion] < rank[y.occasion]; });
    if (!raw.ends.empty()) {
      const auto o = r.occasion(raw.ends);
      if (!o) error(t.id, "it ends at \"" + raw.ends + "\", which is not in the record");
      else if (!is_present(r.occasions[*o], t.a) && !is_present(r.occasions[*o], t.b))
        error(t.id, "it ends at \"" + raw.ends + "\", where neither of its beings is present");
      else if (!t.attested.empty() && rank[*o] < rank[t.attested.back().occasion])
        error(t.id, "it ends at \"" + raw.ends + "\", before its last attestation");
      else if (const Occasion& e = r.occasions[*o]; is_present(e, t.a) != is_present(e, t.b) && !ends_there(e, is_present(e, t.a) ? t.a : t.b))
        error(t.id, "it ends at \"" + raw.ends + "\", where only one of its beings is present, and neither dies nor departs there (codebook §8)");
      else
        t.ends = *o;
    }
    if (t.attested.size() == 1 && !t.ends) warning(t.id, "attested at one meeting and never ended: it makes no lens");
  }

  // What is written and never used.
  std::vector<std::size_t> occasions_of(r.beings.size(), 0);
  for (const auto& o : r.occasions)
    for (const auto& p : o.present) occasions_of[p.being]++;
  for (std::size_t b = 0; b < r.beings.size(); b++)
    if (!occasions_of[b]) warning(r.beings[b].id, "a being present at no occasion");
  for (const auto& x : r.sources)
    if (!cited.count(x.id)) warning(x.id, "a source no assertion cites");
  return r;
}

inline Record read_record_text(std::string_view text) {
  std::istringstream in{std::string(text)};
  return read_record(in);
}

// Reads the record at a path; the only file the domain reads, and only when it is passed.
inline Record read_record_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    Record r;
    r.problems.push_back({true, path, "cannot be read"});
    return r;
  }
  return read_record(in);
}

// Joins several records into one, the integrated record (ADR 0041; the theory note, Definition
// 1.3): the union of their assertions, glued along what they share. An id is a name within one
// record, and two records share what an id names only where they agree about it:
//
// - a being, where both give the same Wikidata item, or, where either gives none, the same name;
// - a source, where both give the same citation;
// - an occasion, where both give the same date; a tie, where it joins the same beings.
//
// A shared occasion has the union of the beings present and of the sources; a shared tie the
// union of its attestations. Where two records use one id for different things (Alexander the
// Great and Alexander I of Russia), the later record's is renamed id@canon throughout that record,
// and the join reports it. The result is read as one record, under a canon line of its own.
inline Record join_records(const std::vector<std::string>& paths, const std::string& id = "integrated", const std::string& title = "All canons, joined") {
  using detail::text;
  std::vector<std::string> order;            // "kind|id", in the order first met
  std::map<std::string, Json> lines;
  std::vector<Problem> problems;
  std::vector<std::string> parts;
  for (const auto& path : paths) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { problems.push_back({true, path, "cannot be read"}); continue; }
    std::vector<Json> mine;
    std::string canon = path, s;
    while (std::getline(in, s)) {
      if (!s.empty() && s.back() == '\r') s.pop_back();
      if (s.find_first_not_of(" \t") == std::string::npos) continue;
      Json line;
      try { line = Json::parse(s); } catch (const std::exception& e) { problems.push_back({true, path, std::string("not JSON: ") + e.what()}); continue; }
      if (text(line, "kind") == "canon") {
        canon = text(line, "id");
        parts.push_back(text(line, "title").empty() ? canon : text(line, "title"));
        continue;
      }
      mine.push_back(std::move(line));
    }
    // The ids this record uses for something other than what the records before it name.
    std::map<std::string, std::string> renamed;  // "kind|id" -> its name in the join
    auto rename = [&](const std::string& kind, const std::string& lid, const std::string& what) {
      renamed[kind + "|" + lid] = lid + "@" + canon;
      problems.push_back({false, lid, "names " + what + " in two records; in " + canon + " it is " + lid + "@" + canon});
    };
    auto same_being = [](const Json& x, const Json& y) {
      const std::string qx = text(x, "wikidata"), qy = text(y, "wikidata");
      return !qx.empty() && !qy.empty() ? qx == qy : text(x, "name") == text(y, "name");
    };
    for (const auto& line : mine) {
      const std::string kind = text(line, "kind"), lid = text(line, "id");
      const auto it = lines.find(kind + "|" + lid);
      if (it == lines.end()) continue;
      if (kind == "being" && !same_being(it->second, line)) rename(kind, lid, "two beings");
      if (kind == "source" && text(it->second, "cite") != text(line, "cite")) rename(kind, lid, "two sources");
      if (kind == "occasion" && text(it->second, "date") != text(line, "date")) rename(kind, lid, "two occasions");
    }
    // Carry the new names through every reference of this record.
    auto ref = [&](const char* kind, const Json& x) {
      if (!x.is_string()) return x;
      const auto r = renamed.find(std::string(kind) + "|" + x.as_string());
      return r == renamed.end() ? x : Json(r->second);
    };
    auto refs = [&](Json& line, const char* field, const char* kind, bool first_of_pair) {
      const Json* v = line.find(field);
      if (!v || !v->is_array()) return;
      Json out = Json::array();
      for (const auto& e : v->as_array()) {
        if (!first_of_pair) { out.push(ref(kind, e)); continue; }
        Json pair = e;
        if (pair.is_array() && pair.size()) pair.as_array()[0] = ref(kind, pair.as_array()[0]);
        out.push(std::move(pair));
      }
      line.set(field, std::move(out));
    };
    for (auto& line : mine) {
      const std::string kind = text(line, "kind");
      if (const auto r = renamed.find(kind + "|" + text(line, "id")); r != renamed.end()) line.set("id", Json(r->second));
      refs(line, "sources", "source", false);
      if (kind == "occasion") refs(line, "present", "being", true);
      if (kind == "tie") {
        for (const char* end : {"a", "b"}) if (const Json* v = line.find(end)) line.set(end, ref("being", *v));
        if (const Json* v = line.find("ends")) line.set("ends", ref("occasion", *v));
        refs(line, "attested", "occasion", true);
      }
    }
    // A tie under one id is one tie only if it joins the same beings in both records.
    for (auto& line : mine) {
      if (text(line, "kind") != "tie") continue;
      const std::string lid = text(line, "id");
      const auto it = lines.find("tie|" + lid);
      if (it != lines.end() && (text(it->second, "a") != text(line, "a") || text(it->second, "b") != text(line, "b"))) {
        problems.push_back({false, lid, "names two ties in two records; in " + canon + " it is " + lid + "@" + canon});
        line.set("id", Json(lid + "@" + canon));
      }
    }
    for (auto& line : mine) {
      const std::string kind = text(line, "kind"), lid = text(line, "id");
      const std::string key = kind + "|" + lid;
      auto it = lines.find(key);
      if (it == lines.end()) { order.push_back(key); lines.emplace(key, line); continue; }
      Json& kept = it->second;
      auto union_of = [&](const char* field, auto same) {
        Json merged = kept.find(field) && kept.at(field).is_array() ? kept.at(field) : Json::array();
        if (const Json* more = line.find(field); more && more->is_array())
          for (const auto& x : more->as_array()) {
            bool seen = false;
            for (const auto& y : merged.as_array()) seen = seen || same(x, y);
            if (!seen) merged.push(x);
          }
        kept.set(field, merged);
      };
      auto first_equal = [](const Json& x, const Json& y) {
        return x.is_array() && y.is_array() && x.size() && y.size() && x.as_array()[0].is_string() && y.as_array()[0].is_string() && x.as_array()[0].as_string() == y.as_array()[0].as_string();
      };
      auto string_equal = [](const Json& x, const Json& y) { return x.is_string() && y.is_string() && x.as_string() == y.as_string(); };
      if (kind == "occasion") {
        union_of("present", first_equal);
        union_of("sources", string_equal);
      } else if (kind == "tie") {
        union_of("attested", first_equal);
        union_of("sources", string_equal);
      }
    }
  }
  std::string about = "The union of " + std::to_string(parts.size()) + " canons, glued along the beings, occasions and ties they share:";
  for (std::size_t i = 0; i < parts.size(); i++) about += (i ? "; " : " ") + parts[i];
  about += ".";
  std::string joined = Json(Json::Object{{"kind", Json("canon")}, {"id", Json(id)}, {"title", Json(title)}, {"about", Json(about)}, {"codebook", Json(CODEBOOK)}, {"status", Json("seed")}}).dump() + "\n";
  for (const auto& key : order) joined += lines.at(key).dump() + "\n";
  Record r = read_record_text(joined);
  r.problems.insert(r.problems.begin(), problems.begin(), problems.end());
  return r;
}

// The derivation of codebook §8: when a tie ends at the death or departure of one of its
// beings and the other is not present, the other is present there as bereaved. Idempotent.
inline void add_bereaved(Record& r) {
  for (const auto& t : r.ties) {
    if (!t.ends) continue;
    Occasion& o = r.occasions[*t.ends];
    if (is_present(o, t.a) == is_present(o, t.b)) continue;
    const std::size_t absent = is_present(o, t.a) ? t.b : t.a;
    if (ends_there(o, absent == t.a ? t.b : t.a)) o.present.push_back({absent, "bereaved", true});
  }
}

}  // namespace history
