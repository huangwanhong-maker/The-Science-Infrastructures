"""Checks of the history check (data/products/history_check.py), offline: the calendars against
Meeus's examples and the reform of 1582, Old Style years, Wikidata's time values, the years a
being lived, the comparisons of years and of days of death, and the reading of a canon. Nothing
is fetched.
"""

from __future__ import annotations

import json

from checks import assert_

import history_check as h  # noqa: E402

# ---- the calendars ------------------------------------------------------------------------------
assert_("history check: JDN of 2000-01-01 (Gregorian) is 2451545", h.jdn(2000, 1, 1) == 2451545, str(h.jdn(2000, 1, 1)))
assert_("history check: JDN of 1957-10-04 (Gregorian) is 2436116 (Meeus, ex. 7.a)", h.jdn(1957, 10, 4) == 2436116)
assert_("history check: JDN of 333-01-27 (Julian) is 1842713 (Meeus, ex. 7.b)", h.jdn(333, 1, 27, "julian") == 1842713)
assert_("history check: 1582-10-04 Julian and 1582-10-15 Gregorian are consecutive days",
        h.jdn(1582, 10, 15) - h.jdn(1582, 10, 4, "julian") == 1)
assert_("history check: Newton's birth, 25 Dec 1642 Julian, is 4 Jan 1643 Gregorian", h.jdn(1642, 12, 25, "julian") == h.jdn(1643, 1, 4))
assert_("history check: the October Revolution, 25 Oct 1917 Julian, is 7 Nov 1917 Gregorian", h.jdn(1917, 10, 25, "julian") == h.jdn(1917, 11, 7))
assert_("history check: Old Style puts 20 March 1726 in 1727 (Newton's death, 31 March 1727 Gregorian)",
        h.jdn(1726, 3, 20, "julian-march") == h.jdn(1727, 3, 31))
assert_("history check: Old Style leaves 25 March and after in their own year", h.jdn(1726, 3, 25, "julian-march") == h.jdn(1726, 3, 25, "julian"))

# ---- Wikidata's time values --------------------------------------------------------------------
t = h.parse_wikidata_time({"time": "+1827-03-26T00:00:00Z", "precision": 11, "calendarmodel": "http://www.wikidata.org/entity/Q1985727"})
assert_("history check: a Gregorian day of Wikidata is read", t == {"year": 1827, "month": 3, "day": 26, "precision": 11, "calendar": "gregorian"}, str(t))
t = h.parse_wikidata_time({"time": "+1642-12-25T00:00:00Z", "precision": 11, "calendarmodel": h.JULIAN_MODEL})
assert_("history check: a Julian day of Wikidata is read as Julian", t["calendar"] == "julian" and t["year"] == 1642)
t = h.parse_wikidata_time({"time": "-0470-00-00T00:00:00Z", "precision": 9})
assert_("history check: Wikidata's 470 BC (-0470, historical numbering) is the record's astronomical -469",
        t["year"] == -469 and t["month"] == 0 and t["precision"] == 9, str(t))
entity = {"claims": {"P570": [
    {"rank": "normal", "mainsnak": {"datavalue": {"type": "time", "value": {"time": "+1827-03-26T00:00:00Z", "precision": 11}}}},
    {"rank": "deprecated", "mainsnak": {"datavalue": {"type": "time", "value": {"time": "+1828-01-01T00:00:00Z", "precision": 11}}}},
]}}
assert_("history check: deprecated statements are left out", [x["year"] for x in h.claim_times(entity, "P570")] == [1827])

# ---- the years a being lived --------------------------------------------------------------------
cases = {"1770/1827": (1770, 1827), "1770?/1827~": (1770, 1827), "/1827": (None, 1827), "1770/": (1770, None), "": (None, None), None: (None, None)}
wrong = [f"{k!r}: {h.lived_years(k)}" for k, v in cases.items() if h.lived_years(k) != v]
assert_("history check: lived years read from EDTF intervals", not wrong, "; ".join(wrong))

year = lambda y: {"year": y, "month": 0, "day": 0, "precision": 9, "calendar": "gregorian"}  # noqa: E731
assert_("history check: years that agree", h.compare_years((1770, 1827), [year(1770)], [year(1827)])["status"] == "agrees")
r = h.compare_years((1770, 1828), [year(1770)], [year(1827)])
assert_("history check: a year that disagrees is named", r["status"] == "disagrees" and "1828" in r["detail"], json.dumps(r))
assert_("history check: no years in the record: identified only", h.compare_years((None, None), [year(1770)], [year(1827)])["status"] == "identified")
assert_("history check: either of two values of Wikidata may agree", h.compare_years((1642, 1727), [year(1642), year(1643)], [year(1727)])["status"] == "agrees")

# ---- days of death ----------------------------------------------------------------------------
day = lambda y, m, d, c="gregorian": {"year": y, "month": m, "day": d, "precision": 11, "calendar": c}  # noqa: E731
assert_("history check: a day of death that agrees", h.compare_death("1827-03-26", "gregorian", [day(1827, 3, 26)])["status"] == "agrees")
r = h.compare_death("1827-03-27", "gregorian", [day(1827, 3, 26)])
assert_("history check: a day of death a day off disagrees, and says by how much", r["status"] == "disagrees" and r["detail"] == "1 days apart", json.dumps(r))
assert_("history check: an Old Style day against Wikidata's Julian day of Newton's death",
        h.compare_death("1726-03-20", "julian-march", [day(1727, 3, 20, "julian")])["status"] == "agrees")
assert_("history check: an Old Style day against a Gregorian one", h.compare_death("1726-03-20", "julian-march", [day(1727, 3, 31)])["status"] == "agrees")
assert_("history check: a month is not a day", h.compare_death("1827-03", "gregorian", [day(1827, 3, 26)])["status"] == "not-a-day")

# ---- reading a canon --------------------------------------------------------------------------
text = "\n".join(json.dumps(x) for x in [
    {"kind": "canon", "id": "test", "title": "A test"},
    {"kind": "being", "id": "haydn", "name": "Joseph Haydn", "lived": "1732/1809", "wikidata": "Q7349"},
    {"kind": "being", "id": "mozart", "name": "Wolfgang Amadeus Mozart"},
    {"kind": "occasion", "id": "o1809_haydn_dies", "type": "death", "date": "1809-05-31", "present": [["haydn", "died"]]},
    {"kind": "occasion", "id": "o1785_evening", "type": "meeting", "date": "1785-02-12", "present": [["mozart", "host"], ["haydn", "guest"]]},
]) + "\n\n"
c = h.read_canon(text)
assert_("history check: a canon's beings and deaths are read", c["id"] == "test" and len(c["beings"]) == 2 and c["deaths"] == [
    {"occasion": "o1809_haydn_dies", "being": "haydn", "date": "1809-05-31", "calendar": "gregorian"}], json.dumps(c))
try:
    h.read_canon('{"kind": "being", "id": "x"}')
    refused = False
except ValueError:
    refused = True
assert_("history check: a file without its canon line is refused", refused)
s = h.summary({"a": {"status": "agrees"}, "b": {"status": "proposed"}, "c": {"status": "agrees"}}, {"o": {"status": "agrees"}})
assert_("history check: the summary counts by status", s == {"beings": 3, "deaths": 1, "beings_agrees": 2, "beings_proposed": 1, "deaths_agrees": 1}, json.dumps(s))
