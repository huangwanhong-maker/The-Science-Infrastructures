"""Checks a canon of the history domain against Wikidata (CC0): the identities of its beings,
the years they lived, and the days of the deaths it records.

A canon (``domains/history/canons/*.jsonl``, coded by the domain's codebook) may be drafted by a
language model; its assertions are proposals until checked. This product checks what a
structured source can check, and reports; it never changes the record, which a person corrects
(codebook §9):

- a being with a ``wikidata`` QID: that entity's label and description, and its years of birth
  and death (P569, P570) against the being's ``lived``: ``agrees`` or ``disagrees``, or
  ``identified`` when the being states no years;
- a being without one: the first match of a search by its name, ``proposed`` for a person to
  confirm, or ``unmatched``;
- an occasion of type ``death``: its day against P570 of the being who died there, both read as
  Julian day numbers in their own calendars (Wikidata's proleptic Julian or Gregorian; the
  record's ``gregorian``, ``julian`` or ``julian-march``).

Wikidata is asked politely: a User-Agent naming the project, about 1.5 s between requests, and
backing off when it asks to (HTTP 429). Its answers are kept in ``data/cache/raw/wikidata/``,
and asked again only with ``--refresh``.

Usage::

    python data/products/history_check.py domains/history/canons/europe-1785-1828.jsonl [--refresh]

writes ``data/cache/products/history/<canon id>.check.json``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "cache" / "raw" / "wikidata"
OUT = ROOT / "data" / "cache" / "products" / "history"
USER_AGENT = "TheScienceInfrastructures-history-check/0.1 (https://github.com/huangwanhong-maker/The-Science-Infrastructures)"
PAUSE = 1.5
JULIAN_MODEL = "http://www.wikidata.org/entity/Q1985786"   # Wikidata's proleptic Julian calendar


# ---- the calendars (Meeus, Astronomical Algorithms, 2nd ed., ch. 7) --------------------------------

def jdn(year: int, month: int, day: int, calendar: str = "gregorian") -> int:
    """The Julian day number (the integer day, noon-based) of a calendar date.

    ``julian-march`` is Old Style as England wrote it until 1752: Julian, the legal year
    beginning on 25 March, so that a date from 1 January to 24 March belongs to the next
    historical year."""
    if calendar == "julian-march" and (month < 3 or (month == 3 and day < 25)):
        year += 1
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    if calendar in ("julian", "julian-march"):
        return day + (153 * m + 2) // 5 + 365 * y + y // 4 - 32083
    return day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045


def parse_wikidata_time(value: dict) -> dict | None:
    """A Wikidata time value: {year, month, day, precision, calendar}, month and day 0 when not known.

    Wikidata writes years BC in historical numbering, with no year 0 (Socrates died in -0399, 399
    BC); the record and the domain write astronomical years (-0398 is 399 BC; codebook §4). A
    negative year is therefore moved by one."""
    m = re.match(r"^([+-])(\d+)-(\d\d)-(\d\d)T", value.get("time", ""))
    if not m:
        return None
    year = int(m.group(2)) * (-1 if m.group(1) == "-" else 1)
    if year < 0:
        year += 1
    return {"year": year, "month": int(m.group(3)), "day": int(m.group(4)), "precision": int(value.get("precision", 9)),
            "calendar": "julian" if value.get("calendarmodel") == JULIAN_MODEL else "gregorian"}


def claim_times(entity: dict, prop: str) -> list[dict]:
    """The time values of a property of an entity, deprecated statements left out."""
    out = []
    for c in entity.get("claims", {}).get(prop, []):
        if c.get("rank") == "deprecated":
            continue
        dv = c.get("mainsnak", {}).get("datavalue", {})
        if dv.get("type") == "time":
            t = parse_wikidata_time(dv.get("value", {}))
            if t:
                out.append(t)
    return out


def lived_years(lived: str | None) -> tuple[int | None, int | None]:
    """The years of birth and death in a being's ``lived`` (an EDTF interval of years: 1770/1827,
    1770?/1827, /1827, 1770/); None where not stated."""
    if not lived:
        return None, None
    parts = lived.split("/")
    if len(parts) != 2:
        parts = [parts[0], ""]

    def year(p: str) -> int | None:
        m = re.match(r"^\s*(-?\d{1,4})", p.strip().strip("~?%"))
        return int(m.group(1)) if m else None
    return year(parts[0]), year(parts[1])


def compare_years(record: tuple[int | None, int | None], born: list[dict], died: list[dict]) -> dict:
    """Whether a being's years agree with Wikidata's: a year agrees when any of Wikidata's values
    (at a precision of a year or finer) has it."""
    years = {"wikidata": [sorted({t["year"] for t in born if t["precision"] >= 9}), sorted({t["year"] for t in died if t["precision"] >= 9})],
             "record": list(record)}
    stated = [(r, w) for r, w in zip(record, years["wikidata"]) if r is not None]
    if not stated:
        return {"status": "identified", "years": years}
    disagree = [r for r, w in stated if w and r not in w]
    unknown = [r for r, w in stated if not w]
    if disagree:
        return {"status": "disagrees", "years": years, "detail": "the record's " + " and ".join(str(r) for r in disagree) + " against Wikidata's"}
    if len(unknown) == len(stated):
        return {"status": "identified", "years": years, "detail": "Wikidata states no years to compare"}
    return {"status": "agrees", "years": years}


def edtf_day(date: str) -> tuple[int, int, int] | None:
    """The day of an EDTF date that names one day (1827-03-26, 1827-03-26~), else None."""
    m = re.match(r"^(-?\d{1,4})-(\d\d)-(\d\d)[~?%]?$", date.strip())
    return (int(m.group(1)), int(m.group(2)), int(m.group(3))) if m else None


def compare_death(date: str, calendar: str, died: list[dict]) -> dict:
    """A recorded day of death against Wikidata's P570, as Julian day numbers in their own calendars."""
    day = edtf_day(date)
    wd = [t for t in died if t["precision"] >= 11]
    if day is None:
        return {"status": "not-a-day", "detail": f"{date} names no single day"}
    if not wd:
        return {"status": "unknown", "detail": "Wikidata states no day of death"}
    mine = jdn(*day, calendar or "gregorian")
    theirs = [jdn(t["year"], t["month"], t["day"], t["calendar"]) for t in wd]
    texts = [f"{t['year']:04d}-{t['month']:02d}-{t['day']:02d} ({t['calendar']})" for t in wd]
    if mine in theirs:
        return {"status": "agrees", "wikidata": texts}
    return {"status": "disagrees", "wikidata": texts, "detail": f"{min(abs(mine - t) for t in theirs)} days apart"}


def read_canon(text: str) -> dict:
    """The parts of a canon this check needs: its id, its beings, and its deaths."""
    head, beings, deaths = None, [], []
    for n, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        try:
            a = json.loads(line)
        except ValueError as error:
            raise ValueError(f"line {n}: not JSON ({error})") from None
        kind = a.get("kind")
        if kind == "canon" and head is None:
            head = a
        elif kind == "being":
            beings.append({k: a.get(k) for k in ("id", "name", "lived", "wikidata")})
        elif kind == "occasion" and a.get("type") == "death":
            for being, role in a.get("present", []):
                if role == "died":
                    deaths.append({"occasion": a.get("id"), "being": being, "date": a.get("date", ""), "calendar": a.get("calendar", "gregorian")})
    if head is None:
        raise ValueError("the canon has no first line of kind 'canon'")
    return {"id": head.get("id"), "title": head.get("title"), "beings": beings, "deaths": deaths}


def summary(beings: dict, deaths: dict) -> dict:
    out: dict = {"beings": len(beings), "deaths": len(deaths)}
    for kind, table in (("beings", beings), ("deaths", deaths)):
        for r in table.values():
            key = f"{kind}_{r['status']}"
            out[key] = out.get(key, 0) + 1
    return out


# ---- asking Wikidata ------------------------------------------------------------------------------

class Wikidata:
    """Wikidata's entities and search, asked politely and kept in the cache."""

    def __init__(self, refresh: bool = False):
        self.refresh = refresh
        self.last = 0.0
        RAW.mkdir(parents=True, exist_ok=True)

    def get(self, url: str, cache: Path) -> dict:
        if cache.exists() and not self.refresh:
            return json.loads(cache.read_text(encoding="utf-8"))
        for wait in (5, 10, 20, 40, 80):
            time.sleep(max(0.0, PAUSE - (time.time() - self.last)))
            self.last = time.time()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
                body = json.loads(urllib.request.urlopen(req, timeout=40).read())
                cache.write_text(json.dumps(body, ensure_ascii=False), encoding="utf-8")
                return body
            except urllib.error.HTTPError as error:
                if error.code not in (429, 503):
                    raise
                time.sleep(int(error.headers.get("Retry-After", wait)) if (error.headers.get("Retry-After") or "").isdigit() else wait)
        raise RuntimeError("Wikidata keeps asking us to wait; try again later")

    def entity(self, qid: str) -> dict | None:
        if not re.fullmatch(r"Q\d+", qid or ""):
            return None
        body = self.get(f"https://www.wikidata.org/wiki/Special:EntityData/{qid}.json", RAW / f"{qid}.json")
        entities = body.get("entities", {})
        # A merged item redirects: its entity comes under the new id.
        return entities.get(qid) or next(iter(entities.values()), None)

    def search(self, name: str) -> dict | None:
        slug = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")[:80]
        q = urllib.parse.urlencode({"action": "wbsearchentities", "search": name, "language": "en", "type": "item", "limit": 1, "format": "json"})
        body = self.get(f"https://www.wikidata.org/w/api.php?{q}", RAW / f"search-{slug}.json")
        hits = body.get("search", [])
        return hits[0] if hits else None


def text_of(entity: dict, field: str) -> str | None:
    values = entity.get(field, {})
    return (values.get("en") or next(iter(values.values()), {})).get("value") if values else None


def check(canon: dict, wd: Wikidata, log=print) -> dict:
    beings, born_died = {}, {}
    for b in canon["beings"]:
        qid = b.get("wikidata")
        if qid:
            e = wd.entity(qid)
            if e is None:
                beings[b["id"]] = {"status": "disagrees", "qid": qid, "detail": f"{qid} is not an entity of Wikidata"}
                continue
            born, died = claim_times(e, "P569"), claim_times(e, "P570")
            born_died[b["id"]] = died
            beings[b["id"]] = {"qid": e.get("id", qid), "label": text_of(e, "labels"), "description": text_of(e, "descriptions"),
                               **compare_years(lived_years(b.get("lived")), born, died)}
        else:
            hit = wd.search(b.get("name") or b["id"])
            if hit:
                beings[b["id"]] = {"status": "proposed", "qid": hit.get("id"), "label": hit.get("label"), "description": hit.get("description"),
                                   "detail": "the first match of a search by name: a person confirms it, in the record"}
            else:
                beings[b["id"]] = {"status": "unmatched", "detail": f"no item of Wikidata is called {b.get('name')!r}"}
        log(f"  {b['id']}: {beings[b['id']]['status']}")
    deaths = {}
    for d in canon["deaths"]:
        died = born_died.get(d["being"])
        if died is None:
            deaths[d["occasion"]] = {"status": "unchecked", "being": d["being"], "detail": "the being has no confirmed QID"}
            continue
        deaths[d["occasion"]] = {"being": d["being"], "date": d["date"], "calendar": d["calendar"], **compare_death(d["date"], d["calendar"], died)}
    return {"canon": canon["id"], "title": canon.get("title"), "checked_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "source": {"name": "Wikidata", "licence": "CC0", "url": "https://www.wikidata.org/"},
            "beings": beings, "deaths": deaths, "summary": summary(beings, deaths)}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Check a canon of the history domain against Wikidata.")
    parser.add_argument("canon", type=Path, help="a canon: domains/history/canons/<id>.jsonl")
    parser.add_argument("--refresh", action="store_true", help="ask Wikidata again rather than reading the cache")
    args = parser.parse_args(argv)
    sys.stdout.reconfigure(encoding="utf-8")
    canon = read_canon(args.canon.read_text(encoding="utf-8"))
    print(f"Checking {canon['id']}: {len(canon['beings'])} beings, {len(canon['deaths'])} deaths, against Wikidata.")
    result = check(canon, Wikidata(args.refresh))
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{canon['id']}.check.json"
    out.write_text(json.dumps(result, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    print(f"{json.dumps(result['summary'])}\nWritten: {out.relative_to(ROOT).as_posix()}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
