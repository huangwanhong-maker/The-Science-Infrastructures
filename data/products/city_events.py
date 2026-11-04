"""Build the record of a district's dated facts for the Urban Science Workbench's Event Spacetime (ADR 0067, §2.6):
Melbourne's middle of the Hoddle Grid, 2002-2024, from the City of Melbourne's and Victoria's published records.

The record holds facts, not a complex: the volumes (CLUE's base properties, each joined to a building of the city
model) with their observed states at every census, the establishments that occupy them, and the dated events the
sources record, each with its type (the design note's §1.4), its time label at its own resolution, its match by
identifier and its source. The urban-science domain's C++ builds the couplings, primitives and footprints from it; no
event is inferred here (the note's §9.4-9.5: every event below is "recorded").

1. **The district.** Landscape Dynamics' district, Elizabeth Street to Russell Street and Flinders Street to Lonsdale
   Street (DISTRICT, as applications/urban-science-workbench/areas.json gives it), its corners carried into the city
   model's local frame (east-north metres about the model's origin, by street_records.enu, exactly as the workbench's
   frame does: plugins/landscape-dynamics/shared.js, boxOf), and widened by MARGIN_M on every side.
2. **Properties.** One row a CLUE base property (the land parcel, the identifier that persists across censuses)
   whose point lies in the box. Its point is the most frequent point of the rows whose property_id is the base's own,
   else of all its rows, over both CLUE datasets; its block likewise; its property_ids every one seen under it in any
   year. Its building is the city model's outline (not a building:part) whose footprint holds the point, holes
   respected, the smallest if several do; else a part that holds it; else the nearest outline within JOIN_M; else none.
3. **Censuses.** One row a CLUE building row (a property_id in a census year) of a property in the box: the observed
   state of a volume at a census. has_showers is empty throughout the Hoddle Grid and is not kept.
4. **Establishments.** One row an establishment in a census year at a property in the box: its trading name, as the
   City publishes it, and its ANZSIC 2006 class.
5. **ANZSIC.** Every class that occurs, with its division, subdivision, group and class codes; the division from the
   subdivision by the ranges of ANZSIC 2006 (DIVISIONS), and the description CLUE gives most often. CLUE's 0000 (vacant
   space) has no division.
6. **Events.** Each with a time label [from, to): a year is 1 January to the next 1 January, a day that day to the
   next. Kept when the label begins in 1952 (fifty years before the window) or later and before the window's end;
   detail.before_window marks those before 2002. Recorded:
   - births (kind "constructed"): one a base property, at its most frequent construction_year (the other years
     seen under it in detail.years_seen, with their counts of rows);
   - alterations (kind "refurbished"): one a base property and distinct refurbished_year;
   - the Development Activity Monitor's completions (kind "development completed"), matched to the base property of
     the first of their property_ids that CLUE knows, else to their point; a birth when CLUE records a building at
     that property constructed within a year of the completion (detail.same_as names that birth) or records no
     earlier building there, an alteration when CLUE records a refurbishment within a year of it (same_as) or only
     earlier buildings (detail.type_by says which held);
   - registrations on the Victorian Heritage Register (institution, kind "heritage registered") and their last
     amendment (alteration, kind "heritage registration amended"), each matching the rule "vhr:<number>" and the place:
     the property in the box joined to the building that holds the place's point (the nearest such property to the
     point), else that building, else the point;
   - event permits whose location names a street of the district (regime, kind the permit's category, lightly
     normalised, the raw text in detail), matched by the names of the city model's walkable ways in the box (the
     kinds Landscape Dynamics' network takes), longest first, "St" read as Street and so on.
7. **Report.** Counts by table, type and source; the joins' outcomes; rows dropped and why; street names in the
   permits that are not the district's; and what the data held that a reader should know.

The Victorian Heritage Database's search gives no dates, and fetch.py's follow rules cannot reach its list (nested
under _embedded): each listed place's record is fetched here, from the listed ids, into the dataset's cache
(places/<id>.json, with places/retrieval.json), once (data/registry/victorian-heritage-register.json, `details`).

Output: data/cache/products/city-events-melbourne.json (compact JSON; tables as {"columns": [...], "rows": [[...]]}),
listed in the index under the kind `city-events`.

Usage::

    python data/acquisition/fetch.py melbourne-clue-buildings melbourne-clue-establishments \\
        melbourne-development-activity victorian-heritage-register melbourne-event-permits osm-city
    python data/products/street_records.py --only city          # the city model, if not yet built
    python data/products/city_events.py                         # fetches the heritage records it lacks
    python data/products/city_events.py --offline               # ... or refuses to
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
CACHE = DATA / "cache"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402
from street_records import contains, enu, retrieved, signed_area, write  # noqa: E402

UTC = dt.timezone.utc
AREA = "melbourne"
CITY_MODEL = OUT / f"city-model-{AREA}.json"
BUILDINGS = CACHE / "melbourne-clue-buildings" / "hoddle-grid"
ESTABLISHMENTS = CACHE / "melbourne-clue-establishments" / "all"
DEVELOPMENTS = CACHE / "melbourne-development-activity" / "hoddle-grid"
HERITAGE = CACHE / "victorian-heritage-register" / "hoddle-grid"
PERMITS = CACHE / "melbourne-event-permits" / "all"
REGISTRY = DATA / "registry"

# Landscape Dynamics' district (applications/urban-science-workbench/areas.json, melbourne.district), WGS 84 degrees.
DISTRICT = {"name": "The Hoddle Grid's middle: Elizabeth Street to Russell Street, Flinders Street to Lonsdale Street",
            "south": -37.8173, "west": 144.9617, "north": -37.8117, "east": 144.9697}
MARGIN_M = 60.0                 # about the district's box, so that the streets round it and their frontages are in
JOIN_M = 15.0                   # the farthest a property's point may lie from the footprint it is joined to
WINDOW = (2002, 2025)           # [2002-01-01, 2025-01-01): the censuses 2002-2024
EARLIEST = WINDOW[0] - 50       # events from 1952: what stood or was done within fifty years before the window
SAME_YEAR = 1                   # a completion and a construction or refurbishment year this close are one fact

# ANZSIC 2006 (Revision 2.0): the divisions and the subdivisions (two digits) each spans. Australian Bureau of
# Statistics, Australian and New Zealand Standard Industrial Classification (ANZSIC), 2006 (Revision 2.0), cat. no.
# 1292.0, "Division and subdivision codes and titles",
# https://www.abs.gov.au/statistics/classifications/australian-and-new-zealand-standard-industrial-classification-anzsic/2006-revision-2-0/numbering-system-and-titles/division-and-subdivision-codes-and-titles
# (retrieved 2026-11-04). The numbering leaves gaps (61, 65, 68, 71, 74, 78, 79, 83, 88, 93 are no subdivision).
DIVISIONS = [
    ("A", 1, 5, "Agriculture, Forestry and Fishing"),
    ("B", 6, 10, "Mining"),
    ("C", 11, 25, "Manufacturing"),
    ("D", 26, 29, "Electricity, Gas, Water and Waste Services"),
    ("E", 30, 32, "Construction"),
    ("F", 33, 38, "Wholesale Trade"),
    ("G", 39, 43, "Retail Trade"),
    ("H", 44, 45, "Accommodation and Food Services"),
    ("I", 46, 53, "Transport, Postal and Warehousing"),
    ("J", 54, 60, "Information Media and Telecommunications"),
    ("K", 62, 64, "Financial and Insurance Services"),
    ("L", 66, 67, "Rental, Hiring and Real Estate Services"),
    ("M", 69, 70, "Professional, Scientific and Technical Services"),
    ("N", 72, 73, "Administrative and Support Services"),
    ("O", 75, 77, "Public Administration and Safety"),
    ("P", 80, 82, "Education and Training"),
    ("Q", 84, 87, "Health Care and Social Assistance"),
    ("R", 89, 92, "Arts and Recreation Services"),
    ("S", 94, 96, "Other Services"),
]

# The city model's ways that Landscape Dynamics' walking network takes (plugins/landscape-dynamics/shared.js, WALK).
WALK = {"primary", "secondary", "tertiary", "residential", "unclassified", "living_street", "pedestrian", "service",
        "primary_link", "secondary_link", "tertiary_link"}
# How the booking system abbreviates the words of street names (a word after the first).
ABBREVIATIONS = {"Street": r"(?:Street|St\.?)", "Lane": r"(?:Lane|Ln\.?)", "Place": r"(?:Place|Pl\.?)", "Road": r"(?:Road|Rd\.?)",
                 "Avenue": r"(?:Avenue|Ave\.?)", "Arcade": r"(?:Arcade|Arc\.?)", "Parade": r"(?:Parade|Pde\.?)", "Court": r"(?:Court|Ct\.?)"}
# Words that end a street-like name, for counting the names the permits give that are not the district's.
STREET_WORDS = {"Street": "Street", "St": "Street", "Lane": "Lane", "Ln": "Lane", "Place": "Place", "Pl": "Place", "Road": "Road",
                "Rd": "Road", "Avenue": "Avenue", "Ave": "Avenue", "Alley": "Alley", "Arcade": "Arcade", "Parade": "Parade",
                "Pde": "Parade", "Way": "Way", "Walk": "Walk", "Boulevard": "Boulevard", "Court": "Court", "Mall": "Mall",
                "Promenade": "Promenade", "Esplanade": "Esplanade", "Drive": "Drive", "Highway": "Highway", "Crescent": "Crescent",
                "Terrace": "Terrace", "Grove": "Grove", "Close": "Close"}
NOT_THE_STREET = r"(?!\s+(?:Station|Railway))"      # 'Flinders St Station' names the station, not the street
STREET_LIKE = re.compile(r"\b((?:[A-Z][A-Za-z'’]*\.?\s+){1,3})(" + "|".join(STREET_WORDS) + r")\b\.?" + NOT_THE_STREET)


# ---- time labels ---------------------------------------------------------------------------------------------

def year_label(year: int) -> tuple[str, str]:
    """A year as the interval it resolves: [1 January, the next 1 January)."""
    return f"{year:04d}-01-01", f"{year + 1:04d}-01-01"


def day_label(first: dt.date, last: dt.date | None = None) -> tuple[str, str]:
    """Days as the interval they resolve: [the first day, the day after the last)."""
    return first.isoformat(), ((last or first) + dt.timedelta(days=1)).isoformat()


def in_window(label: tuple[str, str]) -> str:
    """'before' (begins before EARLIEST), 'after' (begins at or after the window's end), else 'in'."""
    if label[0] < f"{EARLIEST:04d}-01-01":
        return "before"
    if label[0] >= f"{WINDOW[1]:04d}-01-01":
        return "after"
    return "in"


# ---- ANZSIC ----------------------------------------------------------------------------------------------------

def anzsic_division(subdivision: int) -> str | None:
    """The ANZSIC 2006 division (a letter) of a subdivision (two digits), None for a number no division spans."""
    return next((d for d, lo, hi, _ in DIVISIONS if lo <= subdivision <= hi), None)


def anzsic_levels(code: str) -> dict:
    """The levels of a four-digit class code: its division, subdivision, group and class."""
    code = code.strip()
    if not re.fullmatch(r"\d{4}", code):
        return {"division": None, "subdivision": None, "group": None, "class": code or None}
    return {"division": anzsic_division(int(code[:2])), "subdivision": code[:2], "group": code[:3], "class": code}


# ---- the frame and the joins -----------------------------------------------------------------------------------

class Frame:
    """The city model's local frame: x east, y north (m) about its origin, by street_records.enu."""

    def __init__(self, origin: dict):
        self.lat0, self.lon0 = origin["lat"], origin["lon"]

    def xy(self, lat, lon):
        e, n, _ = enu(np.asarray(lat, float), np.asarray(lon, float), self.lat0, self.lon0)
        return e, n

    def latlon(self, x: float, y: float) -> tuple[float, float]:
        """The inverse, by Newton's method from the tangent plane (three steps reach well under a millimetre)."""
        k = math.pi / 180 * 6378137.0
        lat, lon = self.lat0 + y / k, self.lon0 + x / (k * math.cos(math.radians(self.lat0)))
        for _ in range(4):
            e, n = self.xy(lat, lon)
            lat -= (float(n) - y) / k
            lon -= (float(e) - x) / (k * math.cos(math.radians(lat)))
        return lat, lon


def district_box(frame: Frame, district: dict = DISTRICT, margin: float = MARGIN_M) -> dict:
    """The district's box in the local frame, from its corners (as the workbench's boxOf takes them), widened by margin."""
    x0, y0 = frame.xy(district["south"], district["west"])
    x1, y1 = frame.xy(district["north"], district["east"])
    inner = [float(min(x0, x1)), float(min(y0, y1)), float(max(x0, x1)), float(max(y0, y1))]
    outer = [inner[0] - margin, inner[1] - margin, inner[2] + margin, inner[3] + margin]
    corners = [frame.latlon(x, y) for x in (outer[0], outer[2]) for y in (outer[1], outer[3])]
    lats, lons = [c[0] for c in corners], [c[1] for c in corners]
    return {"local": [round(v, 2) for v in outer], "lonlat": [round(min(lons), 7), round(min(lats), 7), round(max(lons), 7), round(max(lats), 7)],
            "district": {"local": [round(v, 2) for v in inner], "lonlat": [district["west"], district["south"], district["east"], district["north"]],
                         "name": district["name"]},
            "margin_m": margin}


def inside_box(x, y, box: list[float]):
    return (x >= box[0]) & (x <= box[2]) & (y >= box[1]) & (y <= box[3])


def point_in_footprint(p, footprint, holes=()) -> bool:
    """Whether p lies inside a footprint and outside its holes (each ring by the even-odd rule)."""
    return contains(footprint, p) and not any(contains(h, p) for h in holes)


def distance_to_rings(p, rings) -> float:
    """The least distance (m) from p to the edges of the rings (each a closed polygon, not repeated at its end)."""
    best = math.inf
    px, py = p
    for ring in rings:
        r = np.asarray(ring, float)
        a, b = r, np.roll(r, -1, axis=0)
        d = b - a
        L = np.einsum("ij,ij->i", d, d)
        t = np.clip(np.einsum("ij,ij->i", np.array([px, py]) - a, d) / np.where(L > 0, L, 1.0), 0.0, 1.0)
        q = a + t[:, None] * d
        best = min(best, float(np.min(np.hypot(q[:, 0] - px, q[:, 1] - py))))
    return best


class Footprints:
    """The city model's buildings near a box, for joining points to them."""

    def __init__(self, buildings: list[dict], box: list[float], reach: float = JOIN_M):
        self.items = []
        for b in buildings:
            f = b.get("footprint") or []
            if len(f) < 3:
                continue
            a = np.asarray(f, float)
            lo, hi = a.min(axis=0), a.max(axis=0)
            if hi[0] < box[0] - reach or lo[0] > box[2] + reach or hi[1] < box[1] - reach or lo[1] > box[3] + reach:
                continue
            self.items.append({"id": b["id"], "footprint": f, "holes": b.get("holes") or [], "lo": lo, "hi": hi,
                               "part": str(b.get("kind", "")).startswith("part:"), "area": abs(signed_area(f))})

    def join(self, p, reach: float = JOIN_M) -> tuple[str | None, str, float | None]:
        """(building id, 'inside' | 'nearest' | 'none', distance in m for 'nearest')."""
        x, y = p
        holding = [b for b in self.items if b["lo"][0] <= x <= b["hi"][0] and b["lo"][1] <= y <= b["hi"][1]
                   and point_in_footprint(p, b["footprint"], b["holes"])]
        if holding:
            holding.sort(key=lambda b: (b["part"], b["area"]))
            return holding[0]["id"], "inside", None
        best, best_d = None, math.inf
        for b in self.items:
            if b["part"] or x < b["lo"][0] - reach or x > b["hi"][0] + reach or y < b["lo"][1] - reach or y > b["hi"][1] + reach:
                continue
            d = distance_to_rings(p, [b["footprint"], *b["holes"]])
            if d < best_d:
                best, best_d = b["id"], d
        if best is not None and best_d <= reach:
            return best, "nearest", best_d
        return None, "none", None


def join_label(kind: str, d: float | None) -> str:
    return f"nearest:{d:.1f}" if kind == "nearest" else kind


# ---- reading the sources -----------------------------------------------------------------------------------------

def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=";"))


def integer(text) -> int | None:
    text = (text or "").strip()
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def mode(values):
    """The most frequent value (the first seen among the most frequent), None if there are none."""
    counts = Counter(v for v in values if v is not None)
    if not counts:
        return None
    top = max(counts.values())
    return next(v for v in counts if counts[v] == top)


def manifest(ident: str) -> dict:
    return json.loads((REGISTRY / f"{ident}.json").read_text(encoding="utf-8"))


def source_entry(ident: str, folder: Path) -> dict:
    m = manifest(ident)
    return {"id": ident, "title": m["title"], "credit": m["credit"], "licence": m["licence"], "retrieved": retrieved(folder)}


def heritage_records(folder: Path, offline: bool, refresh: bool = False) -> dict[str, dict]:
    """Each listed place's record (/v1/places/{id}), fetched into folder/places/ when not yet held, with a provenance
    record beside them; returns them by id."""
    listing = json.loads((folder / "places.json").read_text(encoding="utf-8"))
    places = listing.get("_embedded", {}).get("places", [])
    if listing.get("page_count", 1) > 1:
        print(f"  warning: the search has {listing['page_count']} pages; only the first is held", file=sys.stderr)
    spec = manifest("victorian-heritage-register")["details"]
    out_dir = folder / "places"
    out_dir.mkdir(exist_ok=True)
    record_path = out_dir / "retrieval.json"
    provenance = json.loads(record_path.read_text(encoding="utf-8")) if record_path.exists() else {"dataset": "victorian-heritage-register", "requests": []}
    held = {r["file"]: r for r in provenance["requests"]}
    fetched = 0
    for p in places:
        ident = str(p["id"])
        file = Path(spec["file"].format(id=ident)).name
        dest = out_dir / file
        if dest.exists() and not refresh:
            continue
        if offline:
            raise FileNotFoundError(f"the heritage record {ident} is not held; run without --offline to fetch it")
        sys.path.insert(0, str(DATA / "acquisition"))
        import fetch                                                    # its download: user agent, retries, .part files
        url = spec["url"].format(id=ident)
        fetch.download(url, dest, timeout=60, attempts=3)
        held[file] = {"url": url, "file": file, "size": dest.stat().st_size, "sha256": fetch.sha256_of(dest), "retrieved_at": fetch.now()}
        fetched += 1
    if fetched:
        provenance.update({"retrieved_at": dt.datetime.now(UTC).isoformat(timespec="seconds"), "requests": sorted(held.values(), key=lambda r: r["file"])})
        record_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
        print(f"  heritage: {fetched} place records fetched into {out_dir.relative_to(DATA.parent)}")
    records = {}
    for p in places:
        path = out_dir / Path(spec["file"].format(id=p["id"])).name
        records[str(p["id"])] = {"listing": p, "record": json.loads(path.read_text(encoding="utf-8"))}
    return records


def heritage_date(text) -> dt.date | None:
    """A date of the Heritage Database ('October 9, 1974'); None if empty or unreadable."""
    text = (text or "").strip()
    for fmt in ("%B %d, %Y", "%d %B %Y", "%Y-%m-%d"):
        try:
            return dt.datetime.strptime(text, fmt).date()
        except ValueError:
            continue
    return None


# ---- streets named in free text ------------------------------------------------------------------------------------

def street_pattern(name: str) -> re.Pattern:
    """A pattern for a street's name as the booking system may write it: its words after the first abbreviated or not,
    spaces of any width, any case; not followed by 'Station'."""
    words = name.split()
    parts = [re.escape(w) if i == 0 or w not in ABBREVIATIONS else ABBREVIATIONS[w] for i, w in enumerate(words)]
    return re.compile(r"(?<![\w'’])" + r"\s+".join(parts) + r"(?![\w'’])" + NOT_THE_STREET, re.IGNORECASE)


class StreetNames:
    """The names of the district's streets, matched in a text longest first, each span taken once."""

    def __init__(self, names):
        self.names = sorted(set(names), key=lambda n: (-len(n), n))
        self.patterns = [(n, street_pattern(n)) for n in self.names]

    def find(self, text: str) -> tuple[list[str], str]:
        """The district's streets a text names, in the order of the text, and the text with their spans blanked."""
        found, rest = [], text
        for name, pattern in self.patterns:
            for m in pattern.finditer(rest):
                found.append((m.start(), name))
                rest = rest[:m.start()] + " " * (m.end() - m.start()) + rest[m.end():]
        names = []
        for _, n in sorted(found):
            if n not in names:
                names.append(n)
        return names, rest


def other_street_names(text: str) -> list[str]:
    """Street-like names in a text (capitalised words ending in a street's word), the word written out."""
    out = []
    for m in STREET_LIKE.finditer(text):
        lead = " ".join(m.group(1).split())
        name = f"{lead} {STREET_WORDS[m.group(2)]}"
        if name not in out:
            out.append(name)
    return out


CATEGORY_FIXES = [(re.compile(r"\bPubilc\b", re.I), "Public"), (re.compile(r"\bNon[- ]?ticketed\b", re.I), "Non-ticketed"),
                  (re.compile(r"\bRun[ /]Walk\b", re.I), "Run/Walk"), (re.compile(r"--+"), "-")]


def permit_kind(category: str) -> str:
    """A permit's category with its variants of spelling made one: whitespace collapsed, a trailing comma or dash
    dropped, 'Pubilc' read as Public, 'Non Ticketed' as Non-ticketed, 'Run Walk' as Run/Walk."""
    kind = " ".join((category or "").split())
    for pattern, repl in CATEGORY_FIXES:
        kind = pattern.sub(repl, kind)
    return kind.rstrip(" ,-") or "unknown"


def iso_day(text) -> dt.date | None:
    try:
        return dt.date.fromisoformat((text or "").strip()[:10])
    except ValueError:
        return None


# ---- the record ------------------------------------------------------------------------------------------------------

def table(columns: list[str], rows: list[list]) -> dict:
    return {"columns": columns, "rows": rows}


def build(offline: bool = False, refresh_details: bool = False) -> dict:
    model = json.loads(CITY_MODEL.read_text(encoding="utf-8"))
    frame = Frame(model["origin"])
    box = district_box(frame)
    B = box["local"]
    report = {"dropped": Counter(), "notes": []}

    # -- CLUE's rows, both datasets, joined by base property --------------------------------------------------------
    buildings = read_csv(BUILDINGS / "buildings.csv")
    establishments = read_csv(ESTABLISHMENTS / "establishments.csv")
    if establishments and "base_property_id" not in establishments[0]:
        raise SystemExit("the establishments lack base_property_id: run python data/acquisition/fetch.py melbourne-clue-establishments --refresh")
    # Establishments of the whole municipality: a coarse cut in degrees first (the box is about 0.01 degrees across).
    lon0, lat0, lon1, lat1 = box["lonlat"]
    pad = 0.01
    near = [r for r in establishments if r["latitude"] and r["longitude"]
            and lat0 - pad <= float(r["latitude"]) <= lat1 + pad and lon0 - pad <= float(r["longitude"]) <= lon1 + pad]
    report["dropped"]["establishment rows without base_property_id"] = sum(1 for r in establishments if not r["base_property_id"])

    rows_by_base: dict[str, list[tuple[str, dict]]] = defaultdict(list)
    for r in buildings:
        if r["base_property_id"]:
            rows_by_base[r["base_property_id"]].append(("building", r))
    for r in near:
        if r["base_property_id"]:
            rows_by_base[r["base_property_id"]].append(("establishment", r))

    def point_of(rows):
        own = [(float(r["latitude"]), float(r["longitude"])) for _, r in rows
               if r["latitude"] and r["longitude"] and r["property_id"] == r["base_property_id"]]
        every = [(float(r["latitude"]), float(r["longitude"])) for _, r in rows if r["latitude"] and r["longitude"]]
        return mode(own) or mode(every)

    bases = sorted(rows_by_base, key=int)
    points = {b: point_of(rows_by_base[b]) for b in bases}
    lat = np.array([points[b][0] if points[b] else np.nan for b in bases])
    lon = np.array([points[b][1] if points[b] else np.nan for b in bases])
    px, py = frame.xy(lat, lon)
    keep = inside_box(px, py, B)
    in_box = {b: (round(float(px[i]), 2), round(float(py[i]), 2)) for i, b in enumerate(bases) if keep[i]}

    # -- properties -------------------------------------------------------------------------------------------------
    footprints = Footprints(model["buildings"], B)
    joins = Counter()
    properties, spread_over = [], 0
    for b in sorted(in_box, key=int):
        rows = rows_by_base[b]
        pids = sorted({int(r["property_id"]) for _, r in rows if r["property_id"]})
        block = mode(integer(r["block_id"]) for _, r in rows)
        p = in_box[b]
        building, kind, d = footprints.join(p)
        joins[kind] += 1
        properties.append([f"base:{b}", int(b), pids, block, list(p), building, join_label(kind, d)])
        pts = {(float(r["latitude"]), float(r["longitude"])) for _, r in rows if r["latitude"] and r["longitude"]}
        if len(pts) > 1:
            xs, ys = frame.xy([q[0] for q in pts], [q[1] for q in pts])
            if float(np.max(np.hypot(xs - p[0], ys - p[1]))) > 50:
                spread_over += 1
    building_of = {row[0]: row[5] for row in properties}
    report["notes"].append(f"{spread_over} base properties have rows whose points lie more than 50 m from the property's point "
                           "(several buildings or occupancies on a large parcel); the property keeps one point")

    # -- censuses ---------------------------------------------------------------------------------------------------
    census_rows, seen = [], Counter()
    for b in sorted(in_box, key=int):
        for src, r in rows_by_base[b]:
            if src != "building":
                continue
            seen[(r["census_year"], r["property_id"])] += 1
            census_rows.append([int(r["census_year"][:4]), f"base:{b}", int(r["property_id"]), integer(r["construction_year"]),
                                integer(r["refurbished_year"]), integer(r["number_of_floors_above_ground"]),
                                r["predominant_space_use"] or None, r["building_name"] or None,
                                integer(r["accessibility_rating"]), integer(r["bicycle_spaces"])])
    census_rows.sort(key=lambda r: (r[0], r[2], r[1]))
    duplicated = sum(1 for v in seen.values() if v > 1)
    if duplicated:
        report["notes"].append(f"{duplicated} (census year, property_id) pairs appear in more than one row of the building "
                               "information (different buildings under one property_id, such as Melbourne Central and its "
                               "railway station in 2002-2006); every row is kept")

    # -- establishments and their classes ---------------------------------------------------------------------------
    est_rows, descriptions = [], defaultdict(Counter)
    for b in sorted(in_box, key=int):
        for src, r in rows_by_base[b]:
            if src != "establishment":
                continue
            code = (r["industry_anzsic4_code"] or "").strip()
            descriptions[code][(r["industry_anzsic4_description"] or "").strip()] += 1
            est_rows.append([int(r["census_year"][:4]), f"base:{b}", integer(r["property_id"]), r["trading_name"] or None, code or None])
    est_rows.sort(key=lambda r: (r[0], r[2] or 0, r[3] or "", r[4] or ""))
    anzsic = {}
    for code in sorted(c for c in descriptions if c):
        levels = anzsic_levels(code)
        levels["description"] = descriptions[code].most_common(1)[0][0] or None
        if len(descriptions[code]) > 1:
            levels["descriptions_seen"] = len(descriptions[code])
        anzsic[code] = levels
    no_division = sorted(c for c, v in anzsic.items() if v["division"] is None)
    if no_division:
        report["notes"].append(f"ANZSIC codes without a division: {', '.join(no_division)} (CLUE's own: 0000 is vacant space)")

    # -- events -------------------------------------------------------------------------------------------------------
    events = []

    def emit(ident, type_, kind, label, match, source, detail, why_dropped):
        where = in_window(label)
        if where != "in":
            report["dropped"][f"{why_dropped}: before {EARLIEST}" if where == "before" else f"{why_dropped}: from {WINDOW[1]}"] += 1
            return False
        if label[0] < f"{WINDOW[0]:04d}-01-01":
            detail = {**(detail or {}), "before_window": True}
        events.append([ident, type_, kind, label[0], label[1], match, source, "recorded", detail or None])
        return True

    # CLUE: one birth a base property at its most frequent construction year; an alteration a distinct refurbished year.
    births, refurbs = {}, defaultdict(dict)
    for b in sorted(in_box, key=int):
        rows = [r for src, r in rows_by_base[b] if src == "building"]
        years = Counter(integer(r["construction_year"]) for r in rows if integer(r["construction_year"]))
        if not years:
            if rows:
                report["dropped"]["CLUE base properties with no construction_year"] += 1
            continue
        year = mode(integer(r["construction_year"]) for r in rows)
        detail = {"rows": years[year]}
        if len(years) > 1:
            detail["years_seen"] = {str(y): n for y, n in sorted(years.items())}
        ident = f"clue-constructed:{b}"
        if emit(ident, "birth", "constructed", year_label(year), {"property": f"base:{b}"}, "melbourne-clue-buildings", detail,
                "CLUE construction years"):
            births[b] = (year, ident, set(years))
        by_year = defaultdict(lambda: {"property_ids": set(), "censuses": []})
        for r in rows:
            y = integer(r["refurbished_year"])
            if y:
                by_year[y]["property_ids"].add(int(r["property_id"]))
                by_year[y]["censuses"].append(int(r["census_year"][:4]))
        for y, v in sorted(by_year.items()):
            detail = {"property_ids": sorted(v["property_ids"]), "censuses": [min(v["censuses"]), max(v["censuses"])]}
            if y < min(years):
                detail["before_construction"] = True
            ident = f"clue-refurbished:{b}:{y}"
            if emit(ident, "alteration", "refurbished", year_label(y), {"property": f"base:{b}"}, "melbourne-clue-buildings", detail,
                    "CLUE refurbished years"):
                refurbs[b][y] = ident
    years_by_base = {b: {integer(r["construction_year"]) for src, r in rows_by_base[b] if src == "building" and integer(r["construction_year"])}
                     for b in in_box}
    multi_in_window = sum(1 for b, ys in years_by_base.items() if len(ys) > 1 and b in births
                          and any(WINDOW[0] <= y < WINDOW[1] and y != births[b][0] for y in ys))
    report["notes"].append(f"{sum(1 for ys in years_by_base.values() if len(ys) > 1)} base properties carry more than one construction "
                           f"year over the censuses (several buildings, a rebuilding, or corrections); {multi_in_window} of them have a "
                           "year in the window other than the one taken (see detail.years_seen)")

    # The Development Activity Monitor's completions.
    pid_to_base = {}
    for b, rows in rows_by_base.items():
        for _, r in sorted(rows, key=lambda t: t[1]["census_year"]):
            if r["property_id"]:
                pid_to_base[r["property_id"]] = b               # the latest census's parent
        pid_to_base.setdefault(b, b)
    dam = read_csv(DEVELOPMENTS / "developments.csv")
    floor_uses = ["office", "retail", "industrial", "storage", "education", "hospital", "recreation", "publicdispaly", "community"]
    counts_kept = ["resi_dwellings", "studio_dwe", "one_bdrm_dwe", "two_bdrm_dwe", "three_bdrm_dwe", "student_apartments", "student_beds",
                   "student_accommodation_units", "institutional_accom_beds", "hotel_rooms", "serviced_apartments",
                   "hotels_serviced_apartments", "hostel_beds", "childcare_places", "car_spaces", "bike_spaces"]
    dam_match = Counter()
    for r in dam:
        if r["status"].strip().upper() != "COMPLETED" or not integer(r["year_completed"]):
            report["dropped"][f"developments not completed ({r['status'].strip().lower() or 'no status'})"] += 1
            continue
        year = integer(r["year_completed"])
        pids = [r[k].strip() for k in ("property_id", "property_id_2", "property_id_3", "property_id_4", "property_id_5") if r.get(k, "").strip()]
        known = [pid_to_base[p] for p in pids if p in pid_to_base]
        base = known[0] if known else None
        x, y = (float(v) for v in frame.xy(float(r["latitude"]), float(r["longitude"]))) if r["latitude"] and r["longitude"] else (None, None)
        if base is not None and base not in in_box:
            report["dropped"]["developments at a property outside the box"] += 1
            continue
        if base is None and (x is None or not inside_box(x, y, B)):
            report["dropped"]["developments with no known property, their point outside the box"] += 1
            continue
        detail = {"development_key": r["development_key"], "data_format": r["data_format"], "street_address": r["street_address"] or None,
                  "property_ids": [int(p) for p in pids if p.isdigit()], "clue_block": integer(r["clue_block"]),
                  "floors_above": integer(r["floors_above"])}
        detail.update({k: integer(r[k]) for k in counts_kept if integer(r.get(k))})
        space = {u: integer(r[f"{u}_flr"]) for u in floor_uses if integer(r.get(f"{u}_flr"))}
        if space:
            detail["floor_space_m2"] = space
        if r.get("town_planning_application", "").strip():
            detail["town_planning_application"] = r["town_planning_application"].strip()
        if len(set(known)) > 1:
            detail["also_properties"] = [f"base:{k}" for k in dict.fromkeys(known[1:]) if k != base]
        type_ = "birth"
        if base is not None:
            match = {"property": f"base:{base}"}
            dam_match["property"] += 1
            built = years_by_base.get(base, set())
            if any(abs(c - year) <= SAME_YEAR for c in built):
                detail["type_by"] = "CLUE records a building here constructed within a year of the completion"
                if base in births and abs(births[base][0] - year) <= SAME_YEAR:
                    detail["same_as"] = births[base][1]
            elif any(abs(ry - year) <= SAME_YEAR for ry in refurbs.get(base, {})):
                type_ = "alteration"
                ry = min(refurbs[base], key=lambda ry: abs(ry - year))
                detail["type_by"] = "CLUE records a refurbishment here within a year of the completion"
                detail["same_as"] = refurbs[base][ry]
            elif built and max(built) < year - SAME_YEAR:
                type_ = "alteration"
                detail["type_by"] = "CLUE records only earlier buildings here, and none constructed at the completion"
            else:
                detail["type_by"] = "the monitor's own: a major new development (CLUE records no earlier building here)"
        else:
            match = {"point": [round(x, 2), round(y, 2)]}
            dam_match["point"] += 1
            bid, kind, d = footprints.join((x, y))
            detail.update({"building": bid, "join": join_label(kind, d)})
            detail["type_by"] = "the monitor's own: a major new development (no CLUE property known)"
        emit(f"dam:{r['development_key']}", type_, "development completed", year_label(year), match, "melbourne-development-activity",
             detail, "development completions")

    # The Victorian Heritage Register.
    heritage = heritage_records(HERITAGE, offline, refresh_details)
    props_by_building = defaultdict(list)
    for row in properties:
        if row[5]:
            props_by_building[row[5]].append(row)
    vhr_match = Counter()
    for pid_, h in heritage.items():
        rec = h["record"]
        lat_, lon_ = rec.get("latitude"), rec.get("longitude")
        if not lat_ or not lon_:
            lat_, lon_ = (h["listing"].get("latlon") or ",").split(",")
        try:
            x, y = (float(v) for v in frame.xy(float(lat_), float(lon_)))
        except ValueError:
            report["dropped"]["heritage places without a position"] += 1
            continue
        if not inside_box(x, y, B):
            report["dropped"]["heritage places outside the box"] += 1
            continue
        vhr = (rec.get("vhr_number") or h["listing"].get("vhr_number") or "").strip()
        rule = f"vhr:{vhr}" if vhr else f"vhd:{pid_}"
        bid, kind, d = footprints.join((x, y))
        match = {"rule": rule}
        on = props_by_building.get(bid, []) if bid else []
        if on:
            nearest = min(on, key=lambda row: math.hypot(row[4][0] - x, row[4][1] - y))
            match.update({"property": nearest[0], "building": bid})
            vhr_match["property"] += 1
        elif bid:
            match.update({"building": bid, "point": [round(x, 2), round(y, 2)]})
            vhr_match["building"] += 1
        else:
            match["point"] = [round(x, 2), round(y, 2)]
            vhr_match["point"] += 1
        categories = sorted(c.get("name", "") for c in (rec.get("heritage_act_categories") or {}).values()) if isinstance(rec.get("heritage_act_categories"), dict) else []
        detail = {"vhr_number": vhr or None, "place_id": int(pid_), "name": rec.get("name"), "location": rec.get("location"),
                  "status": rec.get("status"), "categories": categories, "join": join_label(kind, d),
                  "url": h["listing"].get("url")}
        registered = heritage_date(rec.get("registered_on"))
        amended = heritage_date(rec.get("registration_last_amended_on"))
        if registered is None:
            report["dropped"]["heritage places without a registration date"] += 1
        else:
            emit(f"vhr-registered:{vhr or pid_}", "institution", "heritage registered", day_label(registered), match,
                 "victorian-heritage-register", detail, "heritage registrations")
        if amended is not None and amended != registered:
            emit(f"vhr-amended:{vhr or pid_}", "alteration", "heritage registration amended", day_label(amended), match,
                 "victorian-heritage-register", {"vhr_number": vhr or None, "place_id": int(pid_), "registered_on": registered.isoformat() if registered else None},
                 "heritage amendments")

    # Event permits that name a street of the district.
    names = sorted({r["name"] for r in model["roads"] if r.get("name") and r.get("kind") in WALK
                    and any(B[0] <= q[0] <= B[2] and B[1] <= q[1] <= B[3] for q in r.get("line", []))})
    streets = StreetNames(names)
    unmatched, permit_streets, occurrences = Counter(), Counter(), Counter()
    for r in read_csv(PERMITS / "permits.csv"):
        text = r["location"] or ""
        found, rest = streets.find(text)
        for n in other_street_names(rest):
            unmatched[n] += 1
        if not found:
            report["dropped"]["permits naming no street of the district"] += 1
            continue
        start, end = iso_day(r["event_start"]), iso_day(r["event_end"])
        if start is None:
            report["dropped"]["permits without a start date"] += 1
            continue
        detail = {"location": text, "category": r["category_1"] or None}
        if r.get("category_2"):
            detail["category_2"] = r["category_2"]
        if end is None:
            detail["end_missing"] = True
            end = start
        elif end < start:
            detail["end_before_start"] = end.isoformat()
            end = start
        key = "|".join([r["event_start"], r["event_end"], r["category_1"], r.get("category_2", ""), text])
        digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:10]
        occurrences[digest] += 1
        ident = f"permit:{digest}" + (f"-{occurrences[digest]}" if occurrences[digest] > 1 else "")
        for n in found:
            permit_streets[n] += 1
        emit(ident, "regime", permit_kind(r["category_1"]), day_label(start, end), {"streets": found}, "melbourne-event-permits", detail, "permits")

    events.sort(key=lambda e: (e[3], e[0]))

    # -- the product --------------------------------------------------------------------------------------------------
    by_type = Counter(e[1] for e in events)
    by_source = defaultdict(Counter)
    for e in events:
        by_source[e[6]][e[1]] += 1
    est_years = Counter(r[0] for r in est_rows)
    product = {
        "product": "city-events",
        "id": AREA,
        "built": dt.datetime.now(UTC).isoformat(timespec="seconds"),
        "frame": model["frame"],
        "origin": model["origin"],
        "city_model": {"file": CITY_MODEL.name, "built": model.get("built")},
        "box": box,
        "window": {"from": f"{WINDOW[0]:04d}-01-01", "to": f"{WINDOW[1]:04d}-01-01", "earliest": f"{EARLIEST:04d}-01-01"},
        "sources": [source_entry("melbourne-clue-buildings", BUILDINGS), source_entry("melbourne-clue-establishments", ESTABLISHMENTS),
                    source_entry("melbourne-development-activity", DEVELOPMENTS), source_entry("victorian-heritage-register", HERITAGE),
                    source_entry("melbourne-event-permits", PERMITS),
                    {"id": f"city-model-{AREA}", "title": "The city model's buildings and ways (data/products/street_records.py)",
                     "credit": model.get("credit"), "licence": "ODbL 1.0", "retrieved": model.get("built")}],
        "conventions": {
            "times": "from and to bound a label's resolution, [from, to), civil dates in Melbourne: a year 1 January to the next 1 January, a day that day to the next",
            "match": ("keys, any of which may be present: property ('base:<n>', a row of properties), building (a city-model building id), "
                      "rule (a rule coupling: 'vhr:<number>'), streets (names of the city model's ways), point ([x, y] local metres)"),
            "join": "inside: the footprint holds the point (holes respected); nearest:<m>: the nearest outline within 15 m; none",
            "ids": "censuses and establishments name their property by the id of its row in properties ('base:<n>')",
        },
        "properties": table(["id", "base_property_id", "property_ids", "block", "at", "building", "join"], properties),
        "censuses": table(["year", "property", "property_id", "construction_year", "refurbished_year", "floors_above", "predominant_use",
                           "name", "accessibility_rating", "bicycle_spaces"], census_rows),
        "establishments": table(["year", "property", "property_id", "name", "anzsic4"], est_rows),
        "anzsic": anzsic,
        "anzsic_divisions": {d: {"title": t, "subdivisions": f"{lo:02d}-{hi:02d}"} for d, lo, hi, t in DIVISIONS},
        "events": table(["id", "type", "kind", "from", "to", "match", "source", "provenance", "detail"], events),
        "report": {
            "counts": {"properties": len(properties), "censuses": len(census_rows), "establishments": len(est_rows), "events": len(events)},
            "events_by_type": dict(sorted(by_type.items())),
            "events_by_source": {s: dict(sorted(c.items())) for s, c in sorted(by_source.items())},
            "establishments_by_year": {str(y): n for y, n in sorted(est_years.items())},
            "censuses_by_year": {str(y): n for y, n in sorted(Counter(r[0] for r in census_rows).items())},
            "joins": {"inside": joins["inside"], "nearest": joins["nearest"], "none": joins["none"]},
            "developments_matched_by": dict(dam_match),
            "heritage_matched_by": dict(vhr_match),
            "permits_by_street": dict(permit_streets.most_common()),
            "dropped": dict(sorted(report["dropped"].items())),
            "unmatched_street_names": dict(unmatched.most_common()),
            "district_street_names": names,
            "notes": report["notes"],
        },
    }
    path = OUT / f"city-events-{AREA}.json"
    write(path, product)
    register("city-events", {"name": f"city-events-{AREA}", "file": path.name,
                             "title": "Melbourne's middle of the Hoddle Grid, 2002-2024: properties, censuses, establishments and dated events",
                             "domain": "urban-science", "window": [product["window"]["from"], product["window"]["to"]],
                             "properties": len(properties), "events": len(events),
                             "source": "City of Melbourne (CLUE, Development Activity Monitor, event permits); Heritage Council of Victoria"})
    print(f"  city-events-{AREA}: {len(properties)} properties (joined {dict(joins)}), {len(census_rows)} census rows, "
          f"{len(est_rows)} establishment-years, {len(events)} events {dict(by_type)}; {path.stat().st_size / 1e6:.2f} MB")
    return product


NEEDS = [CITY_MODEL, BUILDINGS / "buildings.csv", ESTABLISHMENTS / "establishments.csv", DEVELOPMENTS / "developments.csv",
         HERITAGE / "places.json", PERMITS / "permits.csv"]


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--offline", action="store_true", help="never fetch the heritage records the cache lacks")
    parser.add_argument("--refresh-details", action="store_true", help="fetch every heritage record again")
    args = parser.parse_args(argv)
    sys.stdout.reconfigure(errors="replace")
    lacking = [p for p in NEEDS if not p.exists()]
    if lacking:
        print("city-events: missing " + ", ".join(str(p.relative_to(DATA.parent)) for p in lacking) + "; run: python data/acquisition/fetch.py "
              "melbourne-clue-buildings melbourne-clue-establishments melbourne-development-activity victorian-heritage-register "
              "melbourne-event-permits (and python data/products/street_records.py --only city)", file=sys.stderr)
        return 1
    build(offline=args.offline, refresh_details=args.refresh_details)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
