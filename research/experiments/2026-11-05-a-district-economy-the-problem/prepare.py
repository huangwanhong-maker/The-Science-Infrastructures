"""Prepares the inputs of the shared problem "A district's economy" from the data library's cache.

The problem is shared by three case studies that differ only in how much structure their model of the city has
(README §1): without a model of flow, with a fluid model of flow and its circulation, and with a field theory of
flow given by an action. This script gathers what all three read; the C++ program `calibrate` then builds the
district, its frontages and their relations with the domains, calibrates the people's rhythm from the counts, and
writes the machine-readable problem (results/problem.json).

1. **The district.** The case study "Designing a district" (research/experiments/2026-11-05-designing-a-district):
   its walking network (Landscape Dynamics' district, as the Urban Science Workbench builds it), its 61 parcels
   (the buildings standing in the blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street),
   the 524 buildings of their context, and the open spaces inside the area. Its preparation is run if its
   results/district-inputs.json is missing, and its output is read, not made again here, so that the two studies
   stand on the same cells.
2. **The sensors.** The City of Melbourne's pedestrian counters in the Hoddle Grid's box (the data library's
   street-sensors product), placed in the city model's local frame by the data library's exact east-north on WGS 84
   (street_records.enu), each with the kind and width of the nearest way of the city model and its coverage of the
   reference year: the 52 weeks from Monday 27 October 2025 to Sunday 25 October 2026 (civil time).
3. **The businesses.** The City's Census of Land Use and Employment (CLUE): the establishments of its latest census
   (2024) in the Hoddle Grid's box, each with its ANZSIC class mapped to the problem's kinds (KIND_OF below), placed in
   the local frame; and the establishments of each parcel (by the base property the data library's city-events product
   joins to the building), from which its ground floor's kind as built is read (README §3.2).
4. **Typical hours.** The opening hours declared in OpenStreetMap for the holdings of five Melbourne streets (the data
   library's street-holdings product): the median opening and closing hour of each kind on a weekday (Wednesday) and
   on a Saturday, which stand for the hours as built of the parcels (README §3.3).

Usage (from the repository's root, with the repository's Python):

    .venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-the-problem/prepare.py

It writes results/problem-inputs.json. That file holds extracts of OpenStreetMap (ODbL 1.0): it is not committed
(.gitignore), and this script makes it again from the cache.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGNING = ROOT / "research" / "experiments" / "2026-11-05-designing-a-district"
CACHE = ROOT / "data" / "cache"
PRODUCTS = CACHE / "products"
sys.path.insert(0, str(ROOT / "data" / "products"))
from street_records import enu  # noqa: E402  (the data library's exact local frame)

# The kinds of a ground-floor frontage (README §3.1). Their order is the problem's encoding.
KINDS = ["vacant", "retail", "food and drink", "services", "office lobby", "culture"]
VACANT, RETAIL, FOOD, SERVICES, LOBBY, CULTURE = range(6)

# The Hoddle Grid's box (the data library's melbourne-clue-buildings profile hoddle-grid).
HODDLE = {"south": -37.8205, "west": 144.9535, "north": -37.8075, "east": 144.9735}

# The reference year: 52 weeks of civil time from Monday 27 October 2025 00:00 AEDT (UTC+11).
REFERENCE_START_UTC = dt.datetime(2025, 10, 26, 13, 0, tzinfo=dt.timezone.utc)
REFERENCE_WEEKS = 52


def kind_of(anzsic4: str) -> int:
    """A CLUE establishment's kind from its ANZSIC 2006 class (four digits).

    retail          division G, Retail Trade (subdivisions 39-43);
    food and drink  subdivision 45, Food and Beverage Services (cafes, restaurants, takeaway, pubs, bars, clubs);
    services        personal and household services that face a street: K 62 (banking: branches), K 641-642 (auxiliary
                    finance: money changers, brokers' shopfronts), L 672 (real estate agents), N 722 (travel agents),
                    Q 85 (medical, dental, allied health), S 94-95 (repair, personal care: hairdressing, beauty,
                    laundries), I 51 (postal: post offices);
    culture         division R, Arts and Recreation Services (89-92: museums, galleries, performing arts, sport and
                    recreation, gambling), and J 551 motion picture exhibition (cinemas);
    vacant          CLUE's 0000, Vacant Space;
    office lobby    every other establishment: offices, education, public administration, accommodation and the rest,
                    whose ground floor faces the street as a building's entrance.
    """
    c = (anzsic4 or "").strip()
    if not c.isdigit():
        return LOBBY
    if c == "0000":
        return VACANT
    sub, grp = int(c[:2]), c[:3]
    if 39 <= sub <= 43:
        return RETAIL
    if sub == 45:
        return FOOD
    if sub in (62, 85, 94, 95, 51) or grp in ("641", "642", "672", "722"):
        return SERVICES
    if 89 <= sub <= 92 or grp == "551":
        return CULTURE
    return LOBBY


# The street-holdings product's uses mapped to the kinds (for the typical hours).
HOLDING_KIND = {
    "food-and-drink": FOOD, "clothing-and-fashion": RETAIL, "jewellery": RETAIL, "bridal-and-wedding": RETAIL,
    "other-retail": RETAIL, "personal-services": SERVICES, "health": SERVICES, "finance": SERVICES,
    "entertainment-and-culture": CULTURE, "offices-and-professional": LOBBY, "education": LOBBY,
    "community-and-public": LOBBY, "transport": LOBBY, "wholesale-and-industry": LOBBY, "accommodation": LOBBY,
    "vacant": VACANT, "other": LOBBY,
}

# The as-built kind of a ground floor from its parcel's establishments (README §3.2): the street-facing kind with
# most establishments (ties in this order of precedence), else an office lobby if the parcel has any establishment
# but vacant space, else vacant.
PRECEDENCE = [FOOD, RETAIL, SERVICES, CULTURE]


def as_built_kind(kinds: list[int]) -> int:
    active = Counter(k for k in kinds if k in PRECEDENCE)
    if active:
        top = max(active.values())
        return next(k for k in PRECEDENCE if active.get(k, 0) == top)
    if any(k == LOBBY for k in kinds):
        return LOBBY
    return VACANT


def contains(ring, p) -> bool:
    """Even-odd rule."""
    x, y = p
    c = False
    for i in range(len(ring)):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % len(ring)]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            c = not c
    return c


def local(lat: float, lon: float, origin: dict) -> tuple[float, float]:
    e, n, _ = enu(lat, lon, origin["lat"], origin["lon"])
    return float(e), float(n)


def nearest_way(p, roads):
    """The nearest way of the city model to a point: (distance, kind, name, width)."""
    best = (math.inf, None, None, None)
    for r in roads:
        line = r["line"]
        for a, b in zip(line, line[1:]):
            ex, ey = b[0] - a[0], b[1] - a[1]
            L2 = ex * ex + ey * ey
            t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / L2))
            d = math.hypot(p[0] - a[0] - t * ex, p[1] - a[1] - t * ey)
            if d < best[0]:
                best = (d, r["kind"], r.get("name") or "", r.get("width_m"))
    return best


# The ways people walk (the designing study's WALK) and their kinds as the network names them.
WALK = {"primary", "secondary", "tertiary", "residential", "unclassified", "living_street", "pedestrian", "service",
        "primary_link", "secondary_link", "tertiary_link"}


def way_kind(osm_kind: str | None) -> str:
    if osm_kind == "service":
        return "laneway"
    if osm_kind in ("pedestrian", "living_street"):
        return "mall"
    return "street"


def utc(s: str) -> dt.datetime:
    return dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)


def coverage(series: dict, start: dt.datetime, hours: int) -> float:
    t0 = utc(series["t0_utc"])
    i0 = int((start - t0).total_seconds() // 3600)
    counts = series["counts"]
    have = sum(1 for i in range(i0, i0 + hours) if 0 <= i < len(counts) and counts[i] is not None)
    return have / hours


def typical_hours(holdings: dict) -> dict:
    """The median opening and closing hour of each kind on a Wednesday and on a Saturday, from OpenStreetMap's
    declared hours on the product's Melbourne streets; closed days counted as such."""
    by = defaultdict(lambda: {"weekday": [], "weekend": [], "closed_weekend": 0, "holdings": 0})
    for st in holdings["streets"]:
        if st["city"] != "Melbourne":
            continue
        for h in st["holdings"]:
            if not h.get("hours") or h.get("hours_status") not in (None, "ok"):
                continue
            k = HOLDING_KIND.get(h["use"], LOBBY)
            days = defaultdict(list)
            for day, a, b in h["hours"]:
                days[day].append((a, b))
            if 2 not in days:   # no Wednesday: not a weekday trader as declared
                continue
            by[k]["holdings"] += 1
            wd = (min(a for a, _ in days[2]), max(b for _, b in days[2]))
            by[k]["weekday"].append(wd)
            if 5 in days:
                by[k]["weekend"].append((min(a for a, _ in days[5]), max(b for _, b in days[5])))
            else:
                by[k]["closed_weekend"] += 1
    out = {}
    for k in by:
        def med(pairs):
            if not pairs:
                return None
            o = statistics.median(a for a, _ in pairs) / 60
            c = min(24.0, statistics.median(b for _, b in pairs) / 60)
            return [int(round(o)), int(round(c))]
        out[KINDS[k]] = {"weekday": med(by[k]["weekday"]), "weekend": med(by[k]["weekend"]), "holdings": by[k]["holdings"],
                         "closed_on_saturday": by[k]["closed_weekend"], "from": "OpenStreetMap's declared hours, median"}
    # Fewer than five holdings declare a kind's hours: its hours are an assumption (README §3.3).
    for k, (wd, we) in ASSUMED_HOURS.items():
        if KINDS[k] not in out or out[KINDS[k]]["holdings"] < 5:
            n = out.get(KINDS[k], {}).get("holdings", 0)
            out[KINDS[k]] = {"weekday": wd, "weekend": we, "holdings": n, "closed_on_saturday": 0, "from": f"an assumption ({n} holdings declare hours)"}
    return out


# The hours of a kind that fewer than five holdings declare: assumptions.
ASSUMED_HOURS = {CULTURE: ([10, 22], [10, 22]), LOBBY: ([8, 18], [0, 0])}


if __name__ == "__main__":
    inputs = DESIGNING / "results" / "district-inputs.json"
    if not inputs.exists():
        print("making the designing study's district-inputs.json")
        subprocess.run([sys.executable, str(DESIGNING / "prepare.py")], check=True)
    D = json.loads(inputs.read_text(encoding="utf-8"))
    origin = D["origin"]
    model = json.loads((PRODUCTS / "city-model-melbourne.json").read_text(encoding="utf-8"))
    walk_roads = [r for r in model["roads"] if r.get("kind") in WALK and r.get("line")]

    # ---- the sensors ----
    S = json.loads((PRODUCTS / "street-sensors.json").read_text(encoding="utf-8"))
    sensors = []
    hours = REFERENCE_WEEKS * 168
    for s in S["sensors"]:
        if s.get("lat") is None or not s.get("series"):
            continue
        if not (HODDLE["south"] <= s["lat"] <= HODDLE["north"] and HODDLE["west"] <= s["lon"] <= HODDLE["east"]):
            continue
        series = json.loads((PRODUCTS / s["series"]).read_text(encoding="utf-8"))
        x, y = local(s["lat"], s["lon"], origin)
        d, kind, name, width = nearest_way((x, y), walk_roads)
        sensors.append({"id": s["id"], "code": s["code"], "name": s["name"], "status": s["status"], "x": round(x, 2), "y": round(y, 2),
                        "series": "data/cache/products/" + s["series"], "coverage": round(coverage(series, REFERENCE_START_UTC, hours), 4),
                        "way": {"kind": way_kind(kind), "osm_kind": kind, "name": name, "width_m": width, "distance_m": round(d, 1)},
                        "directions": s.get("directions")})
    print(f"{len(sensors)} sensors in the Hoddle Grid's box; {sum(1 for s in sensors if s['coverage'] >= 0.8)} cover 80% of the reference year")

    # ---- the businesses: CLUE 2024 in the Hoddle Grid's box ----
    rows = list(csv.DictReader((CACHE / "melbourne-clue-establishments" / "all" / "establishments.csv").open(encoding="utf-8-sig"), delimiter=";"))
    latest = max(int(r["census_year"]) for r in rows)
    est = []
    for r in rows:
        if int(r["census_year"]) != latest or not r["latitude"] or not r["longitude"]:
            continue
        lat, lon = float(r["latitude"]), float(r["longitude"])
        if not (HODDLE["south"] <= lat <= HODDLE["north"] and HODDLE["west"] <= lon <= HODDLE["east"]):
            continue
        x, y = local(lat, lon, origin)
        est.append([round(x, 2), round(y, 2), kind_of(r["industry_anzsic4_code"])])
    print(f"CLUE {latest}: {len(est)} establishments in the Hoddle Grid's box,", dict(Counter(KINDS[k] for *_, k in est)))

    # ---- each parcel's establishments and its kind as built (city-events joins CLUE base properties to buildings) ----
    events = json.loads((PRODUCTS / "city-events-melbourne.json").read_text(encoding="utf-8"))
    cols = events["establishments"]["columns"]
    by_property = defaultdict(list)
    for r in events["establishments"]["rows"]:
        row = dict(zip(cols, r))
        if row["year"] == latest:
            by_property[row["property"]].append(kind_of(row["anzsic4"]))
    pcols = events["properties"]["columns"]
    properties_of = defaultdict(list)
    for r in events["properties"]["rows"]:
        row = dict(zip(pcols, r))
        if row["building"]:
            properties_of[row["building"]].append(row["id"])
    parcels = []
    for p in D["parcels"]:
        kinds, read = [], "joined"
        for prop in properties_of.get(p["id"], []):
            kinds += by_property.get(prop, [])
        if not kinds:   # no establishment by the join: those whose CLUE position lies inside the footprint
            kinds = [k for x, y, k in est if contains(p["footprint"], (x, y))]
            read = "inside" if kinds else "none"
        k = as_built_kind(kinds)
        parcels.append({"id": p["id"], "establishments": dict(Counter(KINDS[x] for x in kinds)), "kind": k, "read": read})
    print("as built:", dict(Counter(KINDS[p["kind"]] for p in parcels)))

    # ---- typical hours ----
    holdings = json.loads((PRODUCTS / "street-holdings.json").read_text(encoding="utf-8"))
    typical = typical_hours(holdings)
    print("typical hours:", json.dumps(typical))

    out = {
        "problem": "a-district-economy",
        "district": {k: D[k] for k in ("frame", "origin", "projection", "district_box", "context_box", "redesign_quad", "corners", "network",
                                        "parcels", "context", "open_spaces", "spaces_near", "storey_m")},
        "kinds": KINDS,
        "parcels": parcels,
        "sensors": sensors,
        "establishments": {"year": latest, "columns": ["x", "y", "kind"], "rows": est},
        "typical_hours": typical,
        "reference": {"start_utc": REFERENCE_START_UTC.strftime("%Y-%m-%dT%H:%M:%SZ"), "weeks": REFERENCE_WEEKS,
                      "civil": "Monday 2025-10-27 00:00 AEDT to Sunday 2026-10-25 24:00 AEDT"},
        "files": {"holidays": "data/cache/products/street-holidays-victoria.json", "weather": "data/cache/products/street-weather-melbourne.json"},
        "sources": {
            "district": "research/experiments/2026-11-05-designing-a-district/results/district-inputs.json (city model: (c) OpenStreetMap contributors, ODbL 1.0; floors: CLUE, CC BY 4.0)",
            "sensors": {"product": "data/cache/products/street-sensors.json", "built": S["built"], "credit": "City of Melbourne, Pedestrian Counting System",
                        "licence": "the counts' metadata state no licence; the City's open data are CC BY 4.0, and they are used on that understanding",
                        "retrieved": [x.get("retrieved_at") for x in S["sources"]]},
            "establishments": {"dataset": "melbourne-clue-establishments", "credit": "City of Melbourne, CLUE", "licence": "CC BY 4.0",
                               "retrieved": json.loads((CACHE / "melbourne-clue-establishments" / "all" / "retrieval.json").read_text())["retrieved_at"]},
            "joins": {"product": "data/cache/products/city-events-melbourne.json", "built": events["built"]},
            "hours": {"product": "data/cache/products/street-holdings.json", "built": holdings["built"], "credit": "(c) OpenStreetMap contributors, ODbL 1.0"},
        },
    }
    path = HERE / "results" / "problem-inputs.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print("wrote", path, f"{path.stat().st_size / 1e6:.2f} MB")
