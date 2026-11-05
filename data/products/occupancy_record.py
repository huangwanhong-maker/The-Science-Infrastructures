"""Build the record of the uses of the Hoddle Grid's buildings for the urban-science domain's `occupancy` (ADR 0071 and
ADR 0075; the design, development/ideas/2026-11-04-agents-and-decisions-for-urban-planning.md, §10): every CLUE base
property in the Hoddle Grid's box, its use at every census 2002-2024, its premises and their vacancy, the vacant floor
space of its block, its planning zone, and the permits for changes of use that name its address.

The record holds facts; the domain models them. Nothing below is inferred but the joins, each of which is counted.

1. **Properties.** One row a CLUE base property (the land parcel, the identifier that persists across censuses) of the
   building information whose point lies in the box of the CLUE manifests' profile (the Hoddle Grid's box, -37.8205 to
   -37.8075 S, 144.9535 to 144.9735 E). Its point is the most frequent point of the rows whose property_id is the base's
   own, else of all its rows (as city_events.py takes it), in the city model's local frame (east-north metres about the
   model's origin, by street_records.enu); its block likewise; its building the city model's outline that holds the
   point, else the nearest within 15 m (city_events.Footprints); its planning zone the Vicmap Planning polygon that
   holds the point (zone_code, e.g. CCZ1), else none.
2. **Censuses.** One row a base property and a census year: its use, the predominant space use of its smallest
   property_id with a use that year (the rule of the workbench's Event Spacetime engine, so that the memories agree);
   the number of its buildings (property_ids) that year; its floors above ground (the greatest of its buildings'); its
   construction year (the most frequent).
3. **Premises.** One row a base property and a census year: its establishments and, among them, CLUE's vacant space
   (ANZSIC 0000), from the business establishments (the municipality's file, by base_property_id).
4. **Demand.** One row a census year: the establishments of the municipality outside the box, occupied (not 0000) and
   vacant, and those inside it. The occupied establishments outside the box are the demand index the domain's model
   reads (an assumption declared there): the demand for premises in the city about the grid, not counting the grid's
   own.
5. **Blocks.** Every CLUE block whose polygon meets the box: its ring (local metres; the largest ring of a multipolygon)
   and small area; and its floor space at every census (floor-space-by-use-by-block): the total, the five unoccupied
   uses, and the part of the total that the suppressed cells hide (the total less the cells shown). Small cells are
   suppressed, so an unoccupied column is a lower bound.
6. **Zones.** Every zone polygon of the box (local rings), with its code, description and gazettal date.
7. **Permits.** The permits naming a change of use (melbourne-building-permits) whose address names a street and a
   number that CLUE's addresses in the box know: the address's street (abbreviations expanded) and locality, and its
   number or range, matched to the CLUE address of that street and locality whose range holds the number (an exact
   range first, else the narrowest range holding it); several documents of one matter (council_ref) at one property
   are one matter. Every outcome is counted.
8. **Footprints.** The joined buildings' outlines (local metres, to 0.1 m), for maps.

Output: data/cache/products/occupancy-hoddle-grid.json (compact JSON; tables as {"columns": [...], "rows": [[...]]}),
listed in the index under the kind `occupancy`.

Usage::

    python data/acquisition/fetch.py melbourne-clue-buildings melbourne-clue-establishments melbourne-clue-floor-space \\
        melbourne-clue-blocks melbourne-building-permits vicmap-planning-zones osm-city
    python data/products/street_records.py --only city          # the city model, if not yet built
    python data/products/occupancy_record.py
"""

from __future__ import annotations

import argparse
import datetime as dt
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
from city_events import Footprints, Frame, integer, join_label, mode, read_csv, table  # noqa: E402
from product_index import OUT, register  # noqa: E402
from street_records import contains, retrieved, signed_area, write  # noqa: E402

AREA = "hoddle-grid"
CITY_MODEL = OUT / "city-model-melbourne.json"
BUILDINGS = CACHE / "melbourne-clue-buildings" / "hoddle-grid"
ESTABLISHMENTS = CACHE / "melbourne-clue-establishments" / "all"
FLOOR_SPACE = CACHE / "melbourne-clue-floor-space" / "all"
BLOCKS = CACHE / "melbourne-clue-blocks" / "all"
PERMITS = CACHE / "melbourne-building-permits" / "change-of-use"
ZONES = CACHE / "vicmap-planning-zones" / "hoddle-grid"
REGISTRY = DATA / "registry"
BOX = {"south": -37.8205, "west": 144.9535, "north": -37.8075, "east": 144.9735}  # the CLUE manifests' profile
WINDOW = (2002, 2024)
UNOCCUPIED = ["unoccupied_unused", "unoccupied_under_renovation", "unoccupied_under_construction",
              "unoccupied_under_demolition_condemned", "unoccupied_undeveloped_site"]

# Street types as the City's addresses write them, abbreviated or not.
STREET_TYPES = {"street": "street", "st": "street", "lane": "lane", "ln": "lane", "place": "place", "pl": "place",
                "road": "road", "rd": "road", "avenue": "avenue", "ave": "avenue", "av": "avenue", "parade": "parade",
                "pde": "parade", "alley": "alley", "arcade": "arcade", "walk": "walk", "way": "way", "square": "square",
                "boulevard": "boulevard", "bvd": "boulevard", "blvd": "boulevard", "terrace": "terrace", "tce": "terrace",
                "court": "court", "ct": "court", "crescent": "crescent", "cres": "crescent", "drive": "drive", "dr": "drive",
                "highway": "highway", "hwy": "highway", "close": "close", "mews": "mews", "promenade": "promenade",
                "row": "row", "grove": "grove", "esplanade": "esplanade", "mall": "mall", "quay": "quay"}
NUMBERED = re.compile(r"(?<![\w/])(?:\d+[A-Za-z]?/)?(\d+)[A-Za-z]?(?:\s*-\s*(\d+)[A-Za-z]?)?\s+((?:[A-Za-z'.]+\s+){1,4}?)("
                      + "|".join(sorted(STREET_TYPES, key=len, reverse=True)) + r")\b\.?", re.IGNORECASE)
LOCALITY = re.compile(r"([A-Za-z][A-Za-z ]*?)\s+(?:VIC\s+)?(3\d{3})\s*$", re.IGNORECASE)


def street_key(words: str, kind: str) -> str:
    return " ".join(words.lower().replace(".", "").split()) + " " + STREET_TYPES[kind.lower()]


def parse_address(text: str) -> tuple[int, int, str, str] | None:
    """(first number, last number, street, locality) of an address, or None. The locality is the last component's
    suburb and postcode ("melbourne 3000"); the street the last numbered street before it."""
    if not text:
        return None
    parts = [p.strip() for p in text.split(",") if p.strip()]
    tail = parts[-1] if parts else ""
    loc = LOCALITY.search(tail)
    if loc:
        locality = f"{' '.join(loc.group(1).lower().split())} {loc.group(2)}"
        head = ", ".join(parts[:-1]) if len(parts) > 1 else tail[:loc.start()]
    else:
        return None
    found = list(NUMBERED.finditer(head))
    if not found:
        return None
    m = found[-1]
    a = int(m.group(1))
    b = int(m.group(2)) if m.group(2) else a
    if b < a:
        a, b = b, a
    return a, b, street_key(m.group(3), m.group(4)), locality


def parse_clue_address(text: str) -> tuple[int, int, str, str] | None:
    """The same of CLUE's addresses, written "363-397 Flinders Street MELBOURNE 3000" (no commas)."""
    if not text:
        return None
    m = re.match(r"^(.*?\b(?:" + "|".join(sorted(STREET_TYPES, key=len, reverse=True)) + r"))\s+([A-Za-z ]+?)\s+(3\d{3})\s*$", text.strip(), re.IGNORECASE)
    if not m:
        return None
    return parse_address(f"{m.group(1)}, {m.group(2)} VIC {m.group(3)}")


class AddressIndex:
    """CLUE's addresses by street and locality: each (first, last, base property) as written in some census."""

    def __init__(self):
        self.by_street: dict[tuple[str, str], set[tuple[int, int, str]]] = defaultdict(set)

    def add(self, address: str, base: str) -> bool:
        p = parse_clue_address(address)
        if not p:
            return False
        a, b, street, locality = p
        self.by_street[(street, locality)].add((a, b, base))
        return True

    def find(self, a: int, b: int, street: str, locality: str) -> tuple[str | None, str]:
        known = self.by_street.get((street, locality))
        if not known:
            return None, "street not in the box's CLUE addresses"
        exact = {base for lo, hi, base in known if (lo, hi) == (a, b)}
        if len(exact) == 1:
            return exact.pop(), "exact"
        if len(exact) > 1:
            return None, "ambiguous: one range at several properties"
        holding = sorted((hi - lo, base) for lo, hi, base in known if lo <= a and b <= hi)
        if holding:
            narrow = [base for w, base in holding if w == holding[0][0]]
            if len(set(narrow)) == 1:
                return narrow[0], "within a range"
            return None, "ambiguous: several narrowest ranges"
        meeting = sorted((hi - lo, base) for lo, hi, base in known if lo <= b and a <= hi)
        if meeting:
            narrow = [base for w, base in meeting if w == meeting[0][0]]
            if len(set(narrow)) == 1:
                return narrow[0], "overlapping a range"
        return None, "number not in the street's CLUE ranges"


def manifest(ident: str) -> dict:
    return json.loads((REGISTRY / f"{ident}.json").read_text(encoding="utf-8"))


def source_entry(ident: str, folder: Path) -> dict:
    m = manifest(ident)
    return {"id": ident, "title": m["title"], "credit": m["credit"], "licence": m["licence"], "retrieved": retrieved(folder)}


def local_ring(frame: Frame, ring_lonlat) -> list[list[float]]:
    lon = np.array([p[0] for p in ring_lonlat], float)
    lat = np.array([p[1] for p in ring_lonlat], float)
    x, y = frame.xy(lat, lon)
    pts = [[round(float(a), 1), round(float(b), 1)] for a, b in zip(x, y)]
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]
    return pts


def rings_of(geometry: dict) -> list:
    """The outer rings of a Polygon or MultiPolygon (lon, lat)."""
    if not geometry:
        return []
    if geometry["type"] == "Polygon":
        return [geometry["coordinates"][0]]
    if geometry["type"] == "MultiPolygon":
        return [poly[0] for poly in geometry["coordinates"]]
    return []


def build() -> dict:
    model = json.loads(CITY_MODEL.read_text(encoding="utf-8"))
    frame = Frame(model["origin"])
    report = {"dropped": Counter(), "notes": []}
    x0, y0 = frame.xy(BOX["south"], BOX["west"])
    x1, y1 = frame.xy(BOX["north"], BOX["east"])
    box = [float(min(x0, x1)), float(min(y0, y1)), float(max(x0, x1)), float(max(y0, y1))]

    # -- properties ---------------------------------------------------------------------------------------------------
    rows = read_csv(BUILDINGS / "buildings.csv")
    by_base: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        if r["base_property_id"]:
            by_base[r["base_property_id"]].append(r)
        else:
            report["dropped"]["building rows without base_property_id"] += 1

    def point_of(rs):
        own = [(float(r["latitude"]), float(r["longitude"])) for r in rs if r["latitude"] and r["longitude"] and r["property_id"] == r["base_property_id"]]
        every = [(float(r["latitude"]), float(r["longitude"])) for r in rs if r["latitude"] and r["longitude"]]
        return mode(own) or mode(every)

    bases = sorted(by_base, key=int)
    pts = {b: point_of(by_base[b]) for b in bases}
    lat = np.array([pts[b][0] if pts[b] else np.nan for b in bases])
    lon = np.array([pts[b][1] if pts[b] else np.nan for b in bases])
    px, py = frame.xy(lat, lon)
    at = {}
    for i, b in enumerate(bases):
        if not (np.isfinite(px[i]) and np.isfinite(py[i])):
            report["dropped"]["base properties without a point"] += 1
            continue
        if not (box[0] <= px[i] <= box[2] and box[1] <= py[i] <= box[3]):
            report["dropped"]["base properties whose point lies outside the box"] += 1
            continue
        at[b] = (round(float(px[i]), 2), round(float(py[i]), 2))

    zones_raw = json.loads((ZONES / "zones.geojson").read_text(encoding="utf-8"))
    zones = []
    for f in zones_raw["features"]:
        p = f["properties"]
        for ring in rings_of(f["geometry"]):
            zones.append({"code": p.get("zone_code"), "description": p.get("zone_description"), "gazetted": (p.get("gaz_begin_date") or "")[:10] or None,
                          "ring": local_ring(frame, ring)})

    def zone_of(p):
        holding = [z for z in zones if contains(z["ring"], p)]
        if not holding:
            return None
        holding.sort(key=lambda z: abs(signed_area(z["ring"])))
        return holding[0]["code"]

    footprints = Footprints(model["buildings"], box)
    outline = {b["id"]: b for b in footprints.items}
    properties, joins, zone_count, joined = [], Counter(), Counter(), set()
    for b in sorted(at, key=int):
        rs = by_base[b]
        block = mode(integer(r["block_id"]) for r in rs)
        building, kind, d = footprints.join(at[b])
        joins[kind] += 1
        if building:
            joined.add(building)
        z = zone_of(at[b])
        zone_count[z or "none"] += 1
        built = mode(integer(r["construction_year"]) for r in rs)
        properties.append([f"base:{b}", int(b), block, list(at[b]), building, join_label(kind, d), z, built])

    # -- censuses -----------------------------------------------------------------------------------------------------
    censuses = []
    use_counts = Counter()
    for b in sorted(at, key=int):
        per_year = defaultdict(list)
        for r in by_base[b]:
            per_year[int(r["census_year"][:4])].append(r)
        for year in sorted(per_year):
            rs = sorted(per_year[year], key=lambda r: integer(r["property_id"]) or 0)
            with_use = [r for r in rs if r["predominant_space_use"]]
            if not with_use:
                report["dropped"]["census-years of a property with no use"] += 1
                continue
            use = with_use[0]["predominant_space_use"]
            use_counts[use] += 1
            floors = max((integer(r["number_of_floors_above_ground"]) or 0) for r in rs) or None
            censuses.append([year, f"base:{b}", use, len({r["property_id"] for r in rs}), floors])

    # -- premises and demand ----------------------------------------------------------------------------------------
    est = read_csv(ESTABLISHMENTS / "establishments.csv")
    if est and "base_property_id" not in est[0]:
        raise SystemExit("the establishments lack base_property_id: run python data/acquisition/fetch.py melbourne-clue-establishments --refresh")
    premises = defaultdict(lambda: [0, 0])
    demand = defaultdict(lambda: [0, 0, 0, 0])  # outside: occupied, vacant; inside: occupied, vacant
    for r in est:
        year = integer(r["census_year"][:4] if r["census_year"] else None)
        if year is None:
            continue
        vacant = (r["industry_anzsic4_code"] or "").strip() == "0000"
        b = r["base_property_id"]
        if b in at:
            premises[(year, b)][0] += 1
            premises[(year, b)][1] += vacant
            demand[year][2 + vacant] += 1
        else:
            demand[year][vacant] += 1
    premises_rows = [[y, f"base:{b}", n, v] for (y, b), (n, v) in sorted(premises.items(), key=lambda kv: (kv[0][0], int(kv[0][1])))]
    demand_rows = [[y, *demand[y]] for y in sorted(demand)]

    # -- blocks and their floor space ---------------------------------------------------------------------------------
    blocks_raw = json.loads((BLOCKS / "blocks.geojson").read_text(encoding="utf-8"))
    blocks, in_box = [], set()
    for f in blocks_raw["features"]:
        p = f["properties"]
        rings = rings_of(f["geometry"])
        if not rings:
            continue
        local = [local_ring(frame, r) for r in rings]
        ring = max(local, key=lambda r: abs(signed_area(r)))
        xs, ys = [q[0] for q in ring], [q[1] for q in ring]
        if max(xs) < box[0] or min(xs) > box[2] or max(ys) < box[1] or min(ys) > box[3]:
            continue
        blocks.append([integer(str(p.get("block_id") or "")), p.get("area_name") or p.get("region_name"), ring])
        in_box.add(integer(str(p.get("block_id") or "")))
    blocks.sort(key=lambda r: r[0])
    floor = []
    for r in read_csv(FLOOR_SPACE / "floor-space.csv"):
        blk = integer(r["block_id"])
        if blk not in in_box:
            continue
        year = int(r["census_year"][:4])
        total = float(r["total_floor_space_in_block"]) if r["total_floor_space_in_block"] else None
        shown = sum(float(v) for k, v in r.items() if k not in ("census_year", "block_id", "clue_small_area", "total_floor_space_in_block") and v not in ("", None))
        cells = [integer(r[c]) for c in UNOCCUPIED]
        floor.append([year, blk, total, *cells, round(total - shown, 1) if total is not None else None])
    floor.sort(key=lambda r: (r[0], r[1]))

    # -- permits for changes of use -----------------------------------------------------------------------------------
    index = AddressIndex()
    unparsed = 0
    for b in at:
        for r in by_base[b]:
            if r["street_address"] and not index.add(r["street_address"], b):
                unparsed += 1
    report["notes"].append(f"{unparsed} CLUE building rows in the box have an address the index could not read")
    matched, outcome = {}, Counter()
    for r in read_csv(PERMITS / "permits.csv"):
        p = parse_address(r["address"])
        if not p:
            outcome["address without a numbered street and locality"] += 1
            continue
        base, how = index.find(*p)
        outcome[how] += 1
        if base is None:
            continue
        day = (r["issue_date"] or "")[:10]
        key = (r["council_ref"] or r["permit_number"], base)
        cost = integer(r["estimated_cost_of_works"])
        entry = matched.get(key)
        doc = [day, r["permit_certificate_type"] or None]
        if entry is None:
            matched[key] = {"first": day, "last": day, "property": f"base:{base}", "ref": key[0], "cost": cost, "address": r["address"],
                            "works": (r["desc_of_works"] or "")[:200], "documents": [doc], "match": how}
        else:
            entry["first"], entry["last"] = min(entry["first"], day), max(entry["last"], day)
            entry["cost"] = max(entry["cost"] or 0, cost or 0) or None
            entry["documents"].append(doc)
    permit_rows = [[m["ref"], m["property"], m["first"], m["last"], m["cost"], m["match"], m["address"], m["works"], len(m["documents"])]
                   for m in sorted(matched.values(), key=lambda m: (m["first"], m["ref"]))]

    shapes = [[bid, [[round(q[0], 1), round(q[1], 1)] for q in outline[bid]["footprint"]]] for bid in sorted(joined)]

    sources = [source_entry(i, f) for i, f in [("melbourne-clue-buildings", BUILDINGS), ("melbourne-clue-establishments", ESTABLISHMENTS),
                                                 ("melbourne-clue-floor-space", FLOOR_SPACE), ("melbourne-clue-blocks", BLOCKS),
                                                 ("melbourne-building-permits", PERMITS), ("vicmap-planning-zones", ZONES)]]
    sources.append({"id": "osm-city", "title": "The city model's buildings (data/products/street_records.py)", "credit": "© OpenStreetMap contributors, ODbL 1.0",
                    "licence": "ODbL 1.0", "built": model.get("built")})
    years = sorted({c[0] for c in censuses})
    report["notes"].append(f"census years {years[0]}-{years[-1]}; {len(set(c[1] for c in censuses))} properties with a use")
    report["notes"].append("zones of the properties: " + ", ".join(f"{k} {v}" for k, v in zone_count.most_common()))
    return {
        "product": "occupancy",
        "id": AREA,
        "built": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "title": "The uses of the Hoddle Grid's buildings at every census, their premises and vacancy, their blocks' vacant floor space, zones and permits for changes of use",
        "frame": "local ENU metres from origin: x east, y north",
        "origin": model["origin"],
        "box": {"lonlat": [BOX["west"], BOX["south"], BOX["east"], BOX["north"]], "local": [round(v, 2) for v in box]},
        "window": {"from": WINDOW[0], "to": WINDOW[1]},
        "sources": sources,
        "conventions": {
            "use": "a base property's use at a census is the predominant space use of its smallest property_id with a use that year (the rule of the workbench's Event Spacetime engine)",
            "vacant": "CLUE's family Unoccupied (Unused, Under Renovation, Under Construction, Under Demolition/Condemned, Undeveloped Site); a premises is vacant when its ANZSIC code is 0000",
            "floor_space": "m2; an unoccupied cell is blank (null) where it pertains to two or fewer businesses; hidden is the block's total less the cells shown",
            "demand": "establishments outside the box: the demand index of the domain's model (an assumption, declared there)",
        },
        "properties": table(["id", "base_property_id", "block", "at", "building", "join", "zone", "construction_year"], properties),
        "censuses": table(["year", "property", "use", "buildings", "floors"], censuses),
        "premises": table(["year", "property", "establishments", "vacant"], premises_rows),
        "demand": table(["year", "outside_occupied", "outside_vacant", "inside_occupied", "inside_vacant"], demand_rows),
        "blocks": table(["block", "small_area", "ring"], blocks),
        "floor_space": table(["year", "block", "total", *UNOCCUPIED, "hidden"], floor),
        "zones": table(["code", "description", "gazetted", "ring"], [[z["code"], z["description"], z["gazetted"], z["ring"]] for z in zones]),
        "permits": table(["ref", "property", "first", "last", "cost", "match", "address", "works", "documents"], permit_rows),
        "footprints": table(["building", "ring"], shapes),
        "report": {"counts": {"properties": len(properties), "censuses": len(censuses), "premises": len(premises_rows), "blocks": len(blocks),
                              "floor_space": len(floor), "zones": len(zones), "permits": len(permit_rows), "footprints": len(shapes)},
                   "joins": dict(joins), "zones": dict(zone_count), "uses": dict(use_counts.most_common()),
                   "permits": dict(outcome), "dropped": dict(report["dropped"]), "notes": report["notes"]},
    }


NEEDS = [CITY_MODEL, BUILDINGS / "buildings.csv", ESTABLISHMENTS / "establishments.csv", FLOOR_SPACE / "floor-space.csv",
         BLOCKS / "blocks.geojson", PERMITS / "permits.csv", ZONES / "zones.geojson"]


def main(argv: list[str]) -> int:
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args(argv)
    sys.stdout.reconfigure(errors="replace")
    lacking = [p for p in NEEDS if not p.exists()]
    if lacking:
        print("occupancy: missing " + ", ".join(str(p.relative_to(DATA.parent)) for p in lacking) + "; run: python data/acquisition/fetch.py "
              "melbourne-clue-buildings melbourne-clue-establishments melbourne-clue-floor-space melbourne-clue-blocks "
              "melbourne-building-permits vicmap-planning-zones (and python data/products/street_records.py --only city)", file=sys.stderr)
        return 1
    product = build()
    path = OUT / f"occupancy-{AREA}.json"
    write(path, product)
    register("occupancy", {"name": f"occupancy-{AREA}", "file": path.name, "title": product["title"], "built": product["built"],
                           "sources": [s["id"] for s in product["sources"]]})
    c = product["report"]["counts"]
    print(f"  occupancy-{AREA}: {c['properties']} properties (joined {product['report']['joins']}), {c['censuses']} census rows, "
          f"{c['premises']} premises-years, {c['blocks']} blocks, {c['floor_space']} block-years, {c['zones']} zone polygons, "
          f"{c['permits']} permits matched ({product['report']['permits']}); {path.stat().st_size / 1e6:.2f} MB")
    for n in product["report"]["notes"]:
        print("   ", n)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
