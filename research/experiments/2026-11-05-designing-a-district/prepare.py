"""Prepares the inputs of the case study "Designing a district" from the data library's products.

The district is built as the Urban Science Workbench builds Landscape Dynamics' district, so that the study works on
the same cells the domains already know:

1. **The box.** Landscape Dynamics' district (the Hoddle Grid's middle, Elizabeth to Russell Street, Flinders to
   Lonsdale Street), as the data library's city-events product records it (city_events.py, DISTRICT), its corners
   carried into the city model's local frame by street_records.enu (exact east-north on WGS 84 about the model's
   origin), as the workbench's boxOf does (applications/urban-science-workbench/src/district.js).
2. **The walking network.** The city model's ways where people walk, cut to the box, split at every point two ways
   share and at their ends; a closed piece split at its middle; a piece seen twice kept once. This is networkOf of
   src/district.js, line for line (WALK, the 0.01 m keys, the signatures, the renumbering). The study does not import
   the workbench's code (nothing depends on an application, CLAUDE.md §3.4); it follows it.
3. **The buildings.** The city model's outlines whose footprint's centre lies in the box (footprintsIn of
   src/district.js), without the building:part entries (a part lies inside its outline; a redesign gives a parcel one
   height). Each is joined to the CLUE base property the city-events product joins to it (the City of Melbourne's
   Census of Land Use and Employment): its floors above ground and its predominant use at the latest census.
4. **The redesign area.** The few blocks the study redesigns: the buildings whose centre lies inside the
   quadrilateral of four street corners named below. The rest of the district stands as built, the context.
5. **The Sun.** Its altitude and azimuth at each half hour of the winter solstice, from the astronomy domain's
   subsolar point (domains/astronomy/js/sun.js), as the workbench's clock takes it (sunByHour), by sun.mjs.
6. **The air.** The mean of each civil hour over the Junes of 2009-2025 at Melbourne (Open-Meteo's ERA5 series in
   the data library's street-weather-melbourne product), and the ground at the record's mean.
7. **The bound on storeys.** The 90th percentile of the floors of the developments the City's Development Activity
   Monitor records as completed in the Hoddle Grid's box.

Usage (from the repository's root, with the repository's Python):

    .venv/Scripts/python.exe research/experiments/2026-11-05-designing-a-district/prepare.py

It writes results/district-inputs.json, which the C++ study reads.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "data" / "products"))
from street_records import enu  # noqa: E402  (the data library's exact local frame)

PRODUCTS = ROOT / "data" / "cache" / "products"
CACHE = ROOT / "data" / "cache"

# Landscape Dynamics' district, WGS 84 degrees (applications/urban-science-workbench/areas.json, melbourne.district;
# the same as data/products/city_events.py's DISTRICT).
DISTRICT = {"south": -37.8173, "west": 144.9617, "north": -37.8117, "east": 144.9697}

# The ways of the walking network (src/district.js, WALK).
WALK = {"primary", "secondary", "tertiary", "residential", "unclassified", "living_street", "pedestrian", "service",
        "primary_link", "secondary_link", "tertiary_link"}

STOREY_M = 3.2   # the city model's storey (its convention: building:levels x 3.2 m)


def js_round(x: float) -> int:
    """JavaScript's Math.round: halves toward +infinity."""
    return math.floor(x + 0.5)


def box_of(origin: dict, d: dict) -> dict:
    """The district's box in the local frame (src/district.js, boxOf)."""
    x0, y0, _ = enu(d["south"], d["west"], origin["lat"], origin["lon"])
    x1, y1, _ = enu(d["north"], d["east"], origin["lat"], origin["lon"])
    x0, y0, x1, y1 = float(x0), float(y0), float(x1), float(y1)
    return {"x0": min(x0, x1), "y0": min(y0, y1), "x1": max(x0, x1), "y1": max(y0, y1)}


def inside(p, B) -> bool:
    return B["x0"] <= p[0] <= B["x1"] and B["y0"] <= p[1] <= B["y1"]


def network_of(model: dict, B: dict) -> dict:
    """The walking network in the box: src/district.js, networkOf, line for line."""
    runs = []
    for r in model.get("roads") or []:
        if r.get("kind") not in WALK or not r.get("line") or len(r["line"]) < 2:
            continue
        run = []
        for p in r["line"]:
            if inside(p, B):
                run.append(p)
            else:
                if len(run) >= 2:
                    runs.append({"line": run, "kind": r["kind"], "name": r.get("name") or ""})
                run = []
        if len(run) >= 2:
            runs.append({"line": run, "kind": r["kind"], "name": r.get("name") or ""})

    def key(p):
        return f"{js_round(p[0] * 100)},{js_round(p[1] * 100)}"

    uses = Counter()
    for r in runs:
        for k in {key(p) for p in r["line"]}:
            uses[k] += 1
    junctions, index = [], {}

    def junction(p):
        k = key(p)
        if k not in index:
            index[k] = len(junctions)
            junctions.append([round(p[0], 2), round(p[1], 2)])
        return index[k]

    segments, seen = [], set()

    def push(line, kind, name):
        if len(line) < 2:
            return
        u, v = junction(line[0]), junction(line[-1])
        if u == v:
            if len(line) < 3:
                return
            m = len(line) // 2
            push(line[: m + 1], kind, name)
            push(line[m:], kind, name)
            return
        sig = f"{min(u, v)}-{max(u, v)}-{len(line)}-{key(line[len(line) // 2])}"
        if sig in seen:
            return
        seen.add(sig)
        kind_of = "laneway" if kind == "service" else "mall" if kind in ("pedestrian", "living_street") else "street"
        segments.append({"u": u, "v": v, "line": [[round(p[0], 2), round(p[1], 2)] for p in line], "kind": kind_of, "name": name})

    for r in runs:
        start = 0
        for k in range(1, len(r["line"])):
            last = k == len(r["line"]) - 1
            if last or uses[key(r["line"][k])] > 1:
                push(r["line"][start: k + 1], r["kind"], r["name"])
                start = k
    used = {s["u"] for s in segments} | {s["v"] for s in segments}
    renumber, J = {}, []
    for i, p in enumerate(junctions):
        if i in used:
            renumber[i] = len(J)
            J.append(p)
    for s in segments:
        s["u"], s["v"] = renumber[s["u"]], renumber[s["v"]]
    return {"junctions": J, "segments": segments}


def centre(f):
    return [sum(p[0] for p in f) / len(f), sum(p[1] for p in f) / len(f)]


def contains(ring, p) -> bool:
    """Even-odd rule."""
    x, y = p
    c = False
    n = len(ring)
    for i in range(n):
        x1, y1 = ring[i]
        x2, y2 = ring[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            c = not c
    return c


def signed_area(r) -> float:
    return 0.5 * sum(r[i][0] * r[(i + 1) % len(r)][1] - r[(i + 1) % len(r)][0] * r[i][1] for i in range(len(r)))


# ---- uses ----------------------------------------------------------------------------------------------------

# The study's four uses of a parcel. The first three are floor uses; the fourth gives the parcel to the public.
USES = ["dwellings", "offices", "shops", "open space"]
# CLUE's predominant space use mapped to the three floor uses: homes; the places a resident walks to for daily needs
# and leisure (retail, food, entertainment, performance); and every other workplace, vacant floor included.
DWELLINGS = {"Residential Apartment", "House/Townhouse", "Student Accommodation", "Institutional Accommodation"}
SHOPS = {"Retail - Shop", "Retail - Showroom", "Retail - Stall", "Entertainment/Recreation - Indoor",
         "Performances, Conferences, Ceremonies", "Public Display Area"}


def use_class(clue_use: str | None) -> int:
    if clue_use in DWELLINGS:
        return 0
    if clue_use in SHOPS:
        return 2
    return 1


# The redesign area: the blocks between Swanston and Russell Streets, Flinders Lane and Bourke Street, bounded by
# the junctions where those streets meet (where a street is two carriageways, the junction nearest the area).
CORNERS = [("Swanston Street", "Flinders Lane"), ("Russell Street", "Flinders Lane"), ("Russell Street", "Bourke Street"),
           ("Swanston Street", "Bourke Street")]

# The days and hours (Australian Eastern Standard Time, UTC+10; Melbourne keeps no daylight saving in June).
SOLSTICE = dt.date(2026, 6, 21)
AEST = dt.timezone(dt.timedelta(hours=10))
SUN_WINDOW = [10.25 + 0.5 * k for k in range(8)]   # 10:15 ... 13:45, each standing for half an hour of 10:00-14:00
HEAT_HOURS = [h + 0.5 for h in range(24)]          # the middle of each civil hour
DRAW_HOUR = 12.0


def instants(hours):
    out = []
    for h in hours:
        t = dt.datetime.combine(SOLSTICE, dt.time(0, 0), AEST) + dt.timedelta(hours=h)
        out.append(t.astimezone(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"))
    return out


def sun_at(lat, lon, hours):
    req = json.dumps({"lat": lat, "lon": lon, "instants": instants(hours)})
    out = subprocess.run(["node", str(HERE / "sun.mjs"), req], capture_output=True, text=True, check=True).stdout
    return [{"hour": h, "altitude_deg": s["altitude_deg"], "azimuth_deg": s["azimuth_deg"]} for h, s in zip(hours, json.loads(out))]


def june_air(weather: dict) -> tuple[list[float], float, dict]:
    """The mean air temperature of each AEST hour over the Junes of the record's whole years, and the record's mean."""
    t0 = dt.datetime.strptime(weather["t0_utc"], "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)
    T = weather["temperature_c"]
    sums, counts, years = [0.0] * 24, [0] * 24, set()
    total, n = 0.0, 0
    for i, v in enumerate(T):
        if v is None:
            continue
        local = (t0 + dt.timedelta(seconds=i * weather["step_s"])).astimezone(AEST)
        total += v
        n += 1
        if local.month == 6 and local.year <= 2025:
            sums[local.hour] += v
            counts[local.hour] += 1
            years.add(local.year)
    return [s / c for s, c in zip(sums, counts)], total / n, {"june_years": [min(years), max(years)], "hours_in_june": sum(counts), "hours": n}


def floors_bound(path: Path, q: float = 0.9) -> tuple[int, int]:
    rows = list(csv.DictReader(path.open(encoding="utf-8-sig"), delimiter=";"))
    fl = sorted(int(r["floors_above"]) for r in rows if r["status"] == "COMPLETED" and r["floors_above"])
    return int(np.quantile(fl, q, method="lower")), len(fl)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(HERE / "results" / "district-inputs.json"))
    ap.add_argument("--scratch", default=None, help="also write the whole network here, for a look")
    args = ap.parse_args()
    model = json.loads((PRODUCTS / "city-model-melbourne.json").read_text(encoding="utf-8"))
    events = json.loads((PRODUCTS / "city-events-melbourne.json").read_text(encoding="utf-8"))
    weather = json.loads((PRODUCTS / "street-weather-melbourne.json").read_text(encoding="utf-8"))
    B = box_of(model["origin"], DISTRICT)
    net = network_of(model, B)
    print(f"box {B}; {len(net['junctions'])} junctions, {len(net['segments'])} segments")
    if args.scratch:
        Path(args.scratch).write_text(json.dumps({"box": B, "network": net}), encoding="utf-8")

    # The corners of the redesign area.
    at = defaultdict(set)
    for s in net["segments"]:
        at[s["u"]].add(s["name"]), at[s["v"]].add(s["name"])
    found = [[j for j, names in at.items() if a in names and b in names] for a, b in CORNERS]
    single = [net["junctions"][c[0]] for c in found if len(c) == 1]
    mid = [sum(p[0] for p in single) / len(single), sum(p[1] for p in single) / len(single)]
    quad = [net["junctions"][min(c, key=lambda j: math.dist(net["junctions"][j], mid))] for c in found]
    print("redesign corners", quad)

    # CLUE at the latest census of each property, and the properties joined to each building.
    cen_cols = events["censuses"]["columns"]
    latest = {}
    for r in events["censuses"]["rows"]:
        row = dict(zip(cen_cols, r))
        p = row["property"]
        latest.setdefault(p, []).append(row)
    for p, rows in latest.items():
        y = max(r["year"] for r in rows)
        latest[p] = [r for r in rows if r["year"] == y]
    prop_cols = events["properties"]["columns"]
    of_building = defaultdict(list)
    for r in events["properties"]["rows"]:
        row = dict(zip(prop_cols, r))
        if row["building"]:
            of_building[row["building"]].append(row)

    def clue_of(bid):
        best = None
        for prow in of_building.get(bid, []):
            for r in latest.get(prow["id"], []):
                if r["floors_above"] is None:
                    continue
                if best is None or r["floors_above"] > best["floors_above"]:
                    best = {**r, "base": prow["id"], "join": prow["join"]}
        return best

    # The buildings: in the box widened by the city-events product's margin (where CLUE is joined), outlines only.
    margin = events["box"]["margin_m"]
    Bm = {"x0": B["x0"] - margin, "y0": B["y0"] - margin, "x1": B["x1"] + margin, "y1": B["y1"] + margin}
    parcels, context, skipped_parts = [], [], 0
    for b in model["buildings"]:
        f = b.get("footprint")
        if not f or len(f) < 3:
            continue
        c = centre(f)
        if not inside(c, Bm):
            continue
        if b["kind"].startswith("part:"):
            skipped_parts += 1
            continue
        clue = clue_of(b["id"])
        if clue:
            floors, height, source = int(clue["floors_above"]), STOREY_M * int(clue["floors_above"]), "CLUE floors x 3.2 m"
        else:
            height = float(b.get("height_m") or 9.6)
            floors, source = max(1, round(height / STOREY_M)), f"city model ({b.get('height_source')})"
        entry = {"id": b["id"], "name": b.get("name", ""), "footprint": f, "holes": b.get("holes") or [], "area_m2": round(signed_area(f), 2),
                 "floors": max(1, floors), "height_m": round(height, 2), "height_source": source,
                 "clue_property": clue["base"] if clue else None, "clue_use": clue["predominant_use"] if clue else None,
                 "clue_year": clue["year"] if clue else None, "use": use_class(clue["predominant_use"]) if clue else 1,
                 "in_district_box": inside(c, B)}
        (parcels if contains(quad, c) else context).append(entry)
    print(f"{len(parcels)} parcels, {len(context)} context buildings, {skipped_parts} building:part entries left out")
    print("parcels joined to CLUE:", sum(1 for p in parcels if p["clue_property"]), "uses", Counter(USES[p["use"]] for p in parcels))

    # Open spaces of the city model: those inside the redesign area are kept as public open space; the rest drawn.
    spaces_in, spaces_near = [], []
    for s in model["spaces"]:
        if s["kind"] not in ("park", "plaza", "green"):
            continue
        c = centre(s["ring"])
        if contains(quad, c):
            spaces_in.append({"id": s["id"], "kind": s["kind"], "ring": s["ring"], "area_m2": round(signed_area(s["ring"]), 1)})
        elif inside(c, Bm):
            spaces_near.append({"id": s["id"], "kind": s["kind"], "ring": s["ring"]})
    print("open spaces inside", [(s["id"], s["kind"], s["area_m2"]) for s in spaces_in])

    # The Sun at the solstice, at the area's middle.
    qc = centre(quad)
    # The area's middle in latitude and longitude: the frame's inverse by Newton's method, as city_events.py takes it.
    k = math.pi / 180 * 6378137.0
    lat0, lon0 = model["origin"]["lat"], model["origin"]["lon"]
    lat, lon = lat0 + qc[1] / k, lon0 + qc[0] / (k * math.cos(math.radians(lat0)))
    for _ in range(4):
        e_, n_, _u = enu(lat, lon, lat0, lon0)
        lat -= (float(n_) - qc[1]) / k
        lon -= (float(e_) - qc[0]) / (k * math.cos(math.radians(lat)))
    sun = {"window": sun_at(lat, lon, SUN_WINDOW), "hours": sun_at(lat, lon, HEAT_HOURS), "draw": sun_at(lat, lon, [DRAW_HOUR])[0]}
    air, ground, air_meta = june_air(weather)
    s_max, n_dev = floors_bound(CACHE / "melbourne-development-activity" / "hoddle-grid" / "developments.csv")
    print(f"sun at noon {sun['draw']}; June air {min(air):.1f}..{max(air):.1f} C, ground {ground:.2f} C; storeys bound {s_max} (of {n_dev} completions)")

    out = {
        "study": "designing-a-district",
        "frame": model["frame"], "origin": model["origin"], "projection": model["projection"],
        "district_box": B, "context_box": Bm, "redesign_quad": quad, "redesign_centre_latlon": [lat, lon],
        "corners": [f"{a} / {b}" for a, b in CORNERS],
        "network": net, "parcels": parcels, "context": context, "open_spaces": spaces_in, "spaces_near": spaces_near,
        "uses": USES, "storey_m": STOREY_M, "storeys_bound": s_max, "storeys_bound_from": f"90th percentile of floors_above over the {n_dev} completed developments of the Development Activity Monitor in the Hoddle Grid's box",
        "sun": sun, "sun_date": SOLSTICE.isoformat(), "time_zone": "AEST (UTC+10)",
        "air_by_hour_c": [round(a, 3) for a in air], "ground_c": round(ground, 3), "air_from": air_meta,
        "sources": {
            "city_model": {"file": "data/cache/products/city-model-melbourne.json", "built": model["built"], "credit": "© OpenStreetMap contributors, ODbL 1.0"},
            "clue": {"file": "data/cache/products/city-events-melbourne.json", "built": events["built"], "credit": "CLUE: City of Melbourne, CC BY 4.0"},
            "weather": {"file": "data/cache/products/street-weather-melbourne.json", "built": weather["built"], "credit": weather["source"]["credit"]},
            "developments": {"file": "data/cache/melbourne-development-activity/hoddle-grid/developments.csv",
                             "retrieved": json.loads((CACHE / "melbourne-development-activity" / "hoddle-grid" / "retrieval.json").read_text())["retrieved_at"],
                             "credit": "Development Activity Monitor: City of Melbourne, CC BY 4.0"},
            "sun": "domains/astronomy/js/sun.js (Astronomical Almanac low-precision formulas), via sun.mjs",
        },
        "notes": {"parts_left_out": skipped_parts},
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print("wrote", args.out, f"{Path(args.out).stat().st_size / 1e6:.2f} MB")
