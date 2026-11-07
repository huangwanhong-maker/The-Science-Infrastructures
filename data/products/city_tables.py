"""A city's tables for the toolbox of quantitative economics (ADR 0080): the records the data library already holds of Melbourne and
Freetown, as tables of columns whose rows are places, with what is needed to draw an answer back where it belongs. Kept under
data/cache/products/tables/, in the form of economic_tables.py, with

    "places": {"kind": "cells", "column", "time"?, "rings": {id: [[lon, lat], ...]}, "names": {id: name}}
    "places": {"kind": "points", "column", "time"?, "at": {id: [lon, lat]}, "names": {id: name}}

The tables:

    melbourne-blocks           each CLUE block in each census year, 2002–2024: floor space by use (m²), its businesses by kind, and its
                               dwellings. A cell the census withheld (blank, for confidentiality) is not a zero: a use's floor space
                               is the sum of the cells stated, and null if none is. Businesses are counted from the establishments
                               of each census year, by ANZSIC's two-digit groups; dwellings from the census's dwellings.
    melbourne-counters-month   each pedestrian counter in each month since May 2009: its walkers a day over the complete days of
                               that month. The City publishes nothing from November 2022 to October 2024.
    melbourne-counters-day     each counter on each day since 30 October 2024, with the day's rain and temperature at the centre, and
                               whether it is a weekend or a public holiday.
    freetown-wards             Freetown's wards: the assessment's annual average loss, the ward's area and the loss per km².

A day is a civil day in Melbourne (the offsets the counts carry). It counts for a counter only if it has at least 23 hours counted
(a civil day has 23 at the start of daylight saving) and a total above zero: a counter that counted nobody all day in a city
centre is a counter out of order, not an empty street. A month needs fifteen such days.

    python data/products/city_tables.py
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from economic_tables import KEPT, table  # noqa: E402
from product_index import OUT as PRODUCTS, register  # noqa: E402

CACHE = HERE.parents[0] / "cache"
EPOCH = dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc)

# Floor space by use: the census's columns, gathered.
USES = {
    "office": ["office"],
    "retail": ["retail_shop", "retail_showroom", "retail_stall", "retail_cars"],
    "homes": ["residential_apartment", "house_townhouse", "student_accommodation", "institutional_accommodation"],
    "hotels": ["commercial_accommodation"],
    "education": ["educational_research"],
    "health": ["hospital_clinic"],
    "culture": ["entertainment_recreation_indoor", "performances_conferences_ceremonies", "public_display_area"],
    "community": ["community_use"],
    "workshops": ["manufacturing", "workshop_studio", "wholesale", "storage"],
    "parking": ["parking_commercial_covered", "parking_commercial_uncovered", "parking_private_covered", "parking_private_uncovered"],
    "vacant": ["unoccupied_under_construction", "unoccupied_under_demolition_condemned", "unoccupied_under_renovation", "unoccupied_undeveloped_site", "unoccupied_unused"],
}
# Businesses by ANZSIC 2006's two-digit groups.
BUSINESSES = {
    "shops": range(39, 44),
    "food_and_drink": range(45, 46),
    "finance": range(62, 65),
    "professional": range(69, 71),
    "education_providers": range(80, 83),
    "health_care": range(84, 88),
    "arts_recreation": range(89, 93),
}


def rows_of(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f, delimiter=";")


def number(text) -> float | None:
    if text is None or str(text).strip() == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def outer_ring(geometry: dict) -> list:
    """A polygon's outer ring, or a multipolygon's largest, closed end dropped."""
    polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]

    def area(ring):
        return abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(ring, ring[1:] + ring[:1]))) / 2
    ring = max((p[0] for p in polygons), key=area)
    if ring[0] == ring[-1]:
        ring = ring[:-1]
    return [[round(p[0], 6), round(p[1], 6)] for p in ring]


def spherical_area_km2(ring: list) -> float:
    """A small ring's area by the shoelace formula on a local equirectangular plane, in km²."""
    lat0 = math.radians(sum(p[1] for p in ring) / len(ring))
    k = 111.32
    xy = [(p[0] * k * math.cos(lat0), p[1] * k) for p in ring]
    return abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(xy, xy[1:] + xy[:1]))) / 2


def melbourne_blocks() -> dict:
    blocks = json.loads((CACHE / "melbourne-clue-blocks/all/blocks.geojson").read_text(encoding="utf-8"))
    rings, names = {}, {}
    for f in blocks["features"]:
        p = f["properties"]
        ident = str(p["block_id"])
        rings[ident] = outer_ring(f["geometry"])
        names[ident] = f"Block {ident}, {p.get('clue_small_area') or p.get('area_name') or p.get('region_name') or ''}".rstrip(", ")
    businesses: dict[tuple[str, int], dict] = defaultdict(lambda: defaultdict(float))
    for r in rows_of(CACHE / "melbourne-clue-establishments/all/establishments.csv"):
        key = (str(r["block_id"]), int(r["census_year"]))
        businesses[key]["establishments"] += 1
        code = number(r.get("industry_anzsic4_code"))
        if code is None:
            continue
        group = int(code) // 100
        for name, span in BUSINESSES.items():
            if group in span:
                businesses[key][name] += 1
    dwellings: dict[tuple[str, int], float] = defaultdict(float)
    for r in rows_of(CACHE / "melbourne-clue-dwellings/municipality/dwellings.csv"):
        dwellings[(str(r["block_id"]), int(r["census_year"]))] += number(r.get("dwelling_number")) or 0.0
    rows = []
    for r in rows_of(CACHE / "melbourne-clue-floor-space/all/floor-space.csv"):
        ident, year = str(r["block_id"]), int(r["census_year"])
        if ident not in rings:
            continue
        row = {"block_id": ident, "census_year": year, "small_area": r["clue_small_area"], "floor_total": number(r["total_floor_space_in_block"])}
        for use, cols in USES.items():
            stated = [number(r.get(c)) for c in cols]
            stated = [v for v in stated if v is not None]
            row[use] = sum(stated) if stated else None
        counted = businesses.get((ident, year), {})
        row["establishments"] = counted.get("establishments", 0.0)
        for name in BUSINESSES:
            row[name] = counted.get(name, 0.0)
        row["dwellings"] = dwellings.get((ident, year), 0.0)
        rows.append(row)
    rows.sort(key=lambda x: (x["census_year"], int(x["block_id"]) if x["block_id"].isdigit() else 0))
    columns = [("block_id", "text", "CLUE's block"), ("census_year", "number", ""), ("small_area", "text", "CLUE's small area"),
               ("floor_total", "number", "the block's floor space, m²")]
    columns += [(u, "number", f"{u.replace('_', ' ')}: floor space, m² (the cells the census states)") for u in USES]
    columns += [("establishments", "number", "businesses and other establishments, counted")]
    columns += [(b, "number", f"{b.replace('_', ' ')}: establishments, counted") for b in BUSINESSES]
    columns += [("dwellings", "number", "dwellings, counted")]
    t = table("melbourne-blocks", "Melbourne's blocks by census year: floor space, businesses and dwellings",
              "City of Melbourne, Census of Land Use and Employment (CC BY 4.0)", columns, rows)
    t["places"] = {"kind": "cells", "column": "block_id", "time": "census_year", "time_kind": "year", "rings": rings, "names": names}
    return t


def days_of(hours: np.ndarray, t0: float, offsets: list) -> np.ndarray:
    """Each hour's civil day, as days since 1970, from the offsets that hold from each listed instant on."""
    t = t0 + hours * 3600.0
    instants = np.array([instant(a) for a, _ in offsets])
    hours_ahead = np.array([b for _, b in offsets], dtype=float)
    where = np.clip(np.searchsorted(instants, t, side="right") - 1, 0, len(offsets) - 1)
    return np.floor((t + hours_ahead[where] * 3600.0) / 86400.0).astype(np.int64)


def instant(text: str) -> float:
    return dt.datetime.fromisoformat(text.replace("Z", "+00:00")).timestamp()


def iso_of(day: int) -> str:
    return (EPOCH + dt.timedelta(days=int(day))).date().isoformat()


def counters() -> tuple[dict, dict]:
    sensors = json.loads((PRODUCTS / "street-sensors.json").read_text(encoding="utf-8"))["sensors"]
    weather = json.loads((PRODUCTS / "street-weather-melbourne.json").read_text(encoding="utf-8"))
    holidays = {h["date"] for h in json.loads((PRODUCTS / "street-holidays-victoria.json").read_text(encoding="utf-8"))["holidays"]}
    weather_t0 = instant(weather["t0_utc"])
    temperature = np.array([np.nan if v is None else v for v in weather["temperature_c"]], dtype=float)
    rain = np.array([np.nan if v is None else v for v in weather["precipitation_mm"]], dtype=float)
    at, names, rows_month, rows_day = {}, {}, [], []
    live_from = int(dt.datetime(2024, 10, 30, tzinfo=dt.timezone.utc).timestamp() // 86400)
    day_weather: dict[int, tuple[float, float]] = {}
    for s in sensors:
        if s.get("lat") is None or s.get("lon") is None or not s.get("series") or not (PRODUCTS / s["series"]).exists():
            continue
        c = json.loads((PRODUCTS / s["series"]).read_text(encoding="utf-8"))
        counts = np.array([np.nan if v is None else v for v in c["counts"]], dtype=float)
        if not counts.size:
            continue
        ident = str(s["id"])
        days = days_of(np.arange(counts.size, dtype=float), instant(c["t0_utc"]), c["offsets"])
        if not day_weather:
            # The weather on the same civil days, once.
            w_days = days_of(np.arange(temperature.size, dtype=float), weather_t0, c["offsets"])
            order = np.unique(w_days)
            base = order[0]
            size = order[-1] - base + 1
            rain_sum, rain_n, temp_sum, temp_n = np.zeros(size), np.zeros(size), np.zeros(size), np.zeros(size)
            ok = np.isfinite(rain)
            np.add.at(rain_sum, w_days[ok] - base, rain[ok])
            np.add.at(rain_n, w_days[ok] - base, 1)
            ok = np.isfinite(temperature)
            np.add.at(temp_sum, w_days[ok] - base, temperature[ok])
            np.add.at(temp_n, w_days[ok] - base, 1)
            for k in range(size):
                if rain_n[k] >= 23 and temp_n[k] >= 23:
                    day_weather[int(base + k)] = (float(rain_sum[k]), float(temp_sum[k] / temp_n[k]))
        first = days.min()
        size = days.max() - first + 1
        total, hours = np.zeros(size), np.zeros(size)
        ok = np.isfinite(counts)
        np.add.at(total, days[ok] - first, counts[ok])
        np.add.at(hours, days[ok] - first, 1)
        complete = (hours >= 23) & (total > 0)
        if not complete.any():
            continue
        at[ident] = [round(s["lon"], 6), round(s["lat"], 6)]
        names[ident] = s["name"]
        month_sum: dict[tuple[int, int], list] = defaultdict(lambda: [0.0, 0])
        for k in np.nonzero(complete)[0]:
            day = int(first + k)
            d = (EPOCH + dt.timedelta(days=day)).date()
            m = month_sum[(d.year, d.month)]
            m[0] += float(total[k])
            m[1] += 1
            if day >= live_from and day in day_weather:
                r, tm = day_weather[day]
                iso = d.isoformat()
                rows_day.append({"counter": ident, "name": s["name"], "date": iso, "day": day, "walkers": float(total[k]), "rain_mm": round(r, 2), "temperature_c": round(tm, 2),
                                 "weekend": 1.0 if d.weekday() >= 5 else 0.0, "holiday": 1.0 if iso in holidays else 0.0})
        for (y, mo), (sum_, n) in month_sum.items():
            if n < 15:
                continue
            rows_month.append({"counter": ident, "name": s["name"], "month": f"{y:04d}-{mo:02d}", "year": y, "time": y + (mo - 0.5) / 12, "walkers_a_day": sum_ / n, "days": n,
                               "compare": 1.0 if y in (2019, 2025) else 0.0, "after": 1.0 if y == 2025 else 0.0})
    rows_month.sort(key=lambda r: (r["time"], int(r["counter"])))
    rows_day.sort(key=lambda r: (r["day"], int(r["counter"])))
    credit = "City of Melbourne, Pedestrian Counting System (the City's open data, CC BY 4.0)"
    month = table("melbourne-counters-month", "Melbourne's counters by month since 2009: walkers a day", credit,
                  [("counter", "text", "the counter's id"), ("name", "text", ""), ("month", "text", ""), ("year", "number", ""), ("time", "number", "decimal year"),
                   ("walkers_a_day", "number", "walkers a day, over the month's complete days"), ("days", "number", "complete days counted"),
                   ("compare", "number", "1 in 2019 and in 2025, to compare the years before and after the lockdowns"), ("after", "number", "1 in 2025")], rows_month)
    month["places"] = {"kind": "points", "column": "counter", "time": "month", "time_kind": "month", "at": at, "names": names}
    day = table("melbourne-counters-day", "Melbourne's counters by day since October 2024: walkers and the weather", f"{credit}; Open-Meteo, from ERA5 and the IFS (CC BY 4.0); Nager.Date",
                [("counter", "text", "the counter's id"), ("name", "text", ""), ("date", "text", "a civil day in Melbourne"), ("day", "number", "days since 1970"),
                 ("walkers", "number", "walkers that day"), ("rain_mm", "number", "rain that day at the centre, mm"), ("temperature_c", "number", "mean temperature that day at the centre, °C"),
                 ("weekend", "number", "1 on Saturday and Sunday"), ("holiday", "number", "1 on a public holiday in Victoria")], rows_day)
    day["places"] = {"kind": "points", "column": "counter", "time": "date", "time_kind": "date", "at": at, "names": names}
    return month, day


def freetown_wards() -> dict:
    p = json.loads((PRODUCTS / "freetown-flood.json").read_text(encoding="utf-8"))
    rows, rings, names = [], {}, {}
    for w in p["wards"]:
        polygons = w.get("polygons") or []
        if not polygons:
            continue
        ring = max((q["outer"] for q in polygons), key=spherical_area_km2)
        ring = ring[:-1] if ring[0] == ring[-1] else ring
        rings[w["code"]] = [[round(a, 6), round(b, 6)] for a, b in ring]
        names[w["code"]] = w["name"]
        km2 = sum(spherical_area_km2(q["outer"]) for q in polygons)
        rows.append({"ward": w["code"], "name": w["name"], "annual_loss": w.get("aal"), "area_km2": round(km2, 4), "loss_per_km2": (w["aal"] / km2) if w.get("aal") is not None and km2 > 0 else None})
    t = table("freetown-wards", "Freetown's wards: annual average loss from floods", p.get("credit", ""),
              [("ward", "text", "the ward's code"), ("name", "text", ""), ("annual_loss", "number", "the assessment's annual average loss, US$"), ("area_km2", "number", "the ward's area, km²"),
               ("loss_per_km2", "number", "annual average loss per km², US$")], rows)
    t["places"] = {"kind": "cells", "column": "ward", "rings": rings, "names": names}
    return t


def main() -> int:
    KEPT.mkdir(parents=True, exist_ok=True)
    made = []
    if (CACHE / "melbourne-clue-floor-space/all/floor-space.csv").exists():
        made.append(melbourne_blocks())
    if (PRODUCTS / "street-sensors.json").exists():
        made.extend(counters())
    if (PRODUCTS / "freetown-flood.json").exists():
        made.append(freetown_wards())
    for t in made:
        path = KEPT / f"{t['id']}.json"
        path.write_text(json.dumps(t, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
        print(t["id"], t["rows"], "rows", len(t["columns"]), "columns", path.stat().st_size // 1024, "kB", flush=True)
        register("tables", {"name": t["id"], "title": t["title"], "file": f"tables/{t['id']}.json", "rows": t["rows"], "source": t["source"]})
    return 0 if made else 1


if __name__ == "__main__":
    sys.exit(main())
