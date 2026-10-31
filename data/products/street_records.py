"""Build the street records: Melbourne's hourly pedestrian counts with their weather and holidays, the
uses of the premises about each sensor, and the holdings of chosen streets with their declared opening
hours, for the Urban Building Workbench's street oceanography.

1. **Sensors and counts.** The City of Melbourne's Pedestrian Counting System, its archive (May 2009 to
   October 2022) and its live record (from 30 October 2024) joined by the sensor's id, each sensor a
   series over every UTC hour from its first record to its last, null where nothing was recorded (the
   City publishes nothing from 1 November 2022 to 29 October 2024). The records give Melbourne's civil
   date and hour; an hour becomes the UTC hour that begins at that civil time (§2).
2. **Civil time.** Melbourne keeps AEST (UTC+10) and, from the first Sunday of October at 02:00 AEST to
   the first Sunday of April at 03:00 AEDT, AEDT (UTC+11): the rule of Victoria's Summer Time Act as
   amended in 2007, in force since April 2008 (both changes fall at 16:00 UTC on the Saturday). The rule
   is written here because the library's Python has no time-zone database; the checks compare it with
   the IANA database where one is installed. The hour the clock skips (02:00 on the first Sunday of
   October) refers to no instant: a record for it is dropped (the archive has such records, real counts
   in 2009-2012 and zeros in 2018-2022). The hour the clock repeats (02:00 on the first Sunday of April)
   appears once in the records, with the counts of about one hour (its ratio to its neighbours is that
   of other Sundays in the archive, and 1.10 and 1.26 against 0.92 in the live record, where both hours
   would give about 1.8): it is placed at the first of the two (02:00 AEDT, 15:00 UTC on the Saturday),
   and the second is null.
3. **Repairs.** In September 2010 the archive stamps every hour of every sensor 00:00, 24 rows a day;
   their IDs keep the order of the hours (the profiles so ordered correlate with the weekday profiles of
   August and October 2010 at a median of 0.98), and the hours are restored in that order. Any other
   hour recorded twice is left null.
4. **Weather and holidays.** Open-Meteo's hourly temperature, precipitation and wind (from ERA5 and
   the IFS analyses) at Melbourne's CBD, in UTC; the precipitation is moved to the hour it falls in
   (Open-Meteo stamps a sum at the end of its hour). Victoria's public holidays from Nager.Date, those of
   the whole country and those it lists for AU-VIC, with the discrepancies against Victoria's own
   calendar that are known noted, not corrected.
5. **Terrain.** The CLUE census's business establishments within 200 m of each sensor, census year by
   census year, each in a local tangent plane about the sensor and with a use from a small
   classification of street uses (USES) by its ANZSIC class.
6. **Holdings.** The premises along chosen streets from OpenStreetMap (Overpass): each street's axis is
   the centreline of its carriageways between two cross streets, and each holding within 40 m of it is
   projected onto it (arc length s from the axis's first point, offset positive on the left), with its
   use from its tags and its opening_hours read into weekly intervals (data/formats/opening_hours.py),
   and its latitude and longitude; a street lying in an area of the city model names it.
7. **The city model.** The buildings (with heights), roads (with widths) and open spaces of Melbourne's
   centre, from OpenStreetMap, in exact east-north coordinates about the area's centre (§8), for the
   workbench's three-dimensional view.
8. **Planes.** The city model is projected exactly: geodetic to ECEF to east-north-up at the area's
   centre on WGS 84 (`enu`), the up component (the curvature's drop, under 1 m) dropped. The terrain and
   the streets use the equirectangular plane x = N(φ0) cos φ0 (λ − λ0), y = M(φ0)(φ − φ0) about the
   sensor or the street's centre, whose error grows as |x y| tan φ0 / M (the parallels' radius changing
   across the plane): under 1 cm within 200 m of a sensor, and stated for each street (`plane_error_m`)
   and each area of the city model (`equirectangular_error_m`, against the exact projection).

Output, in data/cache/products/: street-sensors.json, street-counts/<id>.json, street-weather-melbourne.json,
street-holidays-victoria.json, street-terrain-melbourne.json, street-holdings.json, city-model.json; listed in the index
under the kind `streets`.

Usage::

    python data/acquisition/fetch.py melbourne-pedestrian-archive melbourne-pedestrian-counts \\
        melbourne-clue-establishments open-meteo-hourly nager-holidays osm-streets osm-city
    python data/products/street_records.py                      # every product whose inputs are fetched
    python data/products/street_records.py --only counts        # sensors, counts, weather, holidays, terrain, holdings or city
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import heapq
import io
import json
import math
import re
import sys
import zipfile
from array import array
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
CACHE = DATA / "cache"
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(DATA / "formats"))
from product_index import OUT, register  # noqa: E402
import opening_hours  # noqa: E402

UTC = dt.timezone.utc
ARCHIVE = CACHE / "melbourne-pedestrian-archive" / "pedestrian-counts-2009-2022.csv.zip"
LIVE = CACHE / "melbourne-pedestrian-counts" / "current"
CLUE = CACHE / "melbourne-clue-establishments" / "all"
WEATHER = CACHE / "open-meteo-hourly" / "melbourne"
HOLIDAYS = CACHE / "nager-holidays" / "australia"
OSM = CACHE / "osm-streets" / "streets"
TERRAIN_RADIUS = 200.0          # m about a sensor
HOLDING_RADIUS = 40.0           # m from a street's axis
MONTHS = {m: i for i, m in enumerate(("January", "February", "March", "April", "May", "June", "July", "August",
                                      "September", "October", "November", "December"), 1)}


def iso(hour: int) -> str:
    """A UTC hour (hours since 1970-01-01T00Z) as ISO 8601."""
    return dt.datetime.fromtimestamp(hour * 3600, UTC).strftime("%Y-%m-%dT%H:00:00Z")


def hour_of(t: dt.datetime) -> int:
    return int(t.timestamp()) // 3600


def write(path: Path, product: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")


def retrieved(folder: Path) -> str | None:
    p = folder / "retrieval.json"
    return json.loads(p.read_text(encoding="utf-8")).get("retrieved_at") if p.exists() else None


# ---- Melbourne's civil time ----------------------------------------------------------------------------

def first_sunday(year: int, month: int) -> dt.date:
    d = dt.date(year, month, 1)
    return d + dt.timedelta(days=(6 - d.weekday()) % 7)


def transitions(year: int) -> tuple[dt.datetime, dt.datetime]:
    """The instants (UTC) at which daylight saving ends (April) and begins (October) in a year, by the rule
    in force since April 2008: the first Sunday of April at 03:00 AEDT and the first Sunday of October at
    02:00 AEST, both 16:00 UTC on the Saturday before."""
    if year < 2008:
        raise ValueError(f"Melbourne's rule of 2008 does not hold in {year}")
    end = dt.datetime.combine(first_sunday(year, 4), dt.time(16), UTC) - dt.timedelta(days=1)
    start = dt.datetime.combine(first_sunday(year, 10), dt.time(16), UTC) - dt.timedelta(days=1)
    return end, start


def offset_hours(t: dt.datetime) -> int:
    """Melbourne's offset from UTC (10 or 11 hours) at an instant t (UTC, 2008 onwards)."""
    end, start = transitions(t.year)
    return 11 if t < end or t >= start else 10


def local_to_utc(date: dt.date, hour: int) -> tuple[int | None, str]:
    """The UTC hour that begins at a civil hour of Melbourne, and its kind: 'ok'; 'repeated' (the hour the
    clock shows twice: the first, in AEDT, is returned); 'skipped' (the hour the clock skips: None)."""
    local = dt.datetime.combine(date, dt.time(hour), UTC)
    candidates = [local - dt.timedelta(hours=k) for k in (11, 10) if offset_hours(local - dt.timedelta(hours=k)) == k]
    if not candidates:
        return None, "skipped"
    return hour_of(candidates[0]), ("repeated" if len(candidates) == 2 else "ok")


def offset_changes(first: int, last: int) -> list[list]:
    """The offsets of Melbourne's civil time over the UTC hours first..last: [utc_iso, hours] at first and at
    each change."""
    t0 = dt.datetime.fromtimestamp(first * 3600, UTC)
    out = [[iso(first), offset_hours(t0)]]
    for year in range(t0.year, dt.datetime.fromtimestamp(last * 3600, UTC).year + 1):
        for t in sorted(transitions(year)):
            h = hour_of(t)
            if first < h <= last:
                out.append([iso(h), offset_hours(t)])
    return out


# ---- sensors and their counts ----------------------------------------------------------------------------

class Series:
    """A sensor's records: (UTC hour, count) pairs, and what was done to them."""

    def __init__(self):
        self.hours, self.counts = array("q"), array("q")
        self.sources: set[str] = set()
        self.names: dict[str, str] = {}
        # The live record counts each direction (direction_1, direction_2 of the locations dataset): (UTC hour, d1, d2).
        self.dir_hours, self.d1, self.d2 = array("q"), array("q"), array("q")

    def add(self, hour: int, count: int) -> None:
        self.hours.append(hour)
        self.counts.append(count)

    def add_directions(self, hour: int, d1: int, d2: int) -> None:
        self.dir_hours.append(hour)
        self.d1.append(d1)
        self.d2.append(d2)


def _cached_local(cache: dict, ymd: tuple, hour: int):
    key = (ymd, hour)
    if key not in cache:
        cache[key] = local_to_utc(dt.date(*ymd), hour)
    return cache[key]


def read_archive(path: Path, series: dict[int, Series], tally: dict) -> None:
    """The archive's rows into series; September 2010's hours restored from the rows' order (module §3)."""
    with zipfile.ZipFile(path) as z:
        member = next(n for n in z.namelist() if n.lower().endswith(".csv") and not n.startswith("__MACOSX"))
        with z.open(member) as f:
            reader = csv.reader(io.TextIOWrapper(f, encoding="utf-8-sig", newline=""))
            head = next(reader)
            col = {name: i for i, name in enumerate(head)}
            need = ("ID", "Year", "Month", "Mdate", "Time", "Sensor_ID", "Sensor_Name", "Hourly_Counts")
            missing = [k for k in need if k not in col]
            if missing:
                raise SystemExit(f"the archive lacks the columns {missing}; it has {head}")
            # per sensor: day ordinal * 24 + hour, row ID, count
            rows: dict[int, tuple[array, array, array]] = {}
            names: dict[int, dict[str, int]] = {}
            for r in reader:
                s = int(r[col["Sensor_ID"]])
                day = dt.date(int(r[col["Year"]]), MONTHS[r[col["Month"]]], int(r[col["Mdate"]])).toordinal()
                key, ident, count = day * 24 + int(r[col["Time"]]), int(r[col["ID"]]), int(r[col["Hourly_Counts"]])
                a = rows.setdefault(s, (array("q"), array("q"), array("q")))
                a[0].append(key), a[1].append(ident), a[2].append(count)
                nm = names.setdefault(s, {})
                nm[r[col["Sensor_Name"]]] = max(nm.get(r[col["Sensor_Name"]], 0), day)
    cache: dict = {}
    for s, (keys, ids, counts) in rows.items():
        keys, ids, counts = np.frombuffer(keys, np.int64).copy(), np.frombuffer(ids, np.int64), np.frombuffer(counts, np.int64)
        order = np.lexsort((ids, keys))
        keys, ids, counts = keys[order], ids[order], counts[order]
        # A day whose 24 rows all carry hour 0: the hours in the rows' order.
        days = keys // 24
        for day in np.unique(days[keys % 24 == 0]):
            at = np.nonzero(days == day)[0]
            if len(at) == 24 and np.all(keys[at] % 24 == 0):
                keys[at] = day * 24 + np.arange(24)
                tally["restored_days"] += 1
        order = np.lexsort((ids, keys))
        keys, counts = keys[order], counts[order]
        same = np.zeros(len(keys), bool)
        same[1:] |= keys[1:] == keys[:-1]
        same[:-1] |= keys[:-1] == keys[1:]
        tally["duplicated_hours"] += int(np.count_nonzero(same))
        out = series.setdefault(s, Series())
        out.sources.add("archive")
        for key, count in zip(keys[~same].tolist(), counts[~same].tolist()):
            day, hour = divmod(key, 24)
            d = dt.date.fromordinal(day)
            utc, kind = _cached_local(cache, (d.year, d.month, d.day), hour)
            if kind == "skipped":
                tally["skipped_hour_records"] += 1
                continue
            tally["repeated_hour_records"] += kind == "repeated"
            out.add(utc, count)
        for name, last in names[s].items():
            out.names[name] = dt.date.fromordinal(last).isoformat()


def read_live(folder: Path, series: dict[int, Series], tally: dict) -> dict[int, dict]:
    """The live record's rows into series; returns each sensor's code and location from its rows."""
    seen: dict[int, dict] = {}
    cache: dict = {}
    with (folder / "counts.csv").open(encoding="utf-8-sig", newline="") as f:
        reader = csv.reader(f, delimiter=";")
        head = next(reader)
        col = {name: i for i, name in enumerate(head)}
        has_location = "location" in col
        keys: dict[int, set] = {}
        for r in reader:
            s = int(r[col["location_id"]])
            y, m, d = (int(x) for x in r[col["sensing_date"]].split("-"))
            utc, kind = _cached_local(cache, (y, m, d), int(r[col["hourday"]]))
            if kind == "skipped":
                tally["skipped_hour_records"] += 1
                continue
            tally["repeated_hour_records"] += kind == "repeated"
            out = series.setdefault(s, Series())
            out.sources.add("live")
            if utc in keys.setdefault(s, set()):
                tally["duplicated_hours"] += 1
                continue
            keys[s].add(utc)
            out.add(utc, int(r[col["pedestriancount"]]))
            d1, d2 = r[col["direction_1"]], r[col["direction_2"]]
            if d1 != "" and d2 != "":
                out.add_directions(utc, int(d1), int(d2))
            if s not in seen:
                lat = lon = None
                if has_location and r[col["location"]]:
                    lat, lon = (float(x) for x in r[col["location"]].split(","))
                seen[s] = {"code": r[col["sensor_name"]] or None, "lat": lat, "lon": lon}
    return seen


def read_locations(folder: Path) -> dict[int, dict]:
    out = {}
    with (folder / "sensors.csv").open(encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f, delimiter=";"):
            out[int(r["location_id"])] = r
    return out


STATUS = {"A": "active", "I": "inactive", "R": "removed"}


def sensor_entries(series: dict[int, Series], live: dict[int, dict], locations: dict[int, dict]) -> list[dict]:
    """The sensors' descriptions (without their series' statistics)."""
    out = []
    for s in sorted(set(series) | set(locations)):
        loc = locations.get(s)
        ser = series.get(s)
        archive_names = sorted(ser.names, key=lambda n: ser.names[n]) if ser else []
        if loc:
            e = {"id": s, "code": loc["sensor_name"] or None, "name": loc["sensor_description"],
                 "lat": round(float(loc["latitude"]), 7), "lon": round(float(loc["longitude"]), 7),
                 "installed": loc["installation_date"] or None, "status": STATUS.get(loc["status"], loc["status"]),
                 "location_type": loc["location_type"] or None,
                 "directions": [d for d in (loc["direction_1"].strip(), loc["direction_2"].strip()) if d]}
            note = loc.get("note", "").strip()
            if note and note != "<Null>":
                e["note"] = note
        else:
            lv = live.get(s, {})
            e = {"id": s, "code": lv.get("code"), "name": archive_names[-1] if archive_names else lv.get("code"),
                 "lat": round(lv["lat"], 7) if lv.get("lat") is not None else None,
                 "lon": round(lv["lon"], 7) if lv.get("lon") is not None else None,
                 "installed": None, "status": "unlisted", "location_type": None, "directions": []}
        if archive_names:
            e["archive_names"] = archive_names
        out.append(e)
    return out


def counts_product(s: int, ser: Series) -> tuple[dict, dict]:
    """A sensor's series over every UTC hour from its first record to its last, and its statistics."""
    hours = np.frombuffer(ser.hours, np.int64)
    counts = np.frombuffer(ser.counts, np.int64)
    first, last = int(hours.min()), int(hours.max())
    values = np.full(last - first + 1, -1, np.int64)
    values[hours - first] = counts
    present = values >= 0
    # Runs of a week or more without a record.
    gaps, run = [], None
    for i, ok in enumerate(present.tolist() + [True]):
        if not ok and run is None:
            run = i
        elif ok and run is not None:
            if i - run >= 168:
                gaps.append([iso(first + run), i - run])
            run = None
    product = {"id": s, "t0_utc": iso(first), "step_s": 3600,
               "counts": [int(v) if v >= 0 else None for v in values.tolist()],
               "offsets": offset_changes(first, last)}
    if len(ser.dir_hours):
        # Each direction's counts over every UTC hour of the live record's span (d1 + d2 is the count).
        dh = np.frombuffer(ser.dir_hours, np.int64)
        f0, f1 = int(dh.min()), int(dh.max())
        d = np.full((2, f1 - f0 + 1), -1, np.int64)
        d[0, dh - f0] = np.frombuffer(ser.d1, np.int64)
        d[1, dh - f0] = np.frombuffer(ser.d2, np.int64)
        product["directions"] = {"t0_utc": iso(f0), "step_s": 3600,
                                 "d1": [int(v) if v >= 0 else None for v in d[0].tolist()],
                                 "d2": [int(v) if v >= 0 else None for v in d[1].tolist()]}
    stats = {"first_utc": iso(first), "last_utc": iso(last), "hours": int(len(values)),
             "missing": int(np.count_nonzero(~present)), "gaps": gaps, "datasets": sorted(ser.sources)}
    return product, stats


def build_sensors_and_counts(write_counts: bool = True) -> dict:
    tally = {"skipped_hour_records": 0, "repeated_hour_records": 0, "duplicated_hours": 0, "restored_days": 0}
    series: dict[int, Series] = {}
    print("  reading the archive (4.6 million rows) ...")
    read_archive(ARCHIVE, series, tally)
    print("  reading the live record ...")
    live = read_live(LIVE, series, tally)
    locations = read_locations(LIVE)
    sensors = sensor_entries(series, live, locations)
    folder = OUT / "street-counts"
    for e in sensors:
        ser = series.get(e["id"])
        if ser is None or not len(ser.hours):
            e.update({"series": None, "first_utc": None, "last_utc": None, "hours": 0, "missing": 0, "gaps": [], "datasets": []})
            continue
        product, stats = counts_product(e["id"], ser)
        e["series"] = f"street-counts/{e['id']}.json"
        e.update(stats)
        if write_counts:
            write(folder / f"{e['id']}.json", product)
    product = {
        "product": "street-sensors",
        "built": dt.datetime.now(UTC).isoformat(timespec="seconds"),
        "timezone": "Australia/Melbourne",
        "sources": [
            {"dataset": "melbourne-pedestrian-archive", "credit": "City of Melbourne", "licence": "no licence stated; the City's open data are CC BY 4.0",
             "retrieved_at": retrieved(ARCHIVE.parent)},
            {"dataset": "melbourne-pedestrian-counts", "profile": "current", "credit": "City of Melbourne",
             "licence": "counts: no licence stated; locations: CC BY 4.0", "retrieved_at": retrieved(LIVE)},
        ],
        "time": ("Each series is over every UTC hour from the sensor's first record to its last; a count is that of the "
                 "hour beginning at t0_utc + i step_s. The records give Melbourne's civil hour; offsets gives the "
                 "offset of civil time from UTC (hours) from each listed instant on."),
        "directions": ("Where the live record counts each direction, a series' directions holds d1 and d2, the counts toward "
                       "the sensor's directions[0] and directions[1], over every UTC hour of the live record's span from its own "
                       "t0_utc; d1 + d2 is the count."),
        "gap": "The City publishes no counts from 2022-11-01 to 2024-10-29 (the archive ends on 2022-10-31): those hours are null.",
        "daylight_saving": ("The rule since 2008: AEDT (+11) from the first Sunday of October 02:00 AEST to the first Sunday of "
                            "April 03:00 AEDT, else AEST (+10). A record of the skipped hour (02:00 on the first Sunday of October) "
                            "is dropped; the single record of the repeated hour (02:00 on the first Sunday of April), whose counts "
                            "are those of about one hour, is placed at the first of the two (15:00 UTC on the Saturday) and the second is null."),
        "repairs": ("September 2010: the archive stamps all 24 hours of each day 00:00; the hours are restored in the order of "
                    "the rows' IDs (profiles so ordered correlate with the neighbouring months' at a median of 0.98)."),
        "tally": tally,
        "status": "active, inactive or removed as the locations dataset says; unlisted: absent from it (removed before 2024), "
                  "placed by the live record's rows where it has them, else without a place",
        "sensors": sensors,
    }
    write(OUT / "street-sensors.json", product)
    register("streets", {"name": "street-sensors", "file": "street-sensors.json", "title": "Melbourne's pedestrian sensors and their hourly series",
                         "domain": "urban-science", "sensors": len(sensors), "series": "street-counts/<id>.json",
                         "source": "City of Melbourne, Pedestrian Counting System"})
    print(f"  street-sensors: {len(sensors)} sensors, {sum(1 for e in sensors if e['series'])} series; {tally}")
    return product


# ---- weather -------------------------------------------------------------------------------------------

def build_weather() -> dict:
    raw = json.loads((WEATHER / "hourly.json").read_text(encoding="utf-8"))
    h = raw["hourly"]
    times = [dt.datetime.fromisoformat(t).replace(tzinfo=UTC) for t in h["time"]]
    first = hour_of(times[0])
    n = hour_of(times[-1]) - first + 1
    temp, rain, wind = [None] * n, [None] * n, [None] * n
    for t, a, b, c in zip(times, h["temperature_2m"], h["precipitation"], h["wind_speed_10m"]):
        i = hour_of(t) - first
        temp[i] = None if a is None else round(a, 1)
        wind[i] = None if c is None else round(c, 2)
        # Open-Meteo's precipitation is the sum over the hour ending at its stamp: it belongs to the hour before.
        if i >= 1:
            rain[i - 1] = None if b is None else round(b, 2)
    product = {
        "product": "street-weather-melbourne",
        "built": dt.datetime.now(UTC).isoformat(timespec="seconds"),
        "t0_utc": iso(first), "step_s": 3600,
        "temperature_c": temp, "precipitation_mm": rain, "wind_m_s": wind,
        "conventions": ("Index i is the UTC hour beginning at t0_utc + i step_s. temperature_c (2 m) and wind_m_s (10 m) are "
                        "at that instant; precipitation_mm is the total over the hour that begins there (Open-Meteo's value "
                        "stamped an hour later), so that it lies on the same hour as a pedestrian count."),
        "grid": {"latitude": raw["latitude"], "longitude": raw["longitude"], "elevation_m": raw.get("elevation"),
                 "requested": [-37.8136, 144.9631]},
        "source": {"dataset": "open-meteo-hourly", "profile": "melbourne", "credit": "Weather data by Open-Meteo.com (CC BY 4.0), from Copernicus ERA5",
                   "retrieved_at": retrieved(WEATHER)},
    }
    write(OUT / "street-weather-melbourne.json", product)
    register("streets", {"name": "street-weather-melbourne", "file": "street-weather-melbourne.json",
                         "title": "Hourly weather at Melbourne's CBD (temperature, precipitation, wind)", "domain": "meteorology",
                         "times": [iso(first), iso(first + n - 1)], "source": "Open-Meteo (ERA5)"})
    missing = sum(v is None for v in temp)
    print(f"  street-weather-melbourne: {n} hours from {iso(first)}, {missing} without temperature")
    return product


# ---- holidays ------------------------------------------------------------------------------------------

# Known differences between Nager.Date's Victorian holidays and Victoria's own calendar: flagged, not corrected.
HOLIDAY_NOTES = [
    "2009-04-12 Easter Sunday is listed by Nager.Date as national; Victoria has kept Easter Sunday as a public holiday "
    "only since 2016 (Nager.Date omits it in 2010-2015).",
    "2015-10-02, the Friday before the AFL Grand Final, was Victoria's first such holiday; Nager.Date lists it from 2016.",
    "Victoria's substitute days for Christmas and Boxing Day falling on a weekend follow Nager.Date's rule; where Victoria "
    "also kept the day itself, only the substitute is listed.",
]


def victorian_holidays(years: dict[int, list]) -> list[dict]:
    out = []
    for year in sorted(years):
        for h in years[year]:
            if h.get("global") or "AU-VIC" in (h.get("counties") or []):
                out.append({"date": h["date"], "name": h["name"], "national": bool(h.get("global"))})
    return sorted(out, key=lambda h: (h["date"], h["name"]))


def build_holidays() -> dict:
    years = {int(p.stem): json.loads(p.read_text(encoding="utf-8")) for p in sorted(HOLIDAYS.glob("[0-9][0-9][0-9][0-9].json"))}
    holidays = victorian_holidays(years)
    product = {
        "product": "street-holidays-victoria",
        "built": dt.datetime.now(UTC).isoformat(timespec="seconds"),
        "holidays": holidays, "years": [min(years), max(years)],
        "notes": HOLIDAY_NOTES,
        "source": {"dataset": "nager-holidays", "profile": "australia", "credit": "Nager.Date",
                   "rule": "national holidays (global) and those whose counties include AU-VIC",
                   "licence": "Nager.Date API terms: private and non-profit use", "retrieved_at": retrieved(HOLIDAYS)},
    }
    write(OUT / "street-holidays-victoria.json", product)
    register("streets", {"name": "street-holidays-victoria", "file": "street-holidays-victoria.json",
                         "title": "Public holidays of Victoria", "domain": "urban-science", "years": product["years"],
                         "holidays": len(holidays), "source": "Nager.Date"})
    print(f"  street-holidays-victoria: {len(holidays)} holidays, {min(years)}-{max(years)}")
    return product


# ---- uses of premises ------------------------------------------------------------------------------------

# A small classification of street uses. CLUE's establishments by their ANZSIC 2006 class (CLUE adds
# 0000 vacant space, 4111 convenience stores and gendered clothing and footwear classes 4254-4258), the
# longest prefix deciding; OpenStreetMap's holdings by their tags (osm_use).
USES = [
    {"id": "food-and-drink", "name": "Food and drink", "anzsic": ["451", "452", "453"],
     "about": "Cafes, restaurants, takeaway food, catering, pubs, bars and clubs (ANZSIC 451-453)"},
    {"id": "clothing-and-fashion", "name": "Clothing and fashion", "anzsic": ["4251", "4252", "4254", "4255", "4256", "4257", "4258", "4259"],
     "about": "Clothing, footwear and personal accessories retailing (ANZSIC 425 but 4253)"},
    {"id": "jewellery", "name": "Jewellery", "anzsic": ["4253"],
     "about": "Watch and jewellery retailing (ANZSIC 4253); jewellery manufacturing (2591) and wholesaling (3732) are wholesale-and-industry"},
    {"id": "bridal-and-wedding", "name": "Bridal and wedding", "anzsic": [],
     "about": "A heuristic on trading names: bridal, bride(s), wedding(s) as words (and, in OpenStreetMap, shop=wedding or bridal, "
              "clothes=wedding or bridal, and names with ブライダル, ウェディング, ウエディング, 花嫁, 婚礼); jewellers stay jewellery"},
    {"id": "other-retail", "name": "Other retail", "anzsic": ["39", "40", "41", "42", "43"],
     "about": "Every other retail: motor vehicles and fuel, food and groceries, household goods, books, pharmacies, department stores, "
              "non-store retail (ANZSIC 39-43 not above)"},
    {"id": "personal-services", "name": "Personal services", "anzsic": ["94", "951", "952", "9531", "9532", "9534", "9539", "7220"],
     "about": "Hairdressing and beauty, laundries, funerals and other personal services, repair and maintenance, travel agencies "
              "(ANZSIC 94, 95 but 9533 and 954-955, 7220)"},
    {"id": "accommodation", "name": "Accommodation", "anzsic": ["44"], "about": "Hotels, motels, serviced apartments, hostels (ANZSIC 44)"},
    {"id": "entertainment-and-culture", "name": "Entertainment and culture", "anzsic": ["5513", "89", "90", "91", "92"],
     "about": "Cinemas, museums, galleries, gardens and parks, performing arts, sport and recreation, gambling (ANZSIC 5513, 89-92)"},
    {"id": "education", "name": "Education", "anzsic": ["80", "81", "82"], "about": "Schools, vocational and higher education, other teaching (ANZSIC 80-82)"},
    {"id": "health", "name": "Health", "anzsic": ["84", "85", "86", "6970"],
     "about": "Hospitals, medical, dental, allied and residential care, veterinary services (ANZSIC 84-86, 6970)"},
    {"id": "offices-and-professional", "name": "Offices and professional services",
     "anzsic": ["54", "55", "56", "57", "58", "59", "6020", "67", "69", "70", "72", "73"],
     "about": "Information media and telecommunications, property and real estate, professional, scientific and technical services, "
              "computing, administrative and building services (ANZSIC 54-59 but 5513 and 6010, 67, 69-70 but 6970, 72-73 but 7220)"},
    {"id": "finance", "name": "Finance", "anzsic": ["62", "63", "64"], "about": "Banking, insurance, superannuation, financial services (ANZSIC 62-64)"},
    {"id": "community-and-public", "name": "Community and public", "anzsic": ["6010", "75", "76", "77", "87", "954", "955"],
     "about": "Libraries, government, defence, public order and safety, social assistance and child care, religious services, "
              "associations (ANZSIC 6010, 75-77, 87, 954-955)"},
    {"id": "transport", "name": "Transport", "anzsic": ["46", "47", "48", "49", "50", "51", "52", "661", "9533"],
     "about": "Road, rail, water and air transport, postal and courier services, transport support, vehicle rental, parking "
              "(ANZSIC 46-52, 661, 9533)"},
    {"id": "wholesale-and-industry", "name": "Wholesale and industry",
     "anzsic": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "53"],
     "about": "Agriculture, mining, manufacturing, utilities, construction, wholesale trade, warehousing (ANZSIC 01-38, 53)"},
    {"id": "vacant", "name": "Vacant", "anzsic": ["0000"], "about": "Vacant space (CLUE's 0000): a premises without a use"},
    {"id": "other", "name": "Other", "anzsic": [], "about": "Anything not classified above (other rental and hiring, ANZSIC 662-664)"},
]
_PREFIXES = sorted(((p, u["id"]) for u in USES for p in u["anzsic"]), key=lambda pu: -len(pu[0]))
BRIDAL = re.compile(r"\b(bridal|brides?|weddings?)\b|ブライダル|ウェディング|ウエディング|花嫁|婚礼", re.IGNORECASE)


def anzsic_use(code: str, name: str = "") -> str:
    """The street use of a CLUE establishment: its ANZSIC class by the longest prefix of USES; a name that
    names a wedding or a bride makes it bridal-and-wedding unless it is a jeweller."""
    use = next((u for p, u in _PREFIXES if code.startswith(p)), "other")
    if use not in ("jewellery", "vacant") and name and BRIDAL.search(name):
        return "bridal-and-wedding"
    return use


HOLDING_KEYS = ("shop", "amenity", "craft", "office", "tourism", "leisure", "healthcare")
# Tags that mark furniture and open spaces, not premises: such an element is not a holding.
NOT_PREMISES = {
    "amenity": {"bench", "waste_basket", "waste_disposal", "recycling", "bicycle_parking", "motorcycle_parking", "parking",
                "parking_space", "parking_entrance", "post_box", "telephone", "drinking_water", "vending_machine", "atm",
                "fountain", "clock", "shelter", "toilets", "taxi", "charging_station", "car_sharing", "bicycle_rental",
                "kick-scooter_rental", "grit_bin", "loading_dock", "smoking_area", "photo_booth", "letter_box", "bbq",
                "water_point", "compressed_air", "ticket_validator", "public_bookcase", "dog_toilet", "give_box", "lounger",
                "bicycle_repair_station", "street_lamp", "fire_hydrant", "bus_stop", "luggage_locker", "post_depot", "trolley_bay",
                "device_charging_station", "boat_rental", "mobility_scooter_rental", "payment_terminal", "table", "watering_place"},
    "tourism": {"artwork", "information", "viewpoint", "picnic_site", "attraction", "yes"},
    "leisure": {"park", "garden", "playground", "pitch", "picnic_table", "outdoor_seating", "track", "common", "nature_reserve",
                "dog_park", "swimming_area", "fitness_station", "firepit", "bandstand", "slipway", "bleachers", "schoolyard"},
}
OSM_USE = {
    "amenity": {"restaurant": "food-and-drink", "cafe": "food-and-drink", "fast_food": "food-and-drink", "bar": "food-and-drink",
                "pub": "food-and-drink", "biergarten": "food-and-drink", "food_court": "food-and-drink", "ice_cream": "food-and-drink",
                "nightclub": "entertainment-and-culture", "cinema": "entertainment-and-culture", "theatre": "entertainment-and-culture",
                "arts_centre": "entertainment-and-culture", "casino": "entertainment-and-culture", "gambling": "entertainment-and-culture",
                "karaoke_box": "entertainment-and-culture", "events_venue": "entertainment-and-culture", "music_venue": "entertainment-and-culture",
                "conference_centre": "entertainment-and-culture", "exhibition_centre": "entertainment-and-culture",
                "school": "education", "kindergarten": "education", "college": "education", "university": "education",
                "language_school": "education", "music_school": "education", "driving_school": "education", "training": "education",
                "prep_school": "education", "dancing_school": "education",
                "clinic": "health", "doctors": "health", "dentist": "health", "hospital": "health", "veterinary": "health",
                "pharmacy": "other-retail",
                "bank": "finance", "bureau_de_change": "finance", "money_transfer": "finance",
                "place_of_worship": "community-and-public", "community_centre": "community-and-public", "library": "community-and-public",
                "townhall": "community-and-public", "police": "community-and-public", "fire_station": "community-and-public",
                "courthouse": "community-and-public", "social_facility": "community-and-public", "childcare": "community-and-public",
                "social_centre": "community-and-public", "public_bath": "personal-services", "marketplace": "other-retail",
                "post_office": "transport", "bus_station": "transport", "ferry_terminal": "transport", "car_rental": "transport",
                "fuel": "other-retail", "car_wash": "personal-services", "coworking_space": "offices-and-professional"},
    "tourism": {"hotel": "accommodation", "motel": "accommodation", "hostel": "accommodation", "guest_house": "accommodation",
                "apartment": "accommodation", "chalet": "accommodation", "museum": "entertainment-and-culture",
                "gallery": "entertainment-and-culture", "theme_park": "entertainment-and-culture", "zoo": "entertainment-and-culture",
                "aquarium": "entertainment-and-culture"},
    "shop": {"jewelry": "jewellery", "jewellery": "jewellery", "watches": "jewellery", "wedding": "bridal-and-wedding",
             "bridal": "bridal-and-wedding", "clothes": "clothing-and-fashion", "shoes": "clothing-and-fashion",
             "boutique": "clothing-and-fashion", "fashion": "clothing-and-fashion", "fashion_accessories": "clothing-and-fashion",
             "bag": "clothing-and-fashion", "leather": "clothing-and-fashion", "hairdresser": "personal-services",
             "beauty": "personal-services", "massage": "personal-services", "laundry": "personal-services",
             "dry_cleaning": "personal-services", "tattoo": "personal-services", "travel_agency": "personal-services",
             "funeral_directors": "personal-services", "tailor": "personal-services", "optician": "health",
             "hearing_aids": "health", "medical_supply": "health", "wholesale": "wholesale-and-industry",
             "estate_agent": "offices-and-professional", "money_lender": "finance", "car_repair": "personal-services",
             "repair": "personal-services", "copyshop": "personal-services", "photo": "personal-services",
             "vacant": "vacant"},
    "office": {"government": "community-and-public", "ngo": "community-and-public", "association": "community-and-public",
               "religion": "community-and-public", "diplomatic": "community-and-public", "political_party": "community-and-public",
               "financial": "finance", "insurance": "finance", "financial_advisor": "finance", "tax_advisor": "offices-and-professional",
               "educational_institution": "education", "travel_agent": "personal-services"},
    "leisure": {"fitness_centre": "entertainment-and-culture", "sports_centre": "entertainment-and-culture",
                "bowling_alley": "entertainment-and-culture", "amusement_arcade": "entertainment-and-culture",
                "adult_gaming_centre": "entertainment-and-culture", "dance": "entertainment-and-culture",
                "escape_game": "entertainment-and-culture", "sauna": "personal-services", "hackerspace": "community-and-public",
                "stadium": "entertainment-and-culture", "swimming_pool": "entertainment-and-culture", "golf_course": "entertainment-and-culture"},
    "craft": {"tailor": "personal-services", "dressmaker": "personal-services", "shoemaker": "personal-services",
              "photographer": "personal-services", "watchmaker": "personal-services", "key_cutter": "personal-services",
              "locksmith": "personal-services", "hairdresser": "personal-services"},
}


def is_premises(tags: dict) -> bool:
    """Whether an element is a holding: it carries one of HOLDING_KEYS with a value that is not furniture."""
    return any(k in tags and tags[k] not in NOT_PREMISES.get(k, ()) for k in HOLDING_KEYS)


def osm_use(tags: dict) -> str:
    """The street use of an OpenStreetMap holding from its tags (see USES)."""
    shop = tags.get("shop")
    if shop in ("jewelry", "jewellery", "watches"):
        return "jewellery"
    clothes = tags.get("clothes", "")
    names = " ".join(tags.get(k, "") for k in ("name", "name:en", "name:ja", "brand"))
    if shop in ("wedding", "bridal") or re.search(r"wedding|bridal", clothes) or BRIDAL.search(names):
        return "bridal-and-wedding"
    for key in ("shop", "amenity", "tourism", "healthcare", "office", "leisure", "craft"):
        value = tags.get(key)
        if value is None or value in NOT_PREMISES.get(key, ()):
            continue
        if key == "healthcare":
            return "health"
        if value in OSM_USE.get(key, {}):
            return OSM_USE[key][value]
        if key == "shop":
            return "other-retail"
        if key == "office":
            return "offices-and-professional"
        if key == "craft":
            return "wholesale-and-industry"
    return "other"


# ---- local tangent planes ----------------------------------------------------------------------------------

A_WGS84, F_WGS84 = 6378137.0, 1 / 298.257223563
E2 = F_WGS84 * (2 - F_WGS84)


class Plane:
    """The local tangent plane about (lat0, lon0) on WGS 84: x = N cos φ0 Δλ east, y = M Δφ north (metres)."""

    def __init__(self, lat0: float, lon0: float):
        self.lat0, self.lon0 = lat0, lon0
        s = math.sin(math.radians(lat0))
        w = 1 - E2 * s * s
        self.n = A_WGS84 / math.sqrt(w)
        self.m = A_WGS84 * (1 - E2) / w ** 1.5
        self.kx = self.n * math.cos(math.radians(lat0)) * math.pi / 180
        self.ky = self.m * math.pi / 180

    def xy(self, lat, lon):
        return (np.asarray(lon) - self.lon0) * self.kx, (np.asarray(lat) - self.lat0) * self.ky

    def latlon(self, x, y):
        return self.lat0 + np.asarray(y) / self.ky, self.lon0 + np.asarray(x) / self.kx

    def error(self, x, y) -> float:
        """A first-order bound of the plane's error at (x, y): |x y| tan φ0 / M (metres)."""
        return float(np.max(np.abs(np.asarray(x) * np.asarray(y))) * abs(math.tan(math.radians(self.lat0))) / self.m)


# ---- terrain about the sensors -----------------------------------------------------------------------------

def read_clue(folder: Path) -> dict:
    cols = {"year": array("i"), "lat": array("d"), "lon": array("d"), "code": [], "name": []}
    skipped = 0
    with (folder / "establishments.csv").open(encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f, delimiter=";"):
            if not r["latitude"] or not r["longitude"]:
                skipped += 1
                continue
            cols["year"].append(int(r["census_year"][:4]))
            cols["lat"].append(float(r["latitude"]))
            cols["lon"].append(float(r["longitude"]))
            cols["code"].append(r["industry_anzsic4_code"])
            cols["name"].append(r["trading_name"])
    return {"year": np.frombuffer(cols["year"], np.int32), "lat": np.frombuffer(cols["lat"]), "lon": np.frombuffer(cols["lon"]),
            "code": cols["code"], "name": cols["name"], "skipped": skipped}


def default_terrain_years(present: list[int]) -> list[int]:
    """The latest census year and every fifth before it back to 2009, the first year of the counts: every year
    would make about 100 MB (a million establishment-years), these four about a fifth of it."""
    return [y for y in present if y >= 2009 and (present[-1] - y) % 5 == 0]


def build_terrain(sensors: list[dict], years: list[int] | None) -> dict:
    clue = read_clue(CLUE)
    present = sorted(set(clue["year"].tolist()))
    years = sorted(set(present) & set(years)) if years else default_terrain_years(present)
    keep_year = np.isin(clue["year"], years)
    uses_cache: dict[tuple, str] = {}
    out_sensors, total = {}, 0
    for e in sensors:
        if e.get("lat") is None:
            continue
        plane = Plane(e["lat"], e["lon"])
        x, y = plane.xy(clue["lat"], clue["lon"])
        near = np.nonzero(keep_year & (x * x + y * y <= TERRAIN_RADIUS ** 2))[0]
        rows = []
        for i in near[np.lexsort((clue["lon"][near], clue["lat"][near], clue["year"][near]))].tolist():
            code, name = clue["code"][i], clue["name"][i]
            key = (code, name)
            if key not in uses_cache:
                uses_cache[key] = anzsic_use(code, name)
            rows.append({"year": int(clue["year"][i]), "x_m": round(float(x[i]), 1), "y_m": round(float(y[i]), 1),
                         "anzsic4": code, "use": uses_cache[key], "name": name})
        out_sensors[str(e["id"])] = {"radius_m": TERRAIN_RADIUS, "establishments": rows}
        total += len(rows)
    product = {
        "product": "street-terrain-melbourne",
        "built": dt.datetime.now(UTC).isoformat(timespec="seconds"),
        "uses": USES,
        "years": years,
        "census_years_available": present,
        "positions": ("x_m east and y_m north of the sensor in its local tangent plane (WGS 84; error under 1 cm within 200 m). "
                      "CLUE places every establishment of a property at the property's point, so that many share a position."),
        "classification": ("use from the ANZSIC class by the longest prefix listed in uses; bridal-and-wedding is a heuristic on the "
                           "trading name (bridal, bride, brides, wedding, weddings as words), applied to every class but jewellery and vacant"),
        "sensors": out_sensors,
        "source": {"dataset": "melbourne-clue-establishments", "profile": "all", "credit": "CLUE: City of Melbourne, CC BY 4.0",
                   "retrieved_at": retrieved(CLUE), "rows_without_position": clue["skipped"]},
    }
    write(OUT / "street-terrain-melbourne.json", product)
    register("streets", {"name": "street-terrain-melbourne", "file": "street-terrain-melbourne.json",
                         "title": "Business establishments within 200 m of Melbourne's pedestrian sensors, by census year and use",
                         "domain": "urban-science", "years": years, "sensors": len(out_sensors), "establishments": total,
                         "source": "City of Melbourne, CLUE"})
    print(f"  street-terrain-melbourne: {total} establishment-years about {len(out_sensors)} sensors, census years {years[0]}-{years[-1]}")
    return product


# ---- streets and their holdings ------------------------------------------------------------------------------

# Each street: the Overpass response's file (data/registry/osm-streets.json), its name, city and time zone, the
# ends of its axis (from OpenStreetMap, 2026-10-30: the cross streets' shared nodes, the chōme boundaries' crossings,
# or a way's ends), the ways that make it (way_ids, or the name), and why it is here.
STREETS = [
    {"id": "tokyo-residential", "file": "tokyo-residential.json", "name": "A residential street of Okusawa 6-chōme, Setagaya",
     "city": "Tokyo", "timezone": "Asia/Tokyo", "way_ids": {48678505},
     "ends": [(35.6053074, 139.6636919), (35.6012358, 139.6624198)],
     "note": "OpenStreetMap way 48678505 (highway=residential, one way, 30 km/h), 466 m from beside Kanpachi-dōri southward, in "
             "Okusawa 6-chōme, Setagaya, beside the Den-en-chōfu garden suburb. Chosen as a street of uniform residential use, the "
             "case the account predicts nothing of (the user's paper, Claim 26.1): of nine long residential streets surveyed "
             "between Oyamadai and Den-en-chōfu, one of the four with no premises within 40 m, and the one in Setagaya; about "
             "170 buildings stand within 40 m of it. Its zoning was not checked against Setagaya's zoning map."},
    {"id": "melbourne-bourke-mall", "file": "melbourne-bourke-mall.json", "name": "Bourke Street Mall", "city": "Melbourne",
     "timezone": "Australia/Melbourne", "ends": [(-37.8142247, 144.9632407), (-37.8135238, 144.9656320)],
     "note": "The pedestrian mall of Bourke Street from Elizabeth Street to Swanston Street; sensors 1 and 2 (Bourke Street Mall North and South)."},
    {"id": "melbourne-swanston", "file": "melbourne-swanston.json", "name": "Swanston Street", "city": "Melbourne",
     "timezone": "Australia/Melbourne", "ends": [(-37.817461, 144.967445), (-37.809587, 144.963829)],
     "note": "Swanston Street through the Hoddle Grid, from Flinders Street to La Trobe Street; its carriageways' centreline."},
    {"id": "melbourne-lygon", "file": "melbourne-lygon.json", "name": "Lygon Street", "city": "Melbourne",
     "timezone": "Australia/Melbourne", "ends": [(-37.806957, 144.965879), (-37.792380, 144.968393)],
     "note": "Lygon Street, Carlton, from Victoria Street to Princes Street; sensors 31, 37 and 50."},
    {"id": "melbourne-errol", "file": "melbourne-errol.json", "name": "Errol Street", "city": "Melbourne",
     "timezone": "Australia/Melbourne", "ends": [(-37.805131, 144.949305), (-37.800987, 144.949996)],
     "note": "Errol Street, North Melbourne, from Victoria Street to Arden Street, the town's shopping street; sensors 70, 86 and 87."},
    {"id": "melbourne-macaulay", "file": "melbourne-macaulay.json", "name": "Macaulay Road", "city": "Melbourne",
     "timezone": "Australia/Melbourne", "ends": [(-37.794116, 144.927706), (-37.800269, 144.943938)],
     "note": "Macaulay Road from Kensington Road and Epsom Road, Kensington, to Arden Street and Dryburgh Street, North Melbourne; sensors 76, 85 and 180."},
]
AXIS_CLASSES = {"motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential", "living_street",
                "pedestrian", "road", "primary_link", "secondary_link", "tertiary_link", "trunk_link"}


def _segments(ways: list[dict], plane: Plane) -> list[tuple]:
    segs = []
    for w in ways:
        g = w["geometry"]
        x, y = plane.xy([p["lat"] for p in g], [p["lon"] for p in g])
        segs += [((x[i], y[i]), (x[i + 1], y[i + 1])) for i in range(len(g) - 1)]
    return segs


def _shortest_path(segs: list[tuple], a: tuple, b: tuple) -> list[tuple]:
    """The shortest path along the segments from the vertex nearest a to the vertex nearest b (Dijkstra).
    Where the street's ways do not join (a stretch mapped under no highway of its name), the ends of ways
    are bridged, first within 25 m, then within 150 m, a bridge counting twice its length."""
    key = lambda p: (round(p[0], 2), round(p[1], 2))  # noqa: E731
    graph: dict[tuple, dict[tuple, float]] = {}
    for p, q in segs:
        d = math.dist(p, q)
        graph.setdefault(key(p), {})[key(q)] = d
        graph.setdefault(key(q), {})[key(p)] = d
    nodes = list(graph)
    start = min(nodes, key=lambda n: math.dist(n, a))
    goal = min(nodes, key=lambda n: math.dist(n, b))
    ends = [n for n in nodes if len(graph[n]) == 1]
    for reach in (0, 25, 150):
        for i, p in enumerate(ends):
            for q in ends[i + 1:]:
                if math.dist(p, q) < reach and q not in graph[p]:
                    graph[p][q] = graph[q][p] = 2 * math.dist(p, q)
        dist, prev, heap = {start: 0.0}, {}, [(0.0, start)]
        while heap:
            d, n = heapq.heappop(heap)
            if n == goal:
                break
            if d > dist.get(n, math.inf):
                continue
            for m, w in graph[n].items():
                if d + w < dist.get(m, math.inf):
                    dist[m], prev[m] = d + w, n
                    heapq.heappush(heap, (d + w, m))
        if goal in dist:
            path = [goal]
            while path[-1] != start:
                path.append(prev[path[-1]])
            return path[::-1]
    raise ValueError("the axis's ways do not join its ends")


def _resample(line: np.ndarray, step: float) -> np.ndarray:
    d = np.r_[0, np.cumsum(np.hypot(*np.diff(line, axis=0).T))]
    s = np.linspace(0, d[-1], max(2, int(math.ceil(d[-1] / step)) + 1))
    return np.c_[np.interp(s, d, line[:, 0]), np.interp(s, d, line[:, 1])]


def _douglas_peucker(points: np.ndarray, tolerance: float) -> np.ndarray:
    keep = np.zeros(len(points), bool)
    keep[[0, -1]] = True
    stack = [(0, len(points) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        p, q = points[i], points[j]
        v = q - p
        n = math.hypot(*v) or 1.0
        w = points[i + 1:j] - p
        dist = np.abs(v[0] * w[:, 1] - v[1] * w[:, 0]) / n
        k = int(np.argmax(dist))
        if dist[k] > tolerance:
            keep[i + 1 + k] = True
            stack += [(i, i + 1 + k), (i + 1 + k, j)]
    return points[keep]


def centreline(segs: list[tuple], a: tuple, b: tuple, half_width: float = 22.0) -> np.ndarray:
    """The street's centreline from a to b: the shortest path along its ways, resampled every 5 m, each point
    moved along the path's normal to the middle of the outermost carriageways crossing that normal within
    half_width, lightly smoothed and simplified to 0.5 m."""
    path = _resample(np.array(_shortest_path(segs, a, b)), 5.0)
    tangent = np.gradient(path, axis=0)
    tangent /= np.hypot(tangent[:, 0], tangent[:, 1])[:, None]
    normal = np.c_[-tangent[:, 1], tangent[:, 0]]
    P = np.array([s[0] for s in segs])
    Q = np.array([s[1] for s in segs])
    centre = path.copy()
    for i, (q, n) in enumerate(zip(path, normal)):
        # Points q + λ n meeting the segment P + μ (Q − P), 0 ≤ μ ≤ 1.
        d = Q - P
        den = n[0] * d[:, 1] - n[1] * d[:, 0]
        ok = np.abs(den) > 1e-9
        w = P - q
        lam = np.where(ok, (w[:, 0] * d[:, 1] - w[:, 1] * d[:, 0]) / np.where(ok, den, 1), np.nan)
        mu = np.where(ok, (w[:, 0] * n[1] - w[:, 1] * n[0]) / np.where(ok, den, 1), np.nan)
        hit = ok & (mu >= 0) & (mu <= 1) & (np.abs(lam) <= half_width)
        if np.any(hit):
            centre[i] = q + n * (np.nanmin(lam[hit]) + np.nanmax(lam[hit])) / 2
    smooth = centre.copy()
    smooth[1:-1] = (centre[:-2] + centre[1:-1] + centre[2:]) / 3
    return _clip(_douglas_peucker(smooth, 0.5), a, b)


def unmapped_length(line: np.ndarray, segs: list[tuple], within: float = 25.0) -> float:
    """The length of the axis farther than `within` from every segment of the street's ways: a stretch the
    axis bridges (metres, at a resolution of 5 m)."""
    pts = _resample(line, 5.0)
    P = np.array([s[0] for s in segs])
    Q = np.array([s[1] for s in segs])
    d = Q - P
    L2 = np.maximum(np.sum(d * d, axis=1), 1e-12)
    far = 0
    for q in pts:
        t = np.clip(np.sum((q - P) * d, axis=1) / L2, 0, 1)
        far += float(np.min(np.hypot(*(P + t[:, None] * d - q).T))) > within
    return 5.0 * far


def project(line: np.ndarray, p: tuple) -> tuple[float, float, bool]:
    """The arc length s of the point of the polyline nearest p, the signed offset (positive on the left of
    the direction of increasing s), and whether p lies beside the polyline (its foot not beyond the ends)."""
    best = (math.inf, 0.0, 0.0, True)
    s0 = 0.0
    for k in range(len(line) - 1):
        u, v = line[k], line[k + 1]
        d = v - u
        L = math.hypot(*d)
        if L == 0:
            continue
        t = ((p[0] - u[0]) * d[0] + (p[1] - u[1]) * d[1]) / (L * L)
        beside = not ((k == 0 and t < 0) or (k == len(line) - 2 and t > 1))
        tc = min(1.0, max(0.0, t))
        foot = u + tc * d
        dist = math.dist(foot, p)
        if dist < best[0]:
            side = d[0] * (p[1] - u[1]) - d[1] * (p[0] - u[0])
            best = (dist, s0 + tc * L, math.copysign(dist, side) if side else 0.0, beside)
        s0 += L
    return best[1], best[2], best[3]


def _clip(line: np.ndarray, a: tuple, b: tuple) -> np.ndarray:
    """The part of the polyline between the feet of a and b, from a's."""
    d = np.r_[0, np.cumsum(np.hypot(*np.diff(line, axis=0).T))]
    sa, _, _ = project(line, a)
    sb, _, _ = project(line, b)
    lo, hi = sorted((sa, sb))
    inner = [line[i] for i in range(len(line)) if lo < d[i] < hi]
    at = lambda s: np.array([np.interp(s, d, line[:, 0]), np.interp(s, d, line[:, 1])])  # noqa: E731
    out = [at(lo)]
    for q in inner + [at(hi)]:
        if math.dist(q, out[-1]) > 1e-6:
            out.append(q)
    out = np.array(out)
    return out if sa <= sb else out[::-1]


def read_overpass(path: Path) -> dict:
    raw = json.loads(path.read_text(encoding="utf-8"))
    remark = raw.get("remark") or ""
    if "error" in remark.lower():
        raise ValueError(f"{path.name}: Overpass answered with an error ({remark.strip()[:200]}); fetch again: "
                         f"python data/acquisition/fetch.py osm-streets --refresh")
    return raw


def street_record(spec: dict, raw: dict) -> dict:
    elements = raw["elements"]
    axis_ways = [e for e in elements if e["type"] == "way" and "geometry" in e and e.get("tags", {}).get("highway") in AXIS_CLASSES
                 and not (len(e["geometry"]) > 3 and e["geometry"][0] == e["geometry"][-1])
                 and (spec.get("way_ids") is None or e["id"] in spec["way_ids"])
                 and (spec.get("name_tag") is None or e["tags"].get("name") == spec["name_tag"])]
    (lat_a, lon_a), (lat_b, lon_b) = spec["ends"]
    plane = Plane((lat_a + lat_b) / 2, (lon_a + lon_b) / 2)
    a, b = plane.xy(lat_a, lon_a), plane.xy(lat_b, lon_b)
    a, b = (float(a[0]), float(a[1])), (float(b[0]), float(b[1]))
    segs = _segments(axis_ways, plane)
    line = centreline(segs, a, b)
    length = float(np.sum(np.hypot(*np.diff(line, axis=0).T)))
    lat, lon = plane.latlon(line[:, 0], line[:, 1])
    holdings, buildings, seen = [], {}, set()
    parsed = {"ok": 0, "partial": 0, "unparsed": 0, "none": 0}
    for e in elements:
        tags = e.get("tags", {})
        ident = f"{e['type']}/{e['id']}"
        if ident in seen or "geometry" in e:
            continue
        where = e.get("center") or ({"lat": e["lat"], "lon": e["lon"]} if "lat" in e else None)
        if where is None:
            continue
        x, y = plane.xy(where["lat"], where["lon"])
        s, offset, beside = project(line, (float(x), float(y)))
        if not beside or abs(offset) > HOLDING_RADIUS:
            continue
        if "building" in tags and not is_premises(tags):
            seen.add(ident)
            buildings[tags["building"]] = buildings.get(tags["building"], 0) + 1
            continue
        if not is_premises(tags):
            continue
        seen.add(ident)
        if "building" in tags:
            buildings[tags["building"]] = buildings.get(tags["building"], 0) + 1
        raw_hours = tags.get("opening_hours")
        hours = opening_hours.parse(raw_hours)
        parsed[hours["status"]] += 1
        h = {"osm": ident, "name": tags.get("name") or tags.get("brand") or tags.get("operator"),
             "tags": {k: tags[k] for k in HOLDING_KEYS if k in tags}, "use": osm_use(tags),
             "lat": round(float(where["lat"]), 7), "lon": round(float(where["lon"]), 7),
             "s_m": round(s, 1), "offset_m": round(offset, 1),
             "hours": hours["hours"], "hours_raw": raw_hours, "hours_status": hours["status"]}
        if hours["ph"] is not None:
            h["hours_ph"] = hours["ph"]
        if tags.get("name:en") and tags.get("name:en") != h["name"]:
            h["name_en"] = tags["name:en"]
        holdings.append(h)
    holdings.sort(key=lambda h: (h["s_m"], h["offset_m"]))
    area = next((a["id"] for a in AREAS if all(a["bbox"][0] <= p <= a["bbox"][2] for p in lat)
                 and all(a["bbox"][1] <= q <= a["bbox"][3] for q in lon)), None)
    return {"id": spec["id"], "name": spec["name"], "city": spec["city"], "timezone": spec["timezone"], "area": area,
            "axis": [[round(float(p), 7), round(float(q), 7)] for p, q in zip(lat, lon)], "length_m": round(length, 1),
            "holdings": holdings, "note": spec["note"],
            "hours_parsed": parsed, "buildings": dict(sorted(buildings.items(), key=lambda kv: -kv[1])),
            "plane_error_m": round(plane.error(line[:, 0], line[:, 1]), 3), "axis_unmapped_m": unmapped_length(line, segs)}


def build_holdings() -> dict:
    streets, missing = [], []
    for spec in STREETS:
        path = OSM / spec["file"]
        if not path.exists() or spec["ends"] is None:
            missing.append(spec["id"])
            continue
        streets.append(street_record(spec, read_overpass(path)))
    if missing:
        print(f"  street-holdings: no response for {', '.join(missing)}")
    product = {
        "product": "street-holdings",
        "built": dt.datetime.now(UTC).isoformat(timespec="seconds"),
        "streets": streets,
        "missing": missing,
        "uses": [{k: u[k] for k in ("id", "name", "about")} for u in USES],
        "classification": ("use from the OpenStreetMap tags (shop, amenity, tourism, healthcare, office, leisure, craft, in that order), "
                           "as street_records.OSM_USE says; shop=jewelry or watches is jewellery; shop=wedding or bridal, clothes=wedding "
                           "or bridal, or a name naming a wedding or bride is bridal-and-wedding; furniture and open spaces "
                           "(benches, parking, artworks, parks) are not holdings"),
        "conventions": ("axis: [lat, lon] (WGS 84) from the first end to the second; s_m: arc length along the axis of the point "
                        "nearest the holding (a way's or relation's centre), offset_m: its distance, positive on the left of the "
                        "axis's direction; holdings within 40 m whose foot lies on the axis. hours: [day, open, close], day 0 Monday, "
                        "minutes from that day's 00:00 local time, close past 1440 after midnight; hours_ph: on public holidays, when "
                        "stated; hours_status: ok, partial, unparsed or none (data/formats/opening_hours.py). buildings: buildings "
                        "within 40 m by their building tag."),
        "sources": [{"dataset": "osm-streets", "profile": "streets", "credit": "© OpenStreetMap contributors, ODbL 1.0",
                     "retrieved_at": retrieved(OSM)}],
    }
    write(OUT / "street-holdings.json", product)
    register("streets", {"name": "street-holdings", "file": "street-holdings.json",
                         "title": "Holdings along chosen streets and their opening hours", "domain": "urban-science",
                         "streets": [s["id"] for s in streets], "holdings": sum(len(s["holdings"]) for s in streets),
                         "source": "OpenStreetMap (ODbL)"})
    for s in streets:
        print(f"  {s['id']}: {s['length_m']:.0f} m, {len(s['holdings'])} holdings, hours {s['hours_parsed']}")
    return product


# ---- the city model -------------------------------------------------------------------------------------------

# The areas of the city model, [lat0, lon0, lat1, lon1] (WGS 84). Melbourne's reaches north to -37.7910 (Lygon Street to
# Princes Street) and west to 144.9260 so that
# Kensington (sensors 76, 85, 180, all of Macaulay Road) is in it.
AREAS = [
    {"id": "melbourne", "name": "Melbourne: the Hoddle Grid and its neighbours, from Southbank to Carlton and from Kensington "
                                "to East Melbourne, where the pedestrian sensors are", "bbox": [-37.8260, 144.9260, -37.7910, 144.9790],
     "file": "melbourne"},
]
CITY = CACHE / "osm-city" / "areas"
LEVEL_M = 3.2                   # a storey
# Heights by building tag when neither height nor levels is tagged (m): a storey of 3.2 m for one- and
# two-storey kinds, about three storeys for the rest.
DEFAULT_HEIGHT = {"house": 6.4, "detached": 6.4, "semidetached_house": 6.4, "bungalow": 3.2, "terrace": 6.4,
                  "residential": 9.6, "apartments": 16.0, "dormitory": 12.8, "hotel": 19.2,
                  "commercial": 12.8, "retail": 6.4, "office": 16.0, "industrial": 8.0, "warehouse": 8.0,
                  "garage": 3.0, "garages": 3.0, "carport": 3.0, "shed": 3.0, "hut": 3.0, "kiosk": 3.0, "toilets": 3.0,
                  "roof": 4.0, "service": 3.0, "transformer_tower": 6.0, "church": 15.0, "cathedral": 25.0, "chapel": 8.0,
                  "temple": 9.0, "shrine": 6.0, "school": 9.6, "university": 12.8, "college": 12.8, "hospital": 16.0,
                  "public": 9.6, "civic": 9.6, "government": 12.8, "train_station": 9.6, "transportation": 8.0,
                  "parking": 9.6, "stadium": 20.0, "grandstand": 10.0, "construction": 9.6}
DEFAULT_HEIGHT_OTHER = 9.6      # yes and any other kind: three storeys
# Road widths (m) by highway class when neither width nor lanes is tagged; a one-way way of the first five
# classes is one carriageway of a road mapped twice, and takes the second figure.
DEFAULT_WIDTH = {"motorway": (14.0, 7.0), "trunk": (12.0, 7.0), "primary": (12.0, 6.5), "secondary": (10.0, 6.0),
                 "tertiary": (8.0, 5.5), "unclassified": (6.0, 6.0), "residential": (5.5, 5.5), "living_street": (4.5, 4.5),
                 "service": (3.5, 3.5), "road": (6.0, 6.0), "busway": (4.0, 4.0), "pedestrian": (6.0, 6.0),
                 "footway": (2.0, 2.0), "path": (1.5, 1.5), "cycleway": (2.0, 2.0), "steps": (2.5, 2.5),
                 "bridleway": (2.0, 2.0), "track": (3.0, 3.0), "motorway_link": (6.0, 5.0), "trunk_link": (6.0, 5.0),
                 "primary_link": (6.0, 5.0), "secondary_link": (6.0, 5.0), "tertiary_link": (5.5, 4.5)}
DEFAULT_WIDTH_OTHER = 3.0
SIDEWALK_M = 1.5                # added for each side a sidewalk tag names, beside lanes x 3.2
NOT_ROADS = {"platform", "proposed", "construction", "corridor", "elevator", "raceway", "bus_stop", "abandoned", "razed",
             "disused", "no", "emergency_bay", "services", "rest_area"}


def ecef(lat, lon, h=0.0):
    """Earth-centred, Earth-fixed coordinates (m) of geodetic latitude, longitude (degrees) and ellipsoidal height
    h on WGS 84: X = (N + h) cos φ cos λ, Y = (N + h) cos φ sin λ, Z = (N (1 − e²) + h) sin φ, N = a / √(1 − e² sin² φ)
    (Hofmann-Wellenhof, Lichtenegger and Wasle 2008, GNSS, eq. 2.10)."""
    phi, lam = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float))
    n = A_WGS84 / np.sqrt(1 - E2 * np.sin(phi) ** 2)
    return (n + h) * np.cos(phi) * np.cos(lam), (n + h) * np.cos(phi) * np.sin(lam), (n * (1 - E2) + h) * np.sin(phi)


def enu(lat, lon, lat0: float, lon0: float, h=0.0, h0: float = 0.0):
    """East, north, up (m) of points from an origin, exactly, through ECEF on WGS 84: the difference of the
    ECEF vectors rotated by R(φ0, λ0) (Hofmann-Wellenhof et al. 2008, eqs. 2.32-2.33). Points on the ellipsoid
    (h = 0) lie below the origin's tangent plane by about d²/2R (up < 0)."""
    x, y, z = ecef(lat, lon, h)
    x0, y0, z0 = ecef(lat0, lon0, h0)
    dx, dy, dz = x - x0, y - y0, z - z0
    sp, cp = math.sin(math.radians(lat0)), math.cos(math.radians(lat0))
    sl, cl = math.sin(math.radians(lon0)), math.cos(math.radians(lon0))
    e = -sl * dx + cl * dy
    n = -sp * cl * dx - sp * sl * dy + cp * dz
    u = cp * cl * dx + cp * sl * dy + sp * dz
    return e, n, u


def meridian_radius(lat: float) -> float:
    """M(φ) = a (1 − e²) / (1 − e² sin² φ)^(3/2), the radius of curvature of the meridian (m)."""
    s = math.sin(math.radians(lat))
    return A_WGS84 * (1 - E2) / (1 - E2 * s * s) ** 1.5


def signed_area(ring) -> float:
    """The shoelace area of a ring of (x, y), positive when counter-clockwise (x east, y north)."""
    r = np.asarray(ring, float)
    x, y = r[:, 0], r[:, 1]
    return float(0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y)))


def oriented(ring, counter_clockwise: bool) -> list:
    """The ring (not closed) turned counter-clockwise or clockwise."""
    ring = list(ring)
    if len(ring) > 1 and ring[0] == ring[-1]:
        ring = ring[:-1]
    return ring if (signed_area(ring) > 0) == counter_clockwise else ring[::-1]


def contains(ring, p) -> bool:
    """Whether point p lies inside the ring (even-odd rule)."""
    inside, (px, py) = False, p
    n = len(ring)
    for i in range(n):
        (ax, ay), (bx, by) = ring[i], ring[(i + 1) % n]
        if (ay > py) != (by > py) and px < ax + (py - ay) * (bx - ax) / (by - ay):
            inside = not inside
    return inside


def assemble_rings(paths: list[list[tuple]]) -> list[list[tuple]] | None:
    """Closed rings from ways joined at shared end nodes (a multipolygon's members of one role), each closed
    (first point repeated last); None when the ways do not close."""
    todo = [list(p) for p in paths if len(p) >= 2]
    rings = []
    while todo:
        ring = todo.pop(0)
        while ring[0] != ring[-1]:
            for k, p in enumerate(todo):
                if p[0] == ring[-1]:
                    ring += p[1:]
                elif p[-1] == ring[-1]:
                    ring += p[::-1][1:]
                elif p[-1] == ring[0]:
                    ring = p + ring[1:]
                elif p[0] == ring[0]:
                    ring = p[::-1] + ring[1:]
                else:
                    continue
                todo.pop(k)
                break
            else:
                return None
        if len(ring) >= 4:
            rings.append(ring)
    return rings


def clip_to_box(ring: list[tuple], box: list[float]) -> list[tuple]:
    """A ring of (lat, lon) clipped to [lat0, lon0, lat1, lon1] (Sutherland and Hodgman 1974)."""
    lat0, lon0, lat1, lon1 = box
    edges = [(0, lat0, True), (0, lat1, False), (1, lon0, True), (1, lon1, False)]
    out = ring[:-1] if ring and ring[0] == ring[-1] else list(ring)
    for axis, value, keep_above in edges:
        if not out:
            break
        inside = lambda p: p[axis] >= value if keep_above else p[axis] <= value  # noqa: E731
        src, out = out, []
        for i, p in enumerate(src):
            q = src[i - 1]
            if inside(p):
                if not inside(q):
                    out.append(_cut(q, p, axis, value))
                out.append(p)
            elif inside(q):
                out.append(_cut(q, p, axis, value))
    return out


def _cut(q, p, axis, value):
    t = (value - q[axis]) / (p[axis] - q[axis])
    return (q[0] + t * (p[0] - q[0]), q[1] + t * (p[1] - q[1]))


def clip_line_to_box(line: list[tuple], box: list[float]) -> list[list[tuple]]:
    """The pieces of a polyline of (lat, lon) inside [lat0, lon0, lat1, lon1], each segment clipped by
    Liang and Barsky's (1984) parametric test."""
    lat0, lon0, lat1, lon1 = box
    pieces, current = [], []
    for p, q in zip(line, line[1:]):
        t0, t1 = 0.0, 1.0
        d = (q[0] - p[0], q[1] - p[1])
        for num, den in ((p[0] - lat0, -d[0]), (lat1 - p[0], d[0]), (p[1] - lon0, -d[1]), (lon1 - p[1], d[1])):
            if den == 0:
                if num < 0:
                    t0, t1 = 1.0, 0.0
                continue
            r = num / den
            if den < 0:
                t0 = max(t0, r)
            else:
                t1 = min(t1, r)
        if t0 > t1:
            if current:
                pieces.append(current)
                current = []
            continue
        a = (p[0] + t0 * d[0], p[1] + t0 * d[1])
        b = (p[0] + t1 * d[0], p[1] + t1 * d[1])
        if current and current[-1] == a:
            current.append(b)
        else:
            if current:
                pieces.append(current)
            current = [a, b]
        if t1 < 1.0:
            pieces.append(current)
            current = []
    if current:
        pieces.append(current)
    return [c for c in pieces if len(c) >= 2]


LENGTH = re.compile(r"^\s*~?\s*(-?\d+(?:[.,]\d+)?)\s*(m|metres?|meters?|ft|feet|')?\s*(?:(\d+(?:\.\d+)?)\s*(?:\"|in))?\s*$", re.IGNORECASE)


def parse_length(value: str | None) -> float | None:
    """A length tag in metres: '12', '12 m', '12.5m', '40 ft', 40'6\"; the first of a list ('12;15'); None if unreadable."""
    if value is None:
        return None
    m = LENGTH.match(str(value).split(";")[0])
    if not m:
        return None
    x = float(m.group(1).replace(",", "."))
    if (m.group(2) or "").lower() in ("ft", "feet", "'"):
        x = x * 0.3048 + (float(m.group(3)) * 0.0254 if m.group(3) else 0.0)
    return x


def _number(value) -> float | None:
    try:
        return float(str(value).split(";")[0].strip().replace(",", "."))
    except (TypeError, ValueError):
        return None


def building_heights(tags: dict, kind: str) -> tuple[float, float, float | None, str]:
    """(height_m, min_height_m, levels, source): the height tag; else building:levels x 3.2 m with roof:levels x 1.6 m;
    else the default of its kind. The base: min_height, else building:min_level x 3.2 m, else 0."""
    levels = _number(tags.get("building:levels"))
    roof = _number(tags.get("roof:levels")) or 0.0
    height = parse_length(tags.get("height"))
    if height is not None and height > 0:
        source = "height"
    elif levels is not None and levels >= 0:
        height, source = levels * LEVEL_M + roof * LEVEL_M * 0.5, "levels"
    else:
        height, source = DEFAULT_HEIGHT.get(kind, DEFAULT_HEIGHT_OTHER), "default"
    base = parse_length(tags.get("min_height"))
    if base is None:
        min_level = _number(tags.get("building:min_level"))
        base = min_level * LEVEL_M if min_level is not None else 0.0
    if height <= base:
        height = base + LEVEL_M
    return height, max(0.0, base), levels, source


def road_width(tags: dict) -> tuple[float, str]:
    """(width_m, source): the width tag; else lanes x 3.2 m with 1.5 m for each side a sidewalk tag names; else
    the default of its class (DEFAULT_WIDTH)."""
    width = parse_length(tags.get("width") or tags.get("est_width"))
    if width is not None and width > 0:
        return width, "width"
    lanes = _number(tags.get("lanes"))
    if lanes is not None and lanes > 0:
        sides = {"both": 2, "left": 1, "right": 1, "yes": 2}.get(tags.get("sidewalk", ""), 0)
        return lanes * LEVEL_M + sides * SIDEWALK_M, "lanes"
    major = tags.get("highway") in ("motorway", "trunk", "primary", "secondary", "tertiary")
    pair = DEFAULT_WIDTH.get(tags.get("highway"), (DEFAULT_WIDTH_OTHER, DEFAULT_WIDTH_OTHER))
    return (pair[1] if major and tags.get("oneway") in ("yes", "-1") else pair[0]), "default"


def space_kind(tags: dict) -> str | None:
    """The kind of an open space: park, water, plaza, green or rail; None if it is none of them."""
    leisure, natural, landuse = tags.get("leisure"), tags.get("natural"), tags.get("landuse")
    if natural == "water" or tags.get("waterway") in ("riverbank", "dock") or landuse in ("basin", "reservoir"):
        return "water"
    if leisure in ("park", "garden", "playground", "common", "nature_reserve", "recreation_ground", "dog_park") \
            or landuse in ("recreation_ground", "village_green"):
        return "park"
    if tags.get("place") == "square" or tags.get("highway") in ("pedestrian", "footway") and tags.get("area") == "yes" \
            or tags.get("area:highway") in ("pedestrian", "footway") or tags.get("type") == "multipolygon" and tags.get("highway") == "pedestrian":
        return "plaza"
    if leisure == "pitch" or natural in ("wood", "scrub", "grassland", "heath", "wetland") \
            or landuse in ("grass", "forest", "meadow", "cemetery"):
        return "green"
    if landuse == "railway":
        return "rail"
    return None


def _paths(member_or_way) -> list[tuple]:
    return [(p["lat"], p["lon"]) for p in member_or_way.get("geometry") or [] if p]


def polygons(element: dict) -> list[tuple[list, list]] | None:
    """The polygons of a closed way or a multipolygon relation, as (outer ring, inner rings) of (lat, lon), each
    inner ring given to the outer containing it; None when they cannot be assembled."""
    if element["type"] == "way":
        ring = _paths(element)
        return [(ring, [])] if len(ring) >= 4 and ring[0] == ring[-1] else None
    members = [m for m in element.get("members", []) if m.get("type") == "way"]
    if any(not m.get("geometry") or any(p is None for p in m["geometry"]) for m in members):
        return None
    outer = assemble_rings([_paths(m) for m in members if m.get("role") in ("outer", "")])
    inner = assemble_rings([_paths(m) for m in members if m.get("role") == "inner"])
    if not outer or inner is None:
        return None
    out = [(o, []) for o in outer]
    for ring in inner:
        host = next((o for o in out if contains(o[0], ring[0])), None)
        if host is not None:
            host[1].append(ring)
    return out


class _Projector:
    """ENU of an area's origin, with the equirectangular plane beside it for the comparison."""

    def __init__(self, lat0: float, lon0: float):
        self.lat0, self.lon0 = lat0, lon0
        self.plane = Plane(lat0, lon0)
        self.worst_plane, self.worst_drop = 0.0, 0.0

    def __call__(self, ring: list[tuple]) -> list[list[float]]:
        lat = np.array([p[0] for p in ring])
        lon = np.array([p[1] for p in ring])
        e, n, u = enu(lat, lon, self.lat0, self.lon0)
        px, py = self.plane.xy(lat, lon)
        self.worst_plane = max(self.worst_plane, float(np.max(np.hypot(px - e, py - n))))
        self.worst_drop = max(self.worst_drop, float(np.max(-u)))
        return [[round(float(a), 2) + 0.0, round(float(b), 2) + 0.0] for a, b in zip(e, n)]      # no -0.0


def city_area(spec: dict) -> tuple[dict, dict]:
    """One area of the city model from its tiles' Overpass responses (each element kept once), and its tally."""
    lat0 = (spec["bbox"][0] + spec["bbox"][2]) / 2
    lon0 = (spec["bbox"][1] + spec["bbox"][3]) / 2
    project = _Projector(lat0, lon0)
    tally = {"buildings": {"height": 0, "levels": 0, "default": 0}, "parts": 0, "buildings_not_assembled": 0, "roads": 0,
             "spaces": 0, "spaces_not_assembled": 0, "spaces_outside": 0, "ways_not_roads": 0}
    elements: dict[tuple, dict] = {}
    tiles = sorted(CITY.glob(f"{spec['file']}-*.json"), key=lambda q: int(q.stem.rsplit("-", 1)[1]))
    if not tiles:
        raise FileNotFoundError(f"no tiles of {spec['id']} in {CITY}; run: python data/acquisition/fetch.py osm-city")
    for path in tiles:
        for e in read_overpass(path)["elements"]:
            elements.setdefault((e["type"], e["id"]), e)
    tally["tiles"] = len(tiles)
    is_building = lambda tags: "building" in tags or "building:part" in tags  # noqa: E731
    in_relations = {m["ref"] for e in elements.values() if e["type"] == "relation" and "building" in e.get("tags", {})
                    for m in e.get("members", []) if m.get("type") == "way"}
    buildings, parts = [], []
    for e in elements.values():
        tags = e.get("tags", {})
        if not is_building(tags):
            continue
        if e["type"] == "way" and e["id"] in in_relations and "building:part" not in tags:
            continue                           # an outline drawn by its multipolygon
        is_part = "building:part" in tags and "building" not in tags
        if is_part and e["type"] != "way":
            continue
        polys = polygons(e)
        if polys is None:
            tally["buildings_not_assembled"] += 1
            continue
        kind = (f"part:{tags['building:part']}" if is_part else tags.get("building", "yes"))
        height, base, levels, source = building_heights(tags, tags.get("building:part") if is_part else tags.get("building", "yes"))
        for k, (outer, holes) in enumerate(polys):
            ident = f"{e['type']}/{e['id']}" + (f"#{k + 1}" if len(polys) > 1 else "")
            b = {"id": ident}
            if tags.get("name"):
                b["name"] = tags["name"]
            b.update({"kind": kind, "footprint": oriented(project(outer), True),
                      "holes": [oriented(project(h), False) for h in holes],
                      "height_m": round(height, 2), "min_height_m": round(base, 2)})
            if levels is not None:
                b["levels"] = levels
            b["height_source"] = source
            (parts if is_part else buildings).append(b)
    # Parts within an outline: the outline is marked, the part names it.
    boxes = [(np.min(np.array(b["footprint"]), axis=0), np.max(np.array(b["footprint"]), axis=0)) for b in buildings]
    for p in parts:
        c = np.mean(np.array(p["footprint"]), axis=0)
        for b, (lo, hi) in zip(buildings, boxes):
            if lo[0] <= c[0] <= hi[0] and lo[1] <= c[1] <= hi[1] and contains(b["footprint"], c):
                p["part_of"] = b["id"]
                b["has_parts"] = True
                break
    for b in buildings:
        tally["buildings"][b["height_source"]] += 1
    tally["parts"] = len(parts)
    roads = []
    for e in elements.values():
        tags = e.get("tags", {})
        if "highway" not in tags or is_building(tags):
            continue
        line = _paths(e)
        kind = tags.get("highway")
        if e["type"] != "way" or len(line) < 2 or kind in NOT_ROADS or tags.get("area") == "yes":
            tally["ways_not_roads"] += 1
            continue
        width, source = road_width(tags)
        pieces = clip_line_to_box(line, spec["bbox"])
        if not pieces:
            tally["ways_not_roads"] += 1
            continue
        for k, piece in enumerate(pieces):
            r = {"id": f"way/{e['id']}" + (f"#{k + 1}" if len(pieces) > 1 else ""), "kind": kind}
            if tags.get("name"):
                r["name"] = tags["name"]
            r.update({"width_m": round(width, 2), "width_source": source, "line": project(piece)})
            roads.append(r)
    spaces = []
    for e in elements.values():
        tags = e.get("tags", {})
        kind = None if is_building(tags) else space_kind(tags)
        polys = polygons(e) if kind else None
        if kind is None:
            continue
        if polys is None:
            tally["spaces_not_assembled"] += 1
            continue
        for k, (outer, _) in enumerate(polys):
            ring = clip_to_box(outer, spec["bbox"])
            if len(ring) < 3:
                tally["spaces_outside"] += 1
                continue
            s = {"id": f"{e['type']}/{e['id']}" + (f"#{k + 1}" if len(polys) > 1 else ""), "kind": kind}
            if tags.get("name"):
                s["name"] = tags["name"]
            s["ring"] = oriented(project(ring), True)
            spaces.append(s)
    tally["roads"], tally["spaces"] = len(roads), len(spaces)
    area = {"id": spec["id"], "name": spec["name"], "origin": {"lat": round(lat0, 7), "lon": round(lon0, 7)},
            "bbox": spec["bbox"], "frame": "local ENU metres from origin: x east, y north",
            "projection": ("exact: geodetic (WGS 84, on the ellipsoid) to ECEF, rotated to east-north-up at the origin; the up "
                           "component, the curvature's drop below the tangent plane, is dropped"),
            "curvature_drop_m": round(project.worst_drop, 2),
            "equirectangular_error_m": round(project.worst_plane, 3),
            "buildings": buildings + parts, "roads": roads, "spaces": spaces}
    return area, tally


def build_city() -> dict:
    import time
    started = time.monotonic()
    areas, tallies = [], {}
    for spec in AREAS:
        area, tally = city_area(spec)
        areas.append(area)
        tallies[spec["id"]] = tally
    product = {
        "product": "city-model",
        "built": dt.datetime.now(UTC).isoformat(timespec="seconds"),
        "sources": [{"dataset": "osm-city", "profile": "areas", "credit": "© OpenStreetMap contributors, ODbL 1.0",
                     "retrieved_at": retrieved(CITY)}],
        "heights": ("height_m: the height tag (units stripped; feet converted), else building:levels x 3.2 m plus roof:levels x 1.6 m, "
                    "else a default of the building's kind (height_source default): " + json.dumps(DEFAULT_HEIGHT) +
                    f", any other {DEFAULT_HEIGHT_OTHER} m. min_height_m: min_height, else building:min_level x 3.2 m, else 0. "
                    "A height not above its base is raised to the base plus a storey. building:part ways are entries of kind "
                    "'part:<value>' with part_of naming the outline that holds their centre, which carries has_parts."),
        "widths": ("width_m: the width tag, else lanes x 3.2 m plus 1.5 m for each side a sidewalk tag names, else a default "
                   "of the highway class, the second figure for a one-way carriageway of a major road: " + json.dumps(DEFAULT_WIDTH) +
                   f", any other {DEFAULT_WIDTH_OTHER} m"),
        "spaces": ("park: leisure park, garden, playground, common, nature reserve, recreation ground, dog park; water: natural=water, "
                   "waterway riverbank or dock, landuse basin or reservoir; plaza: place=square, pedestrian or footway areas; green: "
                   "pitches, woods, scrub, grassland, heath, wetland, landuse grass, forest, meadow, cemetery; rail: landuse=railway. "
                   "Outer rings only, clipped to the area's box (Sutherland-Hodgman in latitude and longitude); roads are clipped "
                   "to it too (Liang-Barsky), a way leaving and re-entering the box in pieces id#k; buildings are whole."),
        "rings": "footprints and rings counter-clockwise, holes clockwise, not closed; coordinates rounded to 0.01 m",
        "tally": tallies,
        "areas": areas,
    }
    path = OUT / "city-model.json"
    write(path, product)
    register("streets", {"name": "city-model", "file": "city-model.json", "title": "Buildings, roads and open spaces of Melbourne's centre in local ENU metres",
                         "domain": "urban-science", "areas": [a["id"] for a in areas], "source": "OpenStreetMap (ODbL)"})
    for a in areas:
        t = tallies[a["id"]]
        print(f"  {a['id']}: {len(a['buildings'])} buildings ({t['buildings']}, {t['parts']} parts), {t['roads']} roads, "
              f"{t['spaces']} spaces; not assembled: {t['buildings_not_assembled']} buildings, {t['spaces_not_assembled']} spaces; equirectangular error {a['equirectangular_error_m']} m")
    print(f"  city-model.json: {path.stat().st_size / 1e6:.1f} MB in {time.monotonic() - started:.1f} s")
    return product


# ---- the command ---------------------------------------------------------------------------------------------

NEEDS = {
    "sensors": [ARCHIVE, LIVE / "counts.csv", LIVE / "sensors.csv"],
    "counts": [ARCHIVE, LIVE / "counts.csv", LIVE / "sensors.csv"],
    "weather": [WEATHER / "hourly.json"],
    "holidays": [HOLIDAYS / "retrieval.json"],
    "terrain": [CLUE / "establishments.csv", LIVE / "sensors.csv"],
    "holdings": [OSM / "retrieval.json"],
    "city": [CITY / "retrieval.json"],
}
FETCH = {"sensors": "melbourne-pedestrian-archive melbourne-pedestrian-counts", "counts": "melbourne-pedestrian-archive melbourne-pedestrian-counts",
         "weather": "open-meteo-hourly", "holidays": "nager-holidays", "terrain": "melbourne-clue-establishments melbourne-pedestrian-counts",
         "holdings": "osm-streets", "city": "osm-city"}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", choices=list(NEEDS), action="append", help="build only this product (repeatable)")
    parser.add_argument("--terrain-years", default="", help="census years of the terrain, comma-separated (default: the latest and every fifth before it back to 2009)")
    args = parser.parse_args(argv)
    sys.stdout.reconfigure(errors="replace")
    chosen = args.only or list(NEEDS)
    ok = True
    sensors = None
    for what in chosen:
        lacking = [p for p in NEEDS[what] if not p.exists()]
        if lacking:
            print(f"{what}: nothing fetched; run: python data/acquisition/fetch.py {FETCH[what]}", file=sys.stderr)
            ok = False
            continue
        print(what)
        if what in ("sensors", "counts"):
            if sensors is None:
                sensors = build_sensors_and_counts(write_counts=True)
        elif what == "weather":
            build_weather()
        elif what == "holidays":
            build_holidays()
        elif what == "terrain":
            if sensors is None:
                path = OUT / "street-sensors.json"
                sensors = json.loads(path.read_text(encoding="utf-8")) if path.exists() else None
            if sensors is None:
                print("terrain: the sensors are needed first; run with --only sensors", file=sys.stderr)
                ok = False
                continue
            years = [int(y) for y in args.terrain_years.split(",") if y.strip()] or None
            build_terrain(sensors["sensors"], years)
        elif what == "holdings":
            build_holdings()
        elif what == "city":
            build_city()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
