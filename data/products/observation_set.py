"""Build observations products: what stations measured, where and when.

An observations product (kind `observations`, schema observations/1) is a
table with one row per station and time:

- the station's identifier and name, its longitude and latitude (degrees,
  WGS 84) and its elevation (m above mean sea level, null where unknown);
- the time of the report (ISO 8601, UTC, as 2026-10-22T09:00:00Z);
- one column per variable, in SI units (K, Pa, m s-1, m, kg m-2; per cent for
  relative humidity; degree, clockwise from true north, for directions), null
  where missing.

Variables are named by their CF standard names where one exists (the CF
Standard Name Table, https://cfconventions.org/), and each carries its long
name, units and the level at which it is measured. Values beyond the extremes
ever observed (LIMITS) are set to null as gross errors, and counted in the
product's `screened`.

Two sources are built here, each from a dataset of data/registry/:

- `metar`: each aerodrome's latest METAR or SPECI, worldwide (NOAA Aviation
  Weather Center), about 5,000 stations;
- `ndbc`: each buoy's, coastal station's and platform's latest observation
  (NOAA National Data Buoy Center), about 850 stations.

A source is a reader, a function from a retrieval folder to rows, with its
variables (SOURCES); another station network is added as another reader. The
writer, write(), serves every source, and our own sensors too
(observation/sensors/store.py).

Each retrieval holds only the latest report of each station. Rebuilt every
hour, after a fresh retrieval, a product keeps a rolling window: the previous
product's rows are merged with the new ones, a station's report at a given
time is kept once (the newer retrieval's copy wins), and rows older than
`hours` (default 48) before the window's end are dropped. The window ends at the
retrieval's time, or at the newest report if that is later. The product thus
becomes a time series per station, sampled as often as it is fetched. Rows are
ordered by station, then time.

Output, in data/cache/products/: observations-<source>.json

Usage::

    python data/acquisition/fetch.py metar ndbc --refresh      # the latest reports
    python data/products/observation_set.py                     # every fetched source
    python data/products/observation_set.py metar --hours 24
    python data/products/observation_set.py ndbc --fresh        # without the previous rows
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import html
import json
import math
import re
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import product_index  # noqa: E402

UTC = dt.timezone.utc
HOURS = 48

# ---- units --------------------------------------------------------------------
# Exact by definition unless stated: the Celsius scale (BIPM, The International
# System of Units, 9th edition, 2019, section 2.3.1); the nautical mile and the
# knot (ibid., table 8); the international statute mile and the inch (NIST
# Special Publication 811, 2008, appendix B.8). The inch of mercury is the
# conventional one (NIST SP 811, B.8: 3.386 389 E+03 Pa), as altimeter settings
# use. A depth of precipitation becomes a mass per area with water's density
# taken as 1000 kg m-3, as WMO's practice is (1 mm = 1 kg m-2; WMO-No. 8, 2018,
# volume I, chapter 6).

ZERO_CELSIUS = 273.15            # K
KNOT = 1852.0 / 3600.0           # m s-1
STATUTE_MILE = 1609.344          # m
NAUTICAL_MILE = 1852.0           # m
INCH = 0.0254                    # m
INCH_OF_MERCURY = 3386.389       # Pa
HECTOPASCAL = 100.0              # Pa
WATER_DENSITY = 1000.0           # kg m-3


def kelvin(celsius: float | None) -> float | None:
    return None if celsius is None else celsius + ZERO_CELSIUS


def saturation_vapour_pressure(celsius: float) -> float:
    """Saturation vapour pressure over a plane surface of pure water, in Pa.

    The Magnus form with the coefficients of Alduchov and Eskridge (1996,
    eq. 21, 'AERK'): e_s = 610.94 exp(17.625 t / (t + 243.04)), t in degrees
    Celsius. Within about 0.4 per cent of the reference formulation from -40 to
    +50 C; over supercooled water below 0 C, as WMO reports humidity
    (WMO-No. 8, 2018, volume I, annex 4.B). Outside -40 to +50 C the error grows.
    Alduchov, O. A. and Eskridge, R. E. (1996), Improved Magnus form
    approximation of saturation vapor pressure, Journal of Applied Meteorology
    35, 601-609.
    """
    return 610.94 * math.exp(17.625 * celsius / (celsius + 243.04))


def relative_humidity(temperature: float | None, dew_point: float | None) -> float | None:
    """Relative humidity with respect to water (per cent) from the temperature
    and the dew point, both in K: 100 e_s(T_d) / e_s(T), by the Magnus form
    above. The error is that of the formula (a few tenths of a per cent of the
    value) and, above all, of the inputs: with whole degrees, as most METARs
    give them, about 3 per cent at 60 per cent."""
    if temperature is None or dew_point is None:
        return None
    return 100.0 * saturation_vapour_pressure(dew_point - ZERO_CELSIUS) / \
        saturation_vapour_pressure(temperature - ZERO_CELSIUS)


# ---- gross errors -------------------------------------------------------------
# Limits beyond the extremes ever observed at the Earth's surface (SI units), so
# that a value outside them is an error of coding or of a sensor and is set to
# null (in the products) or refused (by observation/sensors/store.py). Extremes
# from the WMO Archive of Weather and Climate Extremes (https://wmo.asu.edu/)
# unless stated. They catch gross errors only: a plausible wrong value passes.

LIMITS = {
    "air_temperature": (180.0, 335.0),            # -89.2 C, Vostok, 1983; 56.7 C, Death Valley, 1913
    "dew_point_temperature": (170.0, 315.0),      # the highest reported are near 35 C, by the Persian Gulf
    "relative_humidity": (0.0, 102.0),            # little supersaturation; beyond, the dew point exceeds the temperature
    "air_pressure": (30000.0, 110000.0),          # about 330 hPa on the summit of Everest
    "air_pressure_at_sea_level": (85000.0, 110000.0),   # 870 hPa, Typhoon Tip, 1979; about 1085 hPa, Mongolia, 2001
    "wind_speed": (0.0, 115.0),
    "wind_speed_of_gust": (0.0, 115.0),           # 113.2 m s-1, Barrow Island, 1996 (Tropical Cyclone Olivia)
    "wind_from_direction": (0.0, 360.0),
    "visibility_in_air": (0.0, 250000.0),
    "precipitation_amount": (0.0, 2000.0),
    "sea_surface_temperature": (269.0, 312.0),    # sea water freezes near -1.9 C; the warmest seas reach about 36 C
    "sea_surface_wave_significant_height": (0.0, 25.0),    # 19 m, North Atlantic, 2013, from a buoy
    "sea_surface_wave_period_at_variance_spectral_density_maximum": (0.0, 40.0),
    "sea_surface_wave_from_direction_at_variance_spectral_density_maximum": (0.0, 360.0),
}


def plausible(name: str, value: float) -> bool:
    lo, hi = LIMITS.get(name, (-math.inf, math.inf))
    return lo <= value <= hi


# ---- variables ----------------------------------------------------------------

def variable(name: str, long_name: str, units: str, level: str | None, digits: int,
             standard: bool = True, note: str | None = None) -> dict:
    """A variable's description; `digits` (decimals kept) is used here and not written."""
    v = {"name": name, "long_name": long_name, "standard_name": name if standard else None,
         "units": units, "level": level, "digits": digits}
    if note:
        v["note"] = note
    return v


METAR_VARIABLES = [
    variable("air_temperature", "Air temperature", "K", "screen height, 1.25-2 m", 2,
             note="whole degrees Celsius in the report's body; tenths where the remarks' T group gives them (North America)"),
    variable("dew_point_temperature", "Dew-point temperature", "K", "screen height, 1.25-2 m", 2),
    variable("relative_humidity", "Relative humidity (with respect to water)", "%", "screen height, 1.25-2 m", 1,
             note="computed here from the temperature and the dew point by the Magnus form of Alduchov and Eskridge (1996)"),
    variable("air_pressure_at_sea_level", "Pressure at mean sea level (the SLP group where reported, else the altimeter setting QNH)",
             "Pa", "mean sea level", 0,
             note="QNH reduces the station's pressure to sea level through the ICAO standard atmosphere, the SLP group "
                  "(North America only) through the observed temperature; at high stations in cold weather they differ "
                  "by several hPa"),
    variable("wind_speed", "Wind speed (a 10-minute mean; 2 minutes in North America)", "m s-1", "10 m", 2),
    variable("wind_from_direction", "Direction the wind comes from (true north)", "degree", "10 m", 0,
             note="null for calm or variable wind"),
    variable("wind_speed_of_gust", "Speed of gusts (the peak 3-second wind of the last 10 minutes; 5 s in the United States)",
             "m s-1", "10 m", 2),
    variable("visibility_in_air", "Prevailing visibility", "m", None, 0,
             note="9999 and CAVOK are 10 km or more and given as 10000 m; '10+' statute miles as 16093 m: lower bounds"),
    variable("precipitation_amount", "Precipitation since the last hourly report (the P group, United States)", "kg m-2",
             "surface", 2, note="a trace is given by the source as 0.005 in (0.127 kg m-2)"),
]

NDBC_VARIABLES = [
    variable("air_temperature", "Air temperature", "K", "3-5 m above the sea on buoys (by hull)", 2),
    variable("dew_point_temperature", "Dew-point temperature", "K", "as the air temperature", 2),
    variable("relative_humidity", "Relative humidity (with respect to water)", "%", "as the air temperature", 1,
             note="computed here from the temperature and the dew point by the Magnus form of Alduchov and Eskridge (1996)"),
    variable("air_pressure_at_sea_level", "Pressure at sea level", "Pa", "mean sea level", 0,
             note="reduced to sea level for C-MAN stations and Great Lakes buoys (NWS Technical Procedures Bulletin 291)"),
    variable("wind_speed", "Wind speed (an 8-minute mean on buoys, 2 minutes on land)", "m s-1",
             "3-5 m above the sea on buoys (by hull)", 2),
    variable("wind_from_direction", "Direction the wind comes from (true north)", "degree", "as the wind speed", 0),
    variable("wind_speed_of_gust", "Speed of gusts (the peak 5 or 8 second wind)", "m s-1", "as the wind speed", 2),
    variable("sea_surface_temperature", "Sea-surface temperature", "K", "near the surface (from the hull's waterline on buoys)", 2),
    variable("sea_surface_wave_significant_height", "Significant wave height (the mean of the highest third, 20 minutes)",
             "m", "sea surface", 2),
    variable("sea_surface_wave_period_at_variance_spectral_density_maximum", "Dominant wave period", "s", "sea surface", 1),
    variable("sea_surface_wave_from_direction_at_variance_spectral_density_maximum",
             "Direction the dominant waves come from (true north)", "degree", "sea surface", 0),
    variable("visibility_in_air", "Visibility", "m", None, 0, note="buoys report only 0 to 1.6 nautical miles"),
]


# ---- times --------------------------------------------------------------------

def iso(t: dt.datetime) -> str:
    return t.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_time(text: str) -> dt.datetime:
    t = dt.datetime.fromisoformat(text.strip().replace("Z", "+00:00"))
    return (t if t.tzinfo else t.replace(tzinfo=UTC)).astimezone(UTC)


def number(text: str | None) -> float | None:
    if text is None:
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


# ---- METAR --------------------------------------------------------------------
# Code form FM 15 METAR / FM 16 SPECI (WMO-No. 306, volume I.1). The AWC's
# columns are decoded already; the report's own groups are read where they are
# more exact than the columns: winds reported in m s-1 (which the columns give
# in rounded knots), visibility in metres (which the columns give in rounded
# statute miles, '6+' for 9999), and QNH in whole hPa (given in inches).

WIND_GROUP = re.compile(r"(?<![\w/])(\d{3}|VRB)(\d{2,3})(?:G(\d{2,3}))?(KT|MPS|KMH)(?!\w)")
WIND_UNIT = {"KT": KNOT, "MPS": 1.0, "KMH": 1.0 / 3.6}
QNH_HPA = re.compile(r"(?<!\w)Q(\d{4})(?!\w)")
QNH_INHG = re.compile(r"(?<!\w)A(\d{4})(?!\w)")


def metar_body(raw: str) -> str:
    """The report without its remarks (RMK) and trend (TEMPO, BECMG, NOSIG)."""
    return re.split(r"\s(?:RMK|TEMPO|BECMG|NOSIG)\b", raw, maxsplit=1)[0]


def metar_wind(raw: str) -> tuple[float | None, float, float | None] | None:
    """(direction in degrees or None, speed and gust in m s-1) from the wind group,
    or None without one. Calm (00000) and variable (VRB) winds have no direction;
    a wind from true north is coded 360."""
    m = WIND_GROUP.search(metar_body(raw))
    if not m:
        return None
    direction, speed, gust, unit = m.groups()
    factor = WIND_UNIT[unit]
    speed_ms = int(speed) * factor
    bearing = None if direction == "VRB" or int(speed) == 0 or int(direction) == 0 else float(int(direction))
    return bearing, speed_ms, int(gust) * factor if gust else None


def metar_visibility(raw: str) -> float | None:
    """Visibility in metres from a report that gives it in metres (four digits,
    9999 for 10 km or more) or as CAVOK; None where it is given in statute miles."""
    body = metar_body(raw)
    m = WIND_GROUP.search(body)
    if not m:
        return None
    for token in body[m.end():].split():
        if re.fullmatch(r"\d{3}V\d{3}", token):      # the range of a varying direction
            continue
        if token == "CAVOK":
            return 10000.0
        v = re.fullmatch(r"(\d{4})(?:NDV)?", token)
        if v:
            metres = int(v.group(1))
            return 10000.0 if metres == 9999 else float(metres)
        return None
    return None


def metar_qnh(raw: str) -> float | None:
    """The altimeter setting QNH in Pa, from a Q group (hPa) or an A group (hundredths of inHg)."""
    body = metar_body(raw)
    m = QNH_HPA.search(body)
    if m:
        return int(m.group(1)) * HECTOPASCAL
    m = QNH_INHG.search(body)
    return int(m.group(1)) / 100.0 * INCH_OF_MERCURY if m else None


def statute_miles(text: str | None) -> float | None:
    """Metres from the AWC's statute miles, where '10+' (at least 10) is taken as 10."""
    value = number((text or "").rstrip("+"))
    return None if value is None else value * STATUTE_MILE


def read_metar(folder: Path) -> list[dict]:
    with gzip.open(folder / "metars.cache.csv.gz", "rt", encoding="utf-8", newline="") as handle:
        rows = csv.reader(handle)
        # Older versions of the cache begin with lines of status before the header.
        header = next(r for r in rows if r and r[0] == "raw_text")
        col = {}
        for i, name in enumerate(header):
            col.setdefault(name, i)                   # some names repeat (sky_cover); the first is kept
        body = [r for r in rows if len(r) >= len(header) - 1]
    names = {}
    path = folder / "stations.cache.json.gz"
    if path.exists():
        for s in json.loads(gzip.decompress(path.read_bytes()).decode("utf-8")):
            if s.get("id"):
                names[s["id"]] = (s.get("site") or s["id"], number(str(s.get("elev"))))
    out = []
    for r in body:
        get = lambda name: (r[col[name]].strip() or None) if name in col and col[name] < len(r) else None  # noqa: E731
        ident, raw = get("station_id"), get("raw_text") or ""
        lat, lon = number(get("latitude")), number(get("longitude"))
        if not ident or lat is None or lon is None or not get("observation_time"):
            continue
        t, td = kelvin(number(get("temp_c"))), kelvin(number(get("dewpoint_c")))
        wind = metar_wind(raw)
        if wind is None:
            speed = number(get("wind_speed_kt"))
            gust = number(get("wind_gust_kt"))
            bearing = number(get("wind_dir_degrees"))
            wind = (bearing if bearing and speed else None, None if speed is None else speed * KNOT,
                    None if gust is None else gust * KNOT)
        visibility = metar_visibility(raw)
        if visibility is None:
            visibility = statute_miles(get("visibility_statute_mi"))
        slp = number(get("sea_level_pressure_mb"))
        if slp is not None:
            pressure = slp * HECTOPASCAL
        else:
            pressure = metar_qnh(raw)
            if pressure is None and number(get("altim_in_hg")) is not None:
                pressure = number(get("altim_in_hg")) * INCH_OF_MERCURY
        precipitation = number(get("precip_in"))
        name, elevation = names.get(ident, (ident, None))
        out.append({
            "id": ident, "name": name, "lon": lon, "lat": lat,
            "elevation": number(get("elevation_m")) if get("elevation_m") is not None else elevation,
            "time": iso(parse_time(get("observation_time"))),
            "values": {
                "air_temperature": t, "dew_point_temperature": td, "relative_humidity": relative_humidity(t, td),
                "air_pressure_at_sea_level": pressure, "wind_speed": wind[1], "wind_from_direction": wind[0],
                "wind_speed_of_gust": wind[2], "visibility_in_air": visibility,
                "precipitation_amount": None if precipitation is None else precipitation * INCH * WATER_DENSITY,
            },
        })
    return out


# ---- NDBC ---------------------------------------------------------------------
# latest_obs.txt: definitions at https://www.ndbc.noaa.gov/faq/measdes.shtml.

def ndbc_stations(path: Path) -> dict:
    """{id: (name, platform)} from the station table; names are HTML-escaped there."""
    out = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="latin-1").splitlines():
        if line.startswith("#"):
            continue
        parts = [html.unescape(p).strip() for p in line.split("|")]
        if len(parts) >= 5 and parts[0]:
            out[parts[0].upper()] = (parts[4] or parts[0].upper(), parts[2])
    return out


def read_ndbc(folder: Path) -> list[dict]:
    lines = (folder / "latest_obs.txt").read_text(encoding="latin-1").splitlines()
    header = lines[0].lstrip("#").split()
    stations = ndbc_stations(folder / "station_table.txt")
    out = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) != len(header):
            continue
        cell = dict(zip(header, parts))
        get = lambda name: None if cell.get(name, "MM") == "MM" else number(cell[name])  # noqa: E731
        lat, lon = get("LAT"), get("LON")
        if lat is None or lon is None:
            continue
        try:
            t = dt.datetime(int(cell["YYYY"]), int(cell["MM"]), int(cell["DD"]), int(cell["hh"]), int(cell["mm"]), tzinfo=UTC)
        except (KeyError, ValueError):
            continue
        ident = cell["STN"].upper()
        name, platform = stations.get(ident, (ident, ""))
        ta, td = kelvin(get("ATMP")), kelvin(get("DEWP"))
        pres, vis, wdir = get("PRES"), get("VIS"), get("WDIR")
        out.append({
            "id": ident, "name": name, "lon": lon, "lat": lat,
            # A buoy floats at the sea surface; other platforms' heights are not in the table.
            "elevation": 0.0 if "buoy" in platform.lower() else None,
            "time": iso(t),
            "values": {
                "air_temperature": ta, "dew_point_temperature": td, "relative_humidity": relative_humidity(ta, td),
                "air_pressure_at_sea_level": None if pres is None else pres * HECTOPASCAL,
                "wind_speed": get("WSPD"), "wind_from_direction": wdir, "wind_speed_of_gust": get("GST"),
                "sea_surface_temperature": kelvin(get("WTMP")),
                "sea_surface_wave_significant_height": get("WVHT"),
                "sea_surface_wave_period_at_variance_spectral_density_maximum": get("DPD"),
                "sea_surface_wave_from_direction_at_variance_spectral_density_maximum": get("MWD"),
                "visibility_in_air": None if vis is None else vis * NAUTICAL_MILE,
            },
        })
    return out


SOURCES = {
    "metar": {"dataset": "metar", "profile": "latest", "title": "Weather at aerodromes, METAR",
              "reader": read_metar, "variables": METAR_VARIABLES},
    "ndbc": {"dataset": "ndbc", "profile": "latest", "title": "Buoys and coastal stations, NDBC",
             "reader": read_ndbc, "variables": NDBC_VARIABLES},
}


# ---- the product --------------------------------------------------------------

def rounded(value, digits: int | None):
    """value to `digits` decimals (None: as it is), or None if it is not a finite number."""
    if value is None or isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        return None
    if digits is None:
        return float(value)
    return round(float(value), digits) if digits > 0 else float(round(value))


def to_row(reading: dict, variables: list[dict], screened: dict | None = None) -> list:
    """A reading ({id, name, lon, lat, elevation, time, values}) as a row of the
    product; values beyond the gross-error limits become null, counted by
    variable in `screened`."""
    values = []
    for v in variables:
        x = rounded(reading["values"].get(v["name"]), v.get("digits", 3))
        if x is not None and not plausible(v["name"], x):
            x = None
            if screened is not None:
                screened[v["name"]] = screened.get(v["name"], 0) + 1
        values.append(x)
    return [reading["id"], reading.get("name") or reading["id"], round(float(reading["lon"]), 5),
            round(float(reading["lat"]), 5), rounded(reading.get("elevation"), 1), reading["time"]] + values


def merge(previous: list[list], rows: list[list], start: str, end: str) -> list[list]:
    """Rows of both, once per (station, time), the later list winning, within
    [start, end] (ISO texts in one format, compared as text), by station then time."""
    kept = {}
    for row in previous + rows:
        if start <= row[5] <= end:
            kept[(row[0], row[5])] = row
    return [kept[k] for k in sorted(kept)]


def previous_rows(path: Path, fields: list[str]) -> list[list]:
    """The rows of an earlier product, rearranged to the present fields (null where it lacked one)."""
    if not path.exists():
        return []
    try:
        old = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    if old.get("schema") != "observations/1":
        return []
    where = {f: i for i, f in enumerate(old["fields"])}
    return [[row[where[f]] if f in where else None for f in fields] for row in old["rows"]]


def write(name: str, title: str, domain: str, source: dict, variables: list[dict], readings: list[dict],
          hours: float = HOURS, end: dt.datetime | None = None, keep_previous: bool = True) -> Path:
    """Writes observations-<name>.json (schema observations/1) and registers it.

    readings: [{id, name, lon, lat, elevation, time (ISO, UTC), values: {variable: value in SI}}].
    With keep_previous, the earlier product's rows are merged in (the rolling
    window). The window ends at `end` (default: now) or at the newest reading,
    whichever is later, and begins `hours` before. Written to a temporary file
    and renamed, so that a page reading it never sees half a file.
    """
    out = product_index.OUT
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"observations-{name}.json"
    fields = ["id", "name", "lon", "lat", "elevation", "time"] + [v["name"] for v in variables]
    screened = {}
    # A reading without a true position (a source's fill value, such as -99.99, -99.99) is left out, and counted.
    placed = [r for r in readings if -90 <= float(r["lat"]) <= 90 and -180 <= float(r["lon"]) <= 360 and (float(r["lat"]), float(r["lon"])) != (-99.99, -99.99)]
    if len(placed) < len(readings):
        screened["position"] = len(readings) - len(placed)
    rows = [to_row(r, variables, screened) for r in placed]
    newest = max((r[5] for r in rows), default=None)
    stop = iso(end or dt.datetime.now(UTC))
    stop = max(stop, newest) if newest else stop
    start = iso(parse_time(stop) - dt.timedelta(hours=hours))
    rows = merge(previous_rows(path, fields) if keep_previous else [], rows, start, stop)
    product = {
        "schema": "observations/1", "name": name, "title": title, "domain": domain, "source": source,
        "variables": [{k: v for k, v in var.items() if k != "digits"} for var in variables],
        "fields": fields, "order": "by station, then time",
        "rows": rows, "window": {"start": start, "end": stop},
        # This build's values set to null as gross errors (LIMITS), by variable.
        "screened": {"limits": {v["name"]: LIMITS[v["name"]] for v in variables if v["name"] in LIMITS},
                     "nulled": screened},
    }
    part = path.with_name(path.name + ".part")
    part.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    part.replace(path)
    stations = len({r[0] for r in rows})
    product_index.register("observations", {
        "name": name, "file": path.name, "title": title, "domain": domain, "count": stations, "rows": len(rows),
        "window": product["window"], "latest": max((r[5] for r in rows), default=None),
        "variables": [v["name"] for v in variables],
        "source": {k: source.get(k) for k in ("dataset", "provider", "credit", "retrieved_at")},
    })
    return path


def build(name: str, hours: float = HOURS, keep_previous: bool = True) -> Path:
    spec = SOURCES[name]
    folder = DATA / "cache" / spec["dataset"] / spec["profile"]
    manifest = json.loads((DATA / "registry" / f"{spec['dataset']}.json").read_text(encoding="utf-8"))
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    readings = spec["reader"](folder)
    source = {"dataset": spec["dataset"], "profile": spec["profile"], "provider": manifest["provider"],
              "credit": manifest.get("credit"), "licence": manifest["licence"], "citation": manifest["citation"],
              "retrieved_at": retrieval["retrieved_at"]}
    return write(name, spec["title"], manifest.get("domain", "meteorology"), source, spec["variables"], readings,
                 hours=hours, end=parse_time(retrieval["retrieved_at"]), keep_previous=keep_previous)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sources", nargs="*", help=f"sources to build, of {', '.join(SOURCES)} (default: every fetched one)")
    parser.add_argument("--hours", type=float, default=HOURS, help=f"the rolling window, in hours (default {HOURS})")
    parser.add_argument("--fresh", action="store_true", help="start the window again, without the previous rows")
    args = parser.parse_args(argv)
    unknown = [n for n in args.sources if n not in SOURCES]
    if unknown:
        print(f"unknown source(s): {', '.join(unknown)}; choose from {', '.join(SOURCES)}", file=sys.stderr)
        return 2
    names = args.sources or [n for n, s in SOURCES.items()
                             if (DATA / "cache" / s["dataset"] / s["profile"] / "retrieval.json").exists()]
    if not names:
        print("nothing fetched; run: python data/acquisition/fetch.py metar ndbc", file=sys.stderr)
        return 1
    ok = True
    for name in names:
        spec = SOURCES[name]
        if not (DATA / "cache" / spec["dataset"] / spec["profile"] / "retrieval.json").exists():
            print(f"{name}: not fetched; run: python data/acquisition/fetch.py {spec['dataset']}", file=sys.stderr)
            ok = False
            continue
        path = build(name, args.hours, keep_previous=not args.fresh)
        product = json.loads(path.read_text(encoding="utf-8"))
        stations = {r[0] for r in product["rows"]}
        print(f"{path.name}: {len(product['rows'])} rows from {len(stations)} stations, "
              f"{product['window']['start']} to {product['window']['end']}, {path.stat().st_size / 1e6:.1f} MB")
        if product["screened"]["nulled"]:
            print(f"  set to null as gross errors: {product['screened']['nulled']}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
