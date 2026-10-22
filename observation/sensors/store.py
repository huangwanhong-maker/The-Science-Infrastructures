"""The store of our own sensors' readings: home weather stations, tide poles,
buoys, loggers.

A source is one set of sensors under one description: its variables, with
their units, and its stations, each with a place. Its readings are appended as
they arrive, through this module (from ingest.py, simulate.py or the
application's API), and kept in one directory per source:

    <STORE>/<source>/source.json      the description, its stations, the count and time of the last reading
    <STORE>/<source>/readings.jsonl   one reading a line, in SI units, as appended

A reading is {"station": id, "lon": .., "lat": .., "elevation"?: .., "time":
ISO 8601 or milliseconds since 1970 (UTC), "values": {variable: value}}. Each
is checked before it is kept (validate()): the variables must be the source's,
the values finite numbers within the gross-error limits of
data/products/observation_set.py, the place on the Earth (longitude -180 to
180, latitude -90 to 90 degrees, WGS 84) and the time between 1990-01-01 and a
quarter of an hour from now (an earlier time is a clock not set, or seconds
given where milliseconds are expected). Values are kept in SI units: a source
may declare its variables in other units (degree_C, hPa, km h-1, ...), and its
values are converted on arrival where the conversion is affine and known
(CONVERSIONS); the description records both.

build_product() writes the readings of the last hours as an observations
product (schema observations/1, data/README.md section 4) named
sensor-<source>, which the Earth Atlas draws and its plug-ins read with
context.observations('sensor-<source>'). Positions are published with it.

Every function is safe to call from several threads (one lock). Several
processes appending to one source at once are not coordinated: each append is
one write of whole lines, which is atomic on POSIX for small writes, and the
description is rewritten whole.

STORE is data/cache/observations, or the directory named by the environment
variable SCIENCE_OBSERVATIONS; tests set store.STORE to a directory of their
own.

Limits: every reading of a source is read once into memory, and then only what
is appended (by size), so a source of a few million readings is the practical
limit of this draft; a busier one would be split by month.
"""

from __future__ import annotations

import datetime as dt
import json
import math
import os
import re
import sys
import threading
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "data" / "products"))
import observation_set  # noqa: E402

UTC = dt.timezone.utc
STORE = Path(os.environ.get("SCIENCE_OBSERVATIONS") or ROOT / "data" / "cache" / "observations")
NAME = re.compile(r"[a-z0-9]+(?:-[a-z0-9]+)*")
STATION = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.\-]{0,63}")
VARIABLE = re.compile(r"[a-z][a-z0-9_]{0,99}")
EARLIEST = dt.datetime(1990, 1, 1, tzinfo=UTC)
AHEAD = dt.timedelta(minutes=15)

_lock = threading.RLock()
_subscribers: list = []
_cache: dict = {}               # source -> {"size": bytes read, "readings": [...]}

# The variables a source may declare by their CF standard names, with their SI
# units and long names; others may be declared with units of their own, and are
# kept as they come. The gross-error limits are observation_set.LIMITS.
KNOWN = {
    "air_temperature": ("K", "Air temperature"),
    "dew_point_temperature": ("K", "Dew-point temperature"),
    "relative_humidity": ("%", "Relative humidity (with respect to water)"),
    "air_pressure": ("Pa", "Air pressure at the station"),
    "air_pressure_at_sea_level": ("Pa", "Air pressure at mean sea level"),
    "wind_speed": ("m s-1", "Wind speed"),
    "wind_from_direction": ("degree", "Direction the wind comes from (true north)"),
    "wind_speed_of_gust": ("m s-1", "Speed of gusts"),
    "precipitation_amount": ("kg m-2", "Precipitation since the last reading"),
    "visibility_in_air": ("m", "Visibility"),
    "sea_surface_temperature": ("K", "Sea-surface temperature"),
    "sea_surface_wave_significant_height": ("m", "Significant wave height"),
    "water_surface_height_above_reference_datum": ("m", "Water level above the station's datum"),
    "surface_downwelling_shortwave_flux_in_air": ("W m-2", "Sunlight on a level surface (global irradiance)"),
    "soil_temperature": ("K", "Soil temperature"),
    "mass_concentration_of_pm2p5_ambient_aerosol_particles_in_air": ("kg m-3", "Fine particles (PM2.5)"),
    "mass_concentration_of_pm10_ambient_aerosol_particles_in_air": ("kg m-3", "Particles (PM10)"),
    "mole_fraction_of_carbon_dioxide_in_air": ("1", "Carbon dioxide (mole fraction)"),
}
LIMITS = dict(observation_set.LIMITS) | {
    "water_surface_height_above_reference_datum": (-50.0, 50.0),
    "surface_downwelling_shortwave_flux_in_air": (0.0, 1500.0),   # above the solar constant, 1361 W m-2 (Kopp and Lean 2011)
    "soil_temperature": (200.0, 360.0),
    "mass_concentration_of_pm2p5_ambient_aerosol_particles_in_air": (0.0, 0.01),   # 10 000 ug m-3
    "mass_concentration_of_pm10_ambient_aerosol_particles_in_air": (0.0, 0.02),
    "mole_fraction_of_carbon_dioxide_in_air": (0.0, 0.05),
}

# Affine conversions to the SI unit, si = a * x + b, by (SI unit, declared unit).
# Exact by definition (BIPM SI Brochure, 2019; NIST SP 811, 2008, appendix B.8),
# except the conventional inch and millimetre of mercury (NIST SP 811, B.8).
_F = 5.0 / 9.0
CONVERSIONS = {
    "K": {"K": (1.0, 0.0), "degree_C": (1.0, 273.15), "degC": (1.0, 273.15), "°C": (1.0, 273.15),
          "degree_F": (_F, 459.67 * _F), "degF": (_F, 459.67 * _F), "°F": (_F, 459.67 * _F)},
    "Pa": {"Pa": (1.0, 0.0), "hPa": (100.0, 0.0), "mbar": (100.0, 0.0), "mb": (100.0, 0.0), "kPa": (1000.0, 0.0),
           "bar": (1e5, 0.0), "inHg": (observation_set.INCH_OF_MERCURY, 0.0), "mmHg": (133.322387415, 0.0)},
    "m s-1": {"m s-1": (1.0, 0.0), "m/s": (1.0, 0.0), "km h-1": (1 / 3.6, 0.0), "km/h": (1 / 3.6, 0.0),
              "kt": (observation_set.KNOT, 0.0), "knot": (observation_set.KNOT, 0.0), "mph": (0.44704, 0.0)},
    "degree": {"degree": (1.0, 0.0), "deg": (1.0, 0.0), "degrees": (1.0, 0.0), "°": (1.0, 0.0)},
    "%": {"%": (1.0, 0.0), "percent": (1.0, 0.0), "1": (100.0, 0.0)},
    "kg m-2": {"kg m-2": (1.0, 0.0), "mm": (1.0, 0.0), "in": (25.4, 0.0)},       # depths of water, 1000 kg m-3
    "m": {"m": (1.0, 0.0), "cm": (0.01, 0.0), "mm": (0.001, 0.0), "km": (1000.0, 0.0), "ft": (0.3048, 0.0),
          "mi": (observation_set.STATUTE_MILE, 0.0), "nmi": (observation_set.NAUTICAL_MILE, 0.0)},
    "W m-2": {"W m-2": (1.0, 0.0), "W/m2": (1.0, 0.0)},
    "kg m-3": {"kg m-3": (1.0, 0.0), "ug m-3": (1e-9, 0.0), "µg m-3": (1e-9, 0.0), "ug/m3": (1e-9, 0.0),
               "µg/m3": (1e-9, 0.0)},
    "1": {"1": (1.0, 0.0), "mol mol-1": (1.0, 0.0), "ppm": (1e-6, 0.0), "ppb": (1e-9, 0.0)},
}


# ---- names, times and places -------------------------------------------------

def _root() -> Path:
    return Path(STORE)


def _folder(name: str) -> Path:
    """The source's directory; the name is checked, so that it cannot leave the store."""
    if not isinstance(name, str) or len(name) > 64 or not NAME.fullmatch(name):
        raise ValueError(f"source name {name!r}: lowercase letters, digits and single hyphens, as garden-station")
    return _root() / name


def _iso(t: dt.datetime) -> str:
    return t.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_time(value) -> dt.datetime:
    """An instant (UTC) from an ISO 8601 text (UTC when it has no offset), a
    number of milliseconds since 1970-01-01 UTC, or a datetime."""
    if isinstance(value, dt.datetime):
        return value if value.tzinfo else value.replace(tzinfo=UTC)
    if isinstance(value, bool):
        raise ValueError("a time, not a truth value")
    if isinstance(value, (int, float)):
        if not math.isfinite(value):
            raise ValueError("time is not finite")
        return dt.datetime.fromtimestamp(value / 1000.0, UTC)
    text = str(value).strip()
    if re.fullmatch(r"-?\d+(\.\d*)?", text):
        return dt.datetime.fromtimestamp(float(text) / 1000.0, UTC)
    t = dt.datetime.fromisoformat(text.replace("Z", "+00:00"))
    return (t if t.tzinfo else t.replace(tzinfo=UTC)).astimezone(UTC)


def _finite(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _place(lon, lat, elevation=None) -> str | None:
    """Why a place is not plausible, or None."""
    if not (_finite(lon) and _finite(lat)):
        return "longitude and latitude must be numbers"
    if not (-180.0 <= lon <= 180.0 and -90.0 <= lat <= 90.0):
        return f"({lon}, {lat}) is not on the Earth (longitude -180 to 180, latitude -90 to 90)"
    if elevation is not None and not (_finite(elevation) and -500.0 <= elevation <= 9000.0):
        return f"elevation {elevation!r} is not between -500 and 9000 m"
    return None


def convert(value: float, units: str, si: str) -> float:
    """value in `units` as `si`, by CONVERSIONS; units equal to si pass unchanged."""
    if units == si:
        return value
    a, b = CONVERSIONS[si][units]
    return a * value + b


# ---- descriptions -------------------------------------------------------------

def _variable(spec) -> dict:
    """A variable's description from {name, units?, long_name?} or 'name:units[:long name]'."""
    if isinstance(spec, str):
        parts = spec.split(":", 2)
        spec = {"name": parts[0].strip()}
        if len(parts) > 1 and parts[1].strip():
            spec["units"] = parts[1].strip()
        if len(parts) > 2:
            spec["long_name"] = parts[2].strip()
    name = str(spec.get("name", ""))
    if not VARIABLE.fullmatch(name):
        raise ValueError(f"variable name {name!r}: lowercase letters, digits and underscores")
    if name in KNOWN:
        si, long_name = KNOWN[name]
        units = spec.get("input_units") or spec.get("units") or si
        if units != si and units not in CONVERSIONS.get(si, {}):
            known = ", ".join(CONVERSIONS.get(si, {si: None}))
            raise ValueError(f"{name}: no conversion from {units!r} to {si}; use one of {known}")
        out = {"name": name, "units": si, "long_name": spec.get("long_name") or long_name, "standard_name": name}
        if units != si:
            out["input_units"] = units
        return out
    if not spec.get("units"):
        raise ValueError(f"{name} is not a CF standard name this store knows; give its units, as {name}:V")
    return {"name": name, "units": str(spec["units"]), "long_name": spec.get("long_name") or name.replace("_", " "),
            "standard_name": None}


def _station(spec: dict) -> dict:
    ident = str(spec.get("id", spec.get("station", "")))
    if not STATION.fullmatch(ident):
        raise ValueError(f"station id {ident!r}: letters, digits, '_', '.', '-', at most 64")
    why = _place(spec.get("lon"), spec.get("lat"), spec.get("elevation"))
    if why:
        raise ValueError(f"station {ident}: {why}")
    return {"id": ident, "name": str(spec.get("name") or ident), "lon": float(spec["lon"]), "lat": float(spec["lat"]),
            "elevation": None if spec.get("elevation") is None else float(spec["elevation"])}


def _describe(name: str) -> dict | None:
    path = _folder(name) / "source.json"
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def _save(desc: dict) -> None:
    path = _folder(desc["name"]) / "source.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name("source.json.part")
    part.write_text(json.dumps(desc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    part.replace(path)


def define(name: str, title: str, variables, stations=None, *, credit: str | None = None,
           licence: str | None = None, domain: str | None = None) -> dict:
    """Creates or updates a source's description and returns it.

    variables: [{name, units?, long_name?}] or ['name:units[:long name]']; a CF
    standard name of KNOWN is kept in its SI unit, and `units`, if other, are
    those its readings arrive in. Variables already described and not named
    here are kept, as are stations; stations: [{id, name?, lon, lat, elevation?}].
    """
    with _lock:
        folder = _folder(name)
        desc = _describe(name) or {"name": name, "created": _iso(dt.datetime.now(UTC)), "variables": [],
                                   "stations": [], "count": 0, "last": None}
        given = [_variable(v) for v in variables]
        names = [v["name"] for v in given]
        if len(set(names)) != len(names):
            raise ValueError("a variable is named twice")
        desc["variables"] = [v for v in desc["variables"] if v["name"] not in names] + given
        if not desc["variables"]:
            raise ValueError("a source needs at least one variable")
        by_id = {s["id"]: s for s in desc["stations"]}
        for s in stations or []:
            s = _station(s)
            by_id[s["id"]] = s
        desc["stations"] = sorted(by_id.values(), key=lambda s: s["id"])
        desc["title"] = str(title or name)
        for key, value in (("credit", credit), ("licence", licence), ("domain", domain)):
            if value is not None:
                desc[key] = str(value)
        desc["updated"] = _iso(dt.datetime.now(UTC))
        folder.mkdir(parents=True, exist_ok=True)
        _save(desc)
        return desc


def sources() -> list[dict]:
    """Every source's description: name, title, variables, stations, count, last."""
    with _lock:
        root = _root()
        if not root.exists():
            return []
        out = []
        for folder in sorted(p for p in root.iterdir() if p.is_dir() and NAME.fullmatch(p.name)):
            desc = _describe(folder.name)
            if desc:
                out.append(desc)
        return out


# ---- readings -----------------------------------------------------------------

def validate(desc: dict, reading: dict, now: dt.datetime | None = None) -> tuple[dict | None, str | None]:
    """(the reading as kept, None) or (None, why it is refused)."""
    if not isinstance(reading, dict):
        return None, "a reading is an object"
    ident = str(reading.get("station", ""))
    if not STATION.fullmatch(ident):
        return None, f"station id {ident!r}: letters, digits, '_', '.', '-', at most 64"
    known = next((s for s in desc["stations"] if s["id"] == ident), None)
    lon = reading.get("lon", known and known["lon"])
    lat = reading.get("lat", known and known["lat"])
    elevation = reading.get("elevation", known and known["elevation"])
    why = _place(lon, lat, elevation)
    if why:
        return None, why if known or "lon" in reading else f"station {ident} has no place: give lon and lat"
    try:
        t = parse_time(reading.get("time"))
    except (TypeError, ValueError, OverflowError, OSError):
        return None, f"time {reading.get('time')!r}: ISO 8601, or milliseconds since 1970 (UTC)"
    now = now or dt.datetime.now(UTC)
    if t < EARLIEST:
        return None, f"time {_iso(t)} is before 1990: a clock not set, or seconds where milliseconds are expected"
    if t > now + AHEAD:
        return None, f"time {_iso(t)} is in the future"
    given = reading.get("values")
    if not isinstance(given, dict) or not given:
        return None, "no values"
    variables = {v["name"]: v for v in desc["variables"]}
    values = {}
    for name, x in given.items():
        if name not in variables:
            return None, f"{name} is not a variable of {desc['name']} ({', '.join(variables)})"
        if x is None:
            continue
        if not _finite(x):
            return None, f"{name}: {x!r} is not a finite number"
        v = variables[name]
        x = convert(float(x), v.get("input_units", v["units"]), v["units"])
        lo, hi = LIMITS.get(name, (-math.inf, math.inf))
        if not lo <= x <= hi:
            return None, f"{name}: {x:g} {v['units']} is outside the plausible {lo:g} to {hi:g}"
        values[name] = x
    if not values:
        return None, "every value is missing"
    kept = {"station": ident, "lon": float(lon), "lat": float(lat),
            "elevation": None if elevation is None else float(elevation), "time": _iso(t), "values": values,
            "received": _iso(now)}
    if reading.get("name"):
        kept["name"] = str(reading["name"])[:100]
    return kept, None


def append(name: str, readings: list[dict], rejected: list | None = None) -> int:
    """Checks the readings, appends the good ones, records new or moved stations,
    and tells the subscribers; returns the number accepted. The refused are
    listed in `rejected`, if given, as (index, reason)."""
    with _lock:
        desc = _describe(name)
        if desc is None:
            raise KeyError(f"no source {name!r}; define it first")
        now = dt.datetime.now(UTC)
        accepted = []
        for i, r in enumerate(readings):
            kept, why = validate(desc, r, now)
            if kept is None:
                if rejected is not None:
                    rejected.append((i, why))
                continue
            accepted.append(kept)
        if not accepted:
            return 0
        path = _folder(name) / "readings.jsonl"
        with path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write("".join(json.dumps(r, ensure_ascii=False, separators=(",", ":")) + "\n" for r in accepted))
        stations = {s["id"]: s for s in desc["stations"]}
        for r in accepted:
            s = stations.get(r["station"])
            if s is None or (s["lon"], s["lat"], s["elevation"]) != (r["lon"], r["lat"], r["elevation"]) \
                    or (r.get("name") and r["name"] != s["name"]):
                stations[r["station"]] = {"id": r["station"], "name": r.get("name") or (s and s["name"]) or r["station"],
                                          "lon": r["lon"], "lat": r["lat"], "elevation": r["elevation"]}
        desc["stations"] = sorted(stations.values(), key=lambda s: s["id"])
        desc["count"] = int(desc.get("count", 0)) + len(accepted)
        desc["last"] = max([r["time"] for r in accepted] + ([desc["last"]] if desc.get("last") else []))
        _save(desc)
        subscribers = list(_subscribers)
    for callback in subscribers:
        try:
            callback(name, accepted)
        except Exception as error:   # a subscriber's failure must not lose the readings
            print(f"store: a subscriber failed: {error}", file=sys.stderr)
    return len(accepted)


def _load(name: str) -> list[dict]:
    """Every reading of the source, read once and then only what was appended since."""
    path = _folder(name) / "readings.jsonl"
    if not path.exists():
        return []
    key = str(path)
    size = path.stat().st_size
    entry = _cache.get(key)
    if entry is None or size < entry["size"]:
        entry = _cache[key] = {"size": 0, "readings": []}
    if size > entry["size"]:
        with path.open("rb") as handle:
            handle.seek(entry["size"])
            data = handle.read(size - entry["size"])
        end = data.rfind(b"\n") + 1               # a line still being written waits for the next read
        for line in data[:end].splitlines():
            if line.strip():
                try:
                    entry["readings"].append(json.loads(line))
                except ValueError:
                    continue
        entry["size"] += end
    return entry["readings"]


def _bbox(bbox):
    if bbox is None:
        return None
    if isinstance(bbox, str):
        bbox = [float(x) for x in bbox.split(",")]
    west, south, east, north = (float(x) for x in bbox)
    return west, south, east, north


def _inside(lon: float, lat: float, box) -> bool:
    west, south, east, north = box
    if not south <= lat <= north:
        return False
    return west <= lon <= east if west <= east else (lon >= west or lon <= east)    # across the antimeridian


def read(name: str, since=None, until=None, bbox=None, limit=None) -> list[dict]:
    """The source's readings, in order of time: from `since` to `until` (ISO,
    milliseconds or datetime; inclusive), inside bbox (west, south, east, north
    in degrees, or 'w,s,e,n'; west > east crosses the antimeridian), and the
    last `limit` of them."""
    with _lock:
        if _describe(name) is None:
            raise KeyError(f"no source {name!r}")
        lo = _iso(parse_time(since)) if since is not None else None
        hi = _iso(parse_time(until)) if until is not None else None
        box = _bbox(bbox)
        out = [r for r in _load(name)
               if (lo is None or r["time"] >= lo) and (hi is None or r["time"] <= hi)
               and (box is None or _inside(r["lon"], r["lat"], box))]
    out.sort(key=lambda r: r["time"])
    if limit is not None:
        out = out[-int(limit):] if int(limit) > 0 else []
    return [dict(r) for r in out]


def subscribe(callback) -> callable:
    """Calls callback(name, readings) after each append, in this process; returns
    a function that ends the subscription."""
    with _lock:
        _subscribers.append(callback)

    def unsubscribe() -> None:
        with _lock:
            if callback in _subscribers:
                _subscribers.remove(callback)
    return unsubscribe


# ---- the product --------------------------------------------------------------

def build_product(name: str, hours: float = 48) -> Path:
    """Writes observations-sensor-<name>.json (schema observations/1), the last
    `hours` of readings before the newest (once per station and time, the last
    appended winning), and registers it in the index of products."""
    with _lock:
        desc = _describe(name)
        if desc is None:
            raise KeyError(f"no source {name!r}")
        readings = _load(name)
        last = max((r["time"] for r in readings), default=None)
        end = parse_time(last) if last else dt.datetime.now(UTC)
        start = _iso(end - dt.timedelta(hours=hours))
        names = {s["id"]: s["name"] for s in desc["stations"]}
        rows = [{"id": r["station"], "name": r.get("name") or names.get(r["station"], r["station"]),
                 "lon": r["lon"], "lat": r["lat"], "elevation": r.get("elevation"), "time": r["time"],
                 "values": r["values"]} for r in readings if r["time"] >= start]
        variables = [{**v, "level": v.get("level"), "digits": None} for v in desc["variables"]]   # as measured
        for v in variables:
            v.pop("input_units", None)
        sea = all(v["name"].startswith(("sea_", "water_")) for v in variables)
        source = {"dataset": f"sensor-{name}", "provider": f"{desc['title']}: our own sensors (observation/sensors)",
                  "credit": desc.get("credit") or desc["title"],
                  "licence": desc.get("licence") or "The sensors' owner's; positions are published with the product.",
                  "retrieved_at": _iso(dt.datetime.now(UTC))}
        return observation_set.write(f"sensor-{name}", desc["title"],
                                     desc.get("domain") or ("oceanography" if sea else "meteorology"),
                                     source, variables, rows, hours=hours, end=end, keep_previous=False)
