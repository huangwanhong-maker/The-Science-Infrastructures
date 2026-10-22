"""A stand-in for the sensors' store (observation/sensors/store.py), with its
interface, held in memory: for the API's checks on a computer where the store is
missing. Where it exists, the checks use it, in a directory of their own.

    sources()                                          [{name, title, variables, stations, count, last}]
    define(name, title, variables, stations=None)      variables [{name, units?, long_name?}] or ['name:units']
    append(name, readings, rejected=None) -> int       readings [{station, lon, lat, time, values: {variable: value}}]
    read(name, since=None, until=None, bbox=None, limit=None) -> [reading]
    subscribe(callback) -> unsubscribe                 callback(name, readings) after each append
    build_product(name, hours=48)                      not written here: the checks never build products
"""

from __future__ import annotations

import datetime as dt
import math
import threading

UTC = dt.timezone.utc
_lock = threading.Lock()
_sources: dict[str, dict] = {}
_readings: dict[str, list[dict]] = {}
_listeners: list = []


def _instant(value) -> dt.datetime:
    if isinstance(value, (int, float)):
        return dt.datetime.fromtimestamp(value / 1000, UTC)
    when = dt.datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    return when.replace(tzinfo=UTC) if when.tzinfo is None else when.astimezone(UTC)


def _iso(when: dt.datetime) -> str:
    return when.strftime("%Y-%m-%dT%H:%M:%SZ")


def sources() -> list[dict]:
    with _lock:
        return [dict(s) for s in _sources.values()]


def define(name: str, title: str, variables, stations=None, **_) -> dict:
    if not name or name != name.lower() or " " in name:
        raise ValueError(f"source name {name!r}: lowercase letters, digits and single hyphens")
    given = [{"name": v.split(":")[0], "units": (v.split(":") + [""])[1]} if isinstance(v, str) else dict(v) for v in variables]
    with _lock:
        _sources[name] = {"name": name, "title": title, "variables": given, "stations": list(stations or []), "count": 0, "last": None}
        _readings.setdefault(name, [])
        return dict(_sources[name])


def append(name: str, readings: list[dict], rejected: list | None = None) -> int:
    with _lock:
        if name not in _sources:
            raise KeyError(f"no source {name!r}; define it first")
        names = {v["name"] for v in _sources[name]["variables"]}
    kept = []
    for i, r in enumerate(readings):
        try:
            lon, lat, when = float(r["lon"]), float(r["lat"]), _instant(r["time"])
            values = {k: float(x) for k, x in r["values"].items() if x is not None}
            ok = -180 <= lon <= 180 and -90 <= lat <= 90 and values and set(values) <= names and all(math.isfinite(x) for x in values.values())
        except (KeyError, TypeError, ValueError, AttributeError):
            ok = False
        if not ok:
            if rejected is not None:
                rejected.append((i, "not a reading of this source"))
            continue
        kept.append({"station": str(r.get("station", "")), "lon": lon, "lat": lat, "time": _iso(when), "values": values})
    if not kept:
        return 0
    with _lock:
        _readings[name].extend(kept)
        _sources[name]["count"] += len(kept)
        _sources[name]["last"] = max(r["time"] for r in _readings[name])
        listeners = list(_listeners)
    for callback in listeners:
        callback(name, kept)
    return len(kept)


def read(name: str, since=None, until=None, bbox=None, limit=None) -> list[dict]:
    with _lock:
        if name not in _sources:
            raise KeyError(f"no source {name!r}")
        rows = list(_readings[name])
    if since is not None:
        rows = [r for r in rows if _instant(r["time"]) >= _instant(since)]
    if until is not None:
        rows = [r for r in rows if _instant(r["time"]) <= _instant(until)]
    if bbox:
        w, s, e, n = bbox
        rows = [r for r in rows if s <= r["lat"] <= n and (w <= r["lon"] <= e if w <= e else (r["lon"] >= w or r["lon"] <= e))]
    rows.sort(key=lambda r: r["time"])
    return rows[-int(limit):] if limit else rows


def subscribe(callback):
    with _lock:
        _listeners.append(callback)

    def unsubscribe():
        with _lock:
            if callback in _listeners:
                _listeners.remove(callback)
    return unsubscribe


def build_product(name: str, hours: float = 48):
    raise RuntimeError("the stand-in store builds no products")
