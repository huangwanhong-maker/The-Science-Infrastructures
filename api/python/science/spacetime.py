"""Plug-ins over data in space and time (ADR 0079): scripts in Python and programs in C++ that are given records of what happened where and
when, and answer records of their own, which an application draws through time on its map, with tables, curves and notes in its panel. They
belong to no application: the Urban Science Workbench binds a `Host` to its own directory of scripts, as the Economy Atlas binds networks.py.
The same contract is written for JavaScript in the page (visualization/js/gis/spacetime.js) and for C++ (api/cpp/spacetime_plugin.hpp).

Space is degrees on WGS 84, longitude east and latitude north. Time is milliseconds since 1970-01-01 in UTC, and every record is held to that. A
record is a field on (a part of) the Earth times a line of time, in one of five forms, by what it is sampled on:

    events    points in space and time: a thing that happened at a place (an earthquake, a reading, a transaction)
              {"kind": "events", "title", "time": [ms], "lat": [], "lon": [], "values": {name: []}, "units": {name: unit}, "labels": []?}
    tracks    moving things, each a curve through space and time (a storm, a vessel, a walker)
              {"kind": "tracks", "title", "tracks": [{"id", "name"?, "time": [ms], "lat": [], "lon": [], "values": {name: []}}], "units"}
    cells     places with a series at each step of one line of time (a hexagon, a ward, a sensor)
              {"kind": "cells", "title", "time": [ms], "cells": [{"id", "ring": [[lon, lat], ...], "name"?}], "values": {name: [[value of each cell] at each time]}, "units"}
    grid      a raster at each step of one line of time, rows from the south
              {"kind": "grid", "title", "time": [ms], "box": [south, west, north, east], "width", "height", "values": {name: [[width * height values] at each time]}, "units"}
    flows     what passes between places
              {"kind": "flows", "title", "places": [{"id", "lat", "lon", "name"?}], "flows": [{"from", "to", "time", "until"?, "value"?}], "units"}

Each may carry `show`, how it asks to be drawn ({"mark", "value", "map", "scale", "window_ms"}), and `note`. A missing value is null or NaN.

A plug-in is given the records chosen, its parameters and the context of the view, and answers results: records, and

    table     rows in the panel (table())
    series    curves in the panel (series())
    note      text in the panel (note())

A Python script exports

    def create(api):
        class Mine(api.SpaceTimeAnalysis):
            id = "mine"; title = "..."; params = [api.param("days", "Days", 7, kind="integer")]
            def run(self, records, params, context):
                r = records[0]
                return [api.events("Late ones", r.time[keep], r.lat[keep], r.lon[keep])]
        return Mine()

where `api` is this module; a C++ program derives from science::spacetime::Analysis and ends with SCIENCE_SPACETIME_PLUGIN(Mine). Both run as processes
behind one file contract, networks.py's, so that a script can neither stop the server nor outlive its time: input.json, {"records", "params", "context"}, and
output.json, {"results"}.

    python -m science.spacetime describe <script>
    python -m science.spacetime run <script> <dir>
"""

from __future__ import annotations

import importlib.util
import json
import sys
import traceback
from pathlib import Path

import numpy as np

from .networks import Host as _ProcessHost
from .networks import note, param, plain, say_progress, table, typed

__all__ = ["Context", "Record", "SpaceTimeAnalysis", "Host", "events", "tracks", "cells", "grid", "flows", "series", "table", "note", "param", "bins_of"]

DAY_MS = 86400000.0
KINDS = ("events", "tracks", "cells", "grid", "flows")


def _numbers(values) -> np.ndarray:
    """A column as floats: a missing value (null) is NaN."""
    return np.array([np.nan if v is None else v for v in values], dtype=float) if isinstance(values, list) else np.asarray(values, dtype=float)


class Record:
    """A record as a plug-in reads it, its columns numpy arrays of floats, NaN where a value is missing. By its kind:
        events  time, lat, lon; values {name: array}; labels
        tracks  tracks, each with id, name, time, lat, lon, values
        cells   time; ids, rings, names; values {name: array (times, cells)}
        grid    time, box, width, height; values {name: array (times, height, width)}, rows from the south
        flows   places, flows"""

    def __init__(self, data: dict):
        self.raw = data
        self.kind = data.get("kind")
        if self.kind not in KINDS:
            raise ValueError(f"a record is one of {KINDS}, not {self.kind!r}")
        self.title = data.get("title", "")
        self.units = data.get("units") or {}
        given = data.get("values") or {}
        if self.kind == "events":
            self.time, self.lat, self.lon = _numbers(data["time"]), _numbers(data["lat"]), _numbers(data["lon"])
            self.values = {k: _numbers(v) for k, v in given.items()}
            self.labels = data.get("labels")
        elif self.kind == "tracks":
            self.tracks = [_Track(t) for t in data["tracks"]]
        elif self.kind == "cells":
            self.time = _numbers(data["time"])
            self.ids = [c["id"] for c in data["cells"]]
            self.rings = [c.get("ring") for c in data["cells"]]
            self.names = [c.get("name") for c in data["cells"]]
            self.values = {k: np.array([_numbers(row) for row in v]).reshape(len(self.time), len(self.ids)) for k, v in given.items()}
        elif self.kind == "grid":
            self.time, self.box = _numbers(data["time"]), data["box"]
            self.width, self.height = int(data["width"]), int(data["height"])
            self.values = {k: np.array([_numbers(f) for f in v]).reshape(len(self.time), self.height, self.width) for k, v in given.items()}
        else:
            self.places, self.flows = data["places"], data["flows"]

    def __repr__(self) -> str:
        return f"<{self.kind} {self.title!r}>"


class _Track:
    def __init__(self, t: dict):
        self.id, self.name = t.get("id"), t.get("name")
        self.time, self.lat, self.lon = _numbers(t["time"]), _numbers(t["lat"]), _numbers(t["lon"])
        self.values = {k: _numbers(v) for k, v in (t.get("values") or {}).items()}


class Context:
    """What a run knows of where it was asked: the view's point of regard and range, the time and its window on the application's scale,
    and the place open."""

    def __init__(self, data: dict | None = None, progress=None):
        data = data or {}
        self.lat, self.lon, self.range_m = data.get("lat"), data.get("lon"), data.get("range_m")
        self.now_ms, self.window_ms = data.get("now_ms"), data.get("window_ms")
        self.area = data.get("area")
        self._progress = progress

    def progress(self, fraction: float, text: str = "") -> None:
        if self._progress:
            self._progress(fraction, text)


class SpaceTimeAnalysis:
    """Subclass and write run(); set id, title, description and params (param(), as for a network's analysis)."""
    id, title, description = "", "", ""
    params: list = []

    def describe(self) -> dict:
        return {"id": self.id, "title": self.title, "description": self.description, "kind": "spacetime-analysis", "language": "python", "params": self.params}

    def run(self, records: list, params: dict, context: Context) -> list:
        raise NotImplementedError


# ---- results -------------------------------------------------------------------------------------------------------------------------


def _common(kind: str, title: str, units, show, note_text) -> dict:
    out = {"kind": kind, "title": title, "units": units or {}}
    if show:
        out["show"] = show
    if note_text:
        out["note"] = note_text
    return out


def events(title: str, time, lat, lon, values: dict | None = None, units: dict | None = None, labels=None, show: dict | None = None, note: str = "") -> dict:
    """Points in space and time. Columns of one length; `time` in ms UTC, degrees on WGS 84."""
    out = _common("events", title, units, show, note)
    out.update(time=time, lat=lat, lon=lon, values=values or {})
    if labels is not None:
        out["labels"] = labels
    return out


def tracks(title: str, items: list, units: dict | None = None, show: dict | None = None, note: str = "") -> dict:
    """Moving things: items of {id, name?, time, lat, lon, values?}."""
    out = _common("tracks", title, units, show, note)
    out["tracks"] = items
    return out


def cells(title: str, time, items: list, values: dict, units: dict | None = None, show: dict | None = None, note: str = "") -> dict:
    """Places with a series: items of {id, ring: [[lon, lat], ...], name?}; values {name: array (times, cells)}."""
    out = _common("cells", title, units, show, note)
    out.update(time=time, cells=items, values=values)
    return out


def grid(title: str, time, box, values: dict, units: dict | None = None, show: dict | None = None, note: str = "") -> dict:
    """A raster through time: box [south, west, north, east], values {name: array (times, height, width)}, rows from the south."""
    out = _common("grid", title, units, show, note)
    first = np.asarray(next(iter(values.values())))
    out.update(time=time, box=list(box), height=int(first.shape[1]), width=int(first.shape[2]), values={k: np.asarray(v, dtype=float).reshape(first.shape[0], -1) for k, v in values.items()})
    return out


def flows(title: str, places: list, items: list, units: dict | None = None, show: dict | None = None, note: str = "") -> dict:
    """What passes between places: places of {id, lat, lon, name?}; items of {from, to, time, until?, value?}."""
    out = _common("flows", title, units, show, note)
    out.update(places=places, flows=items)
    return out


def series(title: str, lines: list, xlabel: str = "", ylabel: str = "", time_axis: bool = False, note: str = "") -> dict:
    """Curves in the panel: lines of {name, x, y}; with time_axis, x is in ms UTC."""
    return {"kind": "series", "title": title, "lines": lines, "xlabel": xlabel, "ylabel": ylabel, "time_axis": time_axis, "note": note}


def bins_of(time, step_ms: float) -> tuple:
    """Each time's bin on a line of equal steps from the first, and the bins' starts: (index array, starts array)."""
    t = _numbers(time)
    known = t[np.isfinite(t)]
    if known.size == 0:
        return np.zeros(t.shape, dtype=int), np.array([])
    t0 = known.min()
    n = int((known.max() - t0) // step_ms) + 1
    index = np.where(np.isfinite(t), np.floor((np.nan_to_num(t, nan=t0) - t0) / step_ms), -1).astype(int)
    return index, t0 + step_ms * np.arange(n)


def finish(results) -> list:
    """The results of a run, checked: each of a known kind, its columns of one length, as JSON holds them."""
    out = []
    for r in results or []:
        if not isinstance(r, dict) or r.get("kind") not in (*KINDS, "table", "series", "note"):
            raise ValueError(f"a result is one of {(*KINDS, 'table', 'series', 'note')}: {str(r)[:120]}")
        if r["kind"] == "events":
            n = len(r["time"])
            if len(r["lat"]) != n or len(r["lon"]) != n or any(len(v) != n for v in r["values"].values()):
                raise ValueError(f"the columns of the events {r.get('title')!r} are not of one length")
        out.append(plain(r))
    return out


# ---- the runner ----------------------------------------------------------------------------------------------------------------------


def load(script: Path):
    script = Path(script)
    spec = importlib.util.spec_from_file_location(f"spacetime_plugin_{script.stem.replace('-', '_')}", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "create"):
        raise ValueError(f"{script.name} does not define create(api)")
    from science import spacetime as api
    plugin = module.create(api)
    if not hasattr(plugin, "run") or not hasattr(plugin, "describe"):
        raise ValueError(f"{script.name}: create(api) answers no plug-in (an instance of api.SpaceTimeAnalysis)")
    return plugin


def run_script(script: Path, work: Path) -> None:
    request = json.loads((Path(work) / "input.json").read_text(encoding="utf-8"))
    plugin = load(script)
    params = typed(plugin.describe().get("params", []), request.get("params"))
    records = [Record(r) for r in request.get("records", [])]
    results = finish(plugin.run(records, params, Context(request.get("context"), say_progress)))
    (Path(work) / "output.json").write_text(json.dumps({"results": results}, ensure_ascii=False, allow_nan=False), encoding="utf-8")


class Host(_ProcessHost):
    """The scripts of one directory, listed, described, compiled and run, as networks.py's Host, on records."""
    module = "science.spacetime"
    kind = "spacetime-analysis"

    def request_of(self, subject: dict, chosen: dict) -> dict:
        return {"records": plain(subject.get("records", [])), "params": chosen, "context": plain(subject.get("context") or {})}


def main(argv: list[str]) -> int:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    if len(argv) >= 2 and argv[0] == "describe":
        print(json.dumps(plain(load(Path(argv[1])).describe()), ensure_ascii=False))
        return 0
    if len(argv) >= 3 and argv[0] == "run":
        try:
            run_script(Path(argv[1]), Path(argv[2]))
        except Exception as error:
            traceback.print_exc()
            print(f"{type(error).__name__}: {error}", file=sys.stderr, flush=True)
            return 1
        return 0
    print(__doc__.split("\n\n")[-1], file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
