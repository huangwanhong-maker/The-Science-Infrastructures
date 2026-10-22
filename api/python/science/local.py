"""The library read on this computer, without a server: the products of the data
cache (data/cache/products/ and its index), the registry's manifests, the sensors'
store (observation/sensors/store.py) and the compute bridge to the domains'
JavaScript (bridge.py). `Local(root)` has the read methods of `Client`, and the
API's service (api/service/handler.py) answers its requests through it, so that
the two cannot drift apart.

Resources are named as the index of products names them: a field `gfs-forecast`,
an event set `usgs-earthquakes-week`, a track set `gdacs-month`. Each kind has a
method that lists its entries when called without a name, and returns one when
given its name: `fields()`, `events()`, `events("gcmt")`.

Field products need numpy, and are read with science/fields.py exactly as the
Earth Atlas's page reads them; everything else needs only the standard library.

Times are UTC. An instant is given as ISO 8601 (an instant without an offset is
UTC), as `now`, or relative to the present as `now+6h` or `now-1d`; a datetime is
taken as it is, naive ones as UTC.
"""

from __future__ import annotations

import collections
import datetime as dt
import importlib
import importlib.util
import json
import math
import os
import re
import sys
import threading
from pathlib import Path

UTC = dt.timezone.utc
RELATIVE = re.compile(r"^now(?:\s*([+-])\s*(\d+(?:\.\d+)?)\s*([smhd]))?$")
UNITS = {"s": 1, "m": 60, "h": 3600, "d": 86400}
FIELD_CACHE = 6           # field products kept decoded in memory, most recently used first
SLICE_CACHE = 16          # slices kept per field


class NotFound(LookupError):
    """No such resource: the service answers 404."""


class BadRequest(ValueError):
    """A request that cannot be answered as asked: the service answers 400."""


class Unavailable(RuntimeError):
    """What this computer cannot do (no numpy, no Node, no store): the service answers 503."""


class Unimplemented(Unavailable):
    """What is not written yet (the export of files, until data/formats/export.py exists): 501."""


class Busy(RuntimeError):
    """A run asked for while another runs: the service answers 409."""


def repository(start: Path | None = None) -> Path:
    """The repository holding this package, or the first directory upward from
    `start` (the working directory) that holds data/registry."""
    here = Path(__file__).resolve().parents[3]
    if (here / "data" / "registry").is_dir():
        return here
    for folder in [Path(start or Path.cwd()).resolve(), *Path(start or Path.cwd()).resolve().parents]:
        if (folder / "data" / "registry").is_dir():
            return folder
    raise NotFound("no repository of The Science Infrastructures here: pass its root to Local(root)")


def instant(value, now: dt.datetime | None = None) -> dt.datetime:
    """An instant from ISO 8601 text, `now`, `now±N[smhd]`, or a datetime; None is now."""
    now = now or dt.datetime.now(UTC)
    if value is None or value == "":
        return now
    if isinstance(value, dt.datetime):
        return value.replace(tzinfo=UTC) if value.tzinfo is None else value.astimezone(UTC)
    text = str(value).strip()
    m = RELATIVE.match(text)
    if m:
        if not m.group(1):
            return now
        seconds = float(m.group(2)) * UNITS[m.group(3)]
        return now + dt.timedelta(seconds=seconds if m.group(1) == "+" else -seconds)
    try:
        when = dt.datetime.fromisoformat(text.replace("Z", "+00:00").replace(" ", "T"))
    except ValueError:
        raise BadRequest(f"{text!r} is not an instant: give ISO 8601 (2026-10-23T06:00Z), now, or now+6h") from None
    return when.replace(tzinfo=UTC) if when.tzinfo is None else when.astimezone(UTC)


def iso(when: dt.datetime | None) -> str | None:
    return when.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ") if when else None


def bbox_of(value) -> tuple[float, float, float, float] | None:
    """A box (west, south, east, north) in degrees from 'w,s,e,n' or a sequence; west
    greater than east crosses the antimeridian."""
    if value is None or value == "":
        return None
    parts = value.split(",") if isinstance(value, str) else list(value)
    try:
        w, s, e, n = (float(p) for p in parts)
    except (TypeError, ValueError):
        raise BadRequest("a box is four numbers, west,south,east,north, in degrees") from None
    if not all(math.isfinite(x) for x in (w, s, e, n)) or not (-90 <= s <= n <= 90) or not (-360 <= w <= 360 and -360 <= e <= 360):
        raise BadRequest(f"the box {w},{s},{e},{n} is not west,south,east,north with south <= north, in degrees")
    return w, s, e, n


def in_box(lon: float, lat: float, box) -> bool:
    """Whether a place lies in a box (west, south, east, north); west > east crosses the antimeridian."""
    w, s, e, n = box
    if not (s <= lat <= n):
        return False
    span = 360.0 if e - w >= 360 else (e - w) % 360
    return (lon - w) % 360 <= span + 1e-9


def filter_rows(rows: list, fields: list, since=None, until=None, bbox=None, limit=None) -> list:
    """Rows of an event or observation set within a time span and a box, the latest
    `limit` of them. Rows are lists in the order of `fields` or objects; times are
    `time_ms` (ms UTC) or `time` (ISO 8601)."""
    def getter(name):
        if name in fields:
            i = fields.index(name)
            return lambda r: r[i] if isinstance(r, (list, tuple)) else r.get(name)
        return lambda r: r.get(name) if isinstance(r, dict) else None
    time_ms, time_iso, lon, lat = getter("time_ms"), getter("time"), getter("lon"), getter("lat")

    def ms(r):
        t = time_ms(r)
        if t is not None:
            return float(t)
        t = time_iso(r)
        if isinstance(t, (int, float)):
            return float(t)
        return instant(t).timestamp() * 1000 if t else None
    lo = instant(since).timestamp() * 1000 if since not in (None, "") else None
    hi = instant(until).timestamp() * 1000 if until not in (None, "") else None
    box = bbox_of(bbox)
    out = []
    for r in rows:
        if lo is not None or hi is not None:
            t = ms(r)
            if t is None or (lo is not None and t < lo) or (hi is not None and t > hi):
                continue
        if box:
            x, y = lon(r), lat(r)
            if x is None or y is None or not in_box(float(x), float(y), box):
                continue
        out.append(r)
    n = limit_of(limit)
    if n is not None and len(out) > n:
        timed = [(ms(r), k) for k, r in enumerate(out)]
        if all(t is not None for t, _ in timed):
            keep = sorted(k for _, k in sorted(timed)[-n:])
            out = [out[k] for k in keep]
        else:
            out = out[-n:]
    return out


def limit_of(value) -> int | None:
    if value is None or value == "":
        return None
    try:
        n = int(value)
    except (TypeError, ValueError):
        raise BadRequest(f"limit is a whole number, not {value!r}") from None
    if n < 0:
        raise BadRequest("limit is zero or more")
    return n


def load_store(where):
    """The sensors' store: a module (or object) with its interface, or the path of
    one; None when there is none. A module is imported under its own name (`store`,
    as ingest.py and simulate.py import it), so that every user in one process shares
    one store and one list of subscribers."""
    if where is None or hasattr(where, "read"):
        return where
    path = Path(where).resolve()
    if not path.is_file():
        return None
    loaded = sys.modules.get(path.stem)
    if loaded is not None and Path(getattr(loaded, "__file__", "") or "").resolve() == path:
        return loaded
    if loaded is None:
        sys.path.insert(0, str(path.parent))
        return importlib.import_module(path.stem)
    spec = importlib.util.spec_from_file_location(f"science_store_{path.stem}", path)     # another module has the name
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Local:
    """The library on this computer: the read methods of Client, from the cache."""

    def __init__(self, root: Path | str | None = None, *, cache: Path | str | None = None, store=None, compute_seconds: float = 10.0):
        self.root = Path(root).resolve() if root else repository()
        self.cache = Path(cache or os.environ.get("SCIENCE_CACHE") or self.root / "data" / "cache").resolve()
        self.products = self.cache / "products"
        self.plugin_products = Path(os.environ.get("SCIENCE_PLUGIN_OUTPUT") or self.cache / "plugins")
        self.registry = self.root / "data" / "registry"
        self._store_given = store
        self._store = None
        self._store_loaded = False
        self._index = (None, {})
        self._fields = collections.OrderedDict()     # name -> (path, mtime, LocalField)
        self._lock = threading.RLock()
        self._bridge = None
        self.compute_seconds = compute_seconds

    # ---- the index and the registry ----------------------------------------------------------

    def index(self) -> dict:
        """The index of products, read again when its file changes."""
        path = self.products / "index.json"
        try:
            mtime = path.stat().st_mtime
        except FileNotFoundError:
            return {}
        with self._lock:
            if self._index[0] != mtime:
                self._index = (mtime, json.loads(path.read_text(encoding="utf-8")))
            return self._index[1]

    def list(self, kind: str) -> list[dict]:
        """The index's entries of a kind: 'fields', 'events', 'tracks', 'observations', 'regions', ..."""
        return list(self.index().get(kind, []))

    def entry(self, kind: str, name: str) -> dict:
        found = next((e for e in self.list(kind) if e.get("name") == name), None)
        if not found:
            raise NotFound(f"no {kind} product {name!r}")
        return found

    def product_path(self, kind: str, name: str) -> Path:
        entry = self.entry(kind, name)
        path = self.products / entry["file"] if entry.get("file") else None
        if not path or not path.is_file():
            raise NotFound(f"the {kind} product {name!r} is listed, but its file is missing")
        return path

    def product(self, kind: str, name: str) -> dict:
        """A product's JSON as it is, by kind and name."""
        return json.loads(self.product_path(kind, name).read_text(encoding="utf-8"))

    def geojson(self, kind: str, name: str, *, since=None, until=None, bbox=None, limit=None) -> tuple[bytes, str, str]:
        """An event, track or observation set as a GeoJSON FeatureCollection, written by
        data/formats/export.py, its rows filtered first: (bytes, content type, file name)."""
        import tempfile
        exporter = export_module(self.root)
        path = self.product_path(kind, name)
        if kind == "tracks" or all(x in (None, "") for x in (since, until, bbox, limit)):
            return exporter.features(path, "geojson")
        if kind == "events":
            filtered = self.events(name, since=since, until=until, bbox=bbox, limit=limit)
        else:
            filtered = self.observations(name, since=since, until=until, bbox=bbox, limit=limit)
        with tempfile.TemporaryDirectory(prefix="science-geojson-") as folder:
            temporary = Path(folder) / path.name
            temporary.write_text(json.dumps(filtered, ensure_ascii=False), encoding="utf-8")
            return exporter.features(temporary, "geojson")

    def datasets(self) -> list[dict]:
        """The registry's manifests: what can be fetched, from whom, under which licence, with which parameters."""
        out = []
        for path in sorted(self.registry.glob("*.json")):
            try:
                m = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            out.append({
                "id": m.get("id", path.stem), "kind": m.get("kind", "files"), "title": m.get("title"),
                "domain": m.get("domain"), "model": m.get("model"), "provider": m.get("provider"), "credit": m.get("credit"),
                "licence": m.get("licence"), "citation": m.get("citation"),
                "parameters": {k: {f: v[f] for f in ("type", "default", "min", "max", "description") if f in v}
                               for k, v in (m.get("parameters") or {}).items()},
                "profiles": {k: (v.get("description") if isinstance(v, dict) else v) for k, v in (m.get("profiles") or {}).items()},
                "default_profile": m.get("default_profile"),
                "on_demand": (m.get("on_demand") or {}).get("parameters"),
            })
        return out

    def catalog(self) -> dict:
        """The registry's datasets and the index of products."""
        return {"datasets": self.datasets(), "products": self.index()}

    # ---- fields -------------------------------------------------------------------------------

    def _field_paths(self) -> dict[str, Path]:
        paths = {e["name"]: self.products / e["file"] for e in self.list("fields") if e.get("file")}
        for path in sorted(self.plugin_products.glob("*/fields-*.json")):
            paths.setdefault(path.stem[len("fields-"):], path)
        return paths

    def fields(self) -> list[dict]:
        """Every field product: the index's entries, and those plug-ins have made (computed: true)."""
        out = [{k: v for k, v in e.items() if k != "file"} for e in self.list("fields")]
        names = {e["name"] for e in out}
        for path in sorted(self.plugin_products.glob("*/fields-*.json")):
            try:
                meta = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if meta.get("name") in names:
                continue
            times = meta.get("times") or []
            out.append({"name": meta["name"], "title": meta.get("title"), "domain": meta.get("domain"), "model": meta.get("model"),
                        "computed": True,
                        "variables": [{k: v.get(k) for k in ("name", "long_name", "standard_name", "units", "level")} for v in meta["variables"]],
                        "vectors": meta.get("vectors", []),
                        "times": {"count": len(times), "first": times[0] if times else None, "last": times[-1] if times else None},
                        "grid": meta["grid"], "source": meta.get("source", {})})
        return out

    def field(self, name: str) -> "LocalField":
        """A field product by name, decoded on demand; read again when its file changes."""
        path = self._field_paths().get(name)
        if not path or not path.is_file():
            raise NotFound(f"no field {name!r}")
        mtime = path.stat().st_mtime
        with self._lock:
            hit = self._fields.get(name)
            if hit and hit[0] == path and hit[1] == mtime:
                self._fields.move_to_end(name)
                return hit[2]
        try:
            from . import fields as decoding
        except ImportError as error:
            raise Unavailable(f"reading fields needs numpy: {error}") from None
        made = LocalField(self, decoding.read(path, cache=SLICE_CACHE), path)
        with self._lock:
            self._fields[name] = (path, mtime, made)
            self._fields.move_to_end(name)
            while len(self._fields) > FIELD_CACHE:
                self._fields.popitem(last=False)
        return made

    # ---- events, tracks, regions, verification -----------------------------------------------

    def events(self, name: str | None = None, *, since=None, until=None, bbox=None, limit=None):
        """The event sets (no name), or one, its rows within a span, a box and a limit."""
        if name is None:
            return self.list("events")
        product = self.product("events", name)
        rows = product.get("events", [])
        if any(x not in (None, "") for x in (since, until, bbox, limit)):
            kept = filter_rows(rows, product.get("fields", []), since, until, bbox, limit)
            product = {**product, "events": kept, "count": len(kept), "of": len(rows)}
        return product

    def tracks(self, name: str | None = None):
        return self.list("tracks") if name is None else self.product("tracks", name)

    def regions(self, name: str | None = None):
        return self.list("regions") if name is None else self.product("regions", name)

    def regions_at(self, name: str, lon, lat) -> list[dict | None]:
        """The region of a set at each place: the cell containing it (rounded to the
        nearest cell centre, as the page's class rasters are read), None where there is none."""
        import numpy as np
        product = self.product("regions", name)
        g = product["grid"]
        lons, lats = places(lon, lat)
        raster = np.fromfile(self.products / product["file"], dtype="<u2")
        c = np.round((lons - g["lon0"]) / g["dlon"]).astype(np.int64)
        r = np.round((g["lat0"] - lats) / g["dlat"]).astype(np.int64)
        inside = (r >= 0) & (r < g["nlat"])
        if g.get("periodic"):
            c = np.mod(c, g["nlon"])
        else:
            inside &= (c >= 0) & (c < g["nlon"])
        ids = np.where(inside, raster[np.clip(r, 0, g["nlat"] - 1) * g["nlon"] + np.clip(c, 0, g["nlon"] - 1)], 0)
        regions = product.get("regions", [])
        return [regions[k - 1] if 0 < k <= len(regions) else None for k in ids.tolist()]

    def verification(self, name: str | None = None):
        return self.list("verification") if name is None else self.product("verification", name)

    # ---- observations: products and the sensors' store ------------------------------------------

    @property
    def store(self):
        """The sensors' store, observation/sensors/store.py, or None where there is none."""
        if not self._store_loaded:
            with self._lock:
                if not self._store_loaded:
                    where = self._store_given if self._store_given is not None else self.root / "observation" / "sensors" / "store.py"
                    self._store = load_store(where)
                    self._store_loaded = True
        return self._store

    def sensors(self) -> list[dict]:
        """The sources of the sensors' store: [{name, title, variables, ...}]."""
        store = self.store
        if store is None:
            return []
        return [s if isinstance(s, dict) else {"name": s} for s in store.sources()]

    def observations(self, source: str | None = None, *, since=None, until=None, bbox=None, limit=None):
        """The observation sets and sensor sources (no name), or one: a sensor source's
        readings, or an observations product's rows, within a span, a box and a limit."""
        if source is None:
            return {"products": self.list("observations"), "sensors": self.sensors()}
        if any(s.get("name") == source for s in self.sensors()):
            n = limit_of(limit)
            readings = self.store.read(source, since=iso(instant(since)) if since not in (None, "") else None,
                                       until=iso(instant(until)) if until not in (None, "") else None,
                                       bbox=bbox_of(bbox), limit=n)
            return {"source": source, "kind": "sensor", "readings": list(readings)}
        product = self.product("observations", source)
        key = next((k for k in ("rows", "readings", "observations") if isinstance(product.get(k), list)), None)
        if key and any(x not in (None, "") for x in (since, until, bbox, limit)):
            kept = filter_rows(product[key], product.get("fields", []), since, until, bbox, limit)
            product = {**product, key: kept, "count": len(kept), "of": len(product[key])}
        return product

    def define(self, source: str, title: str, variables, stations=None) -> dict:
        store = self._writable_store()
        store.define(source, title, variables, stations)
        return {"source": source, "defined": True}

    def append(self, source: str, readings: list, rejected: list | None = None) -> int:
        """Appends readings to a sensor source; the refused are listed in `rejected` as (index, reason)."""
        return int(self._writable_store().append(source, readings, rejected=rejected))

    def build(self, source: str, hours: float = 48):
        """Writes the source's last hours as the observations product sensor-<source>."""
        return self._writable_store().build_product(source, hours=hours)

    def _writable_store(self):
        if self.store is None:
            raise Unavailable("there is no sensors' store on this computer (observation/sensors/store.py)")
        return self.store

    # ---- the compute bridge -------------------------------------------------------------------

    def bridge(self):
        if self._bridge is None:
            with self._lock:
                if self._bridge is None:
                    from .bridge import Pool
                    self._bridge = Pool(self.root, seconds=self.compute_seconds)
        return self._bridge

    def units(self) -> list[dict]:
        """The domains and foundations whose functions can be called, with their functions."""
        return self.bridge().units()

    def compute(self, unit: str, function: str, *args):
        """An exported function of a domain or foundation, called with JSON arguments
        in the persistent Node worker; its result as JSON."""
        return self.bridge().call(unit, function, list(args))

    def close(self) -> None:
        if self._bridge is not None:
            self._bridge.close()


# ---- a field, its places, series and grids --------------------------------------------------------

Series = collections.namedtuple("Series", "times values units variable")
Grid = collections.namedtuple("Grid", "values grid time variable units")


class LocalField:
    """A field product read from the cache: its description, and its values at
    places, over time and on (parts of) its grid, as the page computes them."""

    def __init__(self, library: Local, product, path: Path):
        self.library, self.product, self.path = library, product, path
        self.name = product.name
        self.meta = describe(product.meta)
        self.variables = {v["name"]: v for v in self.meta["variables"]}
        self.times = list(product.times)

    def variable(self, var: str) -> dict:
        if not var:
            raise BadRequest(f"name a variable of {self.name}: {', '.join(self.variables)}")
        if var not in self.variables:
            raise NotFound(f"{self.name} has no variable {var!r}; it has {', '.join(self.variables)}")
        return self.variables[var]

    def interpolation(self, var: str) -> str:
        """'direction' for a direction, interpolated through its unit vector; else 'bilinear'."""
        return "direction" if decoding().is_direction(self.variable(var)) else "bilinear"

    def bracket(self, var: str, when: dt.datetime) -> dict:
        """The instant asked, the stored times around it and the weight of the later one."""
        i0, i1, w = self.product.bracket(when, var)
        if i0 < 0:
            return {"time": None, "times": [], "weight": 0.0}
        shown = when if w > 0 else self.times[i0]
        return {"time": iso(shown), "times": [iso(self.times[i0]), iso(self.times[i1])] if i1 != i0 else [iso(self.times[i0])],
                "weight": w}

    def sample(self, var: str, lon, lat, time=None):
        """Values at places (scalars, or sequences of the same length), bilinear in
        space and linear in time, NaN where there is none; a float for one place."""
        import numpy as np
        v = self.variable(var)
        lons, lats = places(lon, lat)
        when = instant(time)
        product = self.product
        values = product.sample_direction(var, when, lons, lats) if decoding().is_direction(v) else product.sample(var, when, lons, lats)
        values = np.asarray(values, dtype=np.float64)
        return float(values.reshape(-1)[0]) if np.ndim(lon) == 0 and np.ndim(lat) == 0 else values

    def series(self, var: str, lon: float, lat: float) -> Series:
        """The value at a place at every time of the product (one value, at no time, for a field without times)."""
        import numpy as np
        v = self.variable(var)
        lons, lats = places(lon, lat)
        if lons.size != 1:
            raise BadRequest("a series is at one place: give one lon and one lat")
        product = self.product
        if not self.times:
            return Series([None], np.array([float(self.sample(var, float(lons[0]), float(lats[0])))]), v.get("units"), var)
        if decoding().is_direction(v):
            values = np.array([float(product.sample_direction(var, t, lons, lats)[0]) for t in self.times])
        else:
            values = np.array([float(decoding().bilinear(product.slice(var, i), product.grid, lons, lats)[0]) for i in range(len(self.times))])
        return Series(list(self.times), values, v.get("units"), var)

    def grid(self, var: str, time=None, bbox=None, stride: int = 1) -> Grid:
        """A variable on its grid, or a box of it, every `stride`-th point: at a stored
        time when `time` is an index, or at an instant, between the two times around
        it. Rows run north to south; float32, NaN where there is none."""
        import numpy as np
        v = self.variable(var)
        stride = stride_of(stride)
        product = self.product
        index = index_of(time)
        if index is not None:
            if self.times and not (-len(self.times) <= index < len(self.times)):
                raise BadRequest(f"{self.name} has {len(self.times)} times; index {index} is beyond them")
            index = index % len(self.times) if self.times else 0
            values = product.slice(var, index)
            shown = iso(self.times[index]) if self.times else None
        else:
            when = instant(time)
            i0, i1, w = product.bracket(when, var)
            if i0 < 0:
                values, shown = product.slice(var, 0), None
            elif w == 0:
                values, shown = product.slice(var, i0), iso(self.times[i0])
            else:
                a, c = product.slice(var, i0).astype(np.float64), product.slice(var, i1).astype(np.float64)
                values, shown = (a + (c - a) * w).astype(np.float32), iso(when)
        rows, cols, out_grid = window(product.grid, bbox_of(bbox), stride)
        return Grid(np.ascontiguousarray(values[np.ix_(rows, cols)], dtype=np.float32), out_grid, shown, var, v.get("units"))

    def export(self, var: str, time=None, format: str = "tif", bbox=None, stride: int = 1, path=None):
        """The variable at a stored time (an index, or the time nearest an instant; 'all'
        for every time, in netCDF or JSON), written by data/formats/export.py as GeoTIFF
        (tif), netCDF (nc), CSV or JSON: (bytes, content type, file name), or the path
        written when `path` is given."""
        self.variable(var)
        fmt = (format or "tif").lower()
        if fmt not in ("tif", "nc", "csv", "json"):
            raise BadRequest(f"export formats are tif, nc, csv and json, not {format!r}")
        exporter = export_module(self.library.root)
        every = isinstance(time, str) and time.strip().lower() == "all"
        if every and fmt not in ("nc", "json"):
            raise BadRequest("every time ('all') goes into netCDF or JSON; a GeoTIFF or CSV holds one")
        index = index_of(time)
        if index is None and not every and self.times:
            index = self.product.nearest_index(instant(time))
        if index is not None and self.times and not (0 <= index < len(self.times)):
            raise BadRequest(f"{self.name} has {len(self.times)} times; index {index} is beyond them")
        box = bbox_of(bbox)
        data, content_type, filename = exporter.field(self.path, var, index if self.times else None, fmt,
                                                      bbox=box, stride=stride_of(stride))
        if path is None:
            return data, content_type, filename
        target = Path(path)
        if target.is_dir():
            target = target / filename
        target.write_bytes(data)
        return target


def index_of(time) -> int | None:
    """A time given as an index of the stored times (an integer, or its digits), or None."""
    if isinstance(time, int) and not isinstance(time, bool):
        return time
    if isinstance(time, str) and re.fullmatch(r"-?\d{1,6}", time.strip()):
        return int(time)
    return None


def decoding():
    from . import fields
    return fields


def describe(meta: dict) -> dict:
    """A field product's description without its byte layout: no file, encoding or
    offsets; each variable's quantum (the step of its codes: values are exact to half of it)."""
    out = {k: v for k, v in meta.items() if k not in ("file", "encoding", "variables", "_path")}
    out["variables"] = []
    for v in meta["variables"]:
        entry = {k: x for k, x in v.items() if k not in ("offsets", "scale", "offset")}
        if "scale" in v:
            entry["quantum"] = v["scale"]
        out["variables"].append(entry)
    return out


def places(lon, lat):
    """Longitudes and latitudes as float arrays of one length (a scalar is repeated)."""
    import numpy as np
    try:
        lons = np.atleast_1d(np.asarray(lon, dtype=np.float64)).reshape(-1)
        lats = np.atleast_1d(np.asarray(lat, dtype=np.float64)).reshape(-1)
    except (TypeError, ValueError):
        raise BadRequest("lon and lat are numbers, or lists of numbers") from None
    if lons.size == 0 or lats.size == 0:
        raise BadRequest("give at least one place: lon and lat")
    if lons.size != lats.size:
        if lons.size == 1:
            lons = np.repeat(lons, lats.size)
        elif lats.size == 1:
            lats = np.repeat(lats, lons.size)
        else:
            raise BadRequest(f"{lons.size} longitudes and {lats.size} latitudes: give as many of each, or one of either")
    if not (np.all(np.isfinite(lons)) and np.all(np.isfinite(lats))) or np.any(np.abs(lats) > 90):
        raise BadRequest("longitudes and latitudes are finite numbers of degrees, latitudes within ±90")
    return lons, lats


def stride_of(value) -> int:
    try:
        n = int(value if value not in (None, "") else 1)
    except (TypeError, ValueError):
        raise BadRequest(f"stride is a whole number, not {value!r}") from None
    if n < 1:
        raise BadRequest("stride is at least 1")
    return n


def window(grid: dict, box, stride: int):
    """The rows and columns of a grid whose centres lie in a box, every `stride`-th,
    and the grid they make. On a grid that wraps, a box from west to east may cross the
    antimeridian (west > east), and its columns run eastward from the west edge."""
    import numpy as np
    nlon, nlat = int(grid["nlon"]), int(grid["nlat"])
    lons = grid["lon0"] + np.arange(nlon) * grid["dlon"]
    lats = grid["lat0"] - np.arange(nlat) * grid["dlat"]
    periodic = bool(grid.get("periodic"))
    if box is None:
        rows, cols = np.arange(nlat), np.arange(nlon)
    else:
        w, s, e, n = box
        rows = np.nonzero((lats >= s - 1e-9) & (lats <= n + 1e-9))[0]
        if periodic:
            span = 360.0 if e - w >= 360 else (e - w) % 360
            offsets = np.mod(lons - w, 360)
            keep = np.nonzero(offsets <= span + 1e-9)[0]
            cols = keep[np.argsort(offsets[keep], kind="stable")]
        else:
            if w > e:
                raise BadRequest("this grid does not go round the Earth: a box's west is less than its east")
            cols = np.nonzero((lons >= w - 1e-9) & (lons <= e + 1e-9))[0]
    rows, cols = rows[::stride], cols[::stride]
    if rows.size == 0 or cols.size == 0:
        raise BadRequest("the box holds no point of the grid")
    lon0 = float(lons[cols[0]])
    lon0 = ((lon0 + 180) % 360) - 180 if periodic else lon0
    whole = box is None and nlon % stride == 0
    return rows, cols, {"lon0": lon0, "dlon": grid["dlon"] * stride, "nlon": int(cols.size),
                        "lat0": float(lats[rows[0]]), "dlat": grid["dlat"] * stride, "nlat": int(rows.size),
                        "periodic": bool(periodic and whole), "crs": grid.get("crs", "WGS 84 geographic coordinates (EPSG:4326), degrees"),
                        "order": "rows north to south, columns eastward"}


def export_module(root: Path):
    """data/formats/export.py, which writes fields as files; Unavailable until it exists."""
    path = Path(root) / "data" / "formats" / "export.py"
    if not path.is_file():
        raise Unimplemented("exporting needs data/formats/export.py, which is not on this computer yet")
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location("science_export", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
