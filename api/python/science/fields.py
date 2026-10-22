"""Field products read and written from Python: for the plug-ins that run in a
server (ADR 0029) and for the API, which samples them as the Earth Atlas's page
does. The format is that of data/products/field_grid.py, which the page reads: a
JSON description and a binary of int16 codes, little-endian, one block per
variable and time, row-major on a grid whose rows run north to south;
value = q * scale + offset, and q = -32768 means no value.

Sampling follows the page exactly (applications/earth-atlas/src/data/fields.js
over foundations/physis/js/grid.js): codes are decoded in double precision and
rounded once to float32; a place is interpolated bilinearly between the four
grid points around it, the weights of those with a value renormalised when they
carry at least half the total, and NaN otherwise or outside the grid; an instant
is interpolated linearly between the two times around it, clamped to the span,
except for a variable whose cell method is a mean, sum, maximum or minimum over
time (the nearest time is taken) and for a climatology (the month is taken).
Directions are interpolated through their unit vectors. The accuracy is that of
the product: half a quantum (its `scale`) in each value, plus the error of
linear interpolation between grid points and times.
"""

from __future__ import annotations

import collections
import datetime as dt
import json
import math
import re
import threading
from pathlib import Path

import numpy as np

MISSING = -32768
UTC = dt.timezone.utc
NEAREST_IN_TIME = re.compile(r"time:\s*(mean|sum|maximum|minimum)")


def parse_time(text: str) -> dt.datetime:
    return dt.datetime.fromisoformat(text.replace("Z", "+00:00")).astimezone(UTC)


def iso(when: dt.datetime) -> str:
    return when.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def ms_of(when: dt.datetime) -> float:
    """Milliseconds since 1970-01-01T00:00:00Z, as the page counts time."""
    return when.timestamp() * 1000


def level_of(variable: dict) -> str:
    """A variable's level as the page reads it: stated, or inferred from its name ('2 m' for 'at 2m')."""
    if variable.get("level"):
        return variable["level"]
    m = re.search(r"\bat\s*(\d+(?:\.\d+)?)\s*m\b", variable.get("long_name") or "", re.I)
    return f"{m.group(1)} m" if m else ""


def is_direction(variable: dict) -> bool:
    """A direction in degrees, interpolated through its unit vector (the page's catalogue rule)."""
    return "direction" in (variable.get("standard_name") or "") or variable.get("units") == "degree"


def time_bracket(times_ms: list[float], ms: float, nearest: bool = False) -> tuple[int, int, float]:
    """The two time indices around an instant and the weight of the later one,
    clamped to the span; with `nearest`, the nearer of the two alone. (-1, -1, 0)
    for no times. As timeBracket in foundations/physis/js/grid.js."""
    n = len(times_ms)
    if n == 0:
        return -1, -1, 0.0
    if ms <= times_ms[0]:
        return 0, 0, 0.0
    if ms >= times_ms[-1]:
        return n - 1, n - 1, 0.0
    lo, hi = 0, n - 1
    while hi - lo > 1:
        mid = (lo + hi) >> 1
        if times_ms[mid] <= ms:
            lo = mid
        else:
            hi = mid
    w = (ms - times_ms[lo]) / (times_ms[hi] - times_ms[lo])
    if nearest:
        i = lo if w < 0.5 else hi
        return i, i, 0.0
    return lo, hi, w


class Field:
    """A field product: its description, and its variables decoded on demand."""

    def __init__(self, meta: dict, codes: np.ndarray | None = None, slices: dict | None = None, cache: int = 0):
        self.meta = meta
        self.codes = codes          # the int16 codes of a product read from disk, or None
        self.slices = slices        # {name: [float32 arrays per time]} of a product made in memory
        self.grid = meta["grid"]
        self.times = [parse_time(t) for t in meta.get("times", [])]
        self.times_ms = [ms_of(t) for t in self.times]
        self.variables = {v["name"]: v for v in meta["variables"]}
        # A climatology describes a typical year: its times are the months of one year,
        # and any instant is shown by its month.
        self.climatology = meta.get("climatology")
        self.cache_size = cache     # decoded slices kept, most recently used last (0: none)
        self._cache = collections.OrderedDict()
        self._lock = threading.Lock()

    @property
    def name(self) -> str:
        return self.meta["name"]

    def has(self, name: str) -> bool:
        return name in self.variables

    def slice(self, name: str, index: int = 0) -> np.ndarray:
        """The values of a variable at a time index: float32 (nlat, nlon), NaN where there is none."""
        if self.slices is not None:
            series = self.slices[name]
            return series[max(0, min(len(series) - 1, index))]
        if self.cache_size:
            with self._lock:
                hit = self._cache.get((name, index))
                if hit is not None:
                    self._cache.move_to_end((name, index))
                    return hit
        v, g = self.variables[name], self.grid
        n = g["nlat"] * g["nlon"]
        start = v["offsets"][max(0, min(len(v["offsets"]) - 1, index))] // 2
        q = self.codes[start:start + n].reshape(g["nlat"], g["nlon"])
        # In double precision, rounded once to float32, as the page decodes.
        out = (q.astype(np.float64) * v["scale"] + v["offset"]).astype(np.float32)
        out[q == MISSING] = np.nan
        if self.cache_size:
            out.flags.writeable = False
            with self._lock:
                self._cache[(name, index)] = out
                while len(self._cache) > self.cache_size:
                    self._cache.popitem(last=False)
        return out

    def index_at_or_before(self, when: dt.datetime) -> int:
        """The index of the last time at or before `when`, or 0."""
        found = 0
        for i, t in enumerate(self.times):
            if t <= when + dt.timedelta(seconds=1):
                found = i
        return found

    def nearest_index(self, when: dt.datetime) -> int:
        """The index of the time nearest to `when` (0 for a field without times)."""
        if not self.times:
            return 0
        ms = ms_of(when)
        return min(range(len(self.times)), key=lambda i: abs(self.times_ms[i] - ms))

    def bracket(self, when: dt.datetime, name: str | None = None) -> tuple[int, int, float]:
        """The two time indices around `when` and the weight of the second, as the page
        brackets an instant for the variable `name`; (-1, -1, 0) for a field without times."""
        if self.climatology and len(self.times) == 12:
            i = when.astimezone(UTC).month - 1
            return i, i, 0.0
        v = self.variables.get(name) if name else None
        nearest = bool(v and NEAREST_IN_TIME.search(v.get("cell_methods") or ""))
        return time_bracket(self.times_ms, ms_of(when), nearest)

    def sample(self, name: str, when: dt.datetime, lons, lats) -> np.ndarray:
        """Values at places, bilinear in space and linear in time; NaN where the data are missing."""
        i0, i1, w = self.bracket(when, name)
        if i0 < 0:
            return bilinear(self.slice(name, 0), self.grid, lons, lats)
        a = bilinear(self.slice(name, i0), self.grid, lons, lats)
        if w == 0:
            return a
        c = bilinear(self.slice(name, i1), self.grid, lons, lats)
        return a + (c - a) * w

    def sample_direction(self, name: str, when: dt.datetime, lons, lats) -> np.ndarray:
        """Directions (degrees) at places, interpolated through their unit vectors,
        since the mean of 350 and 10 degrees is 0 and not 180."""
        i0, i1, w = self.bracket(when, name)

        def at(i):
            d = self.slice(name, i).astype(np.float64) * math.pi / 180
            sn, cs = np.sin(d).astype(np.float32), np.cos(d).astype(np.float32)
            return bilinear(sn, self.grid, lons, lats), bilinear(cs, self.grid, lons, lats)
        x, y = at(max(0, i0))
        if w > 0:
            x1, y1 = at(i1)
            x, y = x + (x1 - x) * w, y + (y1 - y) * w
        with np.errstate(invalid="ignore"):
            out = np.fmod(np.arctan2(x, y) * 180 / math.pi + 360, 360)
        return np.where(np.isfinite(x) & np.isfinite(y) & ((x != 0) | (y != 0)), out, np.nan)

    def series(self, name: str, lon: float, lat: float) -> np.ndarray:
        """A variable's value at a place for each time of the product."""
        return np.array([float(bilinear(self.slice(name, i), self.grid, lon, lat)) for i in range(len(self.times))])


def bilinear(values: np.ndarray, grid: dict, lons, lats) -> np.ndarray:
    """Bilinear interpolation on the grid, exactly as foundations/physis/js/grid.js does
    it: columns from lon0 by dlon (wrapping round when periodic), rows from lat0 southward
    by dlat; a place within half a cell of the first or last row is clamped onto it, and
    one beyond, or beyond half a cell of the edge columns of a grid that does not wrap,
    has no value; a place whose valid neighbours carry at least half the weight takes
    their weighted mean. Returns float64 values, NaN where there is none, in the shape
    of the places."""
    lons, lats = np.broadcast_arrays(np.asarray(lons, dtype=np.float64), np.asarray(lats, dtype=np.float64))
    flat = np.asarray(values).reshape(-1)
    nlon, nlat = int(grid["nlon"]), int(grid["nlat"])
    periodic = bool(grid.get("periodic"))
    x = (lons - grid["lon0"]) / grid["dlon"]
    y = (grid["lat0"] - lats) / grid["dlat"]
    outside = ~np.isfinite(x) | ~np.isfinite(y) | (y < -0.5) | (y > nlat - 0.5)
    if periodic:
        x = np.fmod(np.fmod(x, nlon) + nlon, nlon)          # JavaScript's ((x % n) + n) % n
    else:
        outside |= (x < -0.5) | (x > nlon - 0.5)
    x, y = np.where(outside, 0, x), np.clip(np.where(outside, 0, y), 0, nlat - 1)
    c0, r0 = np.floor(x).astype(np.int64), np.floor(y).astype(np.int64)
    fx, fy = x - c0, y - r0
    c1 = c0 + 1
    r1 = np.minimum(r0 + 1, nlat - 1)
    c1 = np.where(c1 >= nlon, c1 - nlon if periodic else nlon - 1, c1)
    c0 = np.maximum(c0, 0)
    total = np.zeros(x.shape)
    weight = np.zeros(x.shape)
    for r, c, wt in ((r0, c0, (1 - fx) * (1 - fy)), (r0, c1, fx * (1 - fy)), (r1, c0, (1 - fx) * fy), (r1, c1, fx * fy)):
        v = flat[r * nlon + c].astype(np.float64)
        ok = ~np.isnan(v)
        total = total + np.where(ok, wt * np.where(ok, v, 0), 0)
        weight = weight + np.where(ok, wt, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where((weight >= 0.5) & ~outside, total / weight, np.nan)
    return out


def read(path: Path, cache: int = 0) -> Field:
    """A field product from its JSON file; the codes come from the binary beside it.
    `cache` decoded slices are kept for reuse (read-only arrays)."""
    path = Path(path)
    meta = json.loads(path.read_text(encoding="utf-8"))
    codes = np.fromfile(path.parent / meta["file"], dtype="<i2")
    return Field(meta, codes, cache=cache)


def quantise(x: np.ndarray) -> tuple[np.ndarray, float, float]:
    """int16 codes, scale and offset, as field_grid.py quantises."""
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return np.full(x.shape, MISSING, "<i2"), 1.0, 0.0
    lo, hi = float(finite.min()), float(finite.max())
    offset = (lo + hi) / 2
    scale = (hi - lo) / 65000 if hi > lo else 1.0
    q = np.where(np.isfinite(x), np.round((np.nan_to_num(x) - offset) / scale), MISSING)
    return q.astype("<i2"), scale, offset


def write(field: Field, out_dir: Path) -> Path:
    """Writes a field made in memory as a product, fields-<name>.json and .bin, and returns the JSON path."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    variables, blocks, offset_bytes = [], [], 0
    for v in field.meta["variables"]:
        x = np.stack([np.asarray(s, dtype=np.float32) for s in field.slices[v["name"]]])
        q, scale, offset = quantise(x)
        finite = x[np.isfinite(x)]
        offsets = []
        for t in range(q.shape[0]):
            offsets.append(offset_bytes)
            blocks.append(q[t].tobytes())
            offset_bytes += q[t].nbytes
        entry = dict(v)
        entry.update({
            "scale": scale, "offset": offset,
            "min": float(finite.min()) if finite.size else None,
            "max": float(finite.max()) if finite.size else None,
            "p02": float(np.percentile(finite, 2)) if finite.size else None,
            "p98": float(np.percentile(finite, 98)) if finite.size else None,
            "coverage": round(float(finite.size / x.size), 4),
            "offsets": offsets,
        })
        variables.append(entry)
    name = field.meta["name"]
    (out_dir / f"fields-{name}.bin").write_bytes(b"".join(blocks))
    product = {
        **{k: v for k, v in field.meta.items() if k != "variables"},
        "schema": "field-grid/1",
        "times": [iso(t) for t in field.times],
        "time_scale": "UTC",
        "variables": variables,
        "file": f"fields-{name}.bin",
        "encoding": "int16 little-endian; value = q * scale + offset; q = -32768 means no value",
    }
    json_path = out_dir / f"fields-{name}.json"
    json_path.write_text(json.dumps(product, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return json_path
