"""Field products read and written from Python, for plug-ins that run in the
server (ADR 0029). The format is that of data/products/field_grid.py, which the
page reads: a JSON description and a binary of int16 codes, little-endian, one
block per variable and time, row-major on a grid whose rows run north to south;
value = q * scale + offset, and q = -32768 means no value.
"""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import numpy as np

MISSING = -32768
UTC = dt.timezone.utc


def parse_time(text: str) -> dt.datetime:
    return dt.datetime.fromisoformat(text.replace("Z", "+00:00")).astimezone(UTC)


def iso(when: dt.datetime) -> str:
    return when.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def level_of(variable: dict) -> str:
    """A variable's level as the page reads it: stated, or inferred from its name ('2 m' for 'at 2m')."""
    import re
    if variable.get("level"):
        return variable["level"]
    m = re.search(r"\bat\s*(\d+(?:\.\d+)?)\s*m\b", variable.get("long_name") or "", re.I)
    return f"{m.group(1)} m" if m else ""


class Field:
    """A field product: its description, and its variables decoded on demand."""

    def __init__(self, meta: dict, codes: np.ndarray | None = None, slices: dict | None = None):
        self.meta = meta
        self.codes = codes          # the int16 codes of a product read from disk, or None
        self.slices = slices        # {name: [float32 arrays per time]} of a product made in memory
        self.grid = meta["grid"]
        self.times = [parse_time(t) for t in meta.get("times", [])]
        self.variables = {v["name"]: v for v in meta["variables"]}

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
        v, g = self.variables[name], self.grid
        n = g["nlat"] * g["nlon"]
        start = v["offsets"][max(0, min(len(v["offsets"]) - 1, index))] // 2
        q = self.codes[start:start + n].reshape(g["nlat"], g["nlon"])
        out = q.astype(np.float32) * np.float32(v["scale"]) + np.float32(v["offset"])
        out[q == MISSING] = np.nan
        return out

    def index_at_or_before(self, when: dt.datetime) -> int:
        """The index of the last time at or before `when`, or 0."""
        found = 0
        for i, t in enumerate(self.times):
            if t <= when + dt.timedelta(seconds=1):
                found = i
        return found

    def bracket(self, when: dt.datetime) -> tuple[int, int, float]:
        """The two time indices around `when` and the weight of the second, clamped to the span."""
        if len(self.times) < 2 or when <= self.times[0]:
            return 0, 0, 0.0
        if when >= self.times[-1]:
            return len(self.times) - 1, len(self.times) - 1, 0.0
        i = self.index_at_or_before(when)
        span = (self.times[i + 1] - self.times[i]).total_seconds()
        return i, i + 1, (when - self.times[i]).total_seconds() / span if span > 0 else 0.0

    def sample(self, name: str, when: dt.datetime, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
        """Values at places, bilinear in space and linear in time; NaN where the data are missing."""
        i0, i1, w = self.bracket(when)
        a = bilinear(self.slice(name, i0), self.grid, lons, lats)
        if i1 == i0 or w == 0:
            return a
        b = bilinear(self.slice(name, i1), self.grid, lons, lats)
        return a * (1 - w) + b * w


def bilinear(values: np.ndarray, grid: dict, lons, lats) -> np.ndarray:
    """Bilinear interpolation on the grid, as foundations/physis/js/grid.js does it:
    columns from lon0 by dlon (wrapping round when periodic), rows from lat0 southward by
    dlat; a place whose valid neighbours carry at least half the weight takes their mean."""
    lons, lats = np.asarray(lons, dtype=np.float64), np.asarray(lats, dtype=np.float64)
    nlon, nlat = grid["nlon"], grid["nlat"]
    x = (lons - grid["lon0"]) / grid["dlon"]
    if grid.get("periodic"):
        x = np.mod(x, nlon)
    y = (grid["lat0"] - lats) / grid["dlat"]
    x0 = np.floor(x).astype(int)
    y0 = np.clip(np.floor(y).astype(int), 0, nlat - 1)
    fx, fy = x - x0, np.clip(y - y0, 0, 1)
    if grid.get("periodic"):
        x1 = (x0 + 1) % nlon
        x0 = x0 % nlon
    else:
        x0 = np.clip(x0, 0, nlon - 1)
        x1 = np.clip(x0 + 1, 0, nlon - 1)
    y1 = np.clip(y0 + 1, 0, nlat - 1)
    corners = [(y0, x0, (1 - fx) * (1 - fy)), (y0, x1, fx * (1 - fy)), (y1, x0, (1 - fx) * fy), (y1, x1, fx * fy)]
    total = np.zeros_like(x, dtype=np.float64)
    weight = np.zeros_like(x, dtype=np.float64)
    for r, c, wt in corners:
        v = values[r, c].astype(np.float64)
        ok = np.isfinite(v)
        total += np.where(ok, v * wt, 0)
        weight += np.where(ok, wt, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(weight >= 0.5, total / weight, np.nan)
    return out


def read(path: Path) -> Field:
    """A field product from its JSON file; the codes come from the binary beside it."""
    path = Path(path)
    meta = json.loads(path.read_text(encoding="utf-8"))
    codes = np.fromfile(path.parent / meta["file"], dtype="<i2")
    return Field(meta, codes)


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
