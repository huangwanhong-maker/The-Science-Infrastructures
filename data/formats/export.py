"""Products of the library written in standard formats, for other software.

A field product (data/products/field_grid.py) becomes:

- ``tif``: a GeoTIFF of one time, float32 in EPSG:4326 (data/formats/tiff.py),
  which GDAL and QGIS read with its CRS, no-data value, name and units;
- ``nc``: a netCDF classic file (the 64-bit offset format of
  data/formats/netcdf3.py) following the CF Conventions 1.11: coordinates lat,
  lon and time ("seconds since 1970-01-01 00:00:00 UTC"), the variable packed as
  the product stores it (16-bit integers with scale_factor, add_offset and
  _FillValue, CF section 8.1), with its units, standard_name and long_name, and
  a coordinate of length one for its level (pressure or height) where it has
  one; latitude and longitude are on WGS 84, which GDAL assumes of them;
- ``csv``: one line per grid point of one time, "lon,lat,<name> (<units>)", an
  empty value where there is none;
- ``json``: the grid, the times and the values, [time][row][column], with null
  where there is no value.

A GeoTIFF and a CSV hold one time. Asked for every time (time_index None), they
take the product's only time if it has one, and refuse (ValueError) otherwise:
an implicit choice of instant would be a silent error. netCDF and JSON hold
every time.

The values are decoded as the page decodes them: value = q * scale + offset,
q = -32768 meaning no value (field_grid.py). A box (west, south, east, north,
degrees) keeps the grid points within it, across the antimeridian on a grid
round the Earth (west > east), with longitudes increasing through it; a stride
keeps every n-th row and column.

Observations, events (events-*.json) and tracks (tracks-*.json) become GeoJSON
(RFC 7946): a FeatureCollection in WGS 84, positions [longitude, latitude],
longitudes in [-180, 180]; events and observations as Points with their other
fields as properties, times also as ISO 8601 text, and rows without a valid
position left out and counted (the collection's "left_out"); each track as a
LineString of its positions, with their times and whether each was observed
in the properties, and cut into a MultiLineString where it crosses the
antimeridian (RFC 7946, section 3.1.9).

References: Unidata, "NetCDF Classic and 64-bit Offset File Formats"; CF
Conventions 1.11 (2023), sections 3 (units, standard_name), 4 (coordinates; 4.3,
a vertical coordinate of length one) and 8.1 (packed data); RFC 7946, "The GeoJSON
Format" (2016); RFC 4180 (CSV).
"""

from __future__ import annotations

import datetime as dt
import json
import math
import re
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import netcdf3  # noqa: E402
import tiff  # noqa: E402

MISSING = -32768
EPOCH = dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc)
TYPES = {"tif": "image/tiff; application=geotiff", "nc": "application/x-netcdf", "csv": "text/csv; charset=utf-8",
         "json": "application/json", "geojson": "application/geo+json"}


def _time(text: str) -> dt.datetime:
    return dt.datetime.fromisoformat(text.replace("Z", "+00:00")).astimezone(dt.timezone.utc)


def _stamp(text: str) -> str:
    """'2026-10-22T06:00:00Z' -> '20261022T0600Z', for file names."""
    return _time(text).strftime("%Y%m%dT%H%MZ")


def _safe(name: str) -> str:
    """A name of letters, digits and underscores, as CF recommends (section 2.3)."""
    return re.sub(r"[^A-Za-z0-9_]", "_", name)


class Selection:
    """Rows and columns of a product's grid kept by a box and a stride, and the grid they make."""

    def __init__(self, grid: dict, bbox: tuple | None, stride: int):
        if stride < 1:
            raise ValueError("stride must be at least 1")
        nlon, nlat = grid["nlon"], grid["nlat"]
        lons = grid["lon0"] + np.arange(nlon) * grid["dlon"]
        lats = grid["lat0"] - np.arange(nlat) * grid["dlat"]
        cols, rows = np.arange(nlon), np.arange(nlat)
        if bbox is not None:
            west, south, east, north = (float(v) for v in bbox)
            if south > north:
                raise ValueError("the box's south is north of its north")
            rows = rows[(lats >= south) & (lats <= north)]
            if grid.get("periodic"):
                span = (east - west) % 360 or (360.0 if east != west else 0.0)
                offset = (lons - west) % 360
                cols = cols[offset <= span + 1e-9]
                cols = cols[np.argsort(offset[cols], kind="stable")]       # eastward from the box's west
            else:
                if west > east:
                    raise ValueError("the box's west is east of its east on a grid that does not go round the Earth")
                cols = cols[(lons >= west) & (lons <= east)]
            if not len(rows) or not len(cols):
                raise ValueError(f"the box {bbox} holds no point of the grid")
        self.rows, self.cols = rows[::stride], cols[::stride]
        first = float(lons[self.cols[0]])
        if bbox is not None and grid.get("periodic"):
            first = float(bbox[0]) + (first - float(bbox[0])) % 360         # continuous through the antimeridian
        self.lons = first + (np.arange(len(self.cols)) * grid["dlon"] * stride if len(self.cols) > 1 else np.zeros(1))
        if len(self.cols) > 1 and not np.allclose(np.diff((lons[self.cols] - first) % 360), grid["dlon"] * stride):
            raise ValueError("the box makes an irregular selection of columns")
        self.lats = lats[self.rows]
        whole = bbox is None and stride == 1
        self.grid = {"lon0": float(self.lons[0]), "dlon": grid["dlon"] * stride, "nlon": int(len(self.cols)),
                     "lat0": float(self.lats[0]), "dlat": grid["dlat"] * stride, "nlat": int(len(self.rows)),
                     "periodic": bool(grid.get("periodic")) and whole,
                     "crs": "WGS 84 geographic coordinates (EPSG:4326), degrees", "order": "rows north to south, columns eastward"}

    def take(self, x: np.ndarray) -> np.ndarray:
        return x[np.ix_(self.rows, self.cols)]


class Product:
    """A field product on disk: its description and the codes of its variables."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.meta = json.loads(self.path.read_text(encoding="utf-8"))
        if self.meta.get("schema") != "field-grid/1":
            raise ValueError(f"{self.path.name} is not a field product")
        self.codes = np.fromfile(self.path.parent / self.meta["file"], dtype="<i2")
        self.grid = self.meta["grid"]
        self.times = self.meta.get("times", [])

    def variable(self, name: str) -> dict:
        for v in self.meta["variables"]:
            if v["name"] == name:
                return v
        raise ValueError(f"{self.meta['name']} has no variable {name!r}; it has {[v['name'] for v in self.meta['variables']]}")

    def codes_at(self, v: dict, index: int) -> np.ndarray:
        g = self.grid
        n = g["nlat"] * g["nlon"]
        if not 0 <= index < len(v["offsets"]):
            raise ValueError(f"time index {index} is outside 0..{len(v['offsets']) - 1}")
        start = v["offsets"][index] // 2
        return self.codes[start:start + n].reshape(g["nlat"], g["nlon"])

    def values_at(self, v: dict, index: int) -> np.ndarray:
        """Decoded values, float64, NaN where there is none."""
        q = self.codes_at(v, index)
        return np.where(q == MISSING, np.nan, q.astype(np.float64) * v["scale"] + v["offset"])

    def indices(self, v: dict, time_index: int | None, one: bool, fmt: str) -> list[int]:
        count = len(v["offsets"])
        if time_index is not None:
            if not -count <= time_index < count:
                raise ValueError(f"time index {time_index} is outside 0..{count - 1}")
            return [time_index % count]
        if one and count > 1:
            raise ValueError(f"a {fmt} file holds one time: give time_index, 0 to {count - 1} "
                             f"({self.times[0]} to {self.times[-1]})")
        return list(range(count))


def field(product_json: Path, variable: str, time_index: int | None, fmt: str, bbox: tuple | None = None,
          stride: int = 1) -> tuple[bytes, str, str]:
    """A variable of a field product in a standard format: (content, content type, file name).
    fmt: 'tif', 'nc', 'csv' or 'json'; time_index None means every time (see the module's text)."""
    if fmt not in ("tif", "nc", "csv", "json"):
        raise ValueError(f"unknown format {fmt!r}: tif, nc, csv or json")
    product = Product(product_json)
    v = product.variable(variable)
    indices = product.indices(v, time_index, fmt in ("tif", "csv"), {"tif": "GeoTIFF", "csv": "CSV"}.get(fmt, fmt))
    sel = Selection(product.grid, bbox, int(stride))
    times = [product.times[i] for i in indices] if product.times else []
    name = product.meta["name"]
    stem = f"{name}-{_safe(variable)}" + (f"-{_stamp(times[0])}" if len(times) == 1 else "")
    slices = [] if fmt == "nc" else [sel.take(product.values_at(v, i)) for i in indices]
    if fmt == "tif":
        metadata = {"product": name, "variable": variable}
        if times:
            metadata["time"] = times[0]
        if v.get("standard_name"):
            metadata["standard_name"] = v["standard_name"]
        credit = (product.meta.get("source") or {}).get("credit")
        if credit:
            metadata["credit"] = credit
        content = tiff.encode_grid(slices[0].astype(np.float32), sel.grid, name=v.get("long_name") or variable,
                                   units=v.get("units", ""), description=f"{product.meta.get('title', name)}: {v.get('long_name') or variable}",
                                   metadata=metadata)
    elif fmt == "nc":
        content = _netcdf(product, v, indices, sel)
    elif fmt == "csv":
        content = _csv(slices[0], sel, variable, v.get("units", ""))
    else:
        content = _json(product, v, times, slices, sel)
    return content, TYPES[fmt], f"{stem}.{fmt}"


def _number(x: float) -> float | None:
    return None if not math.isfinite(x) else float(f"{x:.7g}")


def _csv(values: np.ndarray, sel: Selection, variable: str, units: str) -> bytes:
    head = f"lon,lat,{variable}" + (f" ({units})" if units else "")
    lines = [head]
    for r, lat in enumerate(sel.lats):
        row = values[r]
        for c, lon in enumerate(sel.lons):
            x = row[c]
            lines.append(f"{lon:.6f},{lat:.6f},{'' if not math.isfinite(x) else format(x, '.7g')}")
    return ("\r\n".join(lines) + "\r\n").encode("utf-8")           # RFC 4180: lines end in CRLF


def _json(product: Product, v: dict, times: list[str], slices: list[np.ndarray], sel: Selection) -> bytes:
    described = {k: v.get(k) for k in ("name", "long_name", "standard_name", "units", "level", "cell_methods") if v.get(k) is not None}
    out = {
        "schema": "field-values/1",
        "product": product.meta["name"], "title": product.meta.get("title"),
        "variable": described,
        "grid": sel.grid,
        "times": times, "time_scale": "UTC",
        "order": "values[time][row][column]; rows north to south, columns eastward; null where there is no value",
        "values": [[[_number(x) for x in row] for row in s] for s in slices],
        "source": product.meta.get("source", {}),
    }
    return json.dumps(out, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _level_coordinate(level: str | None):
    """A CF scalar coordinate for a level: (name, value, attributes), or None."""
    if not level:
        return None
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*hPa\s*", level)
    if m:
        return "plev", np.float64(float(m.group(1)) * 100), {"standard_name": "air_pressure", "long_name": "pressure", "units": "Pa", "positive": "down", "axis": "Z"}
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*m\s*", level)
    if m:
        return "height", np.float64(float(m.group(1))), {"standard_name": "height", "long_name": "height above the surface", "units": "m", "positive": "up", "axis": "Z"}
    return None


def _netcdf(product: Product, v: dict, indices: list[int], sel: Selection) -> bytes:
    name = _safe(v["name"])
    codes = np.stack([sel.take(product.codes_at(v, i)) for i in indices]).astype(np.int16)
    attrs = {"long_name": v.get("long_name") or v["name"], "units": v.get("units", "1")}
    if v.get("standard_name"):
        attrs["standard_name"] = v["standard_name"]
    if v.get("cell_methods"):
        attrs["cell_methods"] = v["cell_methods"]
    attrs.update({"scale_factor": float(v["scale"]), "add_offset": float(v["offset"]), "_FillValue": np.int16(MISSING)})
    dimensions = {"lat": sel.grid["nlat"], "lon": sel.grid["nlon"]}
    variables = {}
    has_time = bool(product.times)
    if has_time:
        dimensions = {"time": len(indices), **dimensions}
        seconds = np.array([(_time(product.times[i]) - EPOCH).total_seconds() for i in indices], dtype=np.float64)
        variables["time"] = (["time"], seconds, {"standard_name": "time", "long_name": "time", "axis": "T",
                                                 "units": "seconds since 1970-01-01 00:00:00 UTC", "calendar": "standard"})
    variables["lat"] = (["lat"], sel.lats.astype(np.float64), {"standard_name": "latitude", "long_name": "latitude",
                                                               "units": "degrees_north", "axis": "Y", "comment": "WGS 84 (EPSG:4326)"})
    variables["lon"] = (["lon"], sel.lons.astype(np.float64), {"standard_name": "longitude", "long_name": "longitude",
                                                               "units": "degrees_east", "axis": "X", "comment": "WGS 84 (EPSG:4326)"})
    # A level is a coordinate of length one (CF section 4.3): netcdf3.write stores no scalar variables.
    level = _level_coordinate(v.get("level"))
    dims = ["time"] if has_time else []
    if level:
        dimensions = {**{k: n for k, n in dimensions.items() if k == "time"}, level[0]: 1,
                      "lat": dimensions["lat"], "lon": dimensions["lon"]}
        variables[level[0]] = ([level[0]], np.array([level[1]]), level[2])
        dims.append(level[0])
        codes = codes[:, None] if has_time else codes[:1]
    # The grid's datum, for software that reads the CRS (CF section 5.6): WGS 84 in longitude and latitude.
    variables["crs"] = ([], np.array(0, dtype=np.int32), {"grid_mapping_name": "latitude_longitude", "semi_major_axis": 6378137.0,
                                                          "inverse_flattening": 298.257223563, "longitude_of_prime_meridian": 0.0})
    attrs["grid_mapping"] = "crs"
    variables[name] = (dims + ["lat", "lon"], codes if (has_time or level) else codes[0], attrs)
    source = product.meta.get("source") or {}
    glob = {"Conventions": "CF-1.11", "title": product.meta.get("title", product.meta["name"]),
            "source": source.get("provider") or source.get("dataset") or "",
            "history": f"{dt.datetime.now(dt.timezone.utc):%Y-%m-%dT%H:%M:%SZ} written by data/formats/export.py "
                       f"from the product {product.meta['name']} of The Science Infrastructures"}
    for key, attr in (("credit", "acknowledgement"), ("licence", "license"), ("citation", "references")):
        if source.get(key):
            glob[attr] = str(source[key])
    if product.meta.get("climatology"):
        glob["comment"] = f"monthly normals of {product.meta['climatology']['period']}, placed in its middle year"
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / "export.nc"
        netcdf3.write(path, dimensions, variables, glob)
        return path.read_bytes()


# ---- features ---------------------------------------------------------------


def _iso_ms(ms: float) -> str:
    return (EPOCH + dt.timedelta(milliseconds=float(ms))).strftime("%Y-%m-%dT%H:%M:%SZ")


def _position(lon, lat) -> list[float] | None:
    """[lon, lat] in WGS 84, longitude in [-180, 180]; None for a missing or impossible
    position, such as the (-99.99, -99.99) some sources give for an unknown one."""
    if not isinstance(lon, (int, float)) or not isinstance(lat, (int, float)):
        return None
    if not (math.isfinite(lon) and math.isfinite(lat)) or abs(lat) > 90 or not -180 <= lon <= 360:
        return None
    return [float(lon - 360 if lon > 180 else lon), float(lat)]


def _point_features(meta: dict, rows_key: str) -> tuple[list[dict], int]:
    """Points with the other fields as properties, and the number of rows left out for want of a position."""
    fields = meta["fields"]
    if "lon" not in fields or "lat" not in fields:
        raise ValueError(f"{meta.get('name')}: no lon and lat among the fields {fields}")
    ilon, ilat = fields.index("lon"), fields.index("lat")
    features, left_out = [], 0
    for row in meta.get(rows_key, []):
        position = _position(row[ilon], row[ilat])
        if position is None:
            left_out += 1
            continue
        properties = {}
        for k, x in zip(fields, row):
            if k in ("lon", "lat"):
                continue
            if isinstance(x, float) and not math.isfinite(x):
                x = None
            properties[k] = x
            if k.endswith("_ms") and isinstance(x, (int, float)):
                properties[k[:-3] or "time"] = _iso_ms(x)
        feature = {"type": "Feature", "geometry": {"type": "Point", "coordinates": position}, "properties": properties}
        if "id" in properties and properties["id"] is not None:
            feature["id"] = properties["id"]
        features.append(feature)
    return features, left_out


def _line(coords: list[list[float]]) -> dict:
    """A LineString, or a MultiLineString cut at the antimeridian where a step
    crosses it (the shorter way round, |Δλ| > 180°)."""
    parts, part = [], [coords[0]]
    for (x0, y0), (x1, y1) in zip(coords, coords[1:]):
        if abs(x1 - x0) > 180:
            edge = 180.0 if x0 > 0 else -180.0
            x1u = x1 + 360 if x0 > 0 else x1 - 360                   # x1 unwrapped beside x0
            y = y0 + (y1 - y0) * (edge - x0) / (x1u - x0)
            part.append([edge, y])
            parts.append(part)
            part = [[-edge, y]]
        part.append([x1, y1])
    parts.append(part)
    if len(parts) == 1:
        return {"type": "LineString", "coordinates": parts[0]}
    return {"type": "MultiLineString", "coordinates": parts}


def _track_features(meta: dict) -> list[dict]:
    fields = meta["fields"]
    it, ilon, ilat = fields.index("time_ms"), fields.index("lon"), fields.index("lat")
    iobs = fields.index("observed") if "observed" in fields else None
    features = []
    for track in meta.get("tracks", []):
        kept = [(p, c) for p, c in ((p, _position(p[ilon], p[ilat])) for p in track.get("points", [])) if c]
        if not kept:
            continue
        points, coords = [p for p, _ in kept], [c for _, c in kept]
        geometry = _line(coords) if len(coords) > 1 else {"type": "Point", "coordinates": coords[0]}
        properties = {k: x for k, x in track.items() if k != "points"}
        properties["times"] = [_iso_ms(p[it]) for p in points]
        properties["times_ms"] = [p[it] for p in points]
        if iobs is not None:
            properties["observed"] = [p[iobs] for p in points]
        if "last_observed_ms" in track:
            properties["last_observed"] = _iso_ms(track["last_observed_ms"])
        features.append({"type": "Feature", "id": track.get("id"), "geometry": geometry, "properties": properties})
    return features


def features(product_json: Path, fmt: str = "geojson") -> tuple[bytes, str, str]:
    """Observations, events or tracks as a GeoJSON FeatureCollection: (content, content type, file name)."""
    if fmt != "geojson":
        raise ValueError(f"unknown format {fmt!r}: features are written as geojson")
    path = Path(product_json)
    meta = json.loads(path.read_text(encoding="utf-8"))
    schema = meta.get("schema", "")
    left_out = 0
    if schema.startswith("observations/"):
        found, left_out = _point_features(meta, "rows")
    elif schema.startswith("events/"):
        found, left_out = _point_features(meta, "events")
    elif schema.startswith("tracks/"):
        found = _track_features(meta)
    else:
        raise ValueError(f"{path.name}: schema {schema!r} is not observations, events or tracks")
    collection = {"type": "FeatureCollection", "name": meta.get("name"), "title": meta.get("title"),
                  "features": found}
    for key in ("variables", "units", "caveat", "window", "source"):
        if meta.get(key):
            collection[key] = meta[key]                    # foreign members (RFC 7946, section 6.1)
    if left_out:
        collection["left_out"] = {"rows": left_out, "why": "no valid WGS 84 position"}
    content = json.dumps(collection, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return content, TYPES["geojson"], f"{path.stem}.geojson"
