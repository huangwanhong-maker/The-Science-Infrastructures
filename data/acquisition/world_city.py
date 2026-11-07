"""Any city in the world, as it is asked for: its buildings, ways and open spaces from Overture Maps' world archives.

Overture Maps publishes each release as PMTiles archives of the whole Earth, read by range requests (formats/pmtiles.py):
buildings (from OpenStreetMap, Google Open Buildings, Microsoft and others, about 2.7 billion footprints, with heights or
floors where known), transportation (roads, paths, rails) and base (water and land use). A tile is asked for when the camera
comes near, reduced to what a city scene draws, and kept under data/cache/world-city/, so that no place need be prepared in
advance and a place once seen stays.

Where a building has no height, it has an estimate: the average height of the buildings on that 100 m cell in GHS-BUILT-H
R2023A (ANBH, 2018), a tile of 10 degrees fetched behind the scenes the first time a region is seen (`heights_of`). A
building keeps which of the three its height is, so that a guess can be drawn as one.

    city_tile(z, x, y)     buildings and ways of a tile at level 14: JSON text
    base_tile(z, x, y)     water and open land of a tile at level 13: JSON text
    heights_of(z, x, y)    the ground (always) and the buildings' average height (when its tile is here) on a grid over a tile

The tiles are of the Web Mercator pyramid; coordinates are longitude and latitude on WGS 84, degrees, to about 10 cm.
Overture's data are ODbL 1.0 for these themes; credit is given where they are drawn.
"""

from __future__ import annotations

import json
import math
import sys
import threading
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent
CACHE = DATA / "cache"
sys.path.insert(0, str(DATA / "formats"))
sys.path.insert(0, str(DATA))
import pmtiles  # noqa: E402
from tiff import GeoTIFF  # noqa: E402
import remotezip  # noqa: E402

RELEASE = "2026-10-21.0"
ARCHIVE = "https://tiles.overturemaps.org/{release}/{theme}.pmtiles"
CREDIT = "© OpenStreetMap contributors, Overture Maps Foundation; building footprints also from Google Open Buildings and Microsoft"
HEIGHT_CREDIT = "estimated heights: GHS-BUILT-H R2023A, European Commission JRC"
OUT = CACHE / "world-city"
_archives: dict[str, pmtiles.Archive] = {}
_lock = threading.Lock()

# A way's width where the source gives none (m): typical carriageways and paths, for drawing; a guess, said to be one.
WIDTHS = {"motorway": 24, "trunk": 20, "primary": 15, "secondary": 12, "tertiary": 10, "residential": 7, "unclassified": 7, "living_street": 5, "service": 4,
          "pedestrian": 5, "footway": 2, "path": 1.6, "steps": 2, "cycleway": 2.4, "track": 3, "bridleway": 2, "standard_gauge": 4, "subway": 3, "light_rail": 3, "tram": 3}
PEOPLE = {"pedestrian", "footway", "path", "steps", "living_street", "cycleway", "corridor"}
LAND_USE = {"park": "park", "recreation": "park", "cemetery": "park", "horticulture": "park", "protected": "park", "pedestrian": "plaza", "agriculture": None}
LAND_CLASS = {"forest": "park", "grass": "park", "wood": "park", "meadow": "park", "scrub": "park", "wetland": "park", "beach": "plaza"}


def archive(theme: str) -> pmtiles.Archive:
    with _lock:
        if theme not in _archives:
            _archives[theme] = pmtiles.Archive(ARCHIVE.format(release=RELEASE, theme=theme))
        return _archives[theme]


def _flat(points, convert) -> list[float]:
    out = []
    for px, py in points:
        lon, lat = convert(px, py)
        out += [round(lon, 6), round(lat, 6)]
    return out


def _polygons(feature) -> list[tuple[list, list]]:
    """A feature's polygons: (outer ring, holes), by the winding of its rings. Closing points are dropped."""
    polygons, sign = [], 0
    for part in feature["geometry"]:
        ring = part[:-1] if len(part) > 1 and part[0] == part[-1] else part
        if len(ring) < 3:
            continue
        area = pmtiles.signed_area(ring)
        if area == 0:
            continue
        s = 1 if area > 0 else -1
        if sign == 0:
            sign = s
        if s == sign or not polygons:
            polygons.append((ring, []))
        else:
            polygons[-1][1].append(ring)
    return polygons


def _number(v) -> float | None:
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) and v == v else None


def _rules_value(text) -> float | None:
    """A width from Overture's rules, when the first rule is plain."""
    try:
        rules = json.loads(text) if isinstance(text, str) else None
        v = rules[0].get("value") if rules else None
        return float(v) if v else None
    except (ValueError, TypeError, AttributeError, IndexError):
        return None


def _kept(path: Path, make) -> str:
    if path.exists():
        return path.read_text(encoding="utf-8")
    text = make()
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(f".{threading.get_ident()}.tmp")
    temp.write_text(text, encoding="utf-8")
    temp.replace(path)
    return text


def city_tile(z: int, x: int, y: int) -> str:
    def make() -> str:
        started = time.monotonic()
        buildings, ways = [], []
        raw = archive("buildings").tile(z, x, y)
        tile = pmtiles.read_tile(raw, {"building", "building_part"}) if raw else {}
        parts = tile.get("building_part", {"features": [], "extent": 4096})
        outlines = tile.get("building", {"features": [], "extent": 4096})
        with_parts = bool(parts["features"])
        for layer, is_part in ((outlines, False), (parts, True)):
            convert = pmtiles.to_lonlat(z, x, y, layer["extent"])
            for f in layer["features"]:
                p = f["properties"]
                if p.get("is_underground") is True:
                    continue
                # An outline whose parts are drawn stands only as its parts.
                if not is_part and p.get("has_parts") is True and with_parts:
                    continue
                height, floors, base = _number(p.get("height")), _number(p.get("num_floors")), _number(p.get("min_height"))
                if base is None and _number(p.get("min_floor")) is not None:
                    base = _number(p.get("min_floor")) * 3.2
                if height is not None and height > 0:
                    h, code = height, 0
                elif floors is not None and floors > 0:
                    h, code = floors * 3.2, 1
                else:
                    h, code = -1, 2
                name = p.get("@name") or None
                kind = p.get("class") or p.get("subtype") or None
                for ring, holes in _polygons(f):
                    b = [round(h, 1), code, round(base or 0, 1), _flat(ring, convert)]
                    if holes or name or kind:
                        b.append([_flat(q, convert) for q in holes])
                    if name or kind:
                        b += [name, kind]
                    buildings.append(b)
        raw = archive("transportation").tile(z, x, y)
        tile = pmtiles.read_tile(raw, {"segment"}) if raw else {}
        layer = tile.get("segment")
        if layer:
            convert = pmtiles.to_lonlat(z, x, y, layer["extent"])
            for f in layer["features"]:
                p = f["properties"]
                if f["type"] != 2 or p.get("subtype") == "water":
                    continue
                cls = str(p.get("class") or "unclassified")
                width = _rules_value(p.get("width_rules")) or WIDTHS.get(cls, 6)
                name = p.get("@name") or None
                for part in f["geometry"]:
                    if len(part) >= 2:
                        ways.append([cls, round(width, 1), name, _flat(part, convert), 1 if p.get("subtype") == "rail" else 0])
        lat0, lon0, lat1, lon1 = tile_box(z, x, y)
        return json.dumps({"tile": [z, x, y], "box": [lat0, lon0, lat1, lon1], "buildings": buildings, "ways": ways, "release": RELEASE, "credit": CREDIT,
                           "seconds": round(time.monotonic() - started, 2)}, ensure_ascii=False, separators=(",", ":"))
    return _kept(OUT / RELEASE / "city" / str(z) / str(x) / f"{y}.json", make)


def base_tile(z: int, x: int, y: int) -> str:
    def make() -> str:
        spaces = []
        raw = archive("base").tile(z, x, y)
        tile = pmtiles.read_tile(raw, {"water", "land_use", "land"}) if raw else {}
        for name, kind_of in (("land", lambda p: LAND_CLASS.get(p.get("class"))), ("land_use", lambda p: LAND_USE.get(p.get("subtype"), None) or LAND_CLASS.get(p.get("class"))), ("water", lambda p: "water")):
            layer = tile.get(name)
            if not layer:
                continue
            convert = pmtiles.to_lonlat(z, x, y, layer["extent"])
            for f in layer["features"]:
                if f["type"] != 3:
                    continue
                kind = kind_of(f["properties"])
                if not kind:
                    continue
                for ring, holes in _polygons(f):
                    spaces.append([kind, _flat(ring, convert), [_flat(q, convert) for q in holes]])
        return json.dumps({"tile": [z, x, y], "box": list(tile_box(z, x, y)), "spaces": spaces, "release": RELEASE, "credit": CREDIT}, ensure_ascii=False, separators=(",", ":"))
    return _kept(OUT / RELEASE / "base" / str(z) / str(x) / f"{y}.json", make)


def tile_box(z: int, x: int, y: int) -> tuple[float, float, float, float]:
    """A Web Mercator tile's south, west, north and east (degrees)."""
    n = 2 ** z
    lat = lambda t: math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * t / n))))  # noqa: E731
    return lat(y + 1), x / n * 360 - 180, lat(y), (x + 1) / n * 360 - 180


# ---- the buildings' average heights, and the ground -----------------------------------------------------------------

BUILT_H = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_H_GLOBE_R2023A/GHS_BUILT_H_ANBH_E2018_GLOBE_R2023A_4326_3ss/V1-0/tiles/GHS_BUILT_H_ANBH_E2018_GLOBE_R2023A_4326_3ss_V1_0_R{r}_C{c}.zip"
HEIGHTS = CACHE / "ghsl-built-h"
_fetching: set[tuple[int, int]] = set()
_absent: set[tuple[int, int]] = set()
_open: dict[tuple[int, int], GeoTIFF] = {}


def _height_tile_of(lat: float, lon: float) -> tuple[int, int]:
    # The release's rows and columns of 10 degrees, from its grid's corner at -180.0079, 89.0996.
    return int((89.0996 - lat) // 10) + 1, int((lon + 180.0079) // 10) + 1


def _fetch_height_tile(rc: tuple[int, int]) -> None:
    r, c = rc
    url = BUILT_H.format(r=r, c=c)
    try:
        members = [m.filename for m in remotezip.members(url) if m.filename.endswith(".tif")]
        data = remotezip.read_members(url, members)[members[0]]
        HEIGHTS.mkdir(parents=True, exist_ok=True)
        temp = HEIGHTS / f"R{r}_C{c}.tmp"
        temp.write_bytes(data)
        temp.replace(HEIGHTS / f"R{r}_C{c}.tif")
    except Exception as error:  # the sea and the empty parts of the world have no tile; a failure leaves the buildings at floors or the default
        _absent.add(rc)
        print(f"GHS-BUILT-H R{r}_C{c}: {error}", file=sys.stderr, flush=True)
    finally:
        _fetching.discard(rc)


def _height_tiff(lat: float, lon: float) -> tuple[GeoTIFF | None, bool]:
    """The tile of heights holding a point, if it is here; and whether one is on its way."""
    rc = _height_tile_of(lat, lon)
    if rc in _open:
        return _open[rc], False
    path = HEIGHTS / f"R{rc[0]}_C{rc[1]}.tif"
    if path.exists():
        with _lock:
            if len(_open) > 6:
                _open.clear()
            _open[rc] = GeoTIFF.open(path)
        return _open[rc], False
    if rc in _absent:
        return None, False
    with _lock:
        if rc not in _fetching:
            _fetching.add(rc)
            threading.Thread(target=_fetch_height_tile, args=(rc,), daemon=True).start()
    return None, True


def heights_of(z: int, x: int, y: int, terrain, tile_key, decode, n: int = 25) -> dict:
    """On a grid of n x n over a tile (rows north to south): the ground's height above WGS 84 (m), and the buildings' average height (m, 0
    where none stand), or null while its tile is on its way. `terrain`, `tile_key` and `decode` are the data library's tiles.terrain."""
    s, w, north, e = tile_box(z, x, y)
    lats = [north + (s - north) * j / (n - 1) for j in range(n)]
    lons = [w + (e - w) * i / (n - 1) for i in range(n)]
    # The ground, from the geographic pyramid's tile at a level whose cells are finer than the grid's.
    level = min(12, max(2, z - 2))
    span = 180 / 2 ** level
    tiles: dict[tuple[int, int], tuple] = {}
    ground = []
    for lat in lats:
        for lon in lons:
            gx, gy = min(2 ** (level + 1) - 1, int((lon + 180) // span)), min(2 ** level - 1, int((90 - lat) // span))
            if (gx, gy) not in tiles:
                header, h, _, _ = decode(terrain.tile("earth", tile_key(level, gx, gy)))
                tiles[(gx, gy)] = (header, h)
            header, h = tiles[(gx, gy)]
            cells, b = header["n"] - 1, header.get("border", 0)
            u, v = (lon + 180 - gx * span) / span * cells, (90 - lat - gy * span) / span * cells
            i0, j0 = int(min(cells - 1, u)), int(min(cells - 1, v))
            fu, fv = u - i0, v - j0
            at = lambda jj, ii: float(h[jj + b, ii + b])  # noqa: E731
            ground.append(round((1 - fv) * ((1 - fu) * at(j0, i0) + fu * at(j0, i0 + 1)) + fv * ((1 - fu) * at(j0 + 1, i0) + fu * at(j0 + 1, i0 + 1)), 1))
    # The buildings' heights.
    tiff, pending = _height_tiff((s + north) / 2, (w + e) / 2)
    built = None
    if tiff is not None:
        try:
            a, (x0, dx, y0, dy) = tiff.read_box(w, s, e, north)
            built = []
            for lat in lats:
                for lon in lons:
                    c = min(a.shape[1] - 1, max(0, int((lon - x0) / dx)))
                    r = min(a.shape[0] - 1, max(0, int((y0 - lat) / dy)))
                    q = a[max(0, r - 1):r + 2, max(0, c - 1):c + 2]
                    v = q[q > 0]
                    built.append(round(float(v.mean()), 1) if v.size else 0)
        except Exception:
            built = None
    return {"n": n, "ground": ground, "built": built, "pending": pending, "credit": HEIGHT_CREDIT}
