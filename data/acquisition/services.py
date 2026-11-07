"""Services of the world's development data, asked on demand and kept in the cache.

Some sources are services, not files: one asks about an area or a name. These are the ones the Urban Science Workbench
draws on; each request's answer is kept under data/cache/<service>/ by a hash of what was asked, with what was asked and
when beside it, so that the same question is answered from the cache and a page can move about freely. Nothing here is
reachable from outside the computer; the application's server calls it.

    space2stats_tile(lat, lon)      the H3 level 6 cells of a 4 x 4 degree tile: population, flood exposure, lights, built area, urbanisation, hazards
    space2stats_fields()            the service's list of fields
    thinkhazard_search(name)        the administrative units of a name (ThinkHazard!)
    thinkhazard_report(code)        each hazard's level in a unit
    data360_search(text)            indicators by what they mean (the service data360-mcp offers to language models)
    data360_values(database, indicator, year)
                                    an indicator for every economy, the latest at or before a year

    raster_window(kind, lat, lon, rp)
                                    a raster about any place, as bytes the page can lay on the ground and the engine can sample:
                                    settlement (WSF Evolution), urbanisation (GHS-SMOD), flood (river depth by return period, GloFAS)

The registry's manifests of the same names say what each is, its source and its terms.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.parse
import urllib.request
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "formats"))
from tiff import GeoTIFF, TIFFError  # noqa: E402

CACHE = Path(__file__).resolve().parents[1] / "cache"
AGENT = "noetarium-data/1 (+https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure)"
SPACE2STATS = "https://space2stats.ds.io"
THINKHAZARD = "https://thinkhazard.org/en"
DATA360 = "https://data360api.worldbank.org/data360"

# The fields of Space2Stats a cell is asked for; the age and sex of 2025 are asked for separately, at a point.
YEARS_LIGHTS = list(range(2012, 2025))
YEARS_BUILT = list(range(1975, 2031, 5))
YEARS_POP = list(range(2015, 2031))
CLASSES = [11, 12, 13, 21, 22, 23, 30]
TILE_FIELDS = (["pop", "pop_flood", "pop_flood_pct", "fires_density_mean", "cy_frequency_mean", "landslide_susceptibility_mean_2023", "drought_spei_1_5_rp100_mean", "ghs_total_pop"]
               + [f"ghs_{c}_pop" for c in CLASSES] + [f"sum_viirs_ntl_{y}" for y in YEARS_LIGHTS] + [f"sum_built_area_m_{y}" for y in YEARS_BUILT] + [f"sum_pop_{y}" for y in YEARS_POP])
TILE_DEGREES = 4
TILE_BATCH = 16


def _kept(service: str, key: dict, ask):
    """The answer to a question: from the cache if asked before, else asked and kept with its question."""
    folder = CACHE / service
    folder.mkdir(parents=True, exist_ok=True)
    name = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:20]
    path = folder / f"{name}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))["answer"]
    answer = ask()
    path.write_text(json.dumps({"asked": key, "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "answer": answer}, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    return answer


def _get(url: str, timeout: int = 60):
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": AGENT, "Accept": "application/json"}), timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def _post(url: str, body: dict, timeout: int = 120):
    request = urllib.request.Request(url, data=json.dumps(body).encode("utf-8"), headers={"User-Agent": AGENT, "Accept": "application/json", "Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))


def space2stats_fields() -> list:
    return _kept("space2stats", {"fields": True}, lambda: _get(f"{SPACE2STATS}/fields"))


def tile_of(lat: float, lon: float) -> tuple[int, int]:
    """The tile holding a point: its south-west corner, on a grid of TILE_DEGREES."""
    return int(lat // TILE_DEGREES) * TILE_DEGREES, int(lon // TILE_DEGREES) * TILE_DEGREES


def space2stats_tile(lat: float, lon: float) -> dict:
    """A tile of cells: {cells: [hex], rings: [[[lon, lat], ...]], fields: {name: [value per cell]}, tile: [south, west, size]}.
    A cell belongs to the tile holding its centre, so that tiles never repeat a cell."""
    s, w = tile_of(lat, lon)
    e, n = min(180, w + TILE_DEGREES), min(85, s + TILE_DEGREES)

    def ask() -> dict:
        aoi = {"type": "Feature", "properties": {}, "geometry": {"type": "Polygon", "coordinates": [[[w, s], [e, s], [e, n], [w, n], [w, s]]]}}
        # The service fails on an answer past its size, about 6 MB: every field of a 4 degree tile comes in batches, joined by cell.
        out = {"tile": [s, w, TILE_DEGREES], "cells": [], "rings": [], "fields": {f: [] for f in TILE_FIELDS}}
        at: dict[str, int] = {}
        rows = _post(f"{SPACE2STATS}/summary", {"aoi": aoi, "spatial_join_method": "centroid", "fields": TILE_FIELDS[:1], "geometry": "polygon"})
        for row in rows:
            ring = json.loads(row["geometry"])["coordinates"][0][:-1]
            at[row["hex_id"]] = len(out["cells"])
            out["cells"].append(row["hex_id"])
            out["rings"].append([[round(p[0], 5), round(p[1], 5)] for p in ring])
        for f in TILE_FIELDS:
            out["fields"][f] = [None] * len(out["cells"])
        for k in range(0, len(TILE_FIELDS), TILE_BATCH):
            batch = TILE_FIELDS[k:k + TILE_BATCH]
            for row in _post(f"{SPACE2STATS}/summary", {"aoi": aoi, "spatial_join_method": "centroid", "fields": batch}):
                c = at.get(row["hex_id"])
                if c is None:
                    continue
                for f in batch:
                    v = row.get(f)
                    out["fields"][f][c] = None if v is None else round(float(v), 3)
        return out
    return _kept("space2stats", {"tile": [s, w], "fields": TILE_FIELDS}, ask)


def thinkhazard_search(name: str) -> list:
    q = urllib.parse.quote(name)
    return _kept("thinkhazard", {"search": name}, lambda: _get(f"{THINKHAZARD}/administrativedivision?q={q}")["data"])


def thinkhazard_report(code: int) -> list:
    return _kept("thinkhazard", {"report": int(code)}, lambda: _get(f"{THINKHAZARD}/report/{int(code)}.json"))


def data360_search(text: str, top: int = 14) -> list:
    def ask() -> list:
        answer = _post(f"{DATA360}/searchv2", {"search": text, "top": top, "count": False, "filter": "type eq 'indicator'"})
        out = []
        for d in answer.get("value", []):
            series = d.get("series_description") or {}
            out.append({"id": series.get("idno") or d.get("id"), "name": series.get("name") or d.get("title"), "database": series.get("database_id"),
                        "definition": (series.get("definition_long") or series.get("definition_short") or "")[:400], "unit": series.get("measurement_unit")})
        return [x for x in out if x["id"] and x["database"]]
    return _kept("data360", {"search": text, "top": top}, ask)


def data360_values(database: str, indicator: str) -> dict:
    """Every economy's series of an indicator: {iso3: {year: value}} (national totals, both sexes where it is so broken down)."""
    def ask() -> dict:
        out: dict = {}
        skip, total = 0, None
        while total is None or skip < total:
            q = urllib.parse.urlencode({"DATABASE_ID": database, "INDICATOR": indicator, "skip": skip})
            page = _get(f"{DATA360}/data?{q}")
            total = page.get("count", 0)
            rows = page.get("value", [])
            if not rows:
                break
            for r in rows:
                if r.get("SEX") not in (None, "_T", "_Z") or r.get("AGE") not in (None, "_T", "_Z") or r.get("URBANISATION") not in (None, "_T", "_Z"):
                    continue
                try:
                    year, value = str(r["TIME_PERIOD"])[:4], float(r["OBS_VALUE"])
                except (TypeError, ValueError, KeyError):
                    continue
                out.setdefault(r["REF_AREA"], {})[year] = value
            skip += len(rows)
            if skip > 60000:
                break
        return out
    return _kept("data360", {"values": indicator, "database": database}, ask)


# ---- rasters about any place ----------------------------------------------------------------------------------------

WSF = "https://download.geoservice.dlr.de/WSF_EVO/files/WSFevolution_v1_{lon}_{lat}.tif"
JRC_FLOOD = "https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard"
RASTERS = CACHE / "world-rasters"
SMOD_CLASSES = {10: 0, 11: 1, 12: 2, 13: 3, 21: 4, 22: 5, 23: 6, 30: 7}
RETURN_PERIODS = (10, 20, 50, 100, 200, 500)
_flood_tiles: list | None = None


def _mosaic(pieces: list[tuple[np.ndarray, tuple]], fill) -> tuple[np.ndarray, tuple]:
    """Windows of one grid, placed on one array: (values, (x0, dx, y0, dy)), rows north to south."""
    if len(pieces) == 1:
        return pieces[0]
    dx, dy = pieces[0][1][1], pieces[0][1][3]
    x0 = min(t[0] for _, t in pieces)
    y0 = max(t[2] for _, t in pieces)
    w = max(round((t[0] - x0) / dx) + a.shape[1] for a, t in pieces)
    h = max(round((y0 - t[2]) / dy) + a.shape[0] for a, t in pieces)
    out = np.full((h, w), fill, dtype=pieces[0][0].dtype)
    for a, t in pieces:
        c, r = round((t[0] - x0) / dx), round((y0 - t[2]) / dy)
        out[r:r + a.shape[0], c:c + a.shape[1]] = a
    return out, (x0, dx, y0, dy)


def _flood_tiles_at(lat: float, lon: float) -> list[str]:
    """The names of the tiles of the river flood maps holding a place, from the release's own index of tiles."""
    global _flood_tiles
    if _flood_tiles is None:
        index = _get(f"{JRC_FLOOD}/tile_extents.geojson")
        _flood_tiles = []
        for f in index["features"]:
            c = f["geometry"]["coordinates"]
            while isinstance(c[0][0], list):
                c = c[0]
            xs, ys = [p[0] for p in c], [p[1] for p in c]
            _flood_tiles.append((min(xs), min(ys), max(xs), max(ys), f"ID{f['properties']['id']}_{f['properties']['name']}"))
    return [t[4] for t in _flood_tiles if t[0] <= lon <= t[2] and t[1] <= lat <= t[3]]


def raster_window(kind: str, lat: float, lon: float, rp: int | None = None) -> dict:
    """A raster about a place, on a grid snapped so that places nearby share one: its meta, with the address of its bytes under /data/
    (one byte a cell, rows south to north). The window is 1.2 degrees across, or 5 for urbanisation.
        settlement   byte = (year first built - 1984) * 8; 0 where never (WSF Evolution, 30 m)
        urbanisation byte = class index * 32 + 16, from SMOD's ladder of 10, 11, 12, 13, 21, 22, 23, 30 (GHS-SMOD R2023A, 1 km)
        flood        byte = depth / 0.05 m, to 12.75 m; 0 where dry (GloFAS river flood maps, 90 m, by return period)"""
    half = 2.5 if kind == "urbanisation" else 0.6
    step = 1.0 if kind == "urbanisation" else 0.25
    c_lat, c_lon = round(lat / step) * step, round(lon / step) * step
    s, n, w, e = c_lat - half, c_lat + half, c_lon - half, c_lon + half
    tag = f"{kind}{'-rp' + str(rp) if rp else ''}_{c_lat:+.2f}_{c_lon:+.2f}"
    meta_path = RASTERS / f"{tag}.json"
    if meta_path.exists():
        return json.loads(meta_path.read_text(encoding="utf-8"))
    RASTERS.mkdir(parents=True, exist_ok=True)
    if kind == "settlement":
        pieces = []
        for la in range(int(s // 2) * 2, int(n // 2) * 2 + 2, 2):
            for lo in range(int(w // 2) * 2, int(e // 2) * 2 + 2, 2):
                try:
                    pieces.append(GeoTIFF.url(WSF.format(lon=lo, lat=la)).read_box(w, s, e, n))
                except (TIFFError, OSError, ValueError):
                    continue     # no tile over the sea
        if not pieces:
            raise ValueError("no settlement is mapped here")
        a, t = _mosaic(pieces, 0)
        codes = np.where(a > 1984, np.clip((a - 1984) * 8, 0, 255), 0).astype(np.uint8)
        source, credit = "WSF Evolution", "World Settlement Footprint Evolution, DLR, CC BY 4.0"
    elif kind == "urbanisation":
        path = CACHE / "ghsl-smod" / "GHS_SMOD_E2025_GLOBE_R2023A_4326_30ss_V2_0.tif"
        if not path.exists():
            raise ValueError("the degree of urbanisation is not fetched: python data/acquisition/fetch.py ghsl-smod")
        a, t = GeoTIFF.open(path).read_box(w, s, e, n)
        lookup = np.zeros(256, np.uint8)
        for v, k in SMOD_CLASSES.items():
            lookup[v] = k * 32 + 16
        codes = lookup[np.clip(a, 0, 255).astype(np.uint8)]
        source, credit = "GHS-SMOD R2023A", "GHS-SMOD R2023A, European Commission JRC, CC BY 4.0"
    elif kind == "flood":
        if rp not in RETURN_PERIODS:
            raise ValueError(f"the flood maps are of return periods {RETURN_PERIODS}")
        pieces = []
        names = set()
        for la, lo in ((s, w), (s, e), (n, w), (n, e)):
            names.update(_flood_tiles_at(la, lo))
        for name in names:
            try:
                a, t = GeoTIFF.url(f"{JRC_FLOOD}/RP{rp}/{name}_RP{rp}_depth.tif").read_box(w, s, e, n)
                pieces.append((np.nan_to_num(a, nan=0, posinf=0, neginf=0), t))
            except (TIFFError, OSError, ValueError):
                continue
        if not pieces:
            raise ValueError("no river flood map covers this place")
        a, t = _mosaic(pieces, 0)
        a[a < 0] = 0
        codes = np.clip(np.round(a / 0.05), 0, 255).astype(np.uint8)
        source, credit = "GloFAS river flood maps", "CEMS GloFAS river flood hazard maps, European Union, CC BY 4.0"
    else:
        raise ValueError(f"no raster {kind}")
    x0, dx, y0, dy = t
    h, wd = codes.shape
    np.flipud(codes).tofile(RASTERS / f"{tag}.bin")
    meta = {"kind": kind, "return_period": rp, "box": [y0 - h * dy, x0, y0, x0 + wd * dx], "width": wd, "height": h, "url": f"/data/world-rasters/{tag}.bin", "source": source, "credit": credit,
            "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    meta_path.write_text(json.dumps(meta), encoding="utf-8")
    return meta
