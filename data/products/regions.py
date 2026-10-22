"""Build regions products: the cells of a global grid, each named by the region containing it.

A regions product partitions the Earth's surface into places with names, such as
countries or provinces, as a raster of region ids on a regular longitude-latitude
grid. Any field can then be summarised over the regions (the areal mean of a
temperature, the sum of a count of people) by domains/geography/js/zonal.js,
and a value per region can be painted on the map. The product holds geometry,
codes and names only, as the geography product does (geography.py); values per
region come from elsewhere and are joined by code.

Two products are built:

countries
    Natural Earth 1:50 million countries (the geography product's source), coded
    by ISO 3166-1 alpha-3, falling back as domains/geography/js/join.js does
    where Natural Earth gives no ISO code (-99): its corrected ISO column
    (ISO_A3_EH, which gives FRA and NOR), then its own ADM0_A3 (KOS for Kosovo).
    Where two features share a code (Australia's external territories carry
    AUS), the feature whose ADM0_A3 is the code keeps it and the others take
    their ADM0_A3.
admin1
    Natural Earth 1:10 million first-order divisions (states, provinces,
    prefectures), coded by ISO 3166-2 (JP-13) where Natural Earth gives a valid
    code held by that division alone, and by Natural Earth's adm1_code (JPN-1860)
    otherwise; each names its country by the code of the countries product.

A cell belongs to the region whose polygon contains the cell's centre, by the
even-odd rule, scanline by scanline (Foley et al. 1990, Computer Graphics:
Principles and Practice, 2nd ed., sec. 3.6), as domains/geography/js/raster.js
draws polygons: neighbouring polygons that share an edge never both claim a
cell. Where Natural Earth's polygons overlap, the later feature takes the cell.
A region smaller than a cell may hold no cell; it is still listed, with no cells.

The outlines, for drawing, are each region's rings cut where they run along the
antimeridian or the South Pole's edge of the map (Natural Earth cuts polygons
there), and simplified by the algorithm of Ramer (1972) and Douglas and Peucker
(1973) to a tolerance in degrees, with longitudes scaled by the cosine of the
latitude so that the tolerance is about the same distance everywhere.

Output, in data/cache/products/:

regions-<name>.bin
    uint16 little-endian, one id per cell, rows north to south from 90 degrees,
    columns east from -180 degrees; 0 is no region, and id k is regions[k - 1].
regions-<name>.json
    The grid, the regions ({id, code, name, country, parent, ...}, with each
    region's number of cells and their area on the sphere of radius 6371.0088 km)
    and the source.
regions-<name>-outlines.json
    {coords: [lon, lat, ...], runs: [[first pair, number of pairs], ...],
    region: [id per run]}, polylines in degrees rounded to 0.001 (about 100 m).

Usage::

    python data/products/regions.py                     # both, on a grid of 0.1 degrees
    python data/products/regions.py --only countries --step 0.25
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

EARTH_RADIUS_KM = 6371.0088          # mean radius R1 of GRS 80 and WGS 84 (Moritz 2000)
ISO_3166_2 = re.compile(r"^[A-Z]{2}-[A-Z0-9]{1,3}$")
TOLERANCE = {"countries": 0.02, "admin1": 0.02}   # degrees, for the outlines


def no_code(value) -> bool:
    return value in (None, "", "-99", -99)


# ---- sources ------------------------------------------------------------------

def countries() -> tuple[list[dict], dict]:
    """The countries of Natural Earth 1:50m, with their codes; and the source."""
    folder = DATA / "cache" / "natural-earth-50m"
    manifest = json.loads((DATA / "registry" / "natural-earth-50m.json").read_text(encoding="utf-8"))
    features = json.loads((folder / "ne_50m_admin_0_countries.geojson").read_text(encoding="utf-8"))["features"]
    claims = []
    for f in features:
        p = f["properties"]
        code = next((p[k] for k in ("ISO_A3_EH", "ISO_A3", "ADM0_A3") if not no_code(p.get(k))), None)
        claims.append(code)
    regions, taken = [], {}
    for f, code in zip(features, claims):
        p = f["properties"]
        shared = claims.count(code) > 1
        if shared and p["ADM0_A3"] != code:
            code = p["ADM0_A3"]
        if code in taken:
            raise SystemExit(f"two countries claim the code {code}: {taken[code]} and {p['NAME']}")
        taken[code] = p["NAME"]
        regions.append({"code": code, "name": p["NAME"], "country": code, "parent": None,
                        "adm0_a3": p["ADM0_A3"], "continent": p["CONTINENT"], "geometry": f["geometry"]})
    return regions, source_of(manifest, folder, "Natural Earth 1:50 million countries")


def admin1(country_codes: dict) -> tuple[list[dict], dict]:
    """The first-order divisions of Natural Earth 1:10m, with their codes; and the source."""
    folder = DATA / "cache" / "natural-earth-10m-admin1"
    manifest = json.loads((DATA / "registry" / "natural-earth-10m-admin1.json").read_text(encoding="utf-8"))
    features = json.loads((folder / "ne_10m_admin_1_states_provinces.geojson").read_text(encoding="utf-8"))["features"]
    iso = [f["properties"].get("iso_3166_2") or "" for f in features]
    held = {}
    for code in iso:
        held[code] = held.get(code, 0) + 1
    regions = []
    for f, code in zip(features, iso):
        p = f["properties"]
        valid = bool(ISO_3166_2.match(code))
        country = country_codes.get(p["adm0_a3"], p["adm0_a3"])
        region = {"code": code if valid and held[code] == 1 else p["adm1_code"], "name": p["name"] or p["adm1_code"],
                  "country": country, "parent": country, "type": p.get("type_en"), "adm1": p["adm1_code"]}
        if valid and held[code] > 1:
            region["iso"] = code                                   # shared by several divisions
        region["geometry"] = f["geometry"]
        regions.append(region)
    codes = [r["code"] for r in regions]
    if len(set(codes)) != len(codes):
        raise SystemExit("admin-1 codes are not unique")
    return regions, source_of(manifest, folder, "Natural Earth 1:10 million states and provinces")


def source_of(manifest: dict, folder: Path, title: str) -> dict:
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    return {"dataset": manifest["id"], "title": title, "provider": manifest["provider"], "credit": manifest.get("credit"),
            "licence": manifest["licence"], "citation": manifest["citation"], "retrieved_at": retrieval["retrieved_at"]}


def rings_of(geometry: dict) -> list[np.ndarray]:
    polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
    return [np.asarray(ring, dtype=float)[:, :2] for polygon in polygons for ring in polygon]


# ---- rasterising --------------------------------------------------------------

def fill_rings(rings: list[np.ndarray], grid: dict, value: int, target: np.ndarray) -> int:
    """Sets the cells whose centres lie inside the rings (even-odd) to value; returns their number.

    The scanline of row j passes through the centres y = j + 0.5 (in cells from
    the north edge); an edge from y_lo to y_hi crosses the rows with
    y_lo <= j + 0.5 < y_hi (half-open, so a vertex shared by two edges counts
    once), and the cells of a row between the crossings 2m and 2m + 1, sorted,
    are inside.
    """
    nlon, nlat = grid["nlon"], grid["nlat"]
    west, north = grid["lon0"] - grid["dlon"] / 2, grid["lat0"] + grid["dlat"] / 2
    a = np.concatenate(rings)
    b = np.concatenate([np.roll(r, -1, axis=0) for r in rings])
    xa, ya = (a[:, 0] - west) / grid["dlon"], (north - a[:, 1]) / grid["dlat"]
    xb, yb = (b[:, 0] - west) / grid["dlon"], (north - b[:, 1]) / grid["dlat"]
    keep = ya != yb
    xa, ya, xb, yb = xa[keep], ya[keep], xb[keep], yb[keep]
    up = ya < yb
    ylo, yhi = np.where(up, ya, yb), np.where(up, yb, ya)
    xlo, xhi = np.where(up, xa, xb), np.where(up, xb, xa)
    ja = np.maximum(np.ceil(ylo - 0.5), 0).astype(np.int64)
    jb = np.minimum(np.ceil(yhi - 0.5) - 1, nlat - 1).astype(np.int64)
    n = np.maximum(jb - ja + 1, 0)
    if not n.sum():
        return 0
    edge = np.repeat(np.arange(len(n)), n)
    j = ja[edge] + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
    x = xlo[edge] + (j + 0.5 - ylo[edge]) * ((xhi - xlo) / (yhi - ylo))[edge]
    order = np.lexsort((x, j))
    j, x = j[order], x[order]
    if len(j) % 2 or (j[0::2] != j[1::2]).any():
        raise ValueError("a ring is not closed")
    rows = j[0::2]
    i0 = np.maximum(np.ceil(x[0::2] - 0.5), 0).astype(np.int64)
    i1 = np.minimum(np.ceil(x[1::2] - 0.5) - 1, nlon - 1).astype(np.int64)
    filled = 0
    for r, c0, c1 in zip(rows[i1 >= i0], i0[i1 >= i0], i1[i1 >= i0]):
        target[r, c0:c1 + 1] = value
        filled += c1 - c0 + 1
    return int(filled)


def cell_areas_km2(grid: dict) -> np.ndarray:
    """Areas of the cells of each row on the sphere of radius R (km^2)."""
    north = grid["lat0"] + grid["dlat"] / 2
    edges = np.radians(np.clip(north - np.arange(grid["nlat"] + 1) * grid["dlat"], -90, 90))
    return EARTH_RADIUS_KM ** 2 * math.radians(grid["dlon"]) * (np.sin(edges[:-1]) - np.sin(edges[1:]))


# ---- outlines -----------------------------------------------------------------

def simplify(points: np.ndarray, tolerance: float) -> np.ndarray:
    """Ramer-Douglas-Peucker: the points of a polyline that keep it within `tolerance` of itself.

    Distances are taken in degrees, with longitude scaled by the cosine of the
    polyline's mean latitude. The ends are kept; a closed ring (first point
    equal to the last) is split first at the point farthest from its start.
    """
    n = len(points)
    if n <= 3:
        return points
    k = math.cos(math.radians(float(points[:, 1].mean())))
    p = points * (k, 1.0)
    keep = np.zeros(n, dtype=bool)
    keep[0] = keep[-1] = True
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, d, q = p[i], p[j] - p[i], p[i + 1:j] - p[i]
        length = math.hypot(d[0], d[1])
        dist = np.abs(d[0] * q[:, 1] - d[1] * q[:, 0]) / length if length > 0 else np.hypot(q[:, 0], q[:, 1])
        m = int(np.argmax(dist))
        if dist[m] > tolerance:
            keep[i + 1 + m] = True
            stack += [(i, i + 1 + m), (i + 1 + m, j)]
    return points[keep]


def cut(ring: np.ndarray) -> list[np.ndarray]:
    """A ring as polylines, without the edges that run along the map's cuts (|lon| = 180, lat = -90)."""
    on_cut = lambda a, b: (abs(a[0]) >= 179.9999 and abs(b[0]) >= 179.9999) or (a[1] <= -89.9999 and b[1] <= -89.9999)
    pieces, current = [], [ring[0]]
    for a, b in zip(ring[:-1], ring[1:]):
        if on_cut(a, b):
            if len(current) > 1:
                pieces.append(np.asarray(current))
            current = [b]
        else:
            current.append(b)
    if len(current) > 1:
        pieces.append(np.asarray(current))
    # A ring cut nowhere but closed stays one piece; one cut at its seam is joined again.
    if len(pieces) > 1 and np.array_equal(pieces[0][0], pieces[-1][-1]) and not on_cut(ring[-1], ring[0]):
        pieces[0] = np.concatenate([pieces[-1], pieces[0][1:]])
        pieces.pop()
    return pieces


def outlines(regions: list[dict], tolerance: float) -> dict:
    coords, runs, ids = [], [], []
    for rid, region in enumerate(regions, start=1):
        for ring in rings_of(region["geometry"]):
            for piece in cut(ring):
                line = simplify(piece, tolerance)
                extent = np.ptp(piece, axis=0).max() if len(piece) else 0
                if len(line) < 2 or (len(line) < 4 and extent < tolerance):
                    continue
                runs.append([len(coords) // 2, len(line)])
                ids.append(rid)
                coords += [round(float(v), 3) for v in line.ravel()]
    return {"coords": coords, "runs": runs, "region": ids}


# ---- products -----------------------------------------------------------------

def build(name: str, title: str, regions: list[dict], source: dict, step: float) -> dict:
    nlon, nlat = round(360 / step), round(180 / step)
    grid = {"lon0": round(-180 + step / 2, 9), "dlon": step, "nlon": nlon, "lat0": round(90 - step / 2, 9),
            "dlat": step, "nlat": nlat, "periodic": True, "crs": "WGS 84 geographic coordinates (EPSG:4326), degrees",
            "order": "rows north to south, columns eastward"}
    raster = np.zeros((nlat, nlon), dtype="<u2")
    for rid, region in enumerate(regions, start=1):
        fill_rings(rings_of(region["geometry"]), grid, rid, raster)
    counts = np.bincount(raster.ravel(), minlength=len(regions) + 1)
    areas = np.bincount(raster.ravel(), weights=np.repeat(cell_areas_km2(grid), nlon), minlength=len(regions) + 1)
    listed = []
    for rid, region in enumerate(regions, start=1):
        entry = {"id": rid} | {k: v for k, v in region.items() if k != "geometry"}
        entry |= {"cells": int(counts[rid]), "area_km2": round(float(areas[rid]), 1)}
        listed.append(entry)
    lines = outlines(regions, TOLERANCE.get(name, 0.02))
    OUT.mkdir(parents=True, exist_ok=True)
    bin_path, lines_path = OUT / f"regions-{name}.bin", OUT / f"regions-{name}-outlines.json"
    bin_path.write_bytes(raster.tobytes())
    lines_path.write_text(json.dumps(lines, separators=(",", ":")) + "\n", encoding="utf-8")
    product = {
        "schema": "regions/1", "name": name, "title": title, "source": source, "grid": grid,
        "file": bin_path.name,
        "encoding": "uint16 little-endian, row-major, rows north to south; 0 = no region; id k is regions[k-1]",
        "assignment": "a cell belongs to the region whose polygon contains its centre (even-odd rule)",
        "radius_km": EARTH_RADIUS_KM,
        "regions": listed,
        "outlines": {"file": lines_path.name, "pairs": len(lines["coords"]) // 2, "runs": len(lines["runs"]),
                     "tolerance_deg": TOLERANCE.get(name, 0.02)},
    }
    json_path = OUT / f"regions-{name}.json"
    json_path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("regions", {"name": name, "file": json_path.name, "title": title, "count": len(listed),
                         "grid": {k: grid[k] for k in ("lon0", "dlon", "nlon", "lat0", "dlat", "nlat", "periodic")},
                         "source": {k: source[k] for k in ("dataset", "credit", "retrieved_at")}}, first="countries")
    return product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step", type=float, default=0.1, help="cell size in degrees, dividing 180 (default 0.1)")
    parser.add_argument("--only", choices=["countries", "admin1"], help="build one product only")
    args = parser.parse_args(argv)
    for dataset in ("natural-earth-50m", "natural-earth-10m-admin1"):
        if not (DATA / "cache" / dataset / "retrieval.json").exists():
            print(f"{dataset} is not in the cache; run: python data/acquisition/fetch.py {dataset}", file=sys.stderr)
            return 1
    nations, source = countries()
    built = []
    if args.only in (None, "countries"):
        built.append(build("countries", "Countries (Natural Earth 1:50m)", nations, source, args.step))
    if args.only in (None, "admin1"):
        by_adm0 = {r["adm0_a3"]: r["code"] for r in nations}
        divisions, source1 = admin1(by_adm0)
        built.append(build("admin1", "States and provinces (Natural Earth 1:10m)", divisions, source1, args.step))
    for p in built:
        empty = sum(1 for r in p["regions"] if not r["cells"])
        print(f"regions-{p['name']}: {len(p['regions'])} regions on {p['grid']['nlon']} x {p['grid']['nlat']} cells "
              f"({empty} too small to hold a cell); outlines of {p['outlines']['pairs']} points")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
