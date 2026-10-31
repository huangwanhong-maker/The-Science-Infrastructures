"""New York in the city model, from Valdrada's compiled city.

Valdrada (M. K. Turkcan, https://github.com/mkturkcan/valdrada) compiles New York City from its public records into
binary tiles; the data library fetches the tiles of Lower Manhattan at a pinned revision (data/registry/valdrada-nyc.json)
and this module reads them into an area of the city model, in the schema street_records.py gives Melbourne's.

1. **Tiles.** A tile (t_<tx>_<tz>.bin) is a CTL1 container: the magic 'CTL1' and a header length (u32, little-endian), a
   JSON header (the format's version, the tile's tx, tz and origin, street names, the terrain's resolution, its sections'
   names, types, offsets and lengths), then the sections, aligned to four bytes. The layout is the compiler's
   (client/tools/pipeline/binio.mjs and compile.mjs, v0.6.2, commit 09fec97); the records read here are below (BLDG,
   HOLE, FURN).
2. **Coordinates.** Valdrada's frame is equirectangular about (40.7831 N, 73.9712 W): x = (lon - LON0) M_LON east and
   z = -(lat - LAT0) M_LAT south, with M_LAT = 111132 m and M_LON = 111320 cos(LAT0) m a degree; a tile's records are
   tile-local, from its origin (512 tx, 512 tz). The compiler made x and z from WGS 84 longitudes and latitudes (NYC Open
   Data's GeoJSON), so inverting the frame gives them back, to the Float32 of the records (under 1 mm at a tile's far
   corner); they are then projected as Melbourne's are (street_records.enu: geodetic to ECEF to east-north-up at the
   area's centre). The frame's own metres are not true metres (0.4% short north-south here), and are not used.
3. **Buildings.** A building is the area's when the mean of its footprint's vertices lies in the box; it is whole, as in
   Melbourne's. height_m is the compiled height above the ground: the roof height NYC's Building Footprints record
   (aerial survey, feet converted), else a low default of 4-8 m; a building the city's file lacks is OpenStreetMap's, at
   its tagged height (height_source "record" for all: the tiles do not say which). levels are PLUTO's floors as the
   compiler kept them, kind the style it gave the building (prewar_apt, glass_tower_blue, ...), holes its light courts.
4. **Streets.** OpenStreetMap's ways, from osm-city's new-york tiles, by Melbourne's rules (street_records.osm_roads),
   without tunnels. The tiles hold the city's street centrelines too, with their widths, but drawn for the renderer:
   their header lists only some of the junctions, and most centrelines stop short of the junction they meet (by
   several metres, up to tens), so that they do not form a network; OpenStreetMap's ways share their junctions' nodes.
5. **Trees.** The tiles' trees in the box: the street-tree census and NYC Parks' tree points with species and trunk
   diameter (in), and trees the compiler scattered in parks and yards (kind says which).
6. **Water.** Traced on each tile's ground: the terrain grid stands at about +3.2 m on land and -3.5 m in the harbour,
   with -9 m along the bulkheads and in cuts inland (tunnel portals). Water is the ground below 0 m that is connected to
   the harbour's level within the tile; its outline is drawn by marching squares between the grid's nodes (a corner's
   water separated from its diagonal neighbour's), the tile's edges closing it, piers and other land within it its
   holes; collinear points are dropped. A space of kind water, clipped to the box.
7. **Open spaces.** OpenStreetMap's, from osm-city's new-york tiles, by Melbourne's rules (street_records.osm_spaces):
   parks, gardens, playgrounds, green land, railway land, squares and pedestrian areas.

*Rationale.* Each source gives what it holds best: Valdrada the city's buildings at their surveyed heights (where
OpenStreetMap's Manhattan is mostly without heights), its trees and its shore; OpenStreetMap a street network whose
ways meet.
"""

from __future__ import annotations

import json
import math
import struct
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
FOLDER = DATA / "cache" / "valdrada-nyc"
REVISION = "4b3485a034be96acd34a6e76bd36f340d03a9110"

# Valdrada's frame (client/src/shared/geo.js).
LAT0, LON0 = 40.7831, -73.9712
M_LAT = 111132.0
M_LON = 111320.0 * math.cos(math.radians(LAT0))
TILE = 512
MAGIC = 0x43544C31                     # 'CTL1' as a little-endian u32

TYPES = {"Float32Array": "<f4", "Uint8Array": "u1", "Uint32Array": "<u4", "Int16Array": "<i2", "Uint16Array": "<u2",
         "Int32Array": "<i4", "Float64Array": "<f8"}
BLDG = np.dtype([("start", "<u4"), ("len", "<u2"), ("landmarkId", "<i2"), ("baseY", "<f4"), ("height", "<f4"),
                 ("floors", "u1"), ("style", "u1"), ("r", "u1"), ("g", "u1"), ("b", "u1"), ("roofKind", "u1"),
                 ("flags", "u1"), ("lit", "u1"), ("floorH", "<f4"), ("winW", "<f4"), ("storeH", "<f4"),
                 ("blind", "<u4"), ("area", "<u2"), ("colorVar", "u1"), ("frontIdx1", "u1")])
HOLE = np.dtype([("parent", "<u4"), ("start", "<u4"), ("len", "<u2"), ("flags", "u1"), ("area4", "u1")])
FURN = np.dtype([("k", "u1"), ("p0", "u1"), ("p1", "u1"), ("p2", "u1"), ("x", "<f4"), ("y", "<f4"), ("z", "<f4"), ("rot", "<f4")])
STYLES = ["tenement", "prewar_apt", "postwar_brick", "modern_glass", "deco_masonry", "rowhouse", "loft_castiron", "industrial",
          "civic_stone", "church", "retail_strip", "glass_tower_blue", "project_brick", "retail_modern", "frame_house", "condo_new"]
SPECIES = {0: "other", 1: "London planetree", 2: "honeylocust", 3: "Callery pear", 4: "ginkgo", 5: "oak", 6: "linden",
           7: "maple", 8: "cherry", 9: "Japanese pagoda tree", 10: "zelkova", 11: "ash", 12: "elm"}
TREE_KINDS = {0: "park or campus (scattered)", 1: "street (census)", 2: "yard (scattered)", 4: "median (scattered)"}
HARBOUR = (-4.0, -3.0)                 # the harbour's ground (m), -3.5 in the tiles


class Tile:
    """One CTL1 tile: its header and its sections as arrays."""

    def __init__(self, path: Path):
        buf = path.read_bytes()
        magic, length = struct.unpack_from("<II", buf, 0)
        if magic != MAGIC:
            raise ValueError(f"{path.name}: not a CTL1 tile (magic {magic:#x})")
        self.header = json.loads(buf[8:8 + length].decode("utf-8"))
        base = 8 + length + (-(8 + length)) % 4
        self.sections = {}
        for s in self.header["sections"]:
            a = np.frombuffer(buf, TYPES[s["type"]], s["length"], base + s["offset"])
            if a.nbytes != s["bytes"]:
                raise ValueError(f"{path.name}: section {s['name']} is {a.nbytes} bytes, not {s['bytes']}")
            self.sections[s["name"]] = a
        self.ox, self.oz = self.header["origin"]
        self.names = self.header.get("names", [])
        self.res = int(self.header.get("res", 32))

    def records(self, name: str, dtype: np.dtype) -> np.ndarray:
        a = self.sections.get(name)
        return a.view(dtype) if a is not None and a.size else np.zeros(0, dtype)

    def lonlat(self, lx, lz) -> tuple[np.ndarray, np.ndarray]:
        """WGS 84 longitude and latitude of tile-local x, z."""
        x = np.asarray(lx, float) + self.ox
        z = np.asarray(lz, float) + self.oz
        return x / M_LON + LON0, -z / M_LAT + LAT0


def tiles_in(bbox: list[float]) -> list[tuple[int, int, Path]]:
    """The fetched tiles whose squares meet the box [south, west, north, east]."""
    path = FOLDER / "manifest.json"
    if not path.exists():
        raise FileNotFoundError(f"no Valdrada tiles in {FOLDER}; run: python data/acquisition/fetch.py valdrada-nyc")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    s, w, n, e = bbox
    tile = lambda lat, lon: (math.floor((lon - LON0) * M_LON / TILE), math.floor(-(lat - LAT0) * M_LAT / TILE))  # noqa: E731
    (tx0, tz1), (tx1, tz0) = tile(s, w), tile(n, e)
    out = []
    for tx in range(tx0, tx1 + 1):
        for tz in range(tz0, tz1 + 1):
            entry = manifest["tiles"].get(f"{tx}_{tz}")
            if not entry:
                continue
            p = FOLDER / "tiles" / entry["f"]
            if not p.exists():
                raise FileNotFoundError(f"tile {entry['f']} not fetched; run: python data/acquisition/fetch.py valdrada-nyc")
            out.append((tx, tz, p))
    return out


# ---- water ----------------------------------------------------------------------------------------------------

def water_mask(g: np.ndarray) -> np.ndarray:
    """The nodes below 0 m connected (4-neighbours) to a node at the harbour's level."""
    low = g < 0
    mask = (g > HARBOUR[0]) & (g < HARBOUR[1])
    while True:
        grown = mask.copy()
        grown[1:, :] |= mask[:-1, :]
        grown[:-1, :] |= mask[1:, :]
        grown[:, 1:] |= mask[:, :-1]
        grown[:, :-1] |= mask[:, 1:]
        grown &= low
        if (grown == mask).all():
            return mask
        mask = grown


def _left(p, q, w) -> bool:
    return (q[0] - p[0]) * (w[1] - p[1]) - (q[1] - p[1]) * (w[0] - p[0]) > 0


def water_loops(mask: np.ndarray, spacing: float) -> list[list[tuple[float, float]]]:
    """Closed outlines of the water nodes by marching squares, in tile-local (x, z) metres, each with the water on its
    left (so that the water's outer outlines are counter-clockwise in x, z and the land within it clockwise). The mask
    is padded with a ring of land a node beyond the tile, and the segments, oriented there, are then moved onto the
    tile's edge, so that outlines reaching the edge close along it."""
    n = mask.shape[0] - 1
    P = np.zeros((n + 3, n + 3), bool)
    P[1:-1, 1:-1] = mask
    at = lambda p: (p - 1) * spacing                        # noqa: E731
    edge = n * spacing
    onto = lambda c: (min(max(c[0], 0.0), edge), min(max(c[1], 0.0), edge))   # noqa: E731
    out: dict[tuple, list] = {}
    for j in range(n + 2):
        row0, row1 = P[j], P[j + 1]
        for i in range(n + 2):
            corners = ((i, j, row0[i]), (i + 1, j, row0[i + 1]), (i + 1, j + 1, row1[i + 1]), (i, j + 1, row1[i]))
            ws = [c[2] for c in corners]
            if all(ws) or not any(ws):
                continue
            xy = [(at(c[0]), at(c[1])) for c in corners]
            # The crossings on the cell's four edges (top, right, bottom, left), at their middles.
            edges = [(k, (k + 1) % 4) for k in range(4)]
            cross = {k: ((xy[a][0] + xy[b][0]) / 2, (xy[a][1] + xy[b][1]) / 2) for k, (a, b) in enumerate(edges) if ws[a] != ws[b]}
            if len(cross) == 2:
                # One segment; the wet corners' centroid lies on the water's side of it (one, two or three of them).
                wet = [xy[k] for k in range(4) if ws[k]]
                pairs = [(*cross, (sum(c[0] for c in wet) / len(wet), sum(c[1] for c in wet) / len(wet)))]
            else:
                # A saddle: each water corner kept apart from the other, by the edges either side of it.
                pairs = [((w - 1) % 4, w, xy[w]) for w in range(4) if ws[w]]
            for a, b, w in pairs:
                p, q = cross[a], cross[b]
                if not _left(p, q, w):
                    p, q = q, p
                p, q = onto(p), onto(q)
                if p == q:
                    continue
                key = (round(p[0], 3), round(p[1], 3))
                out.setdefault(key, []).append((round(q[0], 3), round(q[1], 3)))
    loops = []
    while out:
        start = next(iter(out))
        loop, cur = [start], start
        while True:
            nxt = out[cur].pop()
            if not out[cur]:
                del out[cur]
            if nxt == start or nxt not in out:
                break
            loop.append(nxt)
            cur = nxt
        if len(loop) >= 3:
            loops.append(_simplified(loop))
    return [lp for lp in loops if len(lp) >= 3]


def _simplified(loop: list[tuple]) -> list[tuple]:
    """A closed outline without its collinear points."""
    pts = loop[:]
    changed = True
    while changed and len(pts) > 3:
        changed = False
        keep = []
        for k, p in enumerate(pts):
            a, b = pts[k - 1], pts[(k + 1) % len(pts)]
            if abs((p[0] - a[0]) * (b[1] - a[1]) - (p[1] - a[1]) * (b[0] - a[0])) < 1e-6:
                changed = True
                continue
            keep.append(p)
        pts = keep if len(keep) >= 3 else pts
        if len(keep) < 3:
            break
    return pts


def _area(loop) -> float:
    return 0.5 * sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(loop, loop[1:] + loop[:1]))


def _contains(loop, p) -> bool:
    inside = False
    for (x1, y1), (x2, y2) in zip(loop, loop[1:] + loop[:1]):
        if (y1 > p[1]) != (y2 > p[1]) and p[0] < (x2 - x1) * (p[1] - y1) / (y2 - y1) + x1:
            inside = not inside
    return inside


def tile_water(t: Tile) -> list[tuple[list, list[list]]]:
    """The tile's water: outer outlines with their holes, in tile-local (x, z) metres."""
    g = t.sections["terrain"].reshape(t.res + 1, t.res + 1)       # [j (z), i (x)]
    mask = water_mask(g)
    if not mask.any():
        return []
    loops = water_loops(mask, TILE / t.res)
    outers = [lp for lp in loops if _area(lp) > 0]
    holes = [lp for lp in loops if _area(lp) < 0]
    polys = [(o, []) for o in outers]
    for h in holes:
        within = [k for k, (o, _) in enumerate(polys) if _contains(o, h[0])]
        if within:
            polys[min(within, key=lambda k: _area(polys[k][0]))][1].append(h)
    return polys


# ---- the area --------------------------------------------------------------------------------------------------

def area(spec: dict, sr) -> tuple[dict, dict, dict]:
    """New York's area of the city model: (area, tally, its sources, credit and conventions). sr is street_records."""
    s, w, n, e = spec["bbox"]
    lat0, lon0 = (s + n) / 2, (w + e) / 2
    project = sr._Projector(lat0, lon0)
    inside = lambda lat, lon: s <= lat <= n and w <= lon <= e     # noqa: E731
    tally = {"tiles": 0, "buildings": 0, "courts": 0, "roads": 0, "trees": 0, "water": 0, "spaces": 0, "buildings_outside": 0}
    buildings, trees, water = [], [], []
    tiles = tiles_in(spec["bbox"])
    for tx, tz, path in tiles:
        t = Tile(path)
        tally["tiles"] += 1
        XZ = t.sections["bldgXZ"].reshape(-1, 2)
        courts: dict[int, list] = {}
        for h in t.records("bholes", HOLE):
            courts.setdefault(int(h["parent"]), []).append(XZ[int(h["start"]):int(h["start"]) + int(h["len"])])
        for k, b in enumerate(t.records("bldg", BLDG)):
            ring = XZ[int(b["start"]):int(b["start"]) + int(b["len"])]
            lon, lat = t.lonlat(ring[:, 0], ring[:, 1])
            if not inside(float(lat.mean()), float(lon.mean())):
                tally["buildings_outside"] += 1
                continue
            footprint = sr.oriented(project(list(zip(lat, lon))), True)
            if len(footprint) < 3:
                continue
            holes = []
            for c in courts.get(k, []):
                clon, clat = t.lonlat(c[:, 0], c[:, 1])
                holes.append(sr.oriented(project(list(zip(clat, clon))), False))
            style = int(b["style"])
            rec = {"id": f"valdrada/{tx}_{tz}/{k}", "kind": STYLES[style] if style < len(STYLES) else f"style-{style}",
                   "footprint": footprint, "holes": holes, "height_m": round(float(b["height"]), 2), "min_height_m": 0.0,
                   "height_source": "record"}
            if b["floors"]:
                rec["levels"] = int(b["floors"])
            buildings.append(rec)
            tally["courts"] += len(holes)
        for f in t.records("furn", FURN):
            if int(f["k"]) != 1:
                continue
            lon, lat = t.lonlat(float(f["x"]), float(f["z"]))
            if not inside(float(lat), float(lon)):
                continue
            (x, y), = project([(float(lat), float(lon))])
            trees.append({"at": [x, y], "species": SPECIES.get(int(f["p0"]), "other"), "dbh_in": int(f["p1"]),
                          "kind": TREE_KINDS.get(int(f["p2"]), "other")})
        for k, (outer, holes) in enumerate(tile_water(t)):
            to_ll = lambda loop: list(zip(*reversed(t.lonlat([p[0] for p in loop], [p[1] for p in loop]))))  # noqa: E731
            ring = sr.clip_to_box([(float(a), float(b)) for a, b in to_ll(outer)], spec["bbox"])
            if len(ring) < 3:
                continue
            kept = []
            for h in holes:
                hr = sr.clip_to_box([(float(a), float(b)) for a, b in to_ll(h)], spec["bbox"])
                if len(hr) >= 3:
                    kept.append(sr.oriented(project(hr), False))
            water.append({"id": f"water/{tx}_{tz}/{k + 1}", "kind": "water", "ring": sr.oriented(project(ring), True), "holes": kept})
    osm = {}
    elements, osm["tiles"] = sr.osm_elements(spec)
    roads = sr.osm_roads(elements, spec, project, osm, tunnels=False)
    spaces = sr.osm_spaces(elements, spec, project, osm)
    tally.update({"buildings": len(buildings), "roads": len(roads), "trees": len(trees), "water": len(water), "spaces": len(spaces),
                  "osm": osm})
    area = {"id": spec["id"], "name": spec["name"], "origin": {"lat": round(lat0, 7), "lon": round(lon0, 7)},
            "bbox": spec["bbox"], "frame": "local ENU metres from origin: x east, y north",
            "projection": ("Valdrada's frame inverted to WGS 84 longitude and latitude, then exact: geodetic (WGS 84, on the ellipsoid) "
                           "to ECEF, rotated to east-north-up at the origin; the up component, the curvature's drop below the tangent "
                           "plane, is dropped"),
            "curvature_drop_m": round(project.worst_drop, 2), "equirectangular_error_m": round(project.worst_plane, 3),
            "buildings": buildings, "roads": roads, "spaces": water + spaces, "trees": trees}
    retrieval = FOLDER / "retrieval.json"
    retrieved = json.loads(retrieval.read_text(encoding="utf-8")).get("retrieved_at") if retrieval.exists() else None
    about = {
        "sources": [{"dataset": "valdrada-nyc", "revision": REVISION, "credit": "Valdrada (M. K. Turkcan), from NYC Open Data; ODbL 1.0, © OpenStreetMap contributors",
                     "retrieved_at": retrieved},
                    {"dataset": "osm-city", "profile": "areas", "credit": "© OpenStreetMap contributors, ODbL 1.0", "retrieved_at": sr.retrieved(sr.CITY)}],
        "credit": "Valdrada (M. K. Turkcan) from NYC Open Data · © OpenStreetMap contributors",
        "conventions": {
            "heights": ("height_m: the compiled height above the ground, the roof height of NYC's Building Footprints (aerial survey, feet "
                        "converted), else a low default of 4-8 m; buildings the city's file lacks are OpenStreetMap's at their tagged height. "
                        "min_height_m 0; levels PLUTO's floors; kind the style Valdrada's compiler gave the building; holes its light courts. "
                        "A building is the area's when its footprint's mean vertex lies in the box, and is whole."),
            "widths": ("roads: OpenStreetMap's ways tagged highway, without tunnels, by Melbourne's rules: width_m the width tag, else "
                       "lanes x 3.2 m plus 1.5 m for each side a sidewalk tag names, else a default of the highway class; bridge true on a "
                       "bridge; clipped to the box (Liang-Barsky). The tiles' own street centrelines stop short of their junctions and are "
                       "not used."),
            "spaces": ("water: the tiles' ground below 0 m connected to the harbour's level (-3.5 m), outlined by marching squares between "
                       "the terrain's nodes (4 or 16 m apart), piers and land within it its holes, clipped to the box; the rest "
                       "OpenStreetMap's by Melbourne's rules (parks, gardens, green land, railway land, squares, pedestrian areas)."),
            "trees": "trees: Valdrada's, the street-tree census and NYC Parks' tree points with species and trunk diameter (in), and trees scattered in parks and yards",
            "rings": "footprints and rings counter-clockwise, holes clockwise, not closed; coordinates rounded to 0.01 m",
        },
    }
    return area, tally, about

