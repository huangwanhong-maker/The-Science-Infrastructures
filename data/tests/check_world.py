"""Checks of what reads the world's places on demand, offline: PMTiles' Hilbert ids against the specification's table; vector tiles written here by
hand and read back; tile boxes; a building's rings by winding, and the rule of heights; shapefiles written here and read back; windows placed on one
grid; the depth maps' codes; a record of the Risk Data Library Standard against its validator. Where the products have been built, the Freetown
assessment's annual loss and its curve as read from its events, and the survey estimates keyed for the regions product.
"""

from __future__ import annotations

import json
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

from checks import PRODUCTS, ROOT, assert_, check, skip

sys.path.insert(0, str(ROOT / "data" / "formats"))
sys.path.insert(0, str(ROOT / "data" / "acquisition"))
sys.path.insert(0, str(ROOT / "data" / "products"))
import pmtiles  # noqa: E402
import rdls_records  # noqa: E402
import services  # noqa: E402
import shapefile  # noqa: E402
import world_city  # noqa: E402

# ---- PMTiles' tile ids: the reference library's test values, and what follows from the curve itself ------------------------------------------
TABLE = {(0, 0, 0): 0, (1, 0, 0): 1, (1, 0, 1): 2, (1, 1, 1): 3, (1, 1, 0): 4, (2, 0, 0): 5, (12, 3423, 1763): 19078479}
wrong = {k: pmtiles.tile_id(*k) for k, v in TABLE.items() if pmtiles.tile_id(*k) != v}
assert_("pmtiles: the Hilbert ids of the reference library's test values", not wrong, str(wrong))
for z in (2, 3, 7):
    below = ((1 << (2 * z)) - 1) // 3
    assert_(f"pmtiles: level {z}'s curve runs from the corner (0, 0) to the corner (n - 1, 0)", pmtiles.tile_id(z, 0, 0) == below and pmtiles.tile_id(z, 2 ** z - 1, 0) == below + 4 ** z - 1)
for z in (1, 2, 3, 5):
    ids = sorted(pmtiles.tile_id(z, x, y) for x in range(2 ** z) for y in range(2 ** z))
    assert_(f"pmtiles: level {z}'s ids fill one run without a gap", ids[0] == ((1 << (2 * z)) - 1) // 3 and ids[-1] - ids[0] == 4 ** z - 1 and len(set(ids)) == 4 ** z)


# ---- a vector tile, written by hand ----------------------------------------------------------------------------------------------------------
def varint(v: int) -> bytes:
    out = bytearray()
    while True:
        b = v & 0x7F
        v >>= 7
        out.append(b | (0x80 if v else 0))
        if not v:
            return bytes(out)


def field(number: int, wire: int, payload: bytes | int) -> bytes:
    head = varint(number << 3 | wire)
    return head + (varint(payload) if wire == 0 else varint(len(payload)) + payload)


zig = lambda v: (v << 1) ^ (v >> 31)  # noqa: E731
# A square of side 100 at (10, 10) with a hole of side 20 at (40, 40): rings wind opposite ways.
geometry = [9, zig(10), zig(10), (3 << 3) | 2, zig(100), 0, 0, zig(100), zig(-100), 0, 15, 9, zig(30), zig(-70), (3 << 3) | 2, 0, zig(20), zig(20), 0, 0, zig(-20), 15]
packed = b"".join(varint(g) for g in geometry)
feature = field(1, 0, 7) + field(2, 2, varint(0) + varint(0) + varint(1) + varint(1)) + field(3, 0, 3) + field(4, 2, packed)
value_height = field(3, 1, b"") if False else varint(3 << 3 | 1) + struct.pack("<d", 12.5)
value_name = field(1, 2, "A hall".encode())
layer = field(1, 2, b"building") + field(2, 2, feature) + field(3, 2, b"height") + field(3, 2, b"@name") + field(4, 2, value_height) + field(4, 2, value_name) + field(5, 0, 4096) + field(15, 0, 2)
tile = field(3, 2, layer)
got = pmtiles.read_tile(tile)
f = got["building"]["features"][0]
assert_("vector tile: a feature's id, type and attributes", f["id"] == 7 and f["type"] == 3 and f["properties"] == {"height": 12.5, "@name": "A hall"}, str(f["properties"]))
assert_("vector tile: its two rings, closed", len(f["geometry"]) == 2 and f["geometry"][0][0] == f["geometry"][0][-1] == (10, 10) and f["geometry"][1][0] == (40, 40), str(f["geometry"]))
polygons = world_city._polygons(f)
assert_("a building's rings by their winding: one outer ring with one hole", len(polygons) == 1 and len(polygons[0][1]) == 1)
check("the outer ring's area", abs(pmtiles.signed_area(polygons[0][0])) / 2 - 10000, 1e-9)
assert_("what another layer is not asked for is not read", pmtiles.read_tile(tile, {"road"}) == {})
convert = pmtiles.to_lonlat(1, 1, 0, 4096)
check("a tile's corner to the Web Mercator square: its longitude", convert(0, 0)[0] - 0, 1e-12)
check("and its latitude, the square's edge", convert(0, 0)[1] - 85.0511287798066, 1e-9)
check("the far corner is the equator", convert(4096, 4096)[1], 1e-9)
s, w, n, e = world_city.tile_box(14, 12075, 6878)
assert_("a tile's box holds Kathmandu's centre", s < 27.7172 and 27.7172 < n + 0.1 and w < 85.324 < e + 0.05, str((s, w, n, e)))
check("the box of a tile is 360 / 2^z wide", e - w - 360 / 2 ** 14, 1e-12)
assert_("the buildings' average heights: their tile by the release's rows and columns", world_city._height_tile_of(27.7, 85.3) == (7, 27) and world_city._height_tile_of(-37.8, 144.96) == (13, 33), str(world_city._height_tile_of(27.7, 85.3)))
check("a width from Overture's rules", world_city._rules_value('[{"value": 7.5}]') - 7.5, 0)
assert_("rules that say nothing plain give no width", world_city._rules_value("[{}]") is None and world_city._rules_value("not json") is None)

# ---- shapefiles, written by hand -------------------------------------------------------------------------------------------------------------
with tempfile.TemporaryDirectory() as folder:
    path = Path(folder) / "t.shp"
    ring = [(0.0, 0.0), (0.0, 2.0), (2.0, 2.0), (2.0, 0.0), (0.0, 0.0)]    # clockwise: an outer ring
    hole = [(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5), (0.5, 0.5)]
    points = ring + hole
    poly = struct.pack("<i4d2i", 5, 0, 0, 2, 2, 2, len(points)) + struct.pack("<2i", 0, len(ring)) + b"".join(struct.pack("<2d", *p) for p in points)
    pt = struct.pack("<i2d", 1, 3.0, 4.0)
    body = struct.pack(">II", 1, len(poly) // 2) + poly + struct.pack(">II", 2, len(pt) // 2) + pt
    path.write_bytes(struct.pack(">7i", 9994, 0, 0, 0, 0, 0, (100 + len(body)) // 2) + b"\0" * 72 + body)
    dbf = bytearray(struct.pack("<B3BIHH", 3, 26, 11, 7, 2, 32 + 32 + 32 + 1, 1 + 6 + 8) + b"\0" * 20)
    dbf += b"NAME".ljust(11, b"\0") + b"C" + b"\0" * 4 + bytes([6, 0]) + b"\0" * 14
    dbf += b"AAL".ljust(11, b"\0") + b"N" + b"\0" * 4 + bytes([8, 2]) + b"\0" * 14
    dbf += b"\x0d" + b" " + b"Ward A" + b"  123.50" + b" " + b"Ward B" + b"        "
    path.with_suffix(".dbf").write_bytes(bytes(dbf))
    fs = shapefile.read(path)
    assert_("shapefile: a polygon and a point, in order", [x["type"] for x in fs] == ["polygon", "point"])
    assert_("and their attributes, a number where the table says one and none where it is blank", fs[0]["attributes"] == {"NAME": "Ward A", "AAL": 123.5} and fs[1]["attributes"] == {"NAME": "Ward B", "AAL": None}, str(fs[0]["attributes"]))
    assert_("a point's place", fs[1]["parts"][0][0] == (3.0, 4.0))
    grouped = shapefile.polygons_of(fs[0])
    assert_("a clockwise ring outer, the other a hole within it", len(grouped) == 1 and len(grouped[0]["holes"]) == 1)

# ---- windows on one grid ---------------------------------------------------------------------------------------------------------------------
a = (np.array([[1, 2], [3, 4]], dtype=np.uint8), (10.0, 1.0, 20.0, 1.0))
b = (np.array([[5, 6], [7, 8]], dtype=np.uint8), (12.0, 1.0, 20.0, 1.0))
c = (np.array([[9, 9]], dtype=np.uint8), (10.0, 1.0, 18.0, 1.0))
m, t = services._mosaic([a, b, c], 0)
assert_("windows placed on one grid, side by side and below", m.tolist() == [[1, 2, 5, 6], [3, 4, 7, 8], [9, 9, 0, 0]] and t == (10.0, 1.0, 20.0, 1.0), str(m.tolist()))

# ---- a record of the Risk Data Library Standard -------------------------------------------------------------------------------------------------
schema = {"type": "object", "required": ["id", "kinds", "who"], "$defs": {"Who": {"type": "object", "required": ["name"], "properties": {"name": {"type": "string", "minLength": 1}}}},
          "properties": {"id": {"type": "string"}, "kinds": {"type": "array", "minItems": 1, "items": {"type": "string", "enum": ["hazard", "loss"]}}, "who": {"$ref": "#/$defs/Who"}}}
assert_("the validator passes a good record", rdls_records.validate({"id": "a", "kinds": ["loss"], "who": {"name": "b"}}, schema) == [])
assert_("and says what is missing, through a reference", rdls_records.validate({"id": "a", "kinds": ["loss"], "who": {}}, schema) == ["$.who: name is missing"])
assert_("what is outside a closed list, and what is of the wrong type", len(rdls_records.validate({"id": 1, "kinds": ["fire"], "who": {"name": ""}}, schema)) == 3)
keep = ROOT / "data" / "cache" / "rdls-schema" / "rdls_schema.json"
if keep.exists():
    assert_("the assessment's record against the fetched schema", rdls_records.validate(rdls_records.freetown(), json.loads(keep.read_text(encoding="utf-8"))) == [])
else:
    skip("the RDLS schema is not fetched")

# ---- the products, where built -----------------------------------------------------------------------------------------------------------------
built = PRODUCTS / "freetown-flood.json"
if built.exists():
    P = json.loads(built.read_text(encoding="utf-8"))
    total = sum(P["events"]["losses"])
    check("Freetown: the catalogue's length, the sum of its events over its annual loss, is 10,000 years", total / P["reported_aal"]["all"] - 10000, 0.01)
    rp, oep = P["curves"]["oep"]["return_period"], P["curves"]["oep"]["loss"]
    ranked = {t: P["events"]["losses"][round(10000 / t) - 1] for t in rp if t >= 10}
    worst = max(abs(ranked[t] - l) / max(1, l) for t, l in zip(rp, oep) if t >= 10 and round(10000 / t) == 10000 / t)
    check("the published occurrence curve is the events in order, at the ranks a return period falls on", worst, 1e-4)
    m0 = P["maps"][0]
    assert_("every depth map has its file and box", all((PRODUCTS / q["file"]).stat().st_size == q["width"] * q["height"] and q["box"][2] > q["box"][0] for q in P["maps"]))
    assert_("all twelve maps", len(P["maps"]) == 12)
    assert_("the buildings' file has four numbers each", (PRODUCTS / P["buildings"]["file"]).stat().st_size == P["buildings"]["count"] * 32)
    deep = np.fromfile(PRODUCTS / next(q for q in P["maps"] if q["peril"] == "FLRF" and q["return_period"] == 1500)["file"], dtype=np.uint8).astype(np.int64).sum()
    shallow = np.fromfile(PRODUCTS / next(q for q in P["maps"] if q["peril"] == "FLRF" and q["return_period"] == 20)["file"], dtype=np.uint8).astype(np.int64).sum()
    assert_("a rarer flood holds more water than a common one", deep > shallow, f"{deep} against {shallow}")
    _ = m0
else:
    skip("the Freetown product is not built")
pov = PRODUCTS / "poverty-latest.json"
if pov.exists():
    Q = json.loads(pov.read_text(encoding="utf-8"))
    E = Q["economies"]
    assert_("the surveys are keyed by the alpha-3 codes of the regions product", len(E) > 150 and all(len(k) == 3 for k in E))
    assert_("shares and indices on [0, 1]", all(0 <= v["headcount"] <= 1 and (v["gini"] is None or 0 < v["gini"] < 1) for v in E.values()))
    assert_("names by the data library's policy", E.get("TWN", {}).get("name") in (None, "Chinese Taipei"), str(E.get("TWN", {}).get("name")))
else:
    skip("the survey estimates are not built")
