"""Checks of the geography products: the reader of TIFF that population_grid.py
uses, the population grid from GHS-POP against the release's own total, and
the regions rasters against the polygons they are drawn from.

    python data/tests/check_geography.py        (run by data/tests/run-tests.mjs)

The checks on products skip, and say so, when a product is missing.
"""

from __future__ import annotations

import io
import json
import math
import sys

import numpy as np

from checks import PRODUCTS, ROOT, assert_, check, skip

import population_grid as pg  # noqa: E402
import regions as rg  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")       # names carry ″, ², π and ×, beyond the console's code page
R_KM = 6371.0088
WORLD_2025 = 8_191_988_536      # GHSL Data Package 2023 (JRC133256), Table 20: GHS-POP R2023A, 2025


def pillow_tiff(values: np.ndarray, compression: str) -> bytes:
    """A TIFF of a uint8 image written by libtiff through Pillow, as an independent encoder."""
    from PIL import Image
    buffer = io.BytesIO()
    Image.fromarray(values).save(buffer, format="TIFF", compression=compression)
    return buffer.getvalue()


# ---- the reader of TIFF ------------------------------------------------------------
try:
    import PIL  # noqa: F401
    have_pillow = True
except ImportError:
    have_pillow = False
if have_pillow:
    rng = np.random.default_rng(20261023)
    # Runs and noise, so that the LZW table fills, is cleared and starts again.
    image = np.repeat(rng.integers(0, 256, (300, 37), dtype=np.uint8), 11, axis=1)[:, :400]
    image[100:200] = rng.integers(0, 256, (100, 400), dtype=np.uint8)
    path = ROOT / "data" / "cache" / "products" / ".check-geography.tif"
    for compression, label in (("tiff_lzw", "LZW"), ("tiff_adobe_deflate", "Deflate"), ("raw", "no compression")):
        path.write_bytes(pillow_tiff(image, compression))
        tiff = pg.Tiff(path)
        rows = np.concatenate([tiff.tile(k)[:, :tiff.width] for k in range(tiff.across * tiff.down)])[:tiff.height]
        tiff.file.close()
        assert_(f"tiff: an image written by libtiff with {label} is read back exactly", np.array_equal(rows, image))
    path.unlink()
else:
    skip("tiff: Pillow is not installed")

source = ROOT / "data" / "cache" / "ghsl-population" / "GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif"
if source.exists() and have_pillow:
    from PIL import Image
    tiff = pg.Tiff(source)
    k = int(np.argmax(np.asarray(tiff.counts)))                       # the tile that compresses least
    tiff.file.seek(int(tiff.offsets[k]))
    data = tiff.file.read(int(tiff.counts[k]))
    # The same bytes declared as a strip of 8-bit samples in a TIFF of our own, for libtiff to decode.
    w, h = tiff.tile_width * tiff.dtype.itemsize, tiff.tile_length
    entries = [(256, 4, w), (257, 4, h), (258, 3, 8), (259, 3, 5), (262, 3, 1), (273, 4, 0), (277, 3, 1), (278, 4, h),
               (279, 4, len(data))]
    head = 8 + 2 + 12 * len(entries) + 4
    raw = bytearray(b"II*\0" + (8).to_bytes(4, "little") + len(entries).to_bytes(2, "little"))
    for tag, kind, v in entries:
        v = head if tag == 273 else v
        raw += tag.to_bytes(2, "little") + kind.to_bytes(2, "little") + (1).to_bytes(4, "little") + v.to_bytes(4, "little")
    raw += b"\0\0\0\0" + data
    ours = pg.lzw_decode(data)
    assert_("tiff: our LZW decodes GHS-POP's densest tile as libtiff does", ours == Image.open(io.BytesIO(bytes(raw))).tobytes(),
            f"{len(ours)} bytes")
    g = tiff.grid()
    check("tiff: GHS-POP's cell is 30″", g["dlon"] * 120 - 1, 1e-8)
    east = g["west"] + g["nlon"] * g["dlon"]
    assert_("tiff: GHS-POP's columns reach past each antimeridian, by about one",
            -180 - 1.5 * g["dlon"] < g["west"] < -180 and 180 < east < 180 + 1.5 * g["dlon"], f"{g['west']:.5f} to {east:.5f}")
else:
    skip("tiff: the GHS-POP GeoTIFF is not in the cache (python data/acquisition/fetch.py ghsl-population)")


# ---- the population grid ------------------------------------------------------------
meta_path = PRODUCTS / "fields-ghsl-population.json"
if meta_path.exists():
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    grid = meta["grid"]
    codes = np.fromfile(PRODUCTS / meta["file"], dtype="<i2")
    n = grid["nlon"] * grid["nlat"]
    decoded = {}
    for v in meta["variables"]:
        q = codes[v["offsets"][0] // 2: v["offsets"][0] // 2 + n]
        decoded[v["name"]] = np.where(q == -32768, np.nan, q * v["scale"] + v["offset"]).reshape(grid["nlat"], grid["nlon"])
    people, density = decoded["population"], decoded["population_density"]
    total = float(np.nansum(people))
    check(f"population: the world in {meta['epoch']}, {total / 1e9:.4f} billion as stored, within 2 % of the release's 8.1920",
          total / WORLD_2025 - 1, 0.02)
    check("population: the sum before storing is the release's total, to the person",
          meta["aggregation"]["total"] - WORLD_2025, 1.0, "people")
    check("population: and the 30″ cells summed to it, every one counted once",
          meta["aggregation"]["total"] - meta["aggregation"]["source_total"], 1.0, "people")
    areas = pg.cell_areas_km2(grid["nlat"], grid["nlon"])
    check("population: the cells' areas on the sphere sum to 4πR²", areas.sum() * grid["nlon"] / (4 * math.pi * R_KM ** 2) - 1, 1e-12)
    # Density times area is the count, to within the quanta of the two.
    vs = {v["name"]: v for v in meta["variables"]}
    err = np.nanmax(np.abs(density * areas[:, None] - people) - (vs["population_density"]["scale"] / 2 * areas[:, None]
                                                               + vs["population"]["scale"] / 2))
    check("population: density × area is the count, within the two quanta", max(0.0, err), 1e-6, "people")
    assert_("population: no cell is negative", bool(np.nanmin(people) >= 0))
    assert_("population: no times, as the relief; the epoch recorded", meta["times"] == [] and meta["epoch"] == 2025, str(meta["times"]))
else:
    skip("population: fields-ghsl-population missing (python data/products/population_grid.py)")


# ---- the regions --------------------------------------------------------------------
def polygon_area_km2(geometry) -> float:
    """The area of a (multi)polygon on the sphere, by the shoelace formula in the plane of the
    cylindrical equal-area projection, as domains/geography/js/raster.js does (featureAreaFraction)."""
    total = 0.0
    for ring in rg.rings_of(geometry):
        x, y = np.radians(ring[:, 0]), np.sin(np.radians(ring[:, 1]))
        total += 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
    return abs(total) * R_KM ** 2


# The scanline fill agrees with raster.js: a square with a hole on a 1° grid fills 100 - 16 cells.
square = [np.array([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]], float), np.array([[3, 3], [3, 7], [7, 7], [7, 3], [3, 3]], float)]
one = {"lon0": -179.5, "dlon": 1.0, "nlon": 360, "lat0": 89.5, "dlat": 1.0, "nlat": 180}
target = np.zeros((180, 360), dtype="<u2")
assert_("regions: a square with a hole fills 100 - 16 cells, as fillRings does", rg.fill_rings(square, one, 1, target) == 84)
two = np.zeros((180, 360), dtype="<u2")
rg.fill_rings([np.array([[0, 0], [5, 0], [5, 5], [0, 5]], float)], one, 1, two)
rg.fill_rings([np.array([[5, 0], [10, 0], [10, 5], [5, 5]], float)], one, 2, two)
assert_("regions: squares sharing an edge claim each cell once", (two == 1).sum() == 25 and (two == 2).sum() == 25)

for name in ("countries", "admin1"):
    path = PRODUCTS / f"regions-{name}.json"
    if not path.exists():
        skip(f"regions: regions-{name} missing (python data/products/regions.py)")
        continue
    meta = json.loads(path.read_text(encoding="utf-8"))
    g = meta["grid"]
    ids = np.fromfile(PRODUCTS / meta["file"], dtype="<u2")
    listed = meta["regions"]
    assert_(f"regions-{name}: one id per cell, none beyond the list", ids.size == g["nlon"] * g["nlat"] and int(ids.max()) <= len(listed))
    assert_(f"regions-{name}: ids run 1 to {len(listed)} in order, codes unique",
            [r["id"] for r in listed] == list(range(1, len(listed) + 1)) and len({r["code"] for r in listed}) == len(listed))
    counts = np.bincount(ids, minlength=len(listed) + 1)
    assert_(f"regions-{name}: each region's count of cells is that of the raster",
            all(int(counts[r["id"]]) == r["cells"] for r in listed))
    lines = json.loads((PRODUCTS / meta["outlines"]["file"]).read_text(encoding="utf-8"))
    pairs = len(lines["coords"]) // 2
    assert_(f"regions-{name}: outlines of {len(lines['runs'])} polylines, each within the coordinates and of a region",
            all(0 <= f and f + c <= pairs and c >= 2 for f, c in lines["runs"]) and len(lines["region"]) == len(lines["runs"])
            and 1 <= min(lines["region"]) and max(lines["region"]) <= len(listed))
    if name == "countries":
        source = json.loads((ROOT / "data" / "cache" / "natural-earth-50m" / "ne_50m_admin_0_countries.geojson").read_text(encoding="utf-8"))
        features = source["features"]
        worst, label = 0.0, ""
        for r, f in zip(listed, features):
            exact = polygon_area_km2(f["geometry"])
            if exact > 3e5:
                e = r["area_km2"] / exact - 1
                if abs(e) > abs(worst):
                    worst, label = e, r["code"]
        check(f"regions-countries: the raster's area of each country over 300 000 km² is its polygon's (worst {label})", worst, 0.01)
        by = {r["code"]: r for r in listed}
        assert_("regions-countries: ISO codes, with Natural Earth's fall-backs (FRA, NOR from ISO_A3_EH; Kosovo KOS)",
                {"FRA", "NOR", "KOS", "JPN", "AUS"} <= set(by) and by["AUS"]["adm0_a3"] == "AUS")
    else:
        japan = [r for r in listed if r["country"] == "JPN"]
        assert_("regions-admin1: Japan's 47 prefectures, coded JP-01 to JP-47",
                sorted(r["code"] for r in japan) == [f"JP-{k:02d}" for k in range(1, 48)])
        countries = PRODUCTS / "regions-countries.json"
        if countries.exists():
            known = {r["code"] for r in json.loads(countries.read_text(encoding="utf-8"))["regions"]}
            strays = sorted({r["country"] for r in listed} - known)
            assert_("regions-admin1: every division names a country of regions-countries but the territories 1:50m omits",
                    len(strays) <= 10, " ".join(strays))
