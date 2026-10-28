"""Checks of the data library's tiles (data/tiles): the encoding of a terrain tile, the pyramid, the
sources' rules (the resolution behind a Terrarium sample, land and sea) and the synthetic detail.

What needs the cache (the land mask, built tiles) is skipped when the cache does not hold it; the
rest is checked from closed forms and from tiles made here.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

from checks import ROOT, assert_, check, skip

sys.path.insert(0, str(ROOT / "data"))
from tiles.pyramid import BORDER, GRID, TileKey, bounds, cell_degrees, sample_lonlat  # noqa: E402
from tiles.sources import SEA, LandMask, Terrarium  # noqa: E402
from tiles.terrain import TerrainService, decode, encode  # noqa: E402
from tiles import synthetic  # noqa: E402

CACHE = ROOT / "data" / "cache"
m = GRID + 2 * BORDER

# ---- the encoding ----------------------------------------------------------------------------------
rng = np.random.default_rng(7)
heights = rng.normal(0, 1000, (m, m)).astype("<f4")
provenance = rng.integers(1, 7, (m, m)).astype(np.uint8)
floor = heights - rng.uniform(0, 500, (m, m)).astype("<f4")
for odd in ("", "x"):  # a header whose length is not a multiple of 4 is padded
    header = {"z": 3, "x": 5, "y": 2, "n": GRID, "border": BORDER, "note": odd}
    data = encode(header, heights, provenance, floor)
    h, hh, pp, ff = decode(data)
    assert_(f"tiles: a tile with its floor reads back as written (header {len(odd)})",
            np.array_equal(hh, heights) and np.array_equal(pp, provenance) and np.array_equal(ff, floor) and h["layers"] == ["heights", "provenance", "floor"])
    length = int.from_bytes(data[4:8], "little")
    assert_("tiles: the heights start on a multiple of 4 bytes", (8 + length) % 4 == 0)
h, hh, pp, ff = decode(encode({"z": 0, "x": 0, "y": 0, "n": GRID, "border": BORDER}, heights, provenance))
assert_("tiles: a tile without a floor reads back without one", ff is None and h["layers"] == ["heights", "provenance"])

# ---- the pyramid -----------------------------------------------------------------------------------
assert_("tiles: level 0 is two tiles of 180 degrees", bounds(TileKey(0, 0, 0)) == (-180.0, -90.0, 0.0, 90.0) and bounds(TileKey(0, 1, 0)) == (0.0, -90.0, 180.0, 90.0))
w, s, e, n = bounds(TileKey(5, 33, 7))
check("tiles: a level-5 tile spans 180/32 degrees", abs((e - w) - 180 / 32) + abs((n - s) - 180 / 32), 1e-12, "deg")
check("tiles: the cell of level z is 180 / 2^z / 64 degrees", abs(cell_degrees(9) - 180 / 512 / 64), 1e-15, "deg")
lons, lats = sample_lonlat(TileKey(1, 0, 0), GRID, BORDER)
step = 90 / (GRID - 1)
assert_("tiles: samples run west to east, with a border sample beyond each edge",
        len(lons) == m and abs(lons[0] - (-180 - step)) < 1e-9 and abs(lons[-1] - (-90 + step)) < 1e-9)
assert_("tiles: the border sample beyond the north pole is reflected back over it",
        abs(lats[0] - (90 - step)) < 1e-9 and abs(lats[1] - 90) < 1e-9 and max(lats) <= 90)
# Neighbours share their edge samples.
a = sample_lonlat(TileKey(6, 10, 20), GRID, 0)[0]
b = sample_lonlat(TileKey(6, 11, 20), GRID, 0)[0]
check("tiles: neighbours sample the same points along their shared edge", abs(a[-1] - b[0]), 1e-12, "deg")
for bad in (TileKey(-1, 0, 0), TileKey(2, 8, 0), TileKey(2, 0, 4)):
    try:
        bad.check()
        assert_(f"tiles: {bad} is refused", False)
    except ValueError:
        assert_(f"tiles: {bad} is refused", True)

# ---- the sources' rules ------------------------------------------------------------------------------
land, sea = Terrarium.resolution_of(["srtm/N36W113.tif", "gmted/30N120W_20101117_gmted_mea075.tif", "etopo1/ETOPO1_Bed_g.tif"])
assert_("tiles: Terrarium's land is at its finest land source, its floor at ETOPO1", land == 30 and sea == 1850, f"{land}, {sea}")
land, sea = Terrarium.resolution_of(["etopo1/ETOPO1_Bed_g.tif"])
assert_("tiles: a tile of ETOPO1 alone is ETOPO1's everywhere", land == 1850 and sea == 1850, f"{land}, {sea}")
land, sea = Terrarium.resolution_of(["mystery/x.tif"])
assert_("tiles: an unknown source counts as 250 m", land == 250 and sea == 250)

if (CACHE / "products" / "land-mask" / "ne50m-20ppd.png").exists():
    mask = LandMask(CACHE / "products")
    places = {"the Dead Sea's shore": (35.5, 31.5, True), "the Caspian Sea": (51.0, 42.0, True), "the Qattara Depression": (27.5, 29.5, True),
              "Flevoland": (5.5, 52.5, True), "Death Valley": (-116.9, 36.4, True), "the Atlantic": (-30.0, 30.0, False),
              "the North Sea": (3.0, 55.0, False), "the Black Sea": (34.0, 43.0, False), "the Mediterranean": (18.0, 35.0, False)}
    for name, (lon, lat, inland) in places.items():
        got = bool(mask.inland(np.array([lon]), np.array([lat]))[0])
        assert_(f"tiles: {name} is {'inland' if inland else 'sea'} in the land mask", got == inland)
else:
    skip("tiles: the land mask is not built (python data/products/land_mask.py)")

# ---- built tiles, where the cache has them ----------------------------------------------------------------
service = TerrainService(CACHE)
ready = service.available()
if ready["earth"]["ready"]:
    # Level 0 comes from ETOPO alone (no network): the sea at the geoid, its floor beneath.
    h, hh, pp, ff = decode(service.tile("earth", TileKey(0, 0, 0)))
    sea_samples = pp == SEA
    assert_("tiles: the Earth's western hemisphere has sea and land", sea_samples.any() and (~sea_samples).any())
    check("tiles: over the sea, the surface is the geoid (EGM96 spans −106 to +85 m)", float(np.abs(hh[sea_samples]).max()), 110, "m")
    assert_("tiles: over the sea, the floor lies beneath the surface", bool((ff[sea_samples] <= hh[sea_samples] + 1e-3).all()))
    assert_("tiles: on land, the floor is the surface", bool(np.array_equal(ff[~sea_samples], hh[~sea_samples])))
    check("tiles: the deepest floor of the western hemisphere is plausible (Puerto Rico Trench, about 8.4 km)", abs(float(ff.min()) + 8400) / 2000, 1)
else:
    skip("tiles: the Earth's terrain sources are not in the cache")

# ---- synthetic detail -------------------------------------------------------------------------------------
p = np.random.default_rng(3).uniform(-1e5, 1e5, (200000, 3))
v = synthetic.value_noise(p, 0)
check("synthetic: the value noise has unit RMS (NOISE_RMS holds)", abs(float(np.sqrt((v ** 2).mean())) - 1), 0.02)
check("synthetic: the value noise has zero mean", abs(float(v.mean())), 0.02)
assert_("synthetic: the noise is a function of the point (the same twice)", np.array_equal(v[:100], synthetic.value_noise(p[:100], 0)))
assert_("synthetic: other octaves have other noise", not np.allclose(v[:100], synthetic.value_noise(p[:100], 1)))
# The detail at a point does not depend on which tile asks: the samples along a shared edge agree.
lon = np.full(50, 23.5)
lat = np.linspace(0.6, 0.7, 50)
g = np.full(50, 0.05)
d1, waves1 = synthetic.detail("moon", lon, lat, g, 237.0, 10.0)
d2, waves2 = synthetic.detail("moon", lon.copy(), lat.copy(), g, 237.0, 10.0)
assert_("synthetic: two tiles sharing an edge compute the same detail there", np.array_equal(d1, d2) and waves1 == waves2)
assert_("synthetic: octaves from half the data's resolution down to two cells", waves1[0] == 237.0 / 2 and waves1[-1] >= 20.0 and waves1[-1] / 2 < 20.0, str(waves1))
# Its amplitude follows the self-affine law a(λ) = g R^(1−H) λ^H, each octave weighted as it enters.
pts = np.random.default_rng(5).uniform(-60, 60, (100000, 2))
# Cells of 25 m under data of 400 m: octaves of 200 m (8 cells, whole) and 100 m (4 cells, half);
# 50 m is at the Nyquist limit and left out. Independent noises add in quadrature.
single, waves = synthetic.detail("moon", pts[:, 0], pts[:, 1], np.full(len(pts), 0.1), 400.0, 25.0)
a = lambda lam: 0.1 * synthetic.COARSE_SCALE["moon"] ** (1 - synthetic.HURST["moon"]) * lam ** synthetic.HURST["moon"]
expected = math.sqrt(sum((synthetic.weight(lam, 25.0) * a(lam)) ** 2 for lam in waves))
assert_("synthetic: an octave enters at two cells and is whole from eight", waves == [200.0, 100.0]
        and synthetic.weight(200.0, 25.0) == 1.0 and synthetic.weight(100.0, 25.0) == 0.5 and synthetic.weight(50.0, 25.0) == 0.0, str(waves))
check("synthetic: the octaves' RMS is that of w g R^(1-H) λ^H in quadrature", abs(float(np.sqrt((single ** 2).mean())) / expected - 1), 0.03)
assert_("synthetic: none is made where the data are fine enough", synthetic.detail("earth", lon, lat, g, 30.0, 10.0)[1] == [])
