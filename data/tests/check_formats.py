"""Checks of the formats in and out of the library: GeoTIFF (data/formats/tiff.py),
the exports of products (data/formats/export.py), the STAC client
(data/acquisition/stac.py) and the pieces of the Earth Atlas's NDVI plug-in
(applications/earth-atlas/plugins/ndvi.py).

Reference values, each with its source:

- LZW codes worked out by hand from TIFF 6.0, section 13, and a reference
  encoder written from the same section, as libtiff's LZWEncode behaves;
- the floating-point predictor's bytes worked out by hand from Adobe's
  Photoshop TIFF Technical Note 3 (2005);
- four small GeoTIFFs written by GDAL 3.12.2 (through rasterio 1.5.2) on
  2026-10-23, embedded below, whose pixels follow formulas recomputed here;
- points of Karney's test set for the transverse Mercator projection (Karney,
  "Test data for the transverse Mercator projection", 2009, doi:10.5281/zenodo.32470,
  CC0; lines as numbered in TMcoords.dat, retrieved 2026-10-23), exact to 0.1 pm;
- the numerical example of Snyder (1987), "Map Projections: A Working Manual",
  USGS Professional Paper 1395, for the transverse Mercator projection of the
  Clarke 1866 ellipsoid (x = 127 106.5 m, y = 4 484 124.4 m).

The network parts (one STAC search and one COG header) are skipped without a network.
"""

from __future__ import annotations

import base64
import datetime as dt
import importlib.util
import io
import json
import math
import sys
import tempfile
import zlib
from pathlib import Path

import numpy as np

from checks import PRODUCTS, ROOT, assert_, check, skip

import export  # noqa: E402
import field_grid  # noqa: E402
import netcdf3  # noqa: E402
import stac  # noqa: E402
import tiff  # noqa: E402

UTC = dt.timezone.utc


def b64(text: str) -> bytes:
    return base64.b64decode("".join(text.split()))


# ---- GeoTIFFs written by GDAL 3.12.2 (rasterio 1.5.2), 2026-10-23 ---------------------------------
# u16_lzw_p2_tiled: uint16 50 x 40, LZW, predictor 2, tiles 16 x 16, EPSG:32654, origin (399960, 4100040), 20 m,
#   nodata 0, one overview (factor 2, nearest); pixel (x, y) = (37x + 101y + (xy mod 17)) mod 4000 + 1000.
# f32_deflate_p3_point: float32 9 x 12, Deflate, predictor 3, strips, EPSG:4326, AREA_OR_POINT=Point, corner
#   (139, 37), 0.25 degree, nodata NaN; pixel = 0.25x - 0.5y + 280, NaN at (4, 3) and (0..2, 7).
# f32_lzw_p3_bigendian_tiled: the same formula, 18 x 20, LZW, predictor 3, big-endian, tiles 16 x 16, EPSG:32733.
# i16_deflate_p2_bigtiff: int16 11 x 10, Deflate, predictor 2, BigTIFF, EPSG:32601; pixel = (113x - 71y) mod 9000 - 4500.
GDAL = {
    "u16_lzw_p2_tiled": """
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""",
    "f32_deflate_p3_point": """
SUkqAAoCAAASAAABAwABAAAACQAAAAEBAwABAAAADAAAAAIBAwABAAAAIAAAAAMBAwABAAAACAAAAAYBAwABAAAAAQAAABEBBAAB
AAAAAAAAABUBAwABAAAAAQAAABYBAwABAAAADAAAABcBBAABAAAAAAAAABwBAwABAAAAAQAAAD0BAwABAAAAAwAAAFMBAwABAAAA
AwAAAA6DDAADAAAA5gAAAIKEDAAGAAAA/gAAAK+HAwAgAAAALgEAALCHDAACAAAAbgEAALGHAgAIAAAAfgEAAIGkAgAEAAAAbmFu
AAAAAAAAAAAAAADQPwAAAAAAANA/AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGBhQAAAAAAAgEJAAAAAAAAA
AAABAAEAAAAHAAAEAAABAAIAAQQAAAEAAQAACAAAAQDmEAEIsYcHAAAABggAAAEAjiMJCLCHAQABAAsIsIcBAAAAiG10lh2kckAA
AABAplRYQVdHUyA4NHwAeJxzZoACTyjNWKwABTAZBmcYw4OBEcIwgalxwKKGAaLqC0xNA7IamyNQNaangaq2AKUTHgCJA1jNAbun
BLd73GHuMcXtHneYe75iuKeegeEImFHLwHAKrKqUgWEBWA2me9xh7inF7R43mHvMcLvHDeaeb1jDB64GrGobTA3cPQBvkSMBEgAA
AQMAAQAAAAkAAAABAQMAAQAAAAwAAAACAQMAAQAAACAAAAADAQMAAQAAAAgAAAAGAQMAAQAAAAEAAAARAQQAAQAAAIYBAAAVAQMA
AQAAAAEAAAAWAQMAAQAAAAwAAAAXAQQAAQAAAIQAAAAcAQMAAQAAAAEAAAA9AQMAAQAAAAMAAABTAQMAAQAAAAMAAAAOgwwAAwAA
AEgDAACChAwABgAAAGADAACvhwMAIAAAAOgCAACwhwwAAwAAACgDAACxhwIACAAAAEADAACBpAIABAAAAG5hbgAAAAAAAQABAAAA
BwAABAAAAQACAAEEAAABAAIAAAgAAAEA5hABCLGHBwAAAAYIAAABAI4jCQiwhwEAAgALCLCHAQABAAAAAAAAAAAAiG10lh2kckAA
AABAplRYQVdHUyA4NHwAAAAAAAAA0D8AAAAAAADQPwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABkYUAAAAAA
AHBCQAAAAAAAAAAA""",
    "f32_lzw_p3_bigendian_tiled": """
TU0AKgAAAAgAEQEAAAMAAAABABIAAAEBAAMAAAABABQAAAECAAMAAAABACAAAAEDAAMAAAABAAUAAAEGAAMAAAABAAEAAAEVAAMA
AAABAAEAAAEcAAMAAAABAAEAAAE9AAMAAAABAAMAAAFCAAMAAAABABAAAAFDAAMAAAABABAAAAFEAAQAAAAEAAAA4gFFAAMAAAAE
AAAA2gFTAAMAAAABAAMAAIMOAAwAAAADAAAA8oSCAAwAAAAGAAABCoevAAMAAAAgAAABOoexAAIAAAAeAAABegAAAAAA2gCmAGkA
XQAAAZgAAAJyAAADGAAAA4FAPgAAAAAAAEA+AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBHoSAAAAAAEFW
42AAAAAAAAAAAAAAAAAAAQABAAAABwQAAAAAAQABBAEAAAABAAEEAoexABYAAAgBh7EABwAWCAYAAAABI44MAAAAAAF/3QwEAAAA
ASMpV0dTIDg0IC8gVVRNIHpvbmUgMzNTfFdHUyA4NHwAgBDAEDgkFgxJg0DAMJABzEEPiERiQghkVgUVhJIAELhMcggzichh5gjE
Ji8lgkahUMjzzkUhUEogsnHjEksqAA1ZcegsLWcQMDgl4gcEygknmU4nsMOlDidGgdIlBHjcsgw0p0RklQqUlqkrjsEelZiExqB/
gc2lB9gbKsEGhZ0gagrNFrlQAFfuEMOtkiF4rsYI1VsMFGt+EFbo2BiuDt9Lgb1xFmxd4x0EnmYAC1xF2ytQy+Qgx2xEUu9QIuEv
cFG2IxUyxkM1OPzUDe2TwF42e1va2zt4gICAEMAL0AQWDQeDI4AHKEQ2DCBwQ6JROEQKCRSCo0AHOMQVgCAQR2RQeLSKNRyOoAQG
CRy2Sx2TyIgCBQS2Ry+MTGRRCbSKcRRGAA6SKPyGexifxOg0OUyuj0iByKlzKaU+KUmJVOdxGrRKsQ5FgA60SQV2vVGO2Gx02WWa
G1+G2qqTW3RW0Ri5Vu63aLxhFAA7WSjXuC3CEX/A2zCQbDQfEXPFwWAggBDAEDgkFgxFg0DAMJAB3EEPiERiQghkVgUVhJEAELhM
cgg3ichh5gjEJi8lgkahUMjz3kUhUEogsnmUqgkegsLW8vibgmU/oFBoVDolFo1HpFJpVLplNp1PqFRqVTqlVq1XrFZrVPgIgBDA
C9AEFg0HgyKAB2hENgwgcEOiUThECgkUgqJAB3jEFYAgEEdkUHi0ijUcjqAEBgkctksdk8iIAgUEtm03nE5nU7nk9n0/oFBoVDol
Fo1HpFJpVLplNp1Pp8BA""",
    "i16_deflate_p2_bigtiff": """
SUkrAAgAAAAQAAAAAAAAABAAAAAAAAAAAAEDAAEAAAAAAAAACwAAAAAAAAABAQMAAQAAAAAAAAAKAAAAAAAAAAIBAwABAAAAAAAA
ABAAAAAAAAAAAwEDAAEAAAAAAAAACAAAAAAAAAAGAQMAAQAAAAAAAAABAAAAAAAAABEBEAABAAAAAAAAAAUCAAAAAAAAFQEDAAEA
AAAAAAAAAQAAAAAAAAAWAQMAAQAAAAAAAAAKAAAAAAAAABcBEAABAAAAAAAAADwAAAAAAAAAHAEDAAEAAAAAAAAAAQAAAAAAAAA9
AQMAAQAAAAAAAAACAAAAAAAAAFMBAwABAAAAAAAAAAIAAAAAAAAADoMMAAMAAAAAAAAAYAEAAAAAAACChAwABgAAAAAAAAB4AQAA
AAAAAK+HAwAgAAAAAAAAAKgBAAAAAAAAsYcCAB0AAAAAAAAA6AEAAAAAAAAAAAAAAAAAAAAAAAAAAE5AAAAAAAAATkAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAQwRBAAAAAICELkEAAAAAAAAAAAEAAQAAAAcAAAQAAAEAAQABBAAAAQABAAIE
sYcVAAAAAQixhwcAFQAGCAAAAQCOIwAMAAABAFl/BAwAAAEAKSNXR1MgODQgLyBVVE0gem9uZSAxTnxXR1MgODR8AHicy3lXyIAJ
fQU972KKsgkWMmCK7xfAJlohACLRxQ0FIDSq+Ct+bKKL+WEsZPEYfgQbIS6KVRQAdy84Iw==""",
}


def u16(h, w):
    y, x = np.mgrid[0:h, 0:w]
    return ((x * 37 + y * 101 + (x * y) % 17) % 4000 + 1000).astype("uint16")


def f32(h, w):
    y, x = np.mgrid[0:h, 0:w]
    v = (x * 0.25 - y * 0.5 + 280.0).astype("float32")
    v[3, 4] = np.nan
    v[7, 0:3] = np.nan
    return v


def i16(h, w):
    y, x = np.mgrid[0:h, 0:w]
    return ((x * 113 - y * 71) % 9000 - 4500).astype("int16")


# ---- LZW --------------------------------------------------------------------------------------------


def pack(codes: list[tuple[int, int]]) -> bytes:
    """Codes (value, width) packed most significant bit first, as TIFF LZW stores them."""
    acc, bits, out = 0, 0, bytearray()
    for code, width in codes:
        acc, bits = (acc << width) | code, bits + width
        while bits >= 8:
            bits -= 8
            out.append((acc >> bits) & 0xFF)
    if bits:
        out.append((acc << (8 - bits)) & 0xFF)
    return bytes(out)


def lzw_encode(data: bytes) -> bytes:
    """A reference encoder from TIFF 6.0, section 13, as libtiff's LZWEncode writes:
    the width grows when the next free entry exceeds the largest code of the width,
    and the table is cleared when the next free entry would be 4094."""
    codes, table, free, width = [(256, 9)], {bytes([i]): i for i in range(256)}, 258, 9
    w = b""

    def grow():
        nonlocal table, free, width
        free += 1
        if free == 4094:
            codes.append((256, width))
            table, free, width = {bytes([i]): i for i in range(256)}, 258, 9
        elif free > (1 << width) - 1:
            width += 1

    for byte in data:
        wc = w + bytes([byte])
        if wc in table:
            w = wc
            continue
        codes.append((table[w], width))
        table[wc] = free
        grow()
        w = bytes([byte])
    if w:
        codes.append((table[w], width))
        grow()
    codes.append((257, width))
    return pack(codes)


def check_lzw() -> None:
    # "ABABABA": A, B, then AB (258) and ABA (260), each 9 bits; worked out by hand.
    assert_("tiff: LZW decodes ABABABA from the codes 256, 65, 66, 258, 260, 257",
            tiff.lzw_decode(pack([(c, 9) for c in (256, 65, 66, 258, 260, 257)])) == b"ABABABA")
    # "AAAAAAA": codes 258 and 259 arrive before the decoder has them (the KwKwK case).
    assert_("tiff: LZW decodes a code that is being defined (AAAAAAA from 256, 65, 258, 259, 65, 257)",
            tiff.lzw_decode(pack([(c, 9) for c in (256, 65, 258, 259, 65, 257)])) == b"AAAAAAA")
    rng = np.random.default_rng(11)
    noisy = bytes(rng.integers(0, 256, 30000, dtype=np.uint8))       # codes of 9 to 12 bits, and clears
    smooth = bytes((np.arange(60000) // 7 % 251).astype(np.uint8))
    for name, data in (("random bytes (widths 9 to 12, cleared tables)", noisy), ("a repetitive sequence", smooth)):
        assert_(f"tiff: LZW decodes the reference encoding of {len(data)} {name}", tiff.lzw_decode(lzw_encode(data)) == data)
    data = bytes(range(200)) * 50
    assert_("tiff: Deflate decodes a zlib stream (RFC 1950)", tiff.inflate(zlib.compress(data)) == data)


def check_predictors() -> None:
    stored = np.array([[[1], [2], [3], [4]]], dtype=np.uint16)
    assert_("tiff: horizontal differencing undone: 1, 2, 3, 4 -> 1, 3, 6, 10 (TIFF 6.0, section 14)",
            tiff.undo_horizontal(stored)[0, :, 0].tolist() == [1, 3, 6, 10])
    wrapped = np.array([[[250], [10]]], dtype=np.uint8)
    assert_("tiff: horizontal differencing wraps modulo 2^8: 250, 10 -> 250, 4",
            tiff.undo_horizontal(wrapped)[0, :, 0].tolist() == [250, 4])
    a = np.arange(24, dtype=np.int16).reshape(2, 4, 3) * 1000 - 7000
    assert_("tiff: horizontal differencing of three samples per pixel is undone",
            np.array_equal(tiff.undo_horizontal(tiff.horizontal_difference(a)), a))
    # Floats 1.0 and 2.0: bytes 3F 80 00 00 and 40 00 00 00 in planes 3F 40 | 80 00 | 00 00 | 00 00,
    # each stored as its difference from the byte before it: 3F 01 40 80 00 00 00 00 (Technical Note 3).
    by_hand = bytes([0x3F, 0x01, 0x40, 0x80, 0, 0, 0, 0])
    assert_("tiff: the floating-point predictor's bytes of 1.0, 2.0 are 3F 01 40 80 00 00 00 00",
            tiff.floating_point_difference(np.array([[[1.0], [2.0]]], dtype=np.float32)) == by_hand)
    assert_("tiff: the floating-point predictor undone gives 1.0, 2.0",
            tiff.undo_floating_point(by_hand, 1, 2, 1, 4).ravel().tolist() == [1.0, 2.0])


def check_gdal_files() -> None:
    expected = {"u16_lzw_p2_tiled": (u16(40, 50), 32654, (399960.0, 20.0, 4100040.0, 20.0), "area"),
                "f32_deflate_p3_point": (f32(12, 9), 4326, (139.0, 0.25, 37.0, 0.25), "point"),
                "f32_lzw_p3_bigendian_tiled": (f32(20, 18), 32733, (500000.0, 30.0, 6000000.0, 30.0), "area"),
                "i16_deflate_p2_bigtiff": (i16(10, 11), 32601, (166000.0, 60.0, 1000000.0, 60.0), "area")}
    for name, (values, epsg, transform, raster) in expected.items():
        t = tiff.GeoTIFF(io.BytesIO(b64(GDAL[name])))
        im = t.main
        got = im.read()
        assert_(f"tiff: {name} (GDAL): every pixel as the formula gives", np.array_equal(got, values, equal_nan=True),
                f"compression {im.compression}, predictor {im.predictor}, order {t.order}, BigTIFF {t.big}")
        assert_(f"tiff: {name}: EPSG:{epsg}, pixels as {raster}s, corner and pixel size as GDAL wrote them",
                t.epsg == epsg and t.raster_type == raster and np.allclose(t.transform, transform), f"{t.epsg} {t.raster_type} {t.transform}")
        window = im.window(3, 2, 5, 4)
        assert_(f"tiff: {name}: a window of 5 x 4 pixels", np.array_equal(window, values[2:6, 3:8], equal_nan=True))
    t = tiff.GeoTIFF(io.BytesIO(b64(GDAL["u16_lzw_p2_tiled"])))
    levels = t.levels()
    assert_("tiff: the overview of u16_lzw_p2_tiled is found, 25 x 20 pixels of 40 m, its pixels those of GDAL's nearest",
            len(levels) == 2 and levels[1][1:] == (40.0, 40.0) and np.array_equal(levels[1][0].read(), u16(40, 50)[::2, ::2]))
    assert_("tiff: the level for 40 m is the overview, for 30 m the full image", t.level_for(40) == 1 and t.level_for(30) == 0)
    box, transform = t.read_box(400000, 4099300, 400300, 4099700, level=0)
    assert_("tiff: a box in UTM coordinates reads the pixels that cover it",
            np.array_equal(box, u16(40, 50)[17:37, 2:17]) and transform == (400000.0, 20.0, 4099700.0, 20.0), f"{box.shape} {transform}")


def check_writer() -> None:
    grid = {"lon0": -179.5, "dlon": 1.0, "nlon": 360, "lat0": 89.5, "dlat": 1.0, "nlat": 180}
    lats, lons = grid["lat0"] - np.arange(180), grid["lon0"] + np.arange(360)
    values = (288 + 30 * np.cos(np.radians(lats))[:, None] + 5 * np.sin(np.radians(lons))[None, :]).astype(np.float32)
    values[10:20, 30:40] = np.nan
    raw = tiff.encode_grid(values, grid, name="Air temperature at 2 m", units="K", description="check",
                           metadata={"time": "2026-10-22T00:00:00Z"})
    t = tiff.GeoTIFF(io.BytesIO(raw))
    back = t.main.read()
    assert_("tiff: a written GeoTIFF reads back value for value, NaN where NaN",
            np.array_equal(back, values, equal_nan=True))
    keys = t.geokeys
    assert_("tiff: its GeoKeys say geographic (1024 = 2), pixels as areas (1025 = 1), EPSG:4326, degrees (2054 = 9102)",
            keys.get(1024) == 2 and keys.get(1025) == 1 and keys.get(2048) == 4326 and keys.get(2054) == 9102 and keys.get(2049) == "WGS 84", str(keys))
    tags = t.main.tags
    assert_("tiff: the directory is GeoTIFF 1.1 (key revision 1.1), Deflate, the floating-point predictor, IEEE floats",
            tags[34735][:3] == (1, 1, 1) and tags[259] == (8,) and tags[317] == (3,) and tags[339] == (3,))
    assert_("tiff: no value is NaN in GDAL_NODATA, and the name and units are in GDAL_METADATA",
            tags.get(42113) == "nan" and 'role="unittype">K<' in tags.get(42112, "") and "Air temperature" in tags.get(42112, ""))
    assert_("tiff: the corner is half a cell north-west of the first node, and the grid reads back",
            t.transform == (-180.0, 1.0, 90.0, 1.0) and tiff.grid_of(t) == {**grid}, f"{t.transform} {tiff.grid_of(t)}")
    # Values longer than four bytes lie on word boundaries (TIFF 6.0, section 2), and the entries ascend.
    n = int.from_bytes(raw[8:10], "little")
    entries = [(int.from_bytes(raw[10 + 12 * k:12 + 12 * k], "little"), raw[10 + 12 * k:22 + 12 * k]) for k in range(n)]
    offsets = [int.from_bytes(e[8:12], "little") for _, e in entries
               if tiff.TYPES[int.from_bytes(e[2:4], "little")][1] * int.from_bytes(e[4:8], "little") > 4]
    assert_("tiff: the entries ascend by tag and every value outside its entry starts on a word boundary",
            [t_ for t_, _ in entries] == sorted(t_ for t_, _ in entries) and all(o % 2 == 0 for o in offsets))


# ---- the transverse Mercator projection --------------------------------------------------------------

# Karney's test set: line, latitude, longitude (degrees, exact), easting, northing (m, to 0.1 pm);
# central meridian 0, k0 = 0.9996, no false easting or northing.
KARNEY = [
    (23, 40.552052236608, 3.10605643618, 263004.7709889716198, 4493669.7624508701848),
    (57, 84.986930137199, 2.067177857261, 20163.0752214330146, 9438635.9888499043815),
    (69, 19.376879758731, 0.344604702538, 36187.3424147238117, 2142565.3248023055831),
    (81, 33.028801664026, 0.765842683855, 71519.4504575306802, 3654740.4448518867014),
    (182, 19.753508884702, 1.603846135125, 168047.3369803282036, 2185000.046530304967),
    (223, 31.750313985725, 4.2481669819, 402511.3467124300447, 3520621.5619106356485),
    (348, 59.600611432975, 5.335613037397, 300975.2245464249443, 6619034.1530948852613),
    (392, 56.44822215131, 4.46082797782, 274875.4380120068524, 6264892.1473064022883),
    (2, 10.018893710119, 23.313323816796, 2624150.74092867098, 1204434.0416046227631),
    (4, 21.072464821428, 29.82868439326, 3206390.6919962592145, 2650745.4004060203847),
    (6, 70.175453699652, 22.865350231728, 847598.2665137468857, 7947180.9624403351564),
]


def check_projection() -> None:
    lat, lon, east, north = (np.array([k[i] for k in KARNEY]) for i in (1, 2, 3, 4))
    tm = tiff.TransverseMercator(0.0, 0.9996)
    x, y = tm.forward(lon, lat)
    check("tiff: transverse Mercator against Karney's test set, forward, 11 points up to 30° from the meridian",
          float(np.hypot(x - east, y - north).max()) * 1e9, 10, "nm")
    back_lon, back_lat = tm.inverse(east, north)
    metres = np.hypot((back_lat - lat) * 111.2e3, (back_lon - lon) * 111.3e3 * np.cos(np.radians(lat)))
    check("tiff: transverse Mercator against Karney's test set, inverse", float(metres.max()) * 1e9, 10, "nm")
    # UTM zones: zone 31 has its central meridian at 3° E; the southern zones' northings are 10 000 km less the northern.
    near = lon < 6
    x31, y31 = tiff.utm(32631).forward(3 + lon[near], lat[near])
    xs, ys = tiff.utm(32731).forward(3 - lon[near], -lat[near])
    err = max(float(np.abs(x31 - 500000 - east[near]).max()), float(np.abs(y31 - north[near]).max()),
              float(np.abs(500000 - xs - east[near]).max()), float(np.abs(10000000 - ys - north[near]).max()))
    check("tiff: UTM zones 31 N and 31 S (EPSG:32631, 32731) against the test set: false easting, false northing, mirror", err * 1e9, 10, "nm")
    lo, la = tiff.utm(32731).inverse(500000 - east[near], 10000000 - north[near])
    check("tiff: UTM 31 S inverse", float(np.abs(lo - (3 - lon[near])).max() + np.abs(la + lat[near]).max()) * 111e3 * 1e9, 10, "nm")
    clarke = tiff.TransverseMercator(-75.0, 0.9996, a=6378206.4, f=1 - math.sqrt(1 - 0.00676866))
    sx, sy = clarke.forward(-73.5, 40.5)
    check("tiff: Snyder's (1987) example, Clarke 1866, 40°30' N 73°30' W, central meridian 75° W (given to 0.1 m)",
          max(abs(float(sx) - 127106.5), abs(float(sy) - 4484124.4)), 0.05, "m")
    try:
        tiff.utm(32661)
        bad = False
    except tiff.TIFFError:
        bad = True
    assert_("tiff: EPSG:32661 (polar stereographic) is not taken for a UTM zone", bad)


# ---- exports -------------------------------------------------------------------------------------------


def synthetic_product(folder: Path) -> tuple[Path, dict, np.ndarray]:
    """A field product round the Earth, 10° cells, three times, one variable at 500 hPa, written by field_grid.quantise."""
    grid = {"lon0": -180.0, "dlon": 10.0, "nlon": 36, "lat0": 80.0, "dlat": 10.0, "nlat": 17, "periodic": True,
            "crs": "WGS 84 geographic coordinates (EPSG:4326), degrees", "order": "rows north to south, columns eastward"}
    lats = grid["lat0"] - np.arange(grid["nlat"]) * grid["dlat"]
    lons = grid["lon0"] + np.arange(grid["nlon"]) * grid["dlon"]
    times = ["2026-10-22T00:00:00Z", "2026-10-22T06:00:00Z", "2026-10-22T12:00:00Z"]
    x = np.stack([250 + lats[:, None] / 10 + lons[None, :] / 100 + 2 * k for k in range(3)])
    x[1, 5, 7] = x[1, 6, 0] = np.nan
    q, scale, offset = field_grid.quantise(x)
    (folder / "fields-check.bin").write_bytes(b"".join(q[k].tobytes() for k in range(3)))
    n = grid["nlat"] * grid["nlon"] * 2
    meta = {"schema": "field-grid/1", "name": "check", "title": "A check", "domain": "meteorology", "grid": grid, "times": times,
            "variables": [{"name": "t@500hPa", "long_name": "Temperature at 500 hPa", "standard_name": "air_temperature",
                           "units": "K", "level": "500 hPa", "scale": scale, "offset": offset, "offsets": [0, n, 2 * n]}],
            "file": "fields-check.bin", "source": {"provider": "the checks", "credit": "made up"}}
    path = folder / "fields-check.json"
    path.write_text(json.dumps(meta), encoding="utf-8")
    return path, meta, x


def check_exports() -> None:
    with tempfile.TemporaryDirectory() as folder:
        path, meta, x = synthetic_product(Path(folder))
        scale = meta["variables"][0]["scale"]
        content, ctype, name = export.field(path, "t@500hPa", None, "nc", bbox=(170, -30, -170, 30))
        (Path(folder) / name).write_bytes(content)
        ds = netcdf3.Dataset.open(Path(folder) / name)
        lon, lat, t = ds.array("lon"), ds.array("lat"), ds.array("time")
        values = ds.array("t_500hPa")
        a = ds.variables["t_500hPa"].attributes
        assert_("export: netCDF of a box across the antimeridian: longitudes 170, 180, 190 and latitudes 30 to -30",
                lon.tolist() == [170.0, 180.0, 190.0] and lat.tolist() == [30.0, 20.0, 10.0, 0.0, -10.0, -20.0, -30.0], f"{lon} {lat}")
        expected_t = [(dt.datetime.fromisoformat(s.replace("Z", "+00:00")) - dt.datetime(1970, 1, 1, tzinfo=UTC)).total_seconds()
                      for s in meta["times"]]
        assert_("export: netCDF times in seconds since 1970-01-01 00:00:00 UTC",
                t.tolist() == expected_t and ds.variables["time"].attributes.get("units") == "seconds since 1970-01-01 00:00:00 UTC")
        assert_("export: netCDF units, standard_name, long_name, and the level as plev = 50000 Pa",
                a.get("units") == "K" and a.get("standard_name") == "air_temperature" and a.get("long_name") == "Temperature at 500 hPa"
                and ds.array("plev").tolist() == [50000.0] and ds.variables["plev"].attributes.get("units") == "Pa"
                and ds.attributes.get("Conventions") == "CF-1.11", str(a))
        cols = [35, 0, 1]                                    # 170, 180 (= -180) and 190 (= -170) degrees
        truth = x[:, 5:12][:, :, cols]
        err = np.nanmax(np.abs(values[:, 0] - truth))
        check("export: netCDF values, unpacked by netcdf3.py, against the values quantised (at most half a step)", err / scale, 0.5001, "steps")
        assert_("export: netCDF keeps no value as no value (20° N, 180° at the second time)",
                bool(np.isnan(values[1, 0, 1, 1])) and int(np.isnan(values).sum()) == 1)

        tif, ctype, name = export.field(path, "t@500hPa", 1, "tif")
        g = tiff.GeoTIFF(io.BytesIO(tif))
        got = g.main.read()
        product = export.Product(path)
        decoded = product.values_at(product.variable("t@500hPa"), 1).astype(np.float32)
        assert_("export: GeoTIFF of time 1 holds the product's values in float32, NaN where it has none",
                np.array_equal(got, decoded, equal_nan=True))
        check("export: those values against the values quantised (half a step and the rounding to float32)",
              float(np.nanmax(np.abs(got - x[1]))) / scale, 0.5 + float(np.spacing(np.float32(260))) / scale, "steps")
        assert_("export: GeoTIFF corner at (-185, 85), EPSG:4326, named for its time",
                g.transform == (-185.0, 10.0, 85.0, 10.0) and g.epsg == 4326 and name == "check-t_500hPa-20261022T0600Z.tif"
                and ctype.startswith("image/tiff"), f"{g.transform} {name}")
        try:
            export.field(path, "t@500hPa", None, "tif")
            refused = False
        except ValueError:
            refused = True
        assert_("export: a GeoTIFF of every time of a product of three is refused (ValueError)", refused)
        csv, ctype, _ = export.field(path, "t@500hPa", 1, "csv", stride=2)
        lines = csv.decode("utf-8").splitlines()
        assert_("export: CSV with stride 2: a header lon,lat,t@500hPa (K) and 9 x 18 lines",
                lines[0] == "lon,lat,t@500hPa (K)" and len(lines) == 1 + 9 * 18 and ctype.startswith("text/csv"), f"{lines[0]} {len(lines)}")
        csv1, _, _ = export.field(path, "t@500hPa", 1, "csv")
        row = csv1.decode("utf-8").splitlines()[1 + 5 * 36 + 7]
        assert_("export: CSV leaves the value of a point without one empty", row == "-110.000000,30.000000,", row)
        js, ctype, _ = export.field(path, "t@500hPa", None, "json", bbox=(-20, -10, 20, 10))
        d = json.loads(js)
        assert_("export: JSON of every time: values[3][3][5], null where there is no value, the grid of the box",
                len(d["values"]) == 3 and len(d["values"][0]) == 3 and len(d["values"][0][0]) == 5
                and d["grid"]["lon0"] == -20.0 and d["grid"]["nlon"] == 5 and d["times"] == meta["times"], json.dumps(d["grid"]))
        js1, _, _ = export.field(path, "t@500hPa", 1, "json")
        assert_("export: JSON writes null for no value", json.loads(js1)["values"][0][5][7] is None)

        # Tracks across the antimeridian are cut there (RFC 7946, section 3.1.9).
        tracks = {"schema": "tracks/1", "name": "check", "fields": ["time_ms", "lon", "lat", "observed"],
                  "tracks": [{"id": "x", "name": "X", "points": [[0, 178.0, 10.0, 1], [21600000, -178.0, 12.0, 1], [43200000, -176.0, 13.0, 0]]}]}
        tp = Path(folder) / "tracks-check.json"
        tp.write_text(json.dumps(tracks), encoding="utf-8")
        gj = json.loads(export.features(tp)[0])
        geom = gj["features"][0]["geometry"]
        assert_("export: a track across the antimeridian becomes a MultiLineString cut at ±180° at latitude 11°",
                geom["type"] == "MultiLineString" and geom["coordinates"][0][-1] == [180.0, 11.0] and geom["coordinates"][1][0] == [-180.0, 11.0]
                and gj["features"][0]["properties"]["times"][1] == "1970-01-01T06:00:00Z", json.dumps(geom))


def valid_geojson(collection: dict) -> tuple[bool, str]:
    """The structure RFC 7946 requires of a FeatureCollection of Points, LineStrings and MultiLineStrings."""
    if collection.get("type") != "FeatureCollection" or not isinstance(collection.get("features"), list):
        return False, "not a FeatureCollection"

    def position(p):
        return (isinstance(p, list) and len(p) >= 2 and all(isinstance(c, (int, float)) for c in p)
                and -180 <= p[0] <= 180 and -90 <= p[1] <= 90)

    for f in collection["features"]:
        g = f.get("geometry") or {}
        if f.get("type") != "Feature" or not isinstance(f.get("properties"), dict):
            return False, f"not a Feature: {str(f)[:80]}"
        c = g.get("coordinates")
        ok = {"Point": lambda: position(c), "LineString": lambda: len(c) >= 2 and all(map(position, c)),
              "MultiLineString": lambda: all(len(l) >= 2 and all(map(position, l)) for l in c)}.get(g.get("type"), lambda: False)()
        if not ok:
            return False, f"bad geometry {str(g)[:120]}"
    return True, ""


def check_features() -> None:
    found = sorted(PRODUCTS.glob("events-*.json")) + sorted(PRODUCTS.glob("tracks-*.json")) + sorted(PRODUCTS.glob("observations-*.json"))
    kinds = set()
    for path in found:
        meta = json.loads(path.read_text(encoding="utf-8"))
        content, ctype, name = export.features(path)
        collection = json.loads(content)
        ok, why = valid_geojson(collection)
        rows = len(meta.get("events", meta.get("rows", meta.get("tracks", []))))
        left = (collection.get("left_out") or {}).get("rows", 0)
        assert_(f"export: {name}: valid GeoJSON (RFC 7946), {len(collection['features'])} features of {rows} rows"
                + (f", {left} without a valid position left out" if left else ""),
                ok and len(collection["features"]) + left == rows and left <= 0.01 * rows and ctype == "application/geo+json", why)
        kinds.add(meta["schema"].split("/")[0])
    for kind in ("events", "tracks", "observations"):
        if kind not in kinds:
            skip(f"export: no {kind} product in the cache to write as GeoJSON")


# ---- STAC -------------------------------------------------------------------------------------------------


def check_stac() -> None:
    start, end = dt.datetime(2026, 9, 23, tzinfo=UTC), dt.datetime(2026, 10, 23, tzinfo=UTC)
    body = stac.search_body(["sentinel-2-l2a"], intersects=stac.point(140.0, 36.1), start=start, end=end,
                            query={"eo:cloud_cover": {"lt": 20}}, sortby=[("eo:cloud_cover", "asc")], limit=10)
    assert_("stac: an Item Search body: collections, a GeoJSON point, an RFC 3339 interval, a query and a sort",
            body == {"collections": ["sentinel-2-l2a"], "limit": 10, "intersects": {"type": "Point", "coordinates": [140.0, 36.1]},
                     "datetime": "2026-09-23T00:00:00Z/2026-10-23T00:00:00Z", "query": {"eo:cloud_cover": {"lt": 20}},
                     "sortby": [{"field": "properties.eo:cloud_cover", "direction": "asc"}]}, json.dumps(body))
    assert_("stac: an open interval is written with ..", stac.interval(None, end) == "../2026-10-23T00:00:00Z")
    page = {"links": [{"rel": "next", "method": "POST", "href": "https://x/search", "body": {"next": "abc"}, "merge": True}]}
    method, href, following = stac.next_request(page, body)
    assert_("stac: a next link to merge keeps the search and adds its token", method == "POST" and following == {**body, "next": "abc"})
    page["links"][0]["merge"] = False
    assert_("stac: a next link not to merge replaces the body", stac.next_request(page, body)[2] == {"next": "abc"})
    assert_("stac: GET paging and the last page", stac.next_request({"links": [{"rel": "next", "href": "https://x?page=2"}]}, body)
            == ("GET", "https://x?page=2", None) and stac.next_request({"links": []}, body) is None)
    square = {"geometry": {"type": "Polygon", "coordinates": [[[139, 35], [141, 35], [141, 37], [139, 37], [139, 35]]]}}
    assert_("stac: a footprint contains a place inside it and not one outside",
            stac.contains(square["geometry"], 140, 36) and not stac.contains(square["geometry"], 141.5, 36))
    assert_("stac: a footprint covers a box inside it and not one across its edge",
            stac.covers(square, 139.5, 35.5, 140.5, 36.5) and not stac.covers(square, 140.5, 35.5, 141.5, 36.5))
    try:
        items = stac.search(stac.EARTH_SEARCH, {**body, "limit": 1, "query": {"eo:cloud_cover": {"lt": 100}}}, max_items=1, timeout=30)
    except (OSError, RuntimeError) as error:
        skip(f"stac: no live search ({type(error).__name__}: {str(error)[:80]})")
        return
    assert_("stac: a live search of Earth Search finds a Sentinel-2 L2A scene over 140.0° E, 36.1° N in 30 days",
            len(items) == 1 and stac.contains(items[0]["geometry"], 140.0, 36.1), str([i.get("id") for i in items]))
    if not items:
        return
    item = items[0]
    band = stac.raster_band(item, "red")
    assert_("stac: the scene has red, near-infrared and classification bands, the red with a scale",
            all(k in item["assets"] for k in ("red", "nir", "scl")) and "scale" in band, str(band))
    try:
        scl = tiff.GeoTIFF.url(stac.asset_href(item, "scl"))
    except OSError as error:
        skip(f"tiff: no live COG header ({error})")
        return
    a = item["assets"]["scl"]
    pt = a.get("proj:transform") or []
    assert_("tiff: the live COG's GeoKeys and tie point agree with the item's proj:epsg and proj:transform",
            scl.epsg == item["properties"].get("proj:epsg") and len(pt) >= 6 and np.allclose(scl.transform, (pt[2], pt[0], pt[5], -pt[4]))
            and len(scl.images) > 1, f"{scl.epsg} {scl.transform} {pt} in {scl.source.requests} requests")


# ---- the NDVI plug-in's pieces -----------------------------------------------------------------------------


def check_ndvi_plugin() -> None:
    path = ROOT / "applications" / "earth-atlas" / "plugins" / "ndvi.py"
    spec = importlib.util.spec_from_file_location("ndvi_check", path)
    ndvi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ndvi)
    grid = ndvi.box_grid(140.0, 36.1, 20000, 20)
    a, f = tiff.WGS84_A, tiff.WGS84_F
    e2, phi = f * (2 - f), math.radians(36.1)
    meridian = a * (1 - e2) / (1 - e2 * math.sin(phi) ** 2) ** 1.5        # the meridian's radius of curvature
    normal = a / math.sqrt(1 - e2 * math.sin(phi) ** 2)                    # and the prime vertical's
    check("ndvi: the grid's rows are 20 m apart along the meridian",
          abs(math.radians(grid["dlat"]) * meridian - 20), 0.01, "m")
    check("ndvi: the grid's columns are 20 m apart along the parallel",
          abs(math.radians(grid["dlon"]) * normal * math.cos(phi) - 20), 0.01, "m")
    assert_("ndvi: the box is about 20 km, centred on the place, not round the Earth",
            grid["nlon"] == grid["nlat"] == 1001 and abs(grid["lon0"] + 500 * grid["dlon"] - 140.0) < 1e-9 and not grid["periodic"], str(grid))
    window = np.add.outer(3.0 * np.arange(10), 2.0 * np.arange(12))         # value 2c + 3r at pixel centres
    transform = (1000.0, 10.0, 5000.0, 10.0)
    xs = 1000 + np.array([12.0, 55.0, 101.0, -1.0])
    ys = 5000 - np.array([17.0, 64.0, 33.0, 5.0])
    lin = ndvi.at_nodes(window, transform, xs, ys)
    expected = 2 * ((xs - 1000) / 10 - 0.5) + 3 * ((5000 - ys) / 10 - 0.5)
    assert_("ndvi: bilinear values of a linear field are exact, and NaN outside the window",
            np.allclose(lin[:3], expected[:3]) and math.isnan(lin[3]), f"{lin} {expected}")
    near = ndvi.at_nodes(window, transform, xs[:3], ys[:3], nearest=True)
    assert_("ndvi: the nearest value is that of the pixel containing the place", near.tolist() == [2 * 1 + 3 * 1, 2 * 5 + 3 * 6, 2 * 10 + 3 * 3])

    def item(flag, offset, baseline):
        return {"id": "x", "properties": {"earthsearch:boa_offset_applied": flag, "s2:processing_baseline": baseline},
                "assets": {"red": {"raster:bands": [{"nodata": 0, "scale": 0.0001, "offset": offset}]}}}
    dark = np.linspace(80, 3000, 1000)
    bright = np.linspace(1300, 4000, 1000)
    rho, said, _ = ndvi.reflectance(dark, item(True, -0.1, "05.13"), "red")
    assert_("ndvi: an offset Earth Search says it has removed is not applied again", np.allclose(rho, dark * 1e-4), said)
    rho, said, _ = ndvi.reflectance(dark, item(False, -0.1, "04.00"), "red")
    assert_("ndvi: a declared offset that would make the darkest pixels below -0.05 is not applied, and the note says so",
            np.allclose(rho, dark * 1e-4) and "not in the data" in said, said)
    rho, said, darkest = ndvi.reflectance(bright, item(False, -0.1, "04.00"), "red")
    assert_("ndvi: a declared offset the pixels allow is applied", np.allclose(rho, bright * 1e-4 - 0.1) and abs(darkest - 0.0327) < 1e-3, f"{said}; {darkest}")
    rho, said, _ = ndvi.reflectance(np.array([0.0, 500.0]), item(False, 0.0, "03.01"), "red")
    assert_("ndvi: before baseline 04.00 no offset; DN 0 is no value", math.isnan(rho[0]) and rho[1] == 0.05, said)
    try:
        sys.path.insert(0, str(ROOT / "api" / "python"))
        from science import plugins as api
    except ImportError as error:
        skip(f"ndvi: the plug-in classes are not importable here ({error})")
        return
    plugin = ndvi.create(api)
    results = plugin.run({}, {}, api.Context(time=dt.datetime.now(UTC), pinned=None))
    assert_("ndvi: without a chosen place the plug-in says so in a note and returns nothing else",
            len(results) == 1 and results[0]["kind"] == "note" and "place" in results[0]["text"])


def main() -> int:
    sys.stdout.reconfigure(encoding="utf-8")       # names carry ρ and −; run-tests.mjs reads UTF-8
    check_lzw()
    check_predictors()
    check_gdal_files()
    check_writer()
    check_projection()
    check_exports()
    check_features()
    check_stac()
    check_ndvi_plugin()
    return 0


if __name__ == "__main__":
    sys.exit(main())
