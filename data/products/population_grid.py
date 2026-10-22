"""Build the population product: the people of the Earth on a grid of 0.25 degrees, from GHSL.

GHS-POP R2023A (Schiavina et al. 2023) gives the residential population of each
cell of a 30-arc-second grid on WGS 84, at epochs from 1975 to 2030. Each
country's total is aligned to the UN World Population Prospects 2022, and its
people are shared among the cells by their residential built-up volume seen by
Sentinel-2 and Landsat (Pesaresi et al. 2024; GHSL Data Package 2023, sec. 2.5).

A count per cell is an extensive quantity: it is aggregated by summing, never
by interpolating, so that the people of any union of cells are the sum of
their counts. Each 30-second cell is given whole to the 0.25-degree cell that
contains its centre. The source grid is not aligned with ours (its corner is at
-180.0079, 89.0996 degrees, the trace of its warping from the Mollweide
projection) and runs one column past each antimeridian; a column centred
beyond +-180 degrees is wrapped onto the other side. A boundary of our grid is
therefore placed to within half a source cell (15 seconds, under 0.5 km), and
the total is kept exactly.

Density is the count divided by the cell's area on a sphere of radius
R = 6371.0088 km, the mean radius of the WGS 84 ellipsoid (Moritz 2000):

    A = R^2 (lambda_2 - lambda_1) (sin phi_2 - sin phi_1).

On the ellipsoid, a cell between the same parallels (of geodetic latitude) is
0.45 % smaller at the equator, 0.2 % larger at 45 degrees and 0.7 % larger at
70 degrees: the error of the density so computed.

The source is a BigTIFF of 7.4 GB of float64 values, in 256 x 256 tiles each
compressed by LZW. It is read tile by tile by the reader below, which needs
only the standard library and numpy.

Output, in data/cache/products/: fields-ghsl-population.json and .bin, a field
product (field_grid.py) without times, as the relief is: a map of where people
live changes over years, not over the days of a forecast, and is shown at any
instant. Its epoch, the year of the UN's mid-year estimates to which it is
aligned, is recorded as `epoch` and in the variables' names. Two variables:

population
    people in the cell (CF cell_methods "area: sum"), units 1;
population_density
    people per square kilometre of the cell (cell_methods "area: mean"), units km-2.

Both are stored as 16-bit integers on a linear scale, as every field product
is: the quantum of the count is its maximum / 65 000 (254 people, the maximum
being 16.5 million, in the cell of central Karachi), and that of the density
0.36 per km², so that a cell with fewer than half a quantum reads as empty. The
stored total is compared with the exact one when the product is built
(8 191 211 825 against 8 191 988 536 for 2025: 0.01 % fewer).

Usage::

    python data/products/population_grid.py                  # 0.25 degrees
    python data/products/population_grid.py --step 0.5
"""

from __future__ import annotations

import argparse
import json
import math
import struct
import sys
import time
import zlib
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from field_grid import MISSING, quantise  # noqa: E402
from product_index import OUT, register  # noqa: E402

SOURCE = "ghsl-population"
EARTH_RADIUS_KM = 6371.0088          # mean radius R1 of GRS 80 and WGS 84 (Moritz 2000)


# ---- a minimal reader of TIFF and BigTIFF --------------------------------------
#
# TIFF Revision 6.0 (Adobe, 1992): the header (sec. 2), image file directories
# and field types (sec. 2), tiles (sec. 15), LZW compression (sec. 13) and the
# horizontal predictor (sec. 14); BigTIFF (Aware Systems, 2008), the same with
# 64-bit offsets and counts; Deflate (compression 8 and 32946) as zlib streams
# (RFC 1950), as libtiff writes them. One sample per pixel; images in strips are
# read as tiles the width of the image. GeoTIFF's tags (OGC 19-008r4) give the
# grid: the model tie point and pixel scale, and the raster type.

FIELD_TYPES = {1: "B", 2: "s", 3: "H", 4: "I", 5: "II", 6: "b", 7: "B", 8: "h", 9: "i", 10: "ii",
               11: "f", 12: "d", 16: "Q", 17: "q", 18: "Q"}
SAMPLE_DTYPES = {(1, 8): "u1", (1, 16): "u2", (1, 32): "u4", (2, 8): "i1", (2, 16): "i2", (2, 32): "i4",
                 (3, 32): "f4", (3, 64): "f8"}
WIDTH, LENGTH, BITS, COMPRESSION, SAMPLES, PREDICTOR, SAMPLE_FORMAT = 256, 257, 258, 259, 277, 317, 339
STRIP_OFFSETS, ROWS_PER_STRIP, STRIP_COUNTS = 273, 278, 279
TILE_WIDTH, TILE_LENGTH, TILE_OFFSETS, TILE_COUNTS = 322, 323, 324, 325
PIXEL_SCALE, TIE_POINT, GEO_KEYS, NODATA = 33550, 33922, 34735, 42113
RASTER_TYPE, GEOGRAPHIC_TYPE, PIXEL_IS_AREA, WGS84 = 1025, 2048, 1, 4326


def lzw_decode(data: bytes) -> bytes:
    """Decodes one TIFF LZW strip or tile (TIFF 6.0, sec. 13; Welch 1984).

    Codes are read most significant bit first, from 9 bits wide, growing to 10,
    11 and 12 bits one code before the table would need them ("early change"):
    when the next free entry is 511, 1023 or 2047. Code 256 clears the table and
    257 ends the data. A code not yet in the table (the KwKwK case) stands for
    the previous string followed by its own first byte.
    """
    out = bytearray()
    table = [bytes([i]) for i in range(256)] + [b"", b""]
    width, position, end = 9, 0, len(data) * 8
    data = data + b"\0\0\0"
    previous = None
    while position + width <= end:
        byte = position >> 3
        code = (int.from_bytes(data[byte:byte + 3], "big") >> (24 - (position & 7) - width)) & ((1 << width) - 1)
        position += width
        if code == 256:
            del table[258:]
            width, previous = 9, None
            continue
        if code == 257:
            break
        if previous is None:
            entry = table[code]
        elif code < len(table):
            entry = table[code]
            table.append(previous + entry[:1])
        else:
            entry = previous + previous[:1]
            table.append(entry)
        out += entry
        previous = entry
        if len(table) == (1 << width) - 1 and width < 12:
            width += 1
    return bytes(out)


class Tiff:
    """The first image of a TIFF or BigTIFF file of one sample per pixel, read tile by tile."""

    def __init__(self, path: Path):
        self.file = open(path, "rb")
        head = self.file.read(16)
        self.order = {b"II": "<", b"MM": ">"}[head[:2]]
        magic = struct.unpack(self.order + "H", head[2:4])[0]
        self.big = magic == 43
        if magic not in (42, 43):
            raise ValueError(f"{path}: not a TIFF file")
        first = struct.unpack(self.order + ("Q" if self.big else "I"), head[8:16] if self.big else head[4:8])[0]
        self.tags = self._directory(first)
        t = self.tags
        self.width, self.height = int(t[WIDTH][0]), int(t[LENGTH][0])
        if int(t.get(SAMPLES, [1])[0]) != 1:
            raise ValueError("only images of one sample per pixel are read")
        self.compression = int(t.get(COMPRESSION, [1])[0])
        self.predictor = int(t.get(PREDICTOR, [1])[0])
        key = (int(t.get(SAMPLE_FORMAT, [1])[0]), int(t[BITS][0]))
        self.dtype = np.dtype(self.order + SAMPLE_DTYPES[key])
        if TILE_WIDTH in t:
            self.tile_width, self.tile_length = int(t[TILE_WIDTH][0]), int(t[TILE_LENGTH][0])
            self.offsets, self.counts = t[TILE_OFFSETS], t[TILE_COUNTS]
        else:
            self.tile_width, self.tile_length = self.width, int(t.get(ROWS_PER_STRIP, [self.height])[0])
            self.offsets, self.counts = t[STRIP_OFFSETS], t[STRIP_COUNTS]
        self.across = -(-self.width // self.tile_width)
        self.down = -(-self.height // self.tile_length)

    def _directory(self, offset: int) -> dict:
        f, o = self.file, self.order
        f.seek(offset)
        n = struct.unpack(o + ("Q" if self.big else "H"), f.read(8 if self.big else 2))[0]
        entry, inline = (20, 8) if self.big else (12, 4)
        raw = f.read(n * entry)
        tags = {}
        for k in range(n):
            e = raw[k * entry:(k + 1) * entry]
            tag, kind = struct.unpack(o + "HH", e[:4])
            count = struct.unpack(o + ("Q" if self.big else "I"), e[4:4 + inline])[0]
            code = FIELD_TYPES.get(kind)
            if code is None:
                continue
            size = struct.calcsize(code) * count
            value = e[4 + inline:]
            if size > inline:
                f.seek(struct.unpack(o + ("Q" if self.big else "I"), value)[0])
                value = f.read(size)
            if kind == 2:
                tags[tag] = value[:size].rstrip(b"\0").decode("latin-1")
            else:
                tags[tag] = struct.unpack(o + code * count, value[:size])
        return tags

    def tile(self, k: int) -> np.ndarray:
        """Tile k (row-major over the tiles) as a (tile_length, tile_width) array."""
        self.file.seek(int(self.offsets[k]))
        data = self.file.read(int(self.counts[k]))
        if self.compression == 5:
            data = lzw_decode(data)
        elif self.compression in (8, 32946):
            data = zlib.decompress(data)
        elif self.compression != 1:
            raise ValueError(f"compression {self.compression} is not read")
        size = self.tile_width * self.tile_length * self.dtype.itemsize
        if len(data) < size:                  # the last strip of an image in strips is shorter
            data = data + bytes(size - len(data))
        values = np.frombuffer(data, self.dtype, self.tile_width * self.tile_length)
        values = values.reshape(self.tile_length, self.tile_width)
        if self.predictor == 2:       # horizontal differencing of integers (TIFF 6.0, sec. 14)
            values = np.cumsum(values, axis=1, dtype=self.dtype)
        elif self.predictor != 1:
            raise ValueError(f"predictor {self.predictor} is not read")
        return values

    def grid(self) -> dict:
        """The corner of the image and its cell size, from the GeoTIFF tags, on WGS 84."""
        keys = self.tags.get(GEO_KEYS, ())
        geo = {keys[4 + 4 * i]: keys[4 + 4 * i + 3] for i in range(keys[3] if keys else 0) if keys[4 + 4 * i + 1] == 0}
        if geo.get(GEOGRAPHIC_TYPE) != WGS84:
            raise ValueError(f"expected geographic coordinates on WGS 84 (EPSG:4326), found GeoKeys {geo}")
        sx, sy = self.tags[PIXEL_SCALE][:2]
        i, j, _, x, y, _ = self.tags[TIE_POINT][:6]
        corner = geo.get(RASTER_TYPE, PIXEL_IS_AREA) == PIXEL_IS_AREA
        # The corner of cell (0, 0): a tie point to a cell's centre (PixelIsPoint) is moved by half a cell.
        west = x - i * sx - (0 if corner else sx / 2)
        north = y + j * sy + (0 if corner else sy / 2)
        nodata = self.tags.get(NODATA)
        return {"west": west, "north": north, "dlon": sx, "dlat": sy, "nlon": self.width, "nlat": self.height,
                "nodata": float(nodata) if nodata else None}


# ---- aggregation --------------------------------------------------------------

def aggregate(tiff: Tiff, step: float, progress=None) -> tuple[np.ndarray, dict]:
    """Sums of the source cells over a global grid of `step` degrees, each cell by its centre.

    Returns (sums, report): sums (nlat, nlon) in float64, rows north to south from
    90 degrees, columns east from -180; the report gives the source total and grid.
    """
    g = tiff.grid()
    nlon, nlat = round(360 / step), round(180 / step)
    lon = g["west"] + (np.arange(tiff.width) + 0.5) * g["dlon"]
    lat = g["north"] - (np.arange(tiff.height) + 0.5) * g["dlat"]
    column = np.floor((lon + 180) / step).astype(np.int64)           # unwrapped: -1 and nlon wrap
    row = np.clip(np.floor((90 - lat) / step).astype(np.int64), 0, nlat - 1)
    sums = np.zeros((nlat, nlon))
    total, negative, blank = 0.0, 0, {}
    tw, tl = tiff.tile_width, tiff.tile_length
    started = time.time()
    for k in range(tiff.across * tiff.down):
        tr, tc = divmod(k, tiff.across)
        r0, c0 = tr * tl, tc * tw
        r1, c1 = min(r0 + tl, tiff.height), min(c0 + tw, tiff.width)
        count = int(tiff.counts[k])
        if count < 4096 and count in blank:                           # a tile already seen to be empty
            tiff.file.seek(int(tiff.offsets[k]))
            if tiff.file.read(count) == blank[count]:
                continue
        v = tiff.tile(k)[: r1 - r0, : c1 - c0]
        if g["nodata"] is not None:
            v = np.where(v == g["nodata"], 0.0, v)
        if not v.any():
            if count < 4096:
                tiff.file.seek(int(tiff.offsets[k]))
                blank[count] = tiff.file.read(count)
            continue
        negative += int((v < 0).sum())
        v = np.where(v > 0, v, 0.0).astype(np.float64)
        total += float(v.sum())
        # Sum runs of source columns, then of rows, that fall in one cell of ours.
        cols, rows = column[c0:c1], row[r0:r1]
        cstart = np.flatnonzero(np.r_[True, cols[1:] != cols[:-1]])
        rstart = np.flatnonzero(np.r_[True, rows[1:] != rows[:-1]])
        block = np.add.reduceat(np.add.reduceat(v, cstart, axis=1), rstart, axis=0)
        np.add.at(sums, (rows[rstart][:, None], cols[cstart][None, :] % nlon), block)
        if progress and k % 400 == 0:
            progress(k / (tiff.across * tiff.down), time.time() - started)
    return sums, {"source_total": total, "negative_cells": negative, "source_grid": g, "seconds": time.time() - started}


def cell_areas_km2(nlat: int, nlon: int) -> np.ndarray:
    """Areas of the cells of each row of a global grid on the sphere of radius R (km^2)."""
    edges = np.radians(90 - np.arange(nlat + 1) * (180 / nlat))
    return EARTH_RADIUS_KM ** 2 * (2 * math.pi / nlon) * (np.sin(edges[:-1]) - np.sin(edges[1:]))


def build(step: float, name: str) -> dict:
    manifest = json.loads((DATA / "registry" / f"{SOURCE}.json").read_text(encoding="utf-8"))
    folder = DATA / "cache" / SOURCE
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    entry = manifest["files"][0]
    tiff = Tiff(folder / entry["file"])

    def progress(fraction, seconds):
        print(f"\r  {fraction:6.1%} of the tiles, {seconds:5.0f} s", end="", flush=True)

    sums, report = aggregate(tiff, step, progress)
    print()
    nlat, nlon = sums.shape
    density = sums / cell_areas_km2(nlat, nlon)[:, None]
    grid = {"lon0": round(-180 + step / 2, 9), "dlon": step, "nlon": nlon, "lat0": round(90 - step / 2, 9),
            "dlat": step, "nlat": nlat, "periodic": True, "crs": "WGS 84 geographic coordinates (EPSG:4326), degrees",
            "order": "rows north to south, columns eastward"}
    epoch = int(manifest["epoch"])
    source = {"dataset": manifest["id"], "provider": manifest["provider"], "credit": manifest.get("credit"),
              "licence": manifest["licence"], "citation": manifest["citation"], "release": manifest.get("release"),
              "epoch": epoch, "retrieved_at": retrieval["retrieved_at"], "file": entry["file"],
              "request": entry["url"], "member": entry.get("member")}
    variables, blocks, offset_bytes = [], [], 0
    stored = {}
    for vname, long_name, units, method, x in (
            ("population", f"Residential population of the cell, {epoch}", "1", "area: sum", sums),
            ("population_density", f"Residential population per km² of the cell, {epoch}", "km-2", "area: mean", density)):
        q, scale, offset = quantise(x[None].astype(np.float64))
        decoded = np.where(q == MISSING, np.nan, q.astype(np.float64) * scale + offset)[0]
        stored[vname] = decoded
        finite = x[np.isfinite(x)]
        variables.append({
            "name": vname, "long_name": long_name, "standard_name": vname, "units": units, "cell_methods": method,
            "scale": scale, "offset": offset,
            "min": float(finite.min()), "max": float(finite.max()),
            "p02": float(np.percentile(finite, 2)), "p98": float(np.percentile(finite, 98)),
            "coverage": round(float(finite.size / x.size), 4), "offsets": [offset_bytes],
        })
        blocks.append(q[0].tobytes())
        offset_bytes += q[0].nbytes
    exact, kept = float(sums.sum()), float(np.nansum(stored["population"]))
    summary = {"total": round(exact), "total_stored": round(kept), "source_total": round(report["source_total"]),
               "radius_km": EARTH_RADIUS_KM, "source_grid": {k: report["source_grid"][k] for k in ("west", "north", "dlon", "dlat", "nlon", "nlat")},
               "assignment": "each source cell whole to the cell containing its centre; columns beyond ±180° wrapped"}
    bin_path = OUT / f"fields-{name}.bin"
    OUT.mkdir(parents=True, exist_ok=True)
    bin_path.write_bytes(b"".join(blocks))
    product = {
        "schema": "field-grid/1", "name": name, "title": f"Population, GHS-POP {manifest.get('release')}, {epoch}",
        "domain": "geography", "model": manifest.get("model", "GHSL"), "source": source, "grid": grid,
        "times": [], "time_scale": "UTC", "epoch": epoch, "variables": variables, "vectors": [],
        "aggregation": summary, "file": bin_path.name,
        "encoding": "int16 little-endian; value = q * scale + offset; q = -32768 means no value",
    }
    json_path = OUT / f"fields-{name}.json"
    json_path.write_text(json.dumps(product, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    register("fields", {
        "name": name, "file": json_path.name, "title": product["title"], "domain": "geography", "model": product["model"],
        "climatology": None,
        "variables": [{k: v.get(k) for k in ("name", "long_name", "standard_name", "units", "level")} for v in variables],
        "vectors": [], "times": {"count": 0, "first": None, "last": None}, "epoch": epoch,
        "grid": {k: grid[k] for k in ("dlon", "dlat", "nlon", "nlat")},
        "source": {k: source.get(k) for k in ("dataset", "provider", "credit", "retrieved_at")},
    })
    return product | {"report": report}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step", type=float, default=0.25, help="cell size in degrees, dividing 180 (default 0.25)")
    parser.add_argument("--name", default="ghsl-population", help="product name (default: ghsl-population)")
    args = parser.parse_args(argv)
    if not (DATA / "cache" / SOURCE / "retrieval.json").exists():
        print(f"GHS-POP is not in the cache; run: python data/acquisition/fetch.py {SOURCE}", file=sys.stderr)
        return 1
    if abs(180 / args.step - round(180 / args.step)) > 1e-9:
        parser.error("the step must divide 180 degrees")
    p = build(args.step, args.name)
    a = p["aggregation"]
    print(f"fields-{args.name}: {p['grid']['nlon']} x {p['grid']['nlat']}, epoch {p['epoch']}; "
          f"{a['total'] / 1e9:.4f} billion people ({a['total_stored'] / 1e9:.4f} as stored); "
          f"{p['report']['seconds']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
