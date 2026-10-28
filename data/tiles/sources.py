"""Sources of elevation for the terrain tiles: where each sample comes from, and how good it is.

Every source here is sampled at the points a tile asks for and answers heights with a
provenance code and the resolution of the data behind each sample, so that a tile knows which of
its samples were measured, which were interpolated by the source's makers, and which lie below
the source's resolution (interpolated here, smoothly, and marked so).

- **Terrarium** (Mapzen's terrain tiles on the AWS Open Data registry): global elevation in Web
  Mercator PNG tiles to zoom 15, compiled from SRTM, EU-DEM, NED, GMTED, ETOPO1 and others; each
  tile names its sources in the header x-amz-meta-x-imagery-sources, from which the resolution of
  its samples is known. Heights are above mean sea level.
- **ETOPO 2022** at 5′ (the data library's product etopo-2022-fine): land and sea floor above the
  EGM2008 geoid, for the coarsest levels and the polar caps beyond Mercator's 85.05°.
- **LOLA** (the Lunar Orbiter Laser Altimeter's gridded data records, PDS): the Moon's radius less
  1737.4 km at 4, 16, 64 and 128 pixels per degree, in the mean-Earth frame; with the count of
  laser shots in each pixel, zero where the grid was filled by interpolation (GMT's surface, with
  tension 0.5). The two coarsest grids are fetched whole; the finer ones are read in blocks by HTTP
  multi-range requests and kept in the cache.
- **EGM96** at 15′ (GeographicLib's copy of NGA's grid): the geoid undulation N, by which heights
  above the geoid become heights above the WGS 84 ellipsoid, h = H + N (bilinear: RMS error 4 cm,
  maximum 1.2 m, by GeographicLib's statement in the file).

- **Land** (the data library's product land-mask, Natural Earth's land at 1:50 million eroded by a
  pixel of 3′): whether a sample below sea level is inland (the Dead Sea's shore, a polder) or sea.

Provenance codes (one byte per sample):
  1 measured        a pixel with laser shots (LOLA)
  2 filled          a pixel the source's makers filled by interpolation (LOLA without shots)
  3 compiled        a source that does not say per sample (Terrarium, ETOPO): measured or filled
  4 finer than data the tile's cells are finer than half the source's resolution: interpolated here
  5 synthetic       generated detail (not produced by these sources)
  6 sea surface     the sea, drawn at the geoid (EGM96): mean sea level without the ocean's mean
                    dynamic topography (±1 m) or its tides; the floor beneath is kept in the tile
"""

from __future__ import annotations

import io
import json
import math
import re
import threading
import time
import http.client
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

AGENT = "science-infrastructures-data/1 (+https://github.com/huangwanhong-maker/The-Science-Infrastructures)"

MEASURED, FILLED, COMPILED, FINER, SYNTHETIC, SEA = 1, 2, 3, 4, 5, 6

# The resolution of Terrarium's sources, metres (Mapzen, "Data Sources", joerd/docs/data-sources.md).
TERRARIUM_RESOLUTION = {
    "srtm": 30, "eudem": 25, "ned": 10, "ned13": 10, "ned_topobathy": 10, "cdem": 20, "linz": 8, "austria": 10,
    "gmted": 250, "etopo1": 1850, "arcticdem": 5, "inegi": 15, "nzdem": 8, "geoscience_australia": 5,
    "uk_lidar": 2, "kartverket": 10, "norway": 10,
}


# Fetches go through one pool of threads, each keeping its connections open (HTTP/1.1 keep-alive):
# a tile's sources are a few files from one host, and a new TLS connection to it costs as much as the
# file (two round trips across the Atlantic for Terrarium's bucket).
FETCH = ThreadPoolExecutor(max_workers=16, thread_name_prefix="tiles-fetch")
_connections = threading.local()


def _connection(scheme: str, host: str, timeout: int) -> http.client.HTTPConnection:
    pool = getattr(_connections, "pool", None)
    if pool is None:
        pool = _connections.pool = {}
    c = pool.get((scheme, host))
    if c is None:
        c = (http.client.HTTPSConnection if scheme == "https" else http.client.HTTPConnection)(host, timeout=timeout)
        pool[(scheme, host)] = c
    return c


def _get(url: str, timeout: int = 60, headers: dict | None = None) -> tuple[bytes, dict]:
    """The body and the headers (lower-case names) of a GET, on a kept connection, with retries."""
    u = urllib.parse.urlsplit(url)
    path = u.path + (f"?{u.query}" if u.query else "")
    last = None
    for attempt in range(3):
        c = _connection(u.scheme, u.netloc, timeout)
        try:
            c.request("GET", path, headers={"User-Agent": AGENT, "Connection": "keep-alive", **(headers or {})})
            r = c.getresponse()
            body = r.read()
            if r.status >= 400:
                raise urllib.error.HTTPError(url, r.status, r.reason, r.headers, None)
            return body, {k.lower(): v for k, v in r.getheaders()}
        except urllib.error.HTTPError as e:
            if e.code < 500:
                raise OSError(f"cannot fetch {url}: HTTP {e.code}") from None
            last = e
        except (OSError, http.client.HTTPException) as e:
            last = e
            c.close()
            _connections.pool.pop((u.scheme, u.netloc), None)
        time.sleep(0.5 * (attempt + 1))
    raise OSError(f"cannot fetch {url}: {last}")


def _store(path: Path, data: bytes) -> None:
    """Writes a file of the cache whole: through a temporary file of this thread's own, since two
    threads may fetch the same file at once; the second to finish keeps the first's."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f".{threading.get_ident()}.part")
    tmp.write_bytes(data)
    try:
        tmp.replace(path)
    except OSError:
        if not path.exists():
            raise
        tmp.unlink(missing_ok=True)


def bilinear(grid: np.ndarray, rows: np.ndarray, cols: np.ndarray, periodic: bool) -> np.ndarray:
    """Bilinear samples of a grid at fractional (row, col) positions; columns wrap if periodic."""
    nr, nc = grid.shape
    rows = np.clip(rows, 0, nr - 1)
    r0 = np.clip(np.floor(rows).astype(int), 0, nr - 2)
    fr = rows - r0
    c0f = np.floor(cols)
    fc = cols - c0f
    c0 = c0f.astype(int)
    if periodic:
        c0, c1 = c0 % nc, (c0 + 1) % nc
    else:
        c0 = np.clip(c0, 0, nc - 2)
        fc = np.clip(cols - c0, 0, 1)
        c1 = c0 + 1
    g = grid
    return ((1 - fr) * ((1 - fc) * g[r0, c0] + fc * g[r0, c1]) + fr * ((1 - fc) * g[r0 + 1, c0] + fc * g[r0 + 1, c1]))


# ---- EGM96 -----------------------------------------------------------------------------------------

class Geoid:
    """EGM96 at 15′ from GeographicLib's PGM (1440 x 721, 90 N to 90 S, from 0 E, point-registered)."""

    def __init__(self, pgm: Path):
        data = pgm.read_bytes()
        tokens, comments, pos = [], {}, 0
        while len(tokens) < 4:  # P5, width, height, maxval; comments between them
            end = data.index(b"\n", pos)
            line = data[pos:end].decode("ascii").strip()
            pos = end + 1
            if line.startswith("#"):
                parts = line[1:].split(None, 1)
                if len(parts) == 2:
                    comments[parts[0]] = parts[1]
                continue
            tokens += line.split()
        if tokens[0] != "P5":
            raise ValueError(f"{pgm} is not a binary PGM")
        width, height = int(tokens[1]), int(tokens[2])
        offset, scale = float(comments["Offset"]), float(comments["Scale"])
        raw = np.frombuffer(data[pos:pos + 2 * width * height], dtype=">u2").reshape(height, width)
        self.n = offset + scale * raw.astype(np.float64)
        self.width, self.height = width, height

    def undulation(self, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
        cols = (np.mod(lons, 360.0)) / 360.0 * self.width
        rows = (90.0 - lats) / 180.0 * (self.height - 1)
        return bilinear(self.n, rows, cols, periodic=True)


# ---- Land and sea --------------------------------------------------------------------------------

class LandMask:
    """The product land-mask (data/products/land_mask.py): 255 inland, 0 elsewhere, pixel-registered
    from 180 W and 90 N."""

    def __init__(self, products: Path):
        from PIL import Image
        meta = json.loads((products / "land-mask" / "ne50m-20ppd.json").read_text(encoding="utf-8"))
        with Image.open(products / "land-mask" / meta["file"]) as image:
            self.grid = np.asarray(image.convert("L")) > 127
        self.height, self.width = self.grid.shape
        self.ppd = self.width / 360.0

    def inland(self, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
        cols = np.floor((np.mod(lons + 180.0, 360.0)) * self.ppd).astype(int) % self.width
        rows = np.clip(np.floor((90.0 - lats) * self.ppd).astype(int), 0, self.height - 1)
        return self.grid[rows, cols]


# ---- ETOPO 2022 ----------------------------------------------------------------------------------

class Etopo:
    """The data library's field product etopo-2022-fine: 5′ cells, rows north to south."""

    def __init__(self, products: Path):
        meta = json.loads((products / "fields-etopo-2022-fine.json").read_text(encoding="utf-8"))
        g = meta["grid"]
        var = meta["variables"][0]
        q = np.fromfile(products / meta["file"], dtype="<i2").reshape(g["nlat"], g["nlon"])
        z = q.astype(np.float64) * var["scale"] + var["offset"]
        z[q == -32768] = np.nan
        self.z, self.g = z, g
        self.resolution_m = 9260.0

    def sample(self, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
        g = self.g
        cols = (lons - g["lon0"]) / g["dlon"]
        rows = (g["lat0"] - lats) / g["dlat"]
        return bilinear(self.z, rows, cols, periodic=True)


# ---- Terrarium -----------------------------------------------------------------------------------

class Terrarium:
    URL = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png"
    LIMIT = 85.0511287798066

    def __init__(self, cache: Path):
        self.cache = cache
        self.lock = threading.Lock()

    def _tile(self, z: int, x: int, y: int) -> tuple[np.ndarray, list[str]]:
        from PIL import Image  # the data library's requirement
        path = self.cache / str(z) / str(x) / f"{y}.png"
        meta = path.with_suffix(".json")
        if not path.exists():
            data, headers = _get(self.URL.format(z=z, x=x, y=y))
            # The sources first, so that a tile in the cache always has them.
            _store(meta, json.dumps({"sources": headers.get("x-amz-meta-x-imagery-sources", ""),
                                     "url": self.URL.format(z=z, x=x, y=y),
                                     "retrieved": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}).encode("utf-8"))
            _store(path, data)
        with Image.open(path) as image:
            rgb = np.asarray(image.convert("RGB"), dtype=np.float64)
        height = rgb[..., 0] * 256 + rgb[..., 1] + rgb[..., 2] / 256 - 32768
        sources = []
        if meta.exists():
            sources = [s.strip() for s in json.loads(meta.read_text(encoding="utf-8")).get("sources", "").split(",") if s.strip()]
        return height, sources

    # Sources of the sea floor only; the rest describe land.
    BATHYMETRY = {"etopo1"}

    @staticmethod
    def resolution_of(sources: list[str]) -> tuple[float, float]:
        """The resolution behind a Terrarium tile's land and its sea floor, metres.

        The header names every source used anywhere in the tile. Mapzen composed them by priority,
        finer sources over coarser, with GMTED beneath as the land's fallback and ETOPO1 for the sea
        floor (Mapzen, joerd/docs/data-sources.md); so land is taken at the finest land source named
        (unknown sources: 250 m), and the floor at ETOPO1's 1.85 km where it is named, else as land."""
        land = [float(TERRARIUM_RESOLUTION.get(s.split("/")[0].lower(), 250)) for s in sources
                if s.split("/")[0].lower() not in Terrarium.BATHYMETRY]
        sea = [float(TERRARIUM_RESOLUTION[s.split("/")[0].lower()]) for s in sources
               if s.split("/")[0].lower() in Terrarium.BATHYMETRY]
        land_res = min(land) if land else (max(sea) if sea else 250.0)
        return land_res, (max(sea) if sea else land_res)

    def sample(self, lons: np.ndarray, lats: np.ndarray, zoom: int):
        """Heights, the resolution behind each sample (m) if it is land and if it is sea floor (the
        caller knows which), and the sources used, at points within Mercator's latitudes."""
        n = 2 ** zoom
        lat = np.radians(np.clip(lats, -self.LIMIT, self.LIMIT))
        fx = (lons + 180.0) / 360.0 * n * 256 - 0.5
        fy = (1 - np.log(np.tan(lat) + 1 / np.cos(lat)) / math.pi) / 2 * n * 256 - 0.5
        x0, x1 = int(np.floor(fx.min() / 256)), int(np.floor((fx.max() + 1) / 256))
        y0, y1 = max(0, int(np.floor(fy.min() / 256))), min(n - 1, int(np.floor((fy.max() + 1) / 256)))
        keys = [(tx % n, ty) for tx in range(x0, x1 + 1) for ty in range(y0, y1 + 1)]
        got = dict(zip(keys, FETCH.map(lambda k: self._tile(zoom, k[0], k[1]), keys)))
        width, height = (x1 - x0 + 1) * 256, (y1 - y0 + 1) * 256
        mosaic = np.zeros((height, width))
        land_res = np.zeros((height, width))
        sea_res = np.zeros((height, width))
        sources: set[str] = set()
        for tx in range(x0, x1 + 1):
            for ty in range(y0, y1 + 1):
                h, src = got[(tx % n, ty)]
                r0, c0 = (ty - y0) * 256, (tx - x0) * 256
                mosaic[r0:r0 + 256, c0:c0 + 256] = h
                land_res[r0:r0 + 256, c0:c0 + 256], sea_res[r0:r0 + 256, c0:c0 + 256] = self.resolution_of(src)
                sources.update(src)
        rows, cols = fy - y0 * 256, fx - x0 * 256
        values = bilinear(mosaic, rows, cols, periodic=False)
        at = (np.clip(np.round(rows).astype(int), 0, height - 1), np.clip(np.round(cols).astype(int), 0, width - 1))
        return values, (land_res[at], sea_res[at]), sorted(sources)


# ---- LOLA ----------------------------------------------------------------------------------------

class Lola:
    """LOLA's gridded data records: LDEM (heights, 0.5 m units, relative to 1737.4 km) and LDEC
    (laser shots per pixel), simple cylindrical, pixel-registered, 0 to 360 E, 90 N to 90 S."""

    BASE = "https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/"
    RESOLUTIONS = (4, 16, 64, 128)  # pixels per degree
    BLOCK_ROWS, BLOCK_COLS = 128, 512
    RADIUS = 1737400.0

    def __init__(self, cache: Path):
        self.cache = cache
        self.lock = threading.Lock()
        self.whole: dict[str, np.ndarray] = {}

    def _whole(self, name: str, ppd: int) -> np.ndarray:
        if name not in self.whole:
            path = self.cache / f"{name}.img"
            if not path.exists():
                raise FileNotFoundError(f"{path} is missing: python data/acquisition/fetch.py lola-gdr")
            self.whole[name] = np.fromfile(path, dtype="<i2").reshape(180 * ppd, 360 * ppd)
        return self.whole[name]

    def _block(self, name: str, ppd: int, br: int, bc: int) -> np.ndarray:
        """A block of a remote grid, by one multi-range request, kept in the cache."""
        path = self.cache / "blocks" / name / f"{br}_{bc}.i16"
        rows_total, cols_total = 180 * ppd, 360 * ppd
        r0, c0 = br * self.BLOCK_ROWS, bc * self.BLOCK_COLS
        nr, nc = min(self.BLOCK_ROWS, rows_total - r0), min(self.BLOCK_COLS, cols_total - c0)
        if path.exists():
            return np.fromfile(path, dtype="<i2").reshape(nr, nc)
        ranges = [((r0 + i) * cols_total * 2 + c0 * 2, (r0 + i) * cols_total * 2 + (c0 + nc) * 2 - 1) for i in range(nr)]
        out = np.zeros((nr, nc), dtype="<i2")
        url = self.BASE + f"{name}.img"
        for start in range(0, nr, 64):  # at most 64 ranges a request
            part = ranges[start:start + 64]
            header = "bytes=" + ",".join(f"{a}-{b}" for a, b in part)
            body, headers = _get(url, timeout=120, headers={"Range": header})
            ctype = headers.get("content-type", "")
            pieces = _multipart(body, ctype) if "multipart/byteranges" in ctype else [(part[0][0], body)]
            for offset, data in pieces:
                i = offset // (cols_total * 2) - r0
                out[i, :] = np.frombuffer(data, dtype="<i2")[:nc]
        _store(path, out.tobytes())
        return out

    def _grid_window(self, kind: str, ppd: int, rows: np.ndarray, cols: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The part of a grid around the sample positions, with the positions relative to it."""
        name = f"{kind}_{ppd}"
        rows_total, cols_total = 180 * ppd, 360 * ppd
        ra, rb = int(max(0, np.floor(rows.min()) - 1)), int(min(rows_total - 1, np.ceil(rows.max()) + 1))
        ca, cb = int(np.floor(cols.min()) - 1), int(np.ceil(cols.max()) + 1)
        window = np.zeros((rb - ra + 1, cb - ca + 1), dtype=np.float64)
        wanted = [(br, bcu) for br in range(ra // self.BLOCK_ROWS, rb // self.BLOCK_ROWS + 1)
                  for bcu in range(ca // self.BLOCK_COLS, cb // self.BLOCK_COLS + 1)]
        nblocks = cols_total // self.BLOCK_COLS
        # The blocks not yet in the cache are fetched together, each by its own request.
        blocks = list(FETCH.map(lambda k: self._block(name, ppd, k[0], k[1] % nblocks), wanted))
        for (br, bcu), block in zip(wanted, blocks):
            # Place the block's overlap with the window.
            gr0, gc0 = br * self.BLOCK_ROWS, bcu * self.BLOCK_COLS
            for i in range(block.shape[0]):
                r = gr0 + i
                if r < ra or r > rb:
                    continue
                j0, j1 = max(ca, gc0), min(cb, gc0 + block.shape[1] - 1)
                if j0 > j1:
                    continue
                window[r - ra, j0 - ca:j1 - ca + 1] = block[i, j0 - gc0:j1 - gc0 + 1]
        return window, rows - ra, cols - ca

    def sample(self, lons: np.ndarray, lats: np.ndarray, ppd: int):
        """Radius less 1737.4 km (m) and whether each sample's nearest pixel has laser shots."""
        east = np.mod(lons, 360.0)
        if east.max() - east.min() > 180:  # a tile across 0° E: keep its columns contiguous
            east = np.where(east < 180, east + 360, east)
        cols = east * ppd - 0.5
        rows = (90.0 - lats) * ppd - 0.5
        if ppd <= 16:
            grid, counts = self._whole(f"ldem_{ppd}", ppd), self._whole(f"ldec_{ppd}", ppd)
            heights = bilinear(grid.astype(np.float64), rows, cols, periodic=True) * 0.5
            ri = np.clip(np.round(rows).astype(int), 0, grid.shape[0] - 1)
            ci = np.round(cols).astype(int) % grid.shape[1]
            return heights, counts[ri, ci] > 0
        with ThreadPoolExecutor(max_workers=2) as pool:
            heights_job = pool.submit(self._grid_window, "ldem", ppd, rows, cols)
            counts_job = pool.submit(self._grid_window, "ldec", ppd, rows, cols)
            (grid, r, c), (counts, _, _) = heights_job.result(), counts_job.result()
        heights = bilinear(grid, r, c, periodic=False) * 0.5
        ri = np.clip(np.round(r).astype(int), 0, counts.shape[0] - 1)
        ci = np.clip(np.round(c).astype(int), 0, counts.shape[1] - 1)
        return heights, counts[ri, ci] > 0


def _multipart(body: bytes, content_type: str) -> list[tuple[int, bytes]]:
    """The parts of a multipart/byteranges answer, each as (first byte offset, bytes)."""
    m = re.search(r'boundary="?([^";]+)"?', content_type)
    if not m:
        raise ValueError("multipart answer without a boundary")
    boundary = ("--" + m.group(1)).encode()
    out = []
    for part in body.split(boundary):
        if b"Content-Range" not in part[:400] and b"content-range" not in part[:400]:
            continue
        head, _, data = part.partition(b"\r\n\r\n")
        rng = re.search(rb"bytes (\d+)-(\d+)/", head, re.I)
        if not rng:
            continue
        a, b = int(rng.group(1)), int(rng.group(2))
        out.append((a, data[: b - a + 1]))
    return out
