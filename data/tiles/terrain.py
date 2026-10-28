"""Terrain tiles of the Earth and the Moon, built on demand from their sources and kept in the cache.

A terrain tile (z, x, y) of the geographic pyramid (pyramid.py) is a grid of 65 x 65 heights at
the corners of its 64 x 64 cells, with a provenance code for each (sources.py) and a header that
says what the heights are and where they came from. The encoding, read by the shared
visualization layer (visualization/js/globe/terrain.js):

    'CPT1'                     4 bytes
    header length L            uint32, little-endian
    header                     L bytes of JSON (UTF-8), padded with spaces to a multiple of 4
    heights                    m x m float32, little-endian, rows north to south, metres
    provenance                 m x m uint8, the same order
    floor                      m x m float32 after zero bytes to a multiple of 4, when the header's
                               "layers" names it (the Earth)

with m = n + 2 border: the tile's n x n samples and a ring of `border` samples beyond its edges at
the same spacing, from which the page computes the normals at the edges, so that neighbouring
tiles shade alike across their seam. The header holds the body, the tile, its bounds (degrees),
n, border, the vertical datum, the sources and their resolutions, the smallest and largest height
of the tile proper, and the counts of each provenance code.

- **The Earth.** Heights above the WGS 84 ellipsoid: the orthometric heights of the sources plus
  the EGM96 geoid undulation. The heights drawn are the surface one sees: over the sea, the sea's
  surface at the geoid (provenance 6), where the sources' height is below sea level outside the
  land (the product land-mask); the solid surface beneath, the sea floor, is the layer "floor". The coarsest levels (cells of 18 km or more) and the polar caps
  beyond 85.05° come from ETOPO 2022 at 5′ (above EGM2008; the two geoids differ by a metre or
  two, well below these levels' cells); the rest from Terrarium at the zoom whose pixels are no
  larger than the tile's cells.
- **The Moon.** Heights above the sphere of 1737.4 km in the mean-Earth frame: LOLA's radius less
  1737.4 km, from the grid whose pixels are at most two of the tile's cells (4, 16, 64 or 128
  pixels per degree).
"""

from __future__ import annotations

import json
import math
import threading
import time
from pathlib import Path

import numpy as np

from .pyramid import BORDER, GRID, TileKey, bounds, cell_degrees, sample_lonlat
from .sources import COMPILED, FILLED, FINER, MEASURED, SEA, Etopo, Geoid, LandMask, Lola, Terrarium

EARTH_RADIUS = 6378137.0
MOON_RADIUS = 1737400.0
VERSION = 3


def encode(header: dict, heights: np.ndarray, provenance: np.ndarray, floor: np.ndarray | None = None) -> bytes:
    header = {**header, "layers": ["heights", "provenance"] + (["floor"] if floor is not None else [])}
    text = json.dumps(header, separators=(",", ":")).encode("utf-8")
    text += b" " * ((4 - len(text) % 4) % 4)
    out = b"CPT1" + len(text).to_bytes(4, "little") + text + heights.astype("<f4").tobytes() + provenance.astype(np.uint8).tobytes()
    if floor is not None:
        # Aligned to 4 bytes, for the page's Float32Array.
        out += bytes((4 - len(out) % 4) % 4) + floor.astype("<f4").tobytes()
    return out


def decode(data: bytes) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray | None]:
    """(header, heights, provenance, floor or None)."""
    if data[:4] != b"CPT1":
        raise ValueError("not a terrain tile")
    length = int.from_bytes(data[4:8], "little")
    header = json.loads(data[8:8 + length].decode("utf-8"))
    m = header["n"] + 2 * header.get("border", 0)
    start = 8 + length
    heights = np.frombuffer(data[start:start + 4 * m * m], dtype="<f4").reshape(m, m)
    end = start + 5 * m * m
    provenance = np.frombuffer(data[start + 4 * m * m:end], dtype=np.uint8).reshape(m, m)
    floor = None
    if "floor" in header.get("layers", []):
        at = end + (4 - end % 4) % 4
        floor = np.frombuffer(data[at:at + 4 * m * m], dtype="<f4").reshape(m, m)
    return header, heights, provenance, floor


class TerrainService:
    """Terrain tiles of a body, built once and then read from data/cache/tiles/terrain/."""

    def __init__(self, cache: Path):
        self.cache = cache
        self.out = cache / "tiles" / "terrain"
        self._lock = threading.Lock()
        self._building: dict[tuple, threading.Event] = {}
        self._etopo = self._terrarium = self._geoid = self._lola = self._land = None
        self._opening = threading.Lock()  # a source is opened once, however many threads ask at first

    # ---- sources, opened on first use ----------------------------------------------------------
    def etopo(self) -> Etopo:
        if self._etopo is None:
            with self._opening:
                if self._etopo is None:
                    self._etopo = Etopo(self.cache / "products")
        return self._etopo

    def terrarium(self) -> Terrarium:
        if self._terrarium is None:
            with self._opening:
                if self._terrarium is None:
                    self._terrarium = Terrarium(self.cache / "aws-terrain-tiles" / "terrarium")
        return self._terrarium

    def geoid(self) -> Geoid:
        if self._geoid is None:
            with self._opening:
                if self._geoid is None:
                    self._geoid = Geoid(self.cache / "products" / "geoid" / "egm96-15.pgm")
        return self._geoid

    def land(self) -> LandMask:
        if self._land is None:
            with self._opening:
                if self._land is None:
                    self._land = LandMask(self.cache / "products")
        return self._land

    def lola(self) -> Lola:
        if self._lola is None:
            with self._opening:
                if self._lola is None:
                    self._lola = Lola(self.cache / "lola-gdr")
        return self._lola

    def available(self) -> dict:
        """What each body's tiles need, and whether it is in the cache."""
        need = {
            "earth": {"etopo-2022 (product etopo-2022-fine)": (self.cache / "products" / "fields-etopo-2022-fine.json").exists(),
                      "egm96-geoid (product geoid/egm96-15.pgm)": (self.cache / "products" / "geoid" / "egm96-15.pgm").exists(),
                      "natural-earth-50m (product land-mask)": (self.cache / "products" / "land-mask" / "ne50m-20ppd.png").exists()},
            "moon": {"lola-gdr (ldem_4, ldec_4, ldem_16, ldec_16)": all((self.cache / "lola-gdr" / f).exists() for f in ("ldem_4.img", "ldec_4.img", "ldem_16.img", "ldec_16.img"))},
        }
        return {body: {"ready": all(v.values()), "needs": v} for body, v in need.items()}

    # ---- tiles ----------------------------------------------------------------------------------
    def path(self, body: str, key: TileKey) -> Path:
        return self.out / f"{body}-v{VERSION}" / str(key.z) / str(key.x) / f"{key.y}.cpt"

    def tile(self, body: str, key: TileKey) -> bytes:
        key.check()
        if body not in ("earth", "moon"):
            raise ValueError(f"no terrain for {body}")
        path = self.path(body, key)
        if path.exists():
            return path.read_bytes()
        ident = (body, key)
        with self._lock:
            event = self._building.get(ident)
            owner = event is None
            if owner:
                event = self._building[ident] = threading.Event()
        if not owner:
            event.wait(300)
            if path.exists():
                return path.read_bytes()
        try:
            data = self._build_earth(key) if body == "earth" else self._build_moon(key)
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(f".{threading.get_ident()}.part")
            tmp.write_bytes(data)
            tmp.replace(path)
            return data
        finally:
            if owner:
                with self._lock:
                    self._building.pop(ident, None)
                event.set()

    def _grid(self, key: TileKey):
        lons, lats = sample_lonlat(key, GRID, BORDER)
        lon, lat = np.meshgrid(np.array(lons), np.array(lats))
        return lon, lat

    def _header(self, body: str, key: TileKey, heights: np.ndarray, provenance: np.ndarray, datum: str, sources: list, resolution: tuple) -> dict:
        west, south, east, north = bounds(key)
        inner = slice(BORDER, BORDER + GRID)
        provenance, heights = provenance[inner, inner], heights[inner, inner]
        counts = {str(code): int((provenance == code).sum()) for code in np.unique(provenance)}
        return {"version": VERSION, "body": body, "z": key.z, "x": key.x, "y": key.y, "n": GRID, "border": BORDER,
                "bounds": [west, south, east, north], "datum": datum, "sources": sources,
                "resolution_m": [float(resolution[0]), float(resolution[1])],
                "min": float(np.nanmin(heights)), "max": float(np.nanmax(heights)), "provenance": counts,
                "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}

    def _build_earth(self, key: TileKey) -> bytes:
        lon, lat = self._grid(key)
        cell_m = math.radians(cell_degrees(key.z)) * EARTH_RADIUS
        heights = np.zeros_like(lon)
        land_res = np.full_like(lon, self.etopo().resolution_m)
        sea_res = land_res.copy()
        sources: list[str] = []
        polar = np.abs(lat) > Terrarium.LIMIT
        use_etopo = cell_m >= 18000 or polar.all()
        if use_etopo:
            heights = self.etopo().sample(lon.ravel(), lat.ravel()).reshape(lon.shape)
            sources = ["ETOPO 2022, 5′ (NOAA NCEI), above EGM2008"]
        else:
            zoom = int(min(15, max(0, math.ceil(math.log2(40075016.686 / (256 * cell_m))))))
            inside = ~polar
            values, (res_land, res_sea), src = self.terrarium().sample(lon[inside], lat[inside], zoom)
            heights[inside] = values
            land_res[inside], sea_res[inside] = res_land, res_sea
            sources = [f"Terrarium z{zoom}: " + ", ".join(src)] if src else [f"Terrarium z{zoom}"]
            if polar.any():
                heights[polar] = self.etopo().sample(lon[polar], lat[polar])
                sources.append("ETOPO 2022, 5′, beyond 85.05°")
        heights = np.nan_to_num(heights, nan=0.0)
        # The sea: below sea level and not inland.
        sea = (heights < 0) & ~self.land().inland(lon.ravel(), lat.ravel()).reshape(lon.shape)
        # Orthometric to ellipsoidal heights: h = H + N (EGM96); the sea's surface is the geoid.
        undulation = self.geoid().undulation(lon.ravel(), lat.ravel()).reshape(lon.shape)
        floor = heights + undulation
        surface = np.where(sea, undulation, floor)
        # Below sea level inland is land, and its resolution the land's.
        resolution = np.where((heights < 0) & sea, sea_res, land_res)
        provenance = np.full(lon.shape, COMPILED, dtype=np.uint8)
        provenance[cell_m < resolution / 2] = FINER
        provenance[sea] = SEA
        datum = ("WGS 84 ellipsoidal height (m): the sources' orthometric heights + the EGM96 geoid; over the sea its surface at "
                 "the geoid, and the floor beneath in the layer floor")
        header = self._header("earth", key, surface, provenance, datum, sources, (resolution.min(), resolution.max()))
        inner = slice(BORDER, BORDER + GRID)
        header["floor_min"], header["floor_max"] = float(floor[inner, inner].min()), float(floor[inner, inner].max())
        header["cell_m"] = cell_m
        return encode(header, surface, provenance, floor)

    def _build_moon(self, key: TileKey) -> bytes:
        lon, lat = self._grid(key)
        cell_deg = cell_degrees(key.z)
        cell_m = math.radians(cell_deg) * MOON_RADIUS
        needed = 1.0 / cell_deg  # pixels per degree to match the cells
        ppd = next((p for p in Lola.RESOLUTIONS if p >= 0.5 * needed), Lola.RESOLUTIONS[-1])
        heights, measured = self.lola().sample(lon.ravel(), lat.ravel(), ppd)
        heights, measured = heights.reshape(lon.shape), measured.reshape(lon.shape)
        provenance = np.where(measured, MEASURED, FILLED).astype(np.uint8)
        pixel_m = math.radians(1.0 / ppd) * MOON_RADIUS
        if cell_m < pixel_m / 2:
            provenance[:] = FINER
        datum = "height above the sphere of 1737.4 km (m), Moon mean-Earth frame (MOON_ME)"
        sources = [f"LOLA LDEM_{ppd} and LDEC_{ppd} (NASA/GSFC), {pixel_m:.0f} m per pixel"]
        header = self._header("moon", key, heights, provenance, datum, sources, (pixel_m, pixel_m))
        header["cell_m"] = cell_m
        return encode(header, heights, provenance)
