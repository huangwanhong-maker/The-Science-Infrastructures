"""Synthetic detail below the data's resolution: generated, labelled, and off unless asked for.

Where a tile's cells are finer than half its source's resolution (provenance 4), its heights are
the source's, interpolated; they are smooth where the ground is not. This module can add detail
there, for drawing only, as a self-affine extrapolation of the region's roughness to the scales the
data do not resolve, and marks every sample it changes as synthetic (provenance 5) so that the page
can say so. It is never used for a computation.

The model. Topography is statistically self-affine over a range of scales: the RMS difference of
heights a distance λ apart grows as λ^H, with H the Hurst exponent (Mark and Aronson 1984, Math.
Geol. 16, 671: H ≈ 0.6–0.8 for land between 0.1 and 10 km; Rosenburg et al. 2011, J. Geophys. Res.
116, E02001: the Moon's highlands H ≈ 0.95 and maria H ≈ 0.76 between 17 m and 2.7 km, from LOLA).
The region's RMS gradient g at the coarse scale R is read from a coarse grid (the Earth: ETOPO 2022
at 5′, R = 9.26 km; the Moon: LOLA at 4 pixels per degree, R = 7.58 km). The detail is a sum of
octaves of value noise, of wavelengths λ_k = r / 2^(k+1) from half the source's resolution r down to
the tile's Nyquist limit, twice its cell c, each of RMS amplitude w(λ/c) a(λ) with
    a(λ) = g R^(1−H) λ^H,
the height difference self-affinity gives at λ when the gradient at R is g, and
    w = min(1, log₂(λ / 2c) / 2),
none at the limit and the whole from λ = 8c: an octave enters over two levels of the pyramid, so
that a tile and its parent or its child differ by at most half an octave and the levels do not show
as patches of rougher and smoother ground. The noise is a function
of the point on the body (the lattice is in body-fixed metres on the mean sphere, hashed), and the
amplitude a function of the coarse grid at that point, so that neighbouring tiles agree at their
shared samples.

Caveat. This is a texture with the right statistics in the mean, not terrain: it puts hills where
there are none and knows nothing of rivers, cliffs or craters. Its parameters are in the header of
every tile it touches.
"""

from __future__ import annotations

import math

import numpy as np

from .pyramid import BORDER, GRID, TileKey, sample_lonlat
from .sources import FINER, SYNTHETIC
from .terrain import EARTH_RADIUS, MOON_RADIUS, decode, encode

HURST = {"earth": 0.8, "moon": 0.9}
COARSE_SCALE = {"earth": 9260.0, "moon": math.radians(0.25) * MOON_RADIUS}
G_MIN = 0.002   # the least regional gradient: no ground is perfectly smooth
MAX_OCTAVES = 8

# body -> (|∇z|, longitude and latitude of the first sample's centre, degrees per sample)
_slopes: dict[str, tuple[np.ndarray, float, float, float]] = {}


def _gradient(z: np.ndarray, radius: float, ppd: float) -> np.ndarray:
    """|∇z| of a grid from 90 N, ppd per degree, rows north to south, periodic in longitude."""
    step = math.radians(1.0 / ppd) * radius
    lat = 90.0 - (np.arange(z.shape[0]) + 0.5) / ppd
    coslat = np.maximum(np.cos(np.radians(lat)), 0.01)[:, None]
    zx = (np.roll(z, -1, axis=1) - np.roll(z, 1, axis=1)) / (2 * step * coslat)
    zy = np.zeros_like(z)
    zy[1:-1] = (z[:-2] - z[2:]) / (2 * step)
    return np.hypot(zx, zy)


def slope_grid(service, body: str) -> tuple[np.ndarray, float, float, float]:
    """The regional gradient of a body: the grid, the longitude and latitude of its first sample
    (degrees) and its spacing (degrees)."""
    if body not in _slopes:
        if body == "earth":
            e = service.etopo()
            z = np.nan_to_num(e.z, nan=0.0)
            step = e.g["dlon"]
            _slopes[body] = (_gradient(z, EARTH_RADIUS, 1.0 / step), e.g["lon0"], e.g["lat0"], step)
        else:
            z = service.lola()._whole("ldem_4", 4).astype(np.float64) * 0.5
            _slopes[body] = (_gradient(z, MOON_RADIUS, 4.0), 0.125, 89.875, 0.25)
    return _slopes[body]


def regional_gradient(service, body: str, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
    """The regional gradient at points (degrees), from the nearest sample of the coarse grid."""
    grid, lon0, lat0, step = slope_grid(service, body)
    rows = np.clip(np.round((lat0 - lats) / step).astype(int), 0, grid.shape[0] - 1)
    cols = np.round(np.mod(lons - lon0, 360.0) / step).astype(int) % grid.shape[1]
    return grid[rows, cols]


def _hash(ix: np.ndarray, iy: np.ndarray, iz: np.ndarray, seed: int) -> np.ndarray:
    """A value in [−1, 1] for each lattice point: SplitMix64's finaliser of the packed coordinates."""
    with np.errstate(over="ignore"):
        h = (ix.astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15)) ^ (iy.astype(np.uint64) * np.uint64(0xC2B2AE3D27D4EB4F)) \
            ^ (iz.astype(np.uint64) * np.uint64(0x165667B19E3779F9)) ^ np.uint64(seed * 0x27D4EB2F165667C5 & 0xFFFFFFFFFFFFFFFF)
        h ^= h >> np.uint64(30)
        h *= np.uint64(0xBF58476D1CE4E5B9)
        h ^= h >> np.uint64(27)
        h *= np.uint64(0x94D049BB133111EB)
        h ^= h >> np.uint64(31)
    return (h >> np.uint64(11)).astype(np.float64) / float(1 << 53) * 2.0 - 1.0


def value_noise(p: np.ndarray, seed: int) -> np.ndarray:
    """Value noise on the integer lattice at points p (N x 3), quintic interpolation, unit RMS."""
    i = np.floor(p)
    f = p - i
    i = i.astype(np.int64)
    w = f * f * f * (f * (f * 6 - 15) + 10)
    out = np.zeros(len(p))
    for dx in (0, 1):
        wx = w[:, 0] if dx else 1 - w[:, 0]
        for dy in (0, 1):
            wy = w[:, 1] if dy else 1 - w[:, 1]
            for dz in (0, 1):
                wz = w[:, 2] if dz else 1 - w[:, 2]
                out += wx * wy * wz * _hash(i[:, 0] + dx, i[:, 1] + dy, i[:, 2] + dz, seed)
    # The RMS of quintic value noise with corners uniform on [−1, 1], measured (NOISE_RMS below).
    return out / NOISE_RMS


NOISE_RMS = 0.4006  # measured on 10⁶ random points (data/tests/check_tiles.py holds it)


def weight(lam: float, cell: float) -> float:
    """An octave's weight in a tile of cells `cell`: 0 at its Nyquist limit (λ = 2 cells), 1 from λ = 8 cells."""
    return min(1.0, max(0.0, math.log2(lam / (2 * cell)) / 2))


def detail(body: str, lons: np.ndarray, lats: np.ndarray, gradient: np.ndarray, resolution: float, cell: float) -> tuple[np.ndarray, list]:
    """The synthetic heights (m) at points (degrees), and the octaves' wavelengths (their weights by `weight`)."""
    radius = EARTH_RADIUS if body == "earth" else MOON_RADIUS
    H, R = HURST[body], COARSE_SCALE[body]
    lon, lat = np.radians(lons), np.radians(lats)
    xyz = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=1) * radius
    g = np.maximum(gradient, G_MIN)
    out = np.zeros(len(lons))
    waves = []
    lam = resolution / 2
    k = 0
    while lam > 2 * cell and k < MAX_OCTAVES:
        amplitude = weight(lam, cell) * g * R ** (1 - H) * lam ** H
        out += amplitude * value_noise(xyz / lam, seed=k + (0 if body == "earth" else 1000))
        waves.append(lam)
        lam /= 2
        k += 1
    return out, waves


def with_detail(service, body: str, data: bytes) -> bytes:
    """The tile with synthetic detail where its cells are finer than its data; unchanged elsewhere."""
    header, heights, provenance, floor = decode(data)
    finer = provenance == FINER
    if not finer.any():
        return data
    key = TileKey(header["z"], header["x"], header["y"])
    lons, lats = sample_lonlat(key, GRID, BORDER)
    lon, lat = np.meshgrid(np.array(lons), np.array(lats))
    gradient = regional_gradient(service, body, lon[finer], lat[finer])
    resolution = float(header["resolution_m"][0])
    dh, waves = detail(body, lon[finer], lat[finer], gradient, resolution, float(header["cell_m"]))
    if not waves:
        return data
    heights, provenance = heights.copy(), provenance.copy()
    heights[finer] += dh
    provenance[finer] = SYNTHETIC
    if floor is not None:
        floor = floor.copy()
        floor[finer] += dh
    inner = slice(BORDER, BORDER + GRID)
    p = provenance[inner, inner]
    header = {**header, "min": float(heights[inner, inner].min()), "max": float(heights[inner, inner].max()),
              "provenance": {str(c): int((p == c).sum()) for c in np.unique(p)},
              "synthetic": {"model": "self-affine value noise scaled by the regional gradient (data/tiles/synthetic.py)",
                            "hurst": HURST[body], "coarse_scale_m": COARSE_SCALE[body], "wavelengths_m": [round(w, 2) for w in waves],
                            "weights": [round(weight(w, float(header["cell_m"])), 3) for w in waves]}}
    header.pop("layers", None)
    return encode(header, heights, provenance, floor)
