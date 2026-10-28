"""The geographic tile pyramid of the globes: which part of a body each tile covers.

Level z of the pyramid has 2^(z+1) columns and 2^z rows of tiles, from 180 W and 90 N, in east
longitude and planetocentric (for the Earth, geodetic) latitude: the scheme of NASA's Trek
portals and of the shared visualization layer's tiles.js (visualization/js/globe/tiling.js
holds the same functions for the page, and the checks hold the two equal). A terrain tile is a
grid of n x n samples (n = 65: 64 cells) at the corners of its cells, edges included, so that a
tile and its neighbour sample the same points along their shared edge and no crack opens between
them.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

GRID = 65  # samples per side of a terrain tile
BORDER = 1  # samples beyond each edge, for the normals at the edge (the tile carries GRID + 2 BORDER)


@dataclass(frozen=True)
class TileKey:
    z: int
    x: int
    y: int

    def check(self) -> "TileKey":
        if not (0 <= self.z <= 24):
            raise ValueError(f"level {self.z} is outside 0..24")
        if not (0 <= self.x < 2 ** (self.z + 1) and 0 <= self.y < 2 ** self.z):
            raise ValueError(f"tile {self.z}/{self.x}/{self.y} does not exist")
        return self


def bounds(key: TileKey) -> tuple[float, float, float, float]:
    """(west, south, east, north) of a tile, degrees."""
    span = 180.0 / 2 ** key.z
    west = -180.0 + key.x * span
    north = 90.0 - key.y * span
    return west, north - span, west + span, north


def sample_lonlat(key: TileKey, n: int = GRID, border: int = 0):
    """The longitudes (west to east) and latitudes (north to south) of a tile's samples, degrees,
    with `border` more samples beyond each edge at the same spacing; latitudes beyond a pole are
    reflected back over it (the sample on the other side of the pole, which the longitude then
    names by its antimeridian, is taken at the reflected latitude and the same longitude, which is
    what the normals at the pole need to the order of the cell)."""
    west, south, east, north = bounds(key)
    step = (east - west) / (n - 1)
    lons = [west + step * i for i in range(-border, n + border)]
    lats = []
    for j in range(-border, n + border):
        lat = north - step * j
        lat = 180 - lat if lat > 90 else (-180 - lat if lat < -90 else lat)
        lats.append(lat)
    return lons, lats


def cell_degrees(z: int, n: int = GRID) -> float:
    """The angular size of a cell of a tile of level z, degrees."""
    return 180.0 / 2 ** z / (n - 1)


def cell_metres(z: int, radius: float, latitude_deg: float = 0.0, n: int = GRID) -> float:
    """The east-west size of a cell at a latitude, metres, on a sphere of the given radius."""
    return math.radians(cell_degrees(z, n)) * radius * max(math.cos(math.radians(latitude_deg)), 1e-6)
