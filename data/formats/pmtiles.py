"""PMTiles and Mapbox Vector Tiles: tiles of a whole world in one file, read by HTTP range requests.

A PMTiles archive (version 3; Protomaps, "PMTiles version 3 specification") is one file holding a pyramid of tiles: a
127-byte header, a root directory, leaf directories, and the tiles. A directory is a compressed run of varints that maps a
Hilbert tile id to an offset and a length, so a reader asks for the header, then at most a few directories, then the bytes
of the one tile it wants, and a world of hundreds of gigabytes costs a few requests a tile.

A Mapbox Vector Tile (Mapbox Vector Tile Specification 2.1) is a protocol buffer of named layers, each of features with a
geometry in a command language over integer coordinates in the tile's own square, and with attributes held as indices
into the layer's tables of keys and values.

Only the standard library is needed. Geometries come back as rings and lines of tile coordinates in [0, extent], y down;
`to_lonlat` brings them to WGS 84, through the Web Mercator square a tile is cut from.

References: Protomaps, PMTiles v3 (https://github.com/protomaps/PMTiles/blob/main/spec/v3/spec.md); Mapbox, Vector Tile
Specification 2.1 (https://github.com/mapbox/vector-tile-spec); Hilbert curve ordering after Skilling (2004).
"""

from __future__ import annotations

import gzip
import json
import math
import struct
import threading
import urllib.request
import zlib

AGENT = "noetarium-data/1 (+https://github.com/huangwanhong-maker/Noetarium-The-Science-and-Humanity-Infrastructure)"


def _varint(data: bytes, pos: int) -> tuple[int, int]:
    value, shift = 0, 0
    while True:
        b = data[pos]
        pos += 1
        value |= (b & 0x7F) << shift
        if not b & 0x80:
            return value, pos
        shift += 7


def tile_id(z: int, x: int, y: int) -> int:
    """The tile's position along the archive's Hilbert curve (the specification's ZXY to tile id)."""
    if z == 0:
        return 0
    acc = ((1 << (2 * z)) - 1) // 3
    n = 1 << z
    rx = ry = d = 0
    s = n >> 1
    tx, ty = x, y
    while s > 0:
        rx = 1 if tx & s else 0
        ry = 1 if ty & s else 0
        d += s * s * ((3 * rx) ^ ry)
        # Rotate.
        if ry == 0:
            if rx == 1:
                tx = n - 1 - tx
                ty = n - 1 - ty
            tx, ty = ty, tx
        s >>= 1
    return acc + d


def _decompress(data: bytes, kind: int) -> bytes:
    if kind in (0, 1):      # unknown, none
        return data
    if kind == 2:           # gzip
        return gzip.decompress(data)
    raise ValueError(f"PMTiles: compression {kind} is not read")


def _directory(raw: bytes) -> list[tuple[int, int, int, int]]:
    """A directory's entries: (tile id, run length, offset, length)."""
    n, pos = _varint(raw, 0)
    ids, runs, lengths, offsets = [], [], [], []
    last = 0
    for _ in range(n):
        d, pos = _varint(raw, pos)
        last += d
        ids.append(last)
    for _ in range(n):
        v, pos = _varint(raw, pos)
        runs.append(v)
    for _ in range(n):
        v, pos = _varint(raw, pos)
        lengths.append(v)
    for i in range(n):
        v, pos = _varint(raw, pos)
        offsets.append(offsets[i - 1] + lengths[i - 1] if v == 0 and i > 0 else v - 1)
    return list(zip(ids, runs, offsets, lengths))


class Archive:
    """A PMTiles archive on the web."""

    def __init__(self, url: str, timeout: int = 60):
        self.url, self.timeout = url, timeout
        self._dirs: dict[int, list] = {}
        self._lock = threading.Lock()
        head = self._range(0, 16384)
        if head[:7] != b"PMTiles" or head[7] != 3:
            raise ValueError(f"{url}: not a PMTiles archive of version 3")
        v = struct.unpack("<QQQQQQQQQQQQ", head[8:104])
        self.root_offset, self.root_length, self.meta_offset, self.meta_length = v[0], v[1], v[2], v[3]
        self.leaf_offset, self.leaf_length, self.tile_offset, self.tile_length = v[4], v[5], v[6], v[7]
        self.internal_compression, self.tile_compression, self.tile_type = head[97], head[98], head[99]
        self.min_zoom, self.max_zoom = head[100], head[101]
        self.bounds = [struct.unpack("<i", head[102 + 4 * k:106 + 4 * k])[0] / 1e7 for k in range(4)]
        self.root = _directory(_decompress(head[self.root_offset:self.root_offset + self.root_length], self.internal_compression))

    def _range(self, start: int, length: int) -> bytes:
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={start}-{start + length - 1}", "User-Agent": AGENT})
        with urllib.request.urlopen(request, timeout=self.timeout) as r:
            return r.read()

    def metadata(self) -> dict:
        return json.loads(_decompress(self._range(self.meta_offset, self.meta_length), self.internal_compression))

    def tile(self, z: int, x: int, y: int) -> bytes | None:
        """A tile's bytes as stored, decompressed (a vector tile), or None where the archive holds none."""
        want = tile_id(z, x, y)
        entries = self.root
        for _ in range(4):
            found = None
            lo, hi = 0, len(entries) - 1
            while lo <= hi:
                mid = (lo + hi) // 2
                if entries[mid][0] <= want:
                    found = entries[mid]
                    lo = mid + 1
                else:
                    hi = mid - 1
            if found is None:
                return None
            tid, run, off, length = found
            if run > 0:         # tiles
                if want >= tid + run:
                    return None
                return _decompress(self._range(self.tile_offset + off, length), self.tile_compression)
            # A leaf directory.
            with self._lock:
                entries = self._dirs.get(off)
            if entries is None:
                entries = _directory(_decompress(self._range(self.leaf_offset + off, length), self.internal_compression))
                with self._lock:
                    if len(self._dirs) > 256:
                        self._dirs.clear()
                    self._dirs[off] = entries
        return None


# ---- vector tiles ---------------------------------------------------------------------------------------------------

def _fields(data: bytes, start: int = 0, end: int | None = None):
    """The fields of a protocol buffer message: (number, wire type, value); length-delimited values as (start, end)."""
    pos, stop = start, len(data) if end is None else end
    while pos < stop:
        key, pos = _varint(data, pos)
        number, wire = key >> 3, key & 7
        if wire == 0:
            v, pos = _varint(data, pos)
        elif wire == 2:
            n, pos = _varint(data, pos)
            v = (pos, pos + n)
            pos += n
        elif wire == 5:
            v = data[pos:pos + 4]
            pos += 4
        elif wire == 1:
            v = data[pos:pos + 8]
            pos += 8
        else:
            raise ValueError(f"vector tile: wire type {wire}")
        yield number, wire, v


def _value(data: bytes, start: int, end: int):
    for number, wire, v in _fields(data, start, end):
        if number == 1:
            return data[v[0]:v[1]].decode("utf-8", "replace")
        if number == 2:
            return struct.unpack("<f", v)[0]
        if number == 3:
            return struct.unpack("<d", v)[0]
        if number in (4, 5):
            return v
        if number == 6:
            return (v >> 1) ^ -(v & 1)
        if number == 7:
            return bool(v)
    return None


def _geometry(data: bytes, start: int, end: int) -> list[list[tuple[int, int]]]:
    """A feature's geometry as parts: each a list of points in the tile's coordinates. MoveTo begins a part; ClosePath closes it."""
    parts, part = [], None
    x = y = 0
    pos = start
    while pos < end:
        command, pos = _varint(data, pos)
        op, count = command & 7, command >> 3
        if op == 7:
            if part:
                part.append(part[0])
            continue
        for _ in range(count):
            a, pos = _varint(data, pos)
            b, pos = _varint(data, pos)
            x += (a >> 1) ^ -(a & 1)
            y += (b >> 1) ^ -(b & 1)
            if op == 1:
                part = [(x, y)]
                parts.append(part)
            elif part is not None:
                part.append((x, y))
    return parts


def read_tile(data: bytes, want: set[str] | None = None) -> dict:
    """A vector tile's layers: {name: {extent, features: [{type (1 point, 2 line, 3 polygon), id, geometry, properties}]}}."""
    out = {}
    for number, wire, v in _fields(data):
        if number != 3 or wire != 2:
            continue
        name, extent, keys, values, raw = "", 4096, [], [], []
        for n, w, u in _fields(data, *v):
            if n == 1:
                name = data[u[0]:u[1]].decode("utf-8", "replace")
            elif n == 3:
                keys.append(data[u[0]:u[1]].decode("utf-8", "replace"))
            elif n == 4:
                values.append(_value(data, *u))
            elif n == 5:
                extent = u
            elif n == 2:
                raw.append(u)
        if want is not None and name not in want:
            continue
        features = []
        for s, e in raw:
            f = {"type": 0, "id": None, "geometry": [], "properties": {}}
            tags = None
            for n, w, u in _fields(data, s, e):
                if n == 1:
                    f["id"] = u
                elif n == 3:
                    f["type"] = u
                elif n == 2:
                    tags = u
                elif n == 4:
                    f["geometry"] = _geometry(data, *u)
            if tags:
                pos, stop = tags
                while pos < stop:
                    a, pos = _varint(data, pos)
                    b, pos = _varint(data, pos)
                    f["properties"][keys[a]] = values[b]
            features.append(f)
        out[name] = {"extent": extent, "features": features}
    return out


def to_lonlat(z: int, x: int, y: int, extent: int):
    """A function from a tile's coordinates to (longitude, latitude), degrees on WGS 84 (Web Mercator's sphere for the square)."""
    n = 2 ** z

    def convert(px: float, py: float) -> tuple[float, float]:
        lon = (x + px / extent) / n * 360.0 - 180.0
        lat = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * (y + py / extent) / n))))
        return lon, lat
    return convert


def signed_area(ring: list[tuple[float, float]]) -> float:
    """Twice the signed area by the shoelace formula, in the ring's own coordinates."""
    return sum(ring[i][0] * ring[(i + 1) % len(ring)][1] - ring[(i + 1) % len(ring)][0] * ring[i][1] for i in range(len(ring)))
