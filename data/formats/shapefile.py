"""ESRI shapefiles, as far as the data library meets them: points, polylines and polygons (with their Z and M variants, whose extra
values are passed over) and the attribute table beside them (ESRI Shapefile Technical Description, 1998; dBase III).

    read(path)  ->  [{"type": "point" | "line" | "polygon", "parts": [[(x, y), ...], ...], "attributes": {...}}, ...]

A polygon's rings come in the file's order, outer rings clockwise and holes counter-clockwise by the format's convention; they are
not grouped here (`polygons_of` groups them). Coordinates are as the file holds them; its .prj says what they mean.
"""

from __future__ import annotations

import struct
from pathlib import Path

POINT, LINE, POLYGON = {1, 11, 21}, {3, 13, 23}, {5, 15, 25}


def _table(path: Path) -> list[dict]:
    if not path.exists():
        return []
    data = path.read_bytes()
    count, header, width = struct.unpack("<IHH", data[4:12])
    fields, pos = [], 32
    while data[pos] != 0x0D:
        name = data[pos:pos + 11].split(b"\0")[0].decode("latin-1")
        fields.append((name, chr(data[pos + 11]), data[pos + 16]))
        pos += 32
    rows = []
    for k in range(count):
        at = header + k * width + 1
        row = {}
        for name, kind, size in fields:
            raw = data[at:at + size].decode("utf-8", "replace").strip()
            at += size
            if kind in "NF":
                try:
                    row[name] = float(raw) if raw else None
                except ValueError:
                    row[name] = None
            else:
                row[name] = raw
        rows.append(row)
    return rows


def read(path: str | Path) -> list[dict]:
    path = Path(path)
    data = path.read_bytes()
    table = _table(path.with_suffix(".dbf"))
    out, pos, k = [], 100, 0
    while pos + 8 <= len(data):
        _, length = struct.unpack(">II", data[pos:pos + 8])
        body = data[pos + 8:pos + 8 + 2 * length]
        pos += 8 + 2 * length
        kind = struct.unpack("<i", body[:4])[0]
        feature = {"type": None, "parts": [], "attributes": table[k] if k < len(table) else {}}
        k += 1
        if kind in POINT:
            feature["type"] = "point"
            feature["parts"] = [[struct.unpack("<2d", body[4:20])]]
        elif kind in LINE or kind in POLYGON:
            feature["type"] = "line" if kind in LINE else "polygon"
            parts, points = struct.unpack("<2i", body[36:44])
            starts = list(struct.unpack(f"<{parts}i", body[44:44 + 4 * parts])) + [points]
            at = 44 + 4 * parts
            xy = struct.unpack(f"<{2 * points}d", body[at:at + 16 * points])
            feature["parts"] = [[(xy[2 * i], xy[2 * i + 1]) for i in range(starts[j], starts[j + 1])] for j in range(parts)]
        elif kind != 0:
            raise ValueError(f"{path}: shape type {kind} is not read")
        out.append(feature)
    return out


def _area(ring: list) -> float:
    return sum(ring[i][0] * ring[(i + 1) % len(ring)][1] - ring[(i + 1) % len(ring)][0] * ring[i][1] for i in range(len(ring))) / 2


def _inside(p, ring: list) -> bool:
    x, y = p
    hit = False
    for i in range(len(ring)):
        (x0, y0), (x1, y1) = ring[i], ring[i - 1]
        if (y0 > y) != (y1 > y) and x < (x1 - x0) * (y - y0) / (y1 - y0) + x0:
            hit = not hit
    return hit


def polygons_of(feature: dict) -> list[dict]:
    """A polygon feature's rings as [{"outer": ring, "holes": [ring, ...]}]: clockwise rings outer, the others holes of the outer holding them."""
    polygons = [{"outer": r, "holes": []} for r in feature["parts"] if _area(r) < 0]
    for r in feature["parts"]:
        if _area(r) >= 0:
            for p in polygons:
                if _inside(r[0], p["outer"]):
                    p["holes"].append(r)
                    break
    return polygons
