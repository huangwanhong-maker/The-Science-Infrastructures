"""Build the plates product: the outlines of the tectonic plates and the steps of
their boundaries, with each step's class and the polarity of subduction.

The product carries geometry only. How each plate moves (its rotation about a
pole) is a model, kept with the code that computes velocities from it (the
Earth Atlas's src/core/plates.js, after NNR-MORVEL56), so that the geometry of
one model can be used with the kinematics of another.

Input: data/cache/pb2002/ (data/registry/pb2002.json): Bird's (2003) PB2002
plates as GeoJSON, cut at the antimeridian, and his file of boundary steps.

Output, in data/cache/products/:
  plates-pb2002.json   plates (code, name, rings as runs of points), steps, classes, source
  plates-pb2002.bin    float32 (longitude, latitude) pairs of the rings

A step is [lon1, lat1, lon2, lat2, class, left, right, over, opening,
right_lateral]: its initial and final points (degrees), the index of its class,
the indices of the plates on its left and right (travelling from the initial
point, seen from outside the Earth), which plate overrides where it is a
subduction zone (1 the left, 2 the right, 0 none), and Bird's own components
of the relative velocity (mm/yr, from his PB2002 poles): divergent (negative
for convergence) and right-lateral, kept to check other models against. Bird's identifier "XX/YY" means the right-hand plate YY
subducts under the left-hand XX; "XX\\YY" the opposite.

Usage::

    python data/acquisition/fetch.py pb2002
    python data/products/plate_model.py
"""

from __future__ import annotations

import json
import struct
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

CLASSES = {
    "OSR": "oceanic spreading ridge", "CRB": "continental rift boundary",
    "OTF": "oceanic transform fault", "CTF": "continental transform fault",
    "SUB": "subduction zone", "OCB": "oceanic convergent boundary", "CCB": "continental convergent boundary",
}


def read_steps(path: Path, codes: list[str]) -> list[list]:
    index = {c: i for i, c in enumerate(codes)}
    out = []
    for line in path.read_text(encoding="latin-1").splitlines():
        if len(line) < 95:
            continue
        ident = line[6:11]
        left, mark, right = ident[:2], ident[2], ident[3:5]
        lon1, lat1, lon2, lat2 = (float(line[a:b]) for a, b in ((12, 20), (21, 28), (29, 37), (38, 45)))
        kind = line[92:95]
        if kind not in CLASSES or left not in index or right not in index:
            raise SystemExit(f"unexpected step: {line!r}")
        over = 1 if mark == "/" else 2 if mark == "\\" else 0
        opening, lateral = float(line[66:72]), float(line[73:79])
        if kind == "SUB" and not over:
            raise SystemExit(f"a subduction step without polarity: {line!r}")
        wrap = lambda x: round(((x + 180.0) % 360.0) - 180.0, 3)
        out.append([wrap(lon1), round(lat1, 3), wrap(lon2), round(lat2, 3),
                    list(CLASSES).index(kind), index[left], index[right], over, opening, lateral])
    return out


def build() -> dict:
    folder = DATA / "cache" / "pb2002"
    manifest = json.loads((DATA / "registry" / "pb2002.json").read_text(encoding="utf-8"))
    collection = json.loads((folder / "PB2002_plates.json").read_text(encoding="utf-8"))
    plates, coords = {}, []
    for f in collection["features"]:
        code, name = f["properties"]["Code"], f["properties"]["PlateName"]
        g = f["geometry"]
        polygons = g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]
        entry = plates.setdefault(code, {"code": code, "name": name, "rings": []})
        for polygon in polygons:
            for ring in polygon:
                points = ring[:-1] if ring[0] == ring[-1] else ring
                entry["rings"].append([len(coords), len(points)])
                coords.extend(points)
    codes = sorted(plates)
    steps = read_steps(folder / "PB2002_steps.dat.txt", codes)
    name = "pb2002"
    (OUT / f"plates-{name}.bin").write_bytes(b"".join(struct.pack("<2f", lon, lat) for lon, lat in coords))
    source = {"dataset": "pb2002", "provider": manifest["provider"], "credit": manifest["credit"],
              "licence": manifest["licence"], "citation": manifest["citation"]}
    product = {
        "schema": "plates/1", "name": name, "title": "Tectonic plates and their boundaries (PB2002)",
        "domain": "geophysics", "file": f"plates-{name}.bin",
        "plates": [plates[c] for c in codes],
        "classes": [{"code": k, "name": v} for k, v in CLASSES.items()],
        "step_fields": ["lon1", "lat1", "lon2", "lat2", "class", "left", "right", "over", "bird_opening", "bird_right_lateral"],
        "steps": steps,
        "units": {"coordinates": "degree_east, degree_north (WGS 84)"},
        "source": source,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"plates-{name}.json").write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("plates", {"name": name, "file": f"plates-{name}.json", "title": product["title"], "domain": "geophysics",
                        "count": len(codes), "source": {"dataset": "pb2002", "credit": manifest["credit"]}})
    return product


def main() -> int:
    if not (DATA / "cache" / "pb2002" / "PB2002_steps.dat.txt").exists():
        print("not fetched; run: python data/acquisition/fetch.py pb2002", file=sys.stderr)
        return 1
    product = build()
    by_class = {}
    for s in product["steps"]:
        by_class[product["classes"][s[4]]["code"]] = by_class.get(product["classes"][s[4]]["code"], 0) + 1
    rings = sum(len(p["rings"]) for p in product["plates"])
    print(f"plates-pb2002: {len(product['plates'])} plates in {rings} rings; {len(product['steps'])} steps {by_class}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
