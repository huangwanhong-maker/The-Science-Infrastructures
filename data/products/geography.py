"""Build the base geography product from Natural Earth: geometry, identifiers and names only.

The product is the map on which every other dataset is drawn, so it carries no
thematic values. Natural Earth's countries file also holds population, GDP and
income groups; these are left out on purpose. Thematic data come from their own
datasets and are joined to the countries by key in the application.

Output, in data/cache/products/:

geography-<name>.bin
    float32 little-endian (longitude, latitude) pairs in degrees, WGS 84. Every
    ring and line of the product is a contiguous run of pairs.
geography-<name>.json
    Layers, each a list of features. A feature refers to its runs by
    [first pair, number of pairs]. Polygons are given as rings (outer and inner
    alike: the even-odd rule decides what is inside); lines as polylines.

Usage::

    python data/products/geography.py
"""

from __future__ import annotations

import argparse
import json
import math
import struct
import sys
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

SOURCE = DATA / "cache" / "natural-earth-50m"
LANGUAGES = ["en", "zh", "zht", "ja", "ko", "fr", "de", "es", "pt", "it", "nl", "pl", "sv", "ru", "uk",
             "el", "tr", "ar", "fa", "he", "hi", "bn", "ur", "id", "vi", "hu"]
BOUNDARY_CLASS = {
    "International boundary (verify)": "international",
    "Disputed (please verify)": "disputed",
    "Line of control (please verify)": "control",
    "Indefinite (please verify)": "indefinite",
    "Indeterminant frontier": "indefinite",
}


class Builder:
    def __init__(self):
        self.pairs: list[float] = []

    def run(self, coords) -> list[int]:
        """Appends a ring or line; returns [first pair, number of pairs]."""
        first = len(self.pairs) // 2
        for lon, lat, *_ in coords:
            self.pairs += (lon, lat)
        return [first, len(self.pairs) // 2 - first]

    def rings(self, geometry) -> list[list[int]]:
        polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
        return [self.run(ring) for polygon in polygons for ring in polygon]

    def lines(self, geometry) -> list[list[int]]:
        lines = [geometry["coordinates"]] if geometry["type"] == "LineString" else geometry["coordinates"]
        return [self.run(line) for line in lines]


def read(name: str) -> list[dict]:
    return json.loads((SOURCE / f"ne_50m_{name}.geojson").read_text(encoding="utf-8"))["features"]


def code(value) -> str | None:
    """A country code, or None for Natural Earth's '-99' (no code)."""
    return None if value in (None, "", "-99", -99) else str(value)


def names(props: dict, prefix: str, base: str) -> dict:
    out = {"default": props.get(base)}
    for lang in LANGUAGES:
        value = props.get(f"{prefix}{lang}") or props.get(f"{prefix}{lang}".upper())
        if value:
            out[lang] = value
    return out


# ---- label points -------------------------------------------------------------

def label_point(geometry) -> list[float]:
    """An interior point far from the edges (pole of inaccessibility) of the largest polygon.

    A grid search refined three times around the best point, after Agafonkin's
    polylabel. Distances are measured in longitude and latitude, with longitude
    scaled by the cosine of the mean latitude, which suffices to place a label.
    """
    polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
    polygon = max(polygons, key=lambda p: len(p[0]))
    rings = [np.asarray(ring, dtype=float)[:, :2] for ring in polygon]
    x0, y0 = rings[0].min(axis=0)
    x1, y1 = rings[0].max(axis=0)
    k = math.cos(math.radians((y0 + y1) / 2))
    a = np.concatenate([r for r in rings]) * (k, 1)                       # edge starts
    b = np.concatenate([np.roll(r, -1, axis=0) for r in rings]) * (k, 1)  # edge ends
    d = b - a
    dd = np.maximum((d ** 2).sum(axis=1), 1e-30)
    x0, x1 = x0 * k, x1 * k
    best, best_d = ((x0 + x1) / 2, (y0 + y1) / 2), -1.0
    for _ in range(4):
        gx, gy = np.meshgrid(np.linspace(x0, x1, 33), np.linspace(y0, y1, 33))
        p = np.stack([gx.ravel(), gy.ravel()], axis=1)[:, None, :]         # points x edges
        crosses = ((a[:, 1] > p[..., 1]) != (b[:, 1] > p[..., 1])) & (
            p[..., 0] < d[:, 0] * (p[..., 1] - a[:, 1]) / np.where(d[:, 1] == 0, 1e-30, d[:, 1]) + a[:, 0])
        inside = crosses.sum(axis=1) % 2 == 1
        t = np.clip(((p - a) * d).sum(axis=2) / dd, 0, 1)
        dist = (((p - a - t[..., None] * d) ** 2).sum(axis=2)).min(axis=1)
        dist[~inside] = -1
        i = int(dist.argmax())
        if dist[i] > best_d:
            best, best_d = (float(p[i, 0, 0]), float(p[i, 0, 1])), float(dist[i])
        wx, wy = (x1 - x0) / 8, (y1 - y0) / 8
        x0, x1, y0, y1 = best[0] - wx, best[0] + wx, best[1] - wy, best[1] + wy
    return [round(best[0] / k, 4), round(best[1], 4)]


# ---- layers -------------------------------------------------------------------

def build(name: str) -> dict:
    b = Builder()
    countries = []
    for f in read("admin_0_countries"):
        p = f["properties"]
        countries.append({
            "id": p["ADM0_A3"],
            "keys": {
                "adm0_a3": code(p["ADM0_A3"]), "iso_a3": code(p["ISO_A3"]), "iso_a3_eh": code(p["ISO_A3_EH"]),
                "wb_a3": code(p["WB_A3"]), "iso_a2": code(p["ISO_A2"]), "iso_a2_eh": code(p["ISO_A2_EH"]),
                "wb_a2": code(p["WB_A2"]), "iso_n3": code(p["ISO_N3"]), "sov_a3": code(p["SOV_A3"]),
            },
            "names": names(p, "NAME_", "NAME") | {"long": p["NAME_LONG"], "formal": p.get("FORMAL_EN")},
            "type": p["TYPE"],
            "continent": p["CONTINENT"], "region_un": p["REGION_UN"], "subregion": p["SUBREGION"],
            "label": [p["LABEL_X"], p["LABEL_Y"]],
            "label_rank": p["LABELRANK"], "min_label": p["MIN_LABEL"],
            "rings": b.rings(f["geometry"]),
        })
    land = [{"rings": b.rings(f["geometry"])} for f in read("land")]
    lakes = []
    for f in read("lakes"):
        p = f["properties"]
        lakes.append({"names": names(p, "name_", "name"), "rank": p["scalerank"],
                      "label": label_point(f["geometry"]) if p.get("name") and p["scalerank"] <= 2 else None,
                      "rings": b.rings(f["geometry"])})
    marine = []
    for f in read("geography_marine_polys"):
        p = f["properties"]
        marine.append({"names": names(p, "name_", "name"), "class": p["featurecla"], "rank": p["scalerank"],
                       "min_label": p["min_label"], "label": label_point(f["geometry"]),
                       "rings": b.rings(f["geometry"])})
    glaciers = [{"rings": b.rings(f["geometry"])} for f in read("glaciated_areas")]
    ice_shelves = [{"rings": b.rings(f["geometry"])} for f in read("antarctic_ice_shelves_polys")]
    coastline = [{"lines": b.lines(f["geometry"])} for f in read("coastline")]
    boundaries = [{"class": BOUNDARY_CLASS.get(f["properties"]["FEATURECLA"], "international"),
                   "lines": b.lines(f["geometry"])} for f in read("admin_0_boundary_lines_land")]
    rivers = [{"name": f["properties"].get("name_en") or f["properties"].get("name"),
               "class": "lake centreline" if f["properties"]["featurecla"] == "Lake Centerline" else "river",
               "rank": f["properties"]["scalerank"], "lines": b.lines(f["geometry"])}
              for f in read("rivers_lake_centerlines")]
    places = []
    for f in read("populated_places"):
        p = f["properties"]
        lon, lat = f["geometry"]["coordinates"][:2]
        places.append({"names": names(p, "NAME_", "NAME"), "lon": lon, "lat": lat, "country": code(p["ADM0_A3"]),
                       "class": p["FEATURECLA"], "capital": bool(p["ADM0CAP"]), "rank": p["SCALERANK"],
                       "label_rank": p["LABELRANK"], "timezone": p.get("TIMEZONE")})
    places.sort(key=lambda q: (q["rank"], q["label_rank"]))

    bin_path = OUT / f"geography-{name}.bin"
    OUT.mkdir(parents=True, exist_ok=True)
    bin_path.write_bytes(struct.pack(f"<{len(b.pairs)}f", *b.pairs))
    manifest = json.loads((DATA / "registry" / "natural-earth-50m.json").read_text(encoding="utf-8"))
    retrieval = json.loads((SOURCE / "retrieval.json").read_text(encoding="utf-8"))
    product = {
        "schema": "geography/1",
        "name": name,
        "title": "Natural Earth, 1:50 million",
        "source": {"dataset": manifest["id"], "provider": manifest["provider"], "credit": manifest.get("credit"),
                   "licence": manifest["licence"],
                   "citation": manifest["citation"], "retrieved_at": retrieval["retrieved_at"]},
        "coordinates": {"file": bin_path.name, "pairs": len(b.pairs) // 2,
                        "encoding": "float32 little-endian (longitude, latitude), degrees, WGS 84"},
        "languages": LANGUAGES,
        "boundaries_depiction": "Natural Earth's default (de facto) depiction",
        "layers": {"countries": countries, "land": land, "lakes": lakes, "glaciers": glaciers,
                   "ice_shelves": ice_shelves, "marine": marine,
                   "coastline": coastline, "boundaries": boundaries, "rivers": rivers, "places": places},
    }
    json_path = OUT / f"geography-{name}.json"
    json_path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    register("geography", {
        "name": name, "file": json_path.name, "title": product["title"],
        "layers": {k: len(v) for k, v in product["layers"].items()},
        "source": product["source"],
    }, first=name)
    return product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--name", default="ne50m", help="product name (default: ne50m)")
    args = parser.parse_args(argv)
    if not (SOURCE / "retrieval.json").exists():
        print("Natural Earth is not in the cache; run: python data/acquisition/fetch.py natural-earth-50m",
              file=sys.stderr)
        return 1
    product = build(args.name)
    sizes = {k: len(v) for k, v in product["layers"].items()}
    print(f"geography-{args.name}: {product['coordinates']['pairs']} coordinate pairs; {sizes}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
