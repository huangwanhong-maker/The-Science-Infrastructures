"""Build class products: a category per cell of a regular longitude-latitude grid.

Some gridded data are categories rather than quantities: climate classes, land
cover, soil types. They cannot be interpolated, averaged or quantised like a
field; a cell is in one class or another. A class product holds them as they
are, one byte per cell, with the legend that names each class and the colour
its authors gave it (the application may draw it otherwise).

The source is a set of GeoTIFF maps, one per period (and scenario), each with a
legend file of lines such as

    15: Cfb  Temperate, no dry season, warm summer [100 255 80]

The georeferencing is read from the GeoTIFF tags (model tie point and pixel
scale; the raster type says whether they refer to cell corners or centres),
and checked to be geographic coordinates on WGS 84.

Each product also records the share of the land in each class, weighted by
the area of the cells (cos latitude).

Output, in data/cache/products/: classes-<name>.bin (uint8, rows north to
south, columns eastward from the first longitude at or after -180; 0 is no
value) and classes-<name>.json.

Usage::

    python data/products/class_raster.py --source koppen-geiger
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
from PIL import Image

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

LEGEND_LINE = re.compile(r"^\s*(\d+):\s+(\S+)\s+(.+?)\s+\[(\d+)\s+(\d+)\s+(\d+)\]\s*$")
SCENARIOS = {"ssp119": "SSP1-1.9", "ssp126": "SSP1-2.6", "ssp245": "SSP2-4.5", "ssp370": "SSP3-7.0",
             "ssp434": "SSP4-3.4", "ssp460": "SSP4-6.0", "ssp585": "SSP5-8.5"}
# GeoTIFF tags and keys (GeoTIFF 1.1, OGC 19-008r4)
PIXEL_SCALE, TIE_POINT, GEO_KEYS = 33550, 33922, 34735
RASTER_TYPE, GEOGRAPHIC_TYPE = 1025, 2048
PIXEL_IS_AREA, WGS84 = 1, 4326


def legend(path: Path) -> list[dict]:
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        m = LEGEND_LINE.match(line)
        if m:
            out.append({"value": int(m.group(1)), "code": m.group(2), "name": m.group(3).strip(),
                        "colour": [int(m.group(4)), int(m.group(5)), int(m.group(6))]})
    if not out:
        raise SystemExit(f"no classes found in {path}")
    return out


def geo_keys(tags) -> dict:
    raw = tags.get(GEO_KEYS)
    if not raw:
        return {}
    n = raw[3]
    return {raw[4 + 4 * i]: raw[4 + 4 * i + 3] for i in range(n) if raw[4 + 4 * i + 1] == 0}


def grid_of(image: Image.Image) -> dict:
    """The grid of a GeoTIFF in geographic coordinates: centres of the first column and row, spacing."""
    tags = image.tag_v2
    if PIXEL_SCALE not in tags or TIE_POINT not in tags:
        raise SystemExit("the GeoTIFF has no model tie point and pixel scale")
    keys = geo_keys(tags)
    if keys.get(GEOGRAPHIC_TYPE) != WGS84:
        raise SystemExit(f"expected geographic coordinates on WGS 84 (EPSG:4326), found GeoKeys {keys}")
    sx, sy = tags[PIXEL_SCALE][:2]
    i, j, _, x, y, _ = tags[TIE_POINT][:6]
    corner = keys.get(RASTER_TYPE, PIXEL_IS_AREA) == PIXEL_IS_AREA
    # Centre of cell (0, 0): the tie point refers to the corner of the raster for PixelIsArea.
    lon0 = x + (0.5 - i) * sx if corner else x - i * sx
    lat0 = y - (0.5 - j) * sy if corner else y + j * sy
    nlon, nlat = image.size
    periodic = abs(nlon * sx - 360) < 1e-6 * sx
    return {"lon0": round(lon0, 9), "dlon": sx, "nlon": nlon, "lat0": round(lat0, 9), "dlat": sy, "nlat": nlat,
            "periodic": bool(periodic), "crs": "WGS 84 geographic coordinates (EPSG:4326), degrees",
            "order": "rows north to south, columns eastward"}


def land_shares(values: np.ndarray, grid: dict, classes: list[dict]) -> dict:
    """Share of the classified area in each class, weighting cells by cos latitude."""
    lats = grid["lat0"] - np.arange(grid["nlat"]) * grid["dlat"]
    w = np.cos(np.radians(lats))[:, None] * np.ones((1, grid["nlon"]))
    totals = np.bincount(values.ravel(), weights=w.ravel(), minlength=256)
    land = totals[1:].sum()
    return {c["code"]: round(float(totals[c["value"]] / land), 6) for c in classes}


def build(source: str) -> list[dict]:
    manifest = json.loads((DATA / "registry" / f"{source}.json").read_text(encoding="utf-8"))
    folder = DATA / "cache" / source
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    classes = legend(folder / "legend.txt")
    out = []
    for entry in manifest["files"]:
        if not entry["file"].endswith((".tif", ".tiff")):
            continue
        image = Image.open(folder / entry["file"])
        grid = grid_of(image)
        values = np.asarray(image, dtype=np.uint8)
        if values.max() > max(c["value"] for c in classes):
            raise SystemExit(f"{entry['file']}: values beyond the legend")
        if grid["periodic"]:
            shift = int(round((-180 - (grid["lon0"] - grid["dlon"] / 2)) / grid["dlon"])) % grid["nlon"]
            values = np.roll(values, -shift, axis=1)
            grid["lon0"] = round(grid["lon0"] + shift * grid["dlon"], 9)
        start, end = (int(y) for y in entry["period"].split("-"))
        scenario = SCENARIOS.get(entry.get("scenario") or "", entry.get("scenario"))
        name = f"{source}-{entry['period']}" + (f"-{entry['scenario']}" if entry.get("scenario") else "")
        title = f"{manifest['title'].split(',')[0]}, {start}–{end}" + (f", {scenario}" if scenario else "")
        bin_path = OUT / f"classes-{name}.bin"
        OUT.mkdir(parents=True, exist_ok=True)
        bin_path.write_bytes(np.ascontiguousarray(values).tobytes())
        product = {
            "schema": "class-raster/1", "name": name, "title": title, "domain": manifest.get("domain", "geography"),
            "series": source, "period": [start, end], "scenario": scenario,
            "grid": grid, "classes": classes, "no_value": 0,
            "shares": land_shares(values, grid, classes),
            "file": bin_path.name, "encoding": "uint8, one byte per cell; 0 means no value",
            "source": {"dataset": source, "provider": manifest["provider"], "credit": manifest.get("credit"),
                       "licence": manifest["licence"], "citation": manifest["citation"],
                       "retrieved_at": retrieval["retrieved_at"], "file": entry["file"], "member": entry.get("member")},
        }
        (OUT / f"classes-{name}.json").write_text(json.dumps(product, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        register("classes", {"name": name, "file": f"classes-{name}.json", "title": title, "domain": product["domain"],
                             "series": source, "period": [start, end], "scenario": scenario,
                             "grid": {k: grid[k] for k in ("dlon", "dlat", "nlon", "nlat")},
                             "source": {"dataset": source, "credit": manifest.get("credit"),
                                        "retrieved_at": retrieval["retrieved_at"]}})
        out.append(product)
    return out


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", required=True, help="a registered dataset of class maps, already fetched")
    args = parser.parse_args(argv)
    if not (DATA / "cache" / args.source / "retrieval.json").exists():
        print(f"not fetched; run: python data/acquisition/fetch.py {args.source}", file=sys.stderr)
        return 1
    for p in build(args.source):
        top = sorted(p["shares"].items(), key=lambda kv: -kv[1])[:3]
        print(f"classes-{p['name']}: {p['grid']['nlon']} x {p['grid']['nlat']}; largest: "
              + ", ".join(f"{k} {v:.1%}" for k, v in top))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
