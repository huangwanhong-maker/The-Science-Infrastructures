"""Build imagery products: pictures of the whole Earth in equirectangular projection.

Imagery is a kind of product apart from fields. Its pixels are colours already,
rendered by the provider from measurements (reflectance, radiance), not values
with units, so an application shows it as it is and does not probe it.

- **day**: the surface by day, such as a daily true-colour mosaic from VIIRS, or
  the cloud-free Blue Marble.
- **night-lights**: the Earth's lights at night (Black Marble), to be added on
  the night side.
- **relief-shading**: grey shaded relief of the land (ASTER GDEM), whose detail
  at close range replaces the coarser shading of the relief map. Its grey
  levels are measured: water (and land without data) is one uniform grey, the
  most frequent; flat land is the most frequent grey of the rest. An
  application divides by the flat level to obtain a shading factor, and takes
  the uniform grey, where a whole neighbourhood has it, as water.

A daily mosaic by UTC date leaves gaps where the satellite passed on the
neighbouring day (a wedge near the date line) and in polar night. The gaps of
the chosen day are filled from the day before. Pixels darker than 6 in every
channel are taken as gaps; the true colour of the darkest ocean is well above
that.

Each product also records how to ask the provider for a part of the same
imagery at higher resolution (a WMS request), so that an application can show
detail where the user looks.

Output, in data/cache/products/: imagery-<name>.jpg and imagery-<name>.json

Usage::

    python data/products/imagery.py --source gibs-daily
    python data/products/imagery.py --source gibs-static --profile blue-marble
    python data/products/imagery.py --source gibs-static --profile black-marble
    python data/products/imagery.py --source gibs-relief
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

Image.MAX_IMAGE_PIXELS = 400_000_000
GAP = 6
NAMES = {("gibs-daily", None): ("true-colour", "day"),
         ("gibs-static", "blue-marble"): ("blue-marble", "day"),
         ("gibs-static", "black-marble"): ("night-lights", "night-lights"),
         ("gibs-relief", "preview"): ("relief-shading", "relief-shading")}


def build(source: str, profile: str) -> dict:
    manifest = json.loads((DATA / "registry" / f"{source}.json").read_text(encoding="utf-8"))
    folder = DATA / "cache" / source / profile
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    name, role = NAMES.get((source, None if source == "gibs-daily" else profile), (f"{source}-{profile}", "day"))
    files = [folder / r["file"] for r in retrieval["requests"]]
    if role == "relief-shading":
        return build_shading(manifest, retrieval, files[0], name, role, profile)
    first = np.asarray(Image.open(files[0]).convert("RGB"))
    filled = 0.0
    if len(files) > 1:
        image = first.copy()
        gaps = (image.max(axis=2) < GAP)
        for path in files[1:]:
            other = np.asarray(Image.open(path).convert("RGB"))
            if other.shape != image.shape:
                continue
            use = gaps & (other.max(axis=2) >= GAP)
            image[use] = other[use]
            gaps &= ~use
        filled = float((first.max(axis=2) < GAP).mean() - gaps.mean())
        remaining = float(gaps.mean())
    else:
        image = first
        remaining = float((image.max(axis=2) < GAP).mean()) if role == "day" else 0.0
    height, width = image.shape[:2]
    out = OUT / f"imagery-{name}.jpg"
    OUT.mkdir(parents=True, exist_ok=True)
    Image.fromarray(image).save(out, quality=92, optimize=True)

    params = retrieval["parameters"]
    resolved = retrieval.get("resolved", {})
    wms = retrieval["requests"][0]["url"].split("&BBOX=")[0]
    live = {
        "kind": "wms",
        "url": wms + "&LAYERS=" + params.get("layer", ""),
        "time": (resolved.get("date") or "")[:10] or None,
        "format": "image/jpeg",
        "crs": "EPSG:4326",
        "bbox_order": "lat_min,lon_min,lat_max,lon_max",
    }
    product = {
        "schema": "imagery/1",
        "name": name,
        "title": manifest["title"].split(":")[0],
        "role": role,
        "layer": params.get("layer"),
        "date": live["time"],
        "projection": "equirectangular (EPSG:4326)",
        "extent": [-180, -90, 180, 90],
        "width": width, "height": height,
        "file": out.name,
        "gaps": {"filled": round(filled, 4), "remaining": round(remaining, 4)},
        "live": live,
        "source": {"dataset": manifest["id"], "profile": profile, "provider": manifest["provider"],
                   "credit": manifest.get("credit"), "licence": manifest["licence"],
                   "citation": manifest["citation"], "retrieved_at": retrieval["retrieved_at"],
                   "request": retrieval["requests"][0]["url"]},
    }
    path = OUT / f"imagery-{name}.json"
    path.write_text(json.dumps(product, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    register("imagery", {k: product[k] for k in ("name", "title", "role", "layer", "date", "width", "height")}
             | {"file": path.name, "source": {k: product["source"][k] for k in ("dataset", "profile", "credit", "retrieved_at")}})
    return product


def build_shading(manifest: dict, retrieval: dict, path: Path, name: str, role: str, profile: str) -> dict:
    """A product of grey shaded relief, with its grey levels and its service for detail."""
    grey = np.asarray(Image.open(path).convert("L"))
    counts = np.bincount(grey.ravel(), minlength=256)
    water = int(counts.argmax())
    land = counts.copy()
    land[max(0, water - 2):water + 3] = 0
    flat = int(land.argmax())
    out = OUT / f"imagery-{name}.png"
    OUT.mkdir(parents=True, exist_ok=True)
    Image.fromarray(grey).save(out, optimize=True)
    params = retrieval["parameters"]
    url = retrieval["requests"][0]["url"]
    fmt = url.split("FORMAT=", 1)[1].split("&", 1)[0] if "FORMAT=" in url else "image/png"
    product = {
        "schema": "imagery/1", "name": name, "title": manifest["title"].split(",")[0], "role": role,
        "layer": params.get("layer"), "date": None,
        "projection": "equirectangular (EPSG:4326)", "extent": [-180, -90, 180, 90],
        "width": int(grey.shape[1]), "height": int(grey.shape[0]), "file": out.name,
        "levels": {"water": water, "flat": flat, "land_share_of_pixels": round(float(1 - counts[water] / grey.size), 4)},
        "live": {"kind": "wms", "url": url.split("&BBOX=")[0] + "&LAYERS=" + params.get("layer", ""), "time": None,
                 "format": fmt, "crs": "EPSG:4326", "bbox_order": "lat_min,lon_min,lat_max,lon_max"},
        "source": {"dataset": manifest["id"], "profile": profile, "provider": manifest["provider"],
                   "credit": manifest.get("credit"), "licence": manifest["licence"],
                   "citation": manifest["citation"], "retrieved_at": retrieval["retrieved_at"], "request": url},
    }
    (OUT / f"imagery-{name}.json").write_text(json.dumps(product, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    register("imagery", {k: product[k] for k in ("name", "title", "role", "layer", "date", "width", "height")}
             | {"file": f"imagery-{name}.json", "source": {k: product["source"][k] for k in ("dataset", "profile", "credit", "retrieved_at")}})
    return product


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", required=True, help="gibs-daily, gibs-static or gibs-relief")
    parser.add_argument("--profile", help="the profile (cache folder); default: the dataset's")
    args = parser.parse_args(argv)
    manifest = json.loads((DATA / "registry" / f"{args.source}.json").read_text(encoding="utf-8"))
    profile = args.profile or manifest.get("default_profile")
    if not (DATA / "cache" / args.source / profile / "retrieval.json").exists():
        print(f"not fetched; run: python data/acquisition/fetch.py {args.source} --profile {profile}", file=sys.stderr)
        return 1
    p = build(args.source, profile)
    print(f"imagery-{p['name']}: {p['width']} x {p['height']}, role {p['role']}"
          + (f", {p['date']}" if p["date"] else "")
          + (f"; gaps filled {p['gaps']['filled']:.1%}, remaining {p['gaps']['remaining']:.1%}" if p["role"] == "day" else "")
          + (f"; grey of water {p['levels']['water']}, of flat land {p['levels']['flat']}" if p["role"] == "relief-shading" else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
