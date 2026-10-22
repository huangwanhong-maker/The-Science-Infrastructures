"""Build imagery products of the Sun, the Moon, the planets and the Milky Way.

Each map becomes an imagery product (see imagery.py) with a role:
- `body`: the surface or cloud tops of a body, in east longitude, 0 at the
  centre of the image, north at the top;
- `rings`: Saturn's rings as a radial strip with transparency;
- `sky`: the whole sky in galactic coordinates.

Maps wider than 4096 pixels are resampled to 4096, which a graphics card holds
easily; the sky keeps 8192 (it spans the whole field of view). A body whose
manifest entry names a tile service (`live`: NASA's Trek WMTS for the Moon,
Mars and Mercury) passes it on, so that an application can show detail at
close range.

The sky map's orientation is checked rather than assumed. The brightest
pixels near twelve bright stars of known galactic position are read under the
stated orientation (longitude increasing to the left, the south galactic pole
at the top). The product records the result, and the build fails if the stars
are not where they should be.

Output, in data/cache/products/: imagery-<role>-<name>.jpg|png and .json.

Usage::

    python data/products/body_maps.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

DATASET = "solar-system-maps"
MAX_WIDTH = {"body": 4096, "rings": 4096, "sky": 8192}
# Galactic longitude and latitude (degrees) of bright stars, from Hipparcos (ESA 1997).
BRIGHT = {"Sirius": (227.23, -8.89), "Canopus": (261.21, -25.29), "Arcturus": (15.09, 69.11),
          "Vega": (67.45, 19.24), "Rigel": (209.24, -25.25), "Betelgeuse": (199.79, -8.96),
          "Capella": (162.59, 4.57), "Altair": (47.74, -8.91), "Antares": (351.95, 15.06),
          "Spica": (316.11, 50.84), "Deneb": (84.28, 1.99), "Alpha Centauri": (315.73, -0.68)}


def sky_check(image: np.ndarray) -> dict:
    """Mean peak brightness near the bright stars, against the image's mean."""
    grey = image.astype(float).mean(axis=2)
    h, w = grey.shape
    peaks = []
    for l, b in BRIGHT.values():
        x = int((0.5 - ((l + 180) % 360 - 180) / 360) * w) % w
        y = min(h - 1, int((90 + b) / 180 * h))                 # the south galactic pole at the top
        peaks.append(float(grey[max(0, y - 6):y + 7, max(0, x - 6):x + 7].max()))
    return {"stars": len(peaks), "mean_peak": round(float(np.mean(peaks)), 1),
            "background": round(float(grey.mean()), 1)}


def build() -> list[dict]:
    manifest = json.loads((DATA / "registry" / f"{DATASET}.json").read_text(encoding="utf-8"))
    folder = DATA / "cache" / DATASET
    retrieval = json.loads((folder / "retrieval.json").read_text(encoding="utf-8"))
    out = []
    for entry in manifest["files"]:
        body = entry["body"]
        role = "sky" if body == "milky-way" else "rings" if body.endswith("rings") else "body"
        image = Image.open(folder / entry["file"])
        image = image.convert("RGBA" if role == "rings" else "RGB")
        if image.width > MAX_WIDTH[role]:
            image = image.resize((MAX_WIDTH[role], round(image.height * MAX_WIDTH[role] / image.width)), Image.LANCZOS)
        name = f"{role}-{body}"
        ext = "png" if role == "rings" else "jpg"
        path = OUT / f"imagery-{name}.{ext}"
        OUT.mkdir(parents=True, exist_ok=True)
        if ext == "png":
            image.save(path, optimize=True)
        else:
            image.save(path, quality=90, optimize=True)
        product = {
            "schema": "imagery/1", "name": name, "title": f"{body.replace('-', ' ').capitalize()}, map", "role": role,
            "body": body if role != "sky" else None,
            "width": image.width, "height": image.height, "file": path.name,
            "source": {"dataset": DATASET, "provider": manifest["provider"], "credit": manifest.get("credit"),
                       "licence": manifest["licence"], "citation": manifest["citation"],
                       "retrieved_at": retrieval["retrieved_at"], "request": entry["url"]},
        }
        if role == "body":
            product["projection"] = "equirectangular; east longitude -180 to 180 from left to right, 0 at the centre; north at the top"
            if entry.get("live"):
                product["live"] = entry["live"]
        elif role == "rings":
            product["projection"] = "radial strip: from the inner edge of the C ring (left) to the outer edge of the A ring (right)"
            product["radii_km"] = [74500, 136780]
        else:
            product["projection"] = "equirectangular, galactic coordinates; longitude 0 at the centre, increasing to the left; the south galactic pole at the top"
            product["frame"] = "galactic"
            check = sky_check(np.asarray(image))
            product["orientation_check"] = check
            if check["mean_peak"] < 20 * max(1.0, check["background"]):
                raise SystemExit(f"the sky map's stars are not where its stated orientation puts them: {check}")
        (OUT / f"imagery-{name}.json").write_text(json.dumps(product, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
        register("imagery", {k: product.get(k) for k in ("name", "title", "role", "body", "width", "height")}
                 | {"file": f"imagery-{name}.json", "source": {"dataset": DATASET, "credit": manifest.get("credit")}})
        out.append(product)
    return out


def main(argv: list[str]) -> int:
    if not (DATA / "cache" / DATASET / "retrieval.json").exists():
        print(f"not fetched; run: python data/acquisition/fetch.py {DATASET}", file=sys.stderr)
        return 1
    for p in build():
        extra = f", orientation check {p['orientation_check']}" if "orientation_check" in p else ""
        print(f"imagery-{p['name']}: {p['width']} x {p['height']}{extra}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
