"""A mask of land and sea as a product: Natural Earth's land, rasterised.

The terrain tiles of the Earth (data/tiles/terrain.py) draw the sea at its surface and keep its
floor beneath: a sample below sea level is sea when it lies outside the land. The sources of
heights cannot say so themselves, because land too lies below sea level (the Dead Sea's shore at
−430 m, the Caspian depression, the polders of the Netherlands). This builder rasterises Natural
Earth's land polygons at 1:50 million (ne_50m_land, in the public domain), with the Caspian Sea
(the polygons' one hole) kept as land, onto a grid of 20 pixels per degree from 180 W and 90 N
(cell centres at half-pixel offsets), and erodes the land by one
pixel (about 5.5 km at the equator), so that a sample below sea level within a pixel of the coast
is taken for sea: the coast is then where the heights cross sea level, as the sources draw it,
and the polygons' generalisation at this scale (a kilometre or two) decides only inland.

Output: data/cache/products/land-mask/ne50m-20ppd.png (8-bit greyscale, 7200 x 3600; 255 inland, 0
elsewhere) and its description beside it.

    python data/products/land_mask.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache"
SOURCE = CACHE / "natural-earth-50m" / "ne_50m_land.geojson"
OUT = CACHE / "products" / "land-mask"
PPD = 20


def rings(geometry: dict):
    """The outer rings and holes of a Polygon or MultiPolygon, as (ring, is_hole)."""
    polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
    for polygon in polygons:
        for k, ring in enumerate(polygon):
            yield ring, k > 0


def main() -> int:
    if not SOURCE.exists():
        print(f"{SOURCE} is missing: python data/acquisition/fetch.py natural-earth-50m", file=sys.stderr)
        return 1
    width, height = 360 * PPD, 180 * PPD
    image = Image.new("L", (width, height), 0)
    draw = ImageDraw.Draw(image)
    features = json.loads(SOURCE.read_text(encoding="utf-8"))["features"]
    # Only the outer rings: the one hole Natural Earth cuts out of its land at this scale is the
    # Caspian Sea, a lake whose surface lies 28 m below sea level, kept as land so that its surface
    # is the sources' and not the geoid.
    holes = 0
    for f in features:
        for ring, hole in rings(f["geometry"]):
            xy = [((lon + 180.0) * PPD, (90.0 - lat) * PPD) for lon, lat in ring]
            if hole:
                holes += 1
            elif len(xy) >= 3:
                draw.polygon(xy, fill=255)
    inland = image.filter(ImageFilter.MinFilter(3))
    OUT.mkdir(parents=True, exist_ok=True)
    name = f"ne50m-{PPD}ppd"
    inland.save(OUT / f"{name}.png", optimize=True)
    (OUT / f"{name}.json").write_text(json.dumps({
        "schema": "land-mask/1",
        "name": name,
        "title": "Land, inland of the coast by one pixel (Natural Earth 1:50m, rasterised)",
        "file": f"{name}.png",
        "format": f"8-bit greyscale PNG, {width} x {height}, {PPD} pixels per degree from 180 W and 90 N (pixel-registered); 255 inland, 0 elsewhere",
        "erosion": "one pixel (3 x 3 minimum filter)",
        "holes": f"{holes} hole(s) of the polygons kept as land (the Caspian Sea, a lake below sea level)",
        "source": {"dataset": "natural-earth-50m", "file": SOURCE.name, "sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                   "licence": "public domain (Natural Earth)"},
        "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }, indent=1), encoding="utf-8")
    print(f"wrote {OUT / (name + '.png')} ({width} x {height})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
