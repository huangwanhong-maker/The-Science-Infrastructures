"""Place names of the Earth and the Moon as a product: a gazetteer for searching a globe.

- **The Earth.** Natural Earth at 1:50 million (public domain): its populated places (1251, each
  with its country and largest population) and its countries at their label points. Coordinates
  are geodetic on WGS 84, degrees, east positive.
- **The Moon.** The IAU's Gazetteer of Planetary Nomenclature (USGS Astrogeology, public domain):
  every feature whose name the IAU has adopted, at its centre, with its type, diameter and the
  link to its page. Coordinates are planetocentric on the sphere of 1737.4 km in the mean-Earth
  frame (the gazetteer's "Moon 2000"), east positive; the gazetteer's 0–360 E longitudes are given
  here from −180 to 180.

Output: data/cache/products/places/earth.json and moon.json, each
    {"schema": "places/1", "body", "frame", "source", "fields": [...], "places": [[...], ...]}
with places as rows [name, lon, lat, kind, detail, size_km, link], sorted by importance (a capital
or a large feature first), and the index of products updated.

    python data/products/place_names.py
"""

from __future__ import annotations

import hashlib
import io
import json
import struct
import sys
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "data" / "cache"
OUT = CACHE / "products" / "places"
NE = CACHE / "natural-earth-50m"
MOON = CACHE / "usgs-moon-nomenclature" / "MOON_nomenclature_center_pts.zip"
FIELDS = ["name", "lon", "lat", "kind", "detail", "size_km", "link"]

sys.path.insert(0, str(Path(__file__).resolve().parent))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_dbf(data: bytes) -> list[dict]:
    """The records of a dBASE III table (the attribute table of a shapefile), as strings."""
    count, header_length, record_length = struct.unpack("<IHH", data[4:12])
    fields, pos = [], 32
    while data[pos] != 0x0D:
        f = data[pos:pos + 32]
        fields.append((f[:11].split(b"\0")[0].decode("ascii"), f[16]))
        pos += 32
    out = []
    for k in range(count):
        r = data[header_length + k * record_length: header_length + (k + 1) * record_length]
        if r[:1] == b"*":  # deleted
            continue
        o, row = 1, {}
        for name, length in fields:
            row[name] = r[o:o + length].decode("utf-8", "replace").strip()
            o += length
        out.append(row)
    return out


def earth() -> dict:
    rows = []
    places = json.loads((NE / "ne_50m_populated_places.geojson").read_text(encoding="utf-8"))["features"]
    for f in places:
        p = f["properties"]
        rows.append([p["NAME"], round(float(p["LONGITUDE"]), 5), round(float(p["LATITUDE"]), 5), p["FEATURECLA"],
                     p.get("ADM0NAME") or "", None, None, int(p.get("POP_MAX") or 0), int(p.get("SCALERANK") or 10)])
    countries = json.loads((NE / "ne_50m_admin_0_countries.geojson").read_text(encoding="utf-8"))["features"]
    for f in countries:
        p = f["properties"]
        if p.get("LABEL_X") is None:
            continue
        rows.append([p.get("NAME_LONG") or p["NAME"], round(float(p["LABEL_X"]), 5), round(float(p["LABEL_Y"]), 5), "Country",
                     p.get("CONTINENT") or "", None, None, int(p.get("POP_EST") or 0) * 10, -1])
    # Countries, then capitals and the most populous places.
    rows.sort(key=lambda r: (r[8], -r[7]))
    return {"schema": "places/1", "body": "earth",
            "frame": "geodetic longitude and latitude on WGS 84, degrees, east positive",
            "source": {"dataset": "natural-earth-50m", "files": ["ne_50m_populated_places.geojson", "ne_50m_admin_0_countries.geojson"],
                       "sha256": [sha256(NE / "ne_50m_populated_places.geojson"), sha256(NE / "ne_50m_admin_0_countries.geojson")],
                       "licence": "public domain (Natural Earth)"},
            "fields": FIELDS, "places": [r[:7] for r in rows]}


def moon() -> dict:
    with zipfile.ZipFile(MOON) as z:
        records = read_dbf(z.read("MOON_nomenclature_center_pts.dbf"))
    rows = []
    for r in records:
        if r.get("approval") != "Adopted by IAU":
            continue
        lon = float(r["center_lon"])
        lon = lon - 360.0 if lon > 180.0 else lon
        size = float(r["diameter"]) if r.get("diameter") else None
        rows.append([r["clean_name"] or r["name"], round(lon, 5), round(float(r["center_lat"]), 5), r["type"].split(",")[0],
                     r.get("origin", "")[:160], round(size, 3) if size is not None else None, r.get("link", "").replace("http://", "https://")])
    rows.sort(key=lambda r: -(r[5] or 0))
    return {"schema": "places/1", "body": "moon",
            "frame": "planetocentric longitude and latitude on the sphere of 1737.4 km, mean-Earth frame (MOON_ME), degrees, east positive",
            "source": {"dataset": "usgs-moon-nomenclature", "files": [MOON.name], "sha256": [sha256(MOON)],
                       "licence": "public domain (USGS)", "citation": "Gazetteer of Planetary Nomenclature, IAU WGPSN and USGS Astrogeology"},
            "fields": FIELDS, "places": rows}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    built = []
    for name, make, needs in (("earth", earth, NE / "ne_50m_populated_places.geojson"), ("moon", moon, MOON)):
        if not needs.exists():
            print(f"{needs} is missing: fetch its dataset first", file=sys.stderr)
            continue
        product = make()
        product["built"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        (OUT / f"{name}.json").write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
        print(f"wrote {OUT / (name + '.json')}: {len(product['places'])} places")
        built.append(name)
    try:
        from product_index import register
        for name in built:
            register("places", {"name": f"places/{name}", "title": f"Place names of the {'Earth' if name == 'earth' else 'Moon'}",
                                "body": name, "file": f"places/{name}.json"})
    except Exception as error:  # the index is a convenience; the product stands without it
        print(f"(the index of products was not updated: {error})", file=sys.stderr)
    return 0 if len(built) == 2 else 1


if __name__ == "__main__":
    sys.exit(main())
