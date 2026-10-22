"""A minimal client of the SpatioTemporal Asset Catalog (STAC) API: item search.

A STAC API (STAC API Specification 1.0.0, 2023, "Item Search") answers a POST
to <endpoint>/search whose JSON body names collections, a geometry the items
must intersect (GeoJSON, RFC 7946: longitude first, WGS 84), a datetime
interval (RFC 3339, "start/end"), a filter on properties (the Query extension
1.0.0, as {"eo:cloud_cover": {"lt": 20}}), an order (the Sort extension 1.0.0,
[{"field": "properties.eo:cloud_cover", "direction": "asc"}]) and a page size.
The answer is a GeoJSON FeatureCollection of items, each with its properties
and its assets (files, by key, with their addresses and metadata such as the
raster extension's scale and offset); a link of relation "next" gives the
next page, as a request to make: its method, address and body, to be merged
into the first body or to replace it ("merge").

Requests are built by search_body() alone, so that they can be checked
without a network; search() sends them and follows the pages.

The first catalogue is Element 84's Earth Search, over the Sentinel-2 archive
on Amazon's Registry of Open Data (EARTH_SEARCH, SENTINEL_2_L2A). Copernicus
Sentinel data are free, full and open to everyone, for any use, without warranty
(Commission Delegated Regulation (EU) No 1159/2013, Articles 3 and 9, as the
Sentinel Data Legal Notice states them); what is made from them carries the
notice "Contains modified Copernicus Sentinel data <year>".

Usage::

    python data/acquisition/stac.py 140.0 36.1                      # the scenes of the last 30 days over a place
    python data/acquisition/stac.py 140.0 36.1 --days 60 --cloud 30
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import json
import sys
import time
import urllib.error
import urllib.request

AGENT = "science-infrastructures-data/1 (+https://github.com/huangwanhong-maker/The-Science-Infrastructures)"
UTC = dt.timezone.utc

EARTH_SEARCH = "https://earth-search.aws.element84.com/v1"
SENTINEL_2_L2A = {
    "collection": "sentinel-2-l2a",
    "title": "Sentinel-2 Level-2A, surface reflectance (Earth Search)",
    "provider": "ESA Copernicus; Cloud-Optimised GeoTIFFs on the Registry of Open Data on AWS, catalogued by Element 84 (Earth Search)",
    "licence": ("Copernicus Sentinel data: free, full and open access, without warranty (Commission Delegated "
                "Regulation (EU) No 1159/2013, Articles 3 and 9; Sentinel Data Legal Notice, "
                "https://sentinel.esa.int/documents/247904/690755/Sentinel_Data_Legal_Notice)"),
    "credit": "Contains modified Copernicus Sentinel data {year}",
    # The collection's "cite-as" link, doi:10.5270/S2_-742ikth, resolves to Level-1C; this is Level-2A's.
    "citation": "European Space Agency (2022), Sentinel-2 MSI Level-2A BOA Reflectance, doi:10.5270/S2_-znk9xsj",
}


def rfc3339(when: dt.datetime) -> str:
    return when.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def interval(start: dt.datetime | None, end: dt.datetime | None) -> str:
    """A datetime interval of the Item Search; '..' leaves an end open."""
    return f"{rfc3339(start) if start else '..'}/{rfc3339(end) if end else '..'}"


def search_body(collections: list[str], *, intersects: dict | None = None, bbox: list[float] | None = None,
                start: dt.datetime | None = None, end: dt.datetime | None = None, query: dict | None = None,
                sortby: list[tuple[str, str]] | None = None, limit: int = 50) -> dict:
    """The JSON body of an Item Search. intersects: a GeoJSON geometry; bbox:
    [west, south, east, north] in degrees (one of the two); sortby: [(property,
    'asc' | 'desc')], properties named as in the items ('eo:cloud_cover')."""
    if intersects and bbox:
        raise ValueError("give intersects or bbox, not both (STAC API, Item Search)")
    body: dict = {"collections": list(collections), "limit": int(limit)}
    if intersects:
        body["intersects"] = intersects
    if bbox:
        body["bbox"] = [float(x) for x in bbox]
    if start or end:
        body["datetime"] = interval(start, end)
    if query:
        body["query"] = query
    if sortby:
        body["sortby"] = [{"field": f if f.startswith("properties.") or f in ("id", "collection") else f"properties.{f}",
                           "direction": d} for f, d in sortby]
    return body


def point(lon: float, lat: float) -> dict:
    return {"type": "Point", "coordinates": [float(lon), float(lat)]}


def request_json(url: str, body: dict | None = None, method: str | None = None, timeout: int = 60) -> dict:
    """A JSON answer to a GET, or to a POST of `body`; one more attempt after a
    failure of the server or of the network."""
    data = json.dumps(body).encode("utf-8") if body is not None else None
    headers = {"User-Agent": AGENT, "Accept": "application/geo+json, application/json"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    for attempt in (1, 2):
        request = urllib.request.Request(url, data=data, headers=headers, method=method or ("POST" if data else "GET"))
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            if error.code < 500 or attempt == 2:
                detail = error.read().decode("utf-8", "replace")[:300]
                raise RuntimeError(f"{url}: {error.code} {error.reason}: {detail}") from None
        except urllib.error.URLError:
            if attempt == 2:
                raise
        time.sleep(3)
    raise RuntimeError("unreachable")


def next_request(page: dict, body: dict) -> tuple[str, str, dict | None] | None:
    """The next page's (method, address, body), from the link of relation 'next', or None."""
    link = next((l for l in page.get("links", []) if l.get("rel") == "next"), None)
    if not link:
        return None
    method = link.get("method", "GET").upper()
    if method == "POST":
        following = link.get("body") or {}
        return method, link["href"], ({**copy.deepcopy(body), **following} if link.get("merge") else following)
    return method, link["href"], None


def search(endpoint: str, body: dict, max_items: int = 200, timeout: int = 60) -> list[dict]:
    """The items an Item Search finds, page by page, up to max_items."""
    items: list[dict] = []
    method, url, current = "POST", endpoint.rstrip("/") + "/search", body
    while True:
        page = request_json(url, current if method == "POST" else None, method, timeout)
        items += page.get("features", [])
        following = next_request(page, current or body)
        if len(items) >= max_items or not following or not page.get("features"):
            return items[:max_items]
        method, url, current = following


# ---- items ------------------------------------------------------------------


def asset_href(item: dict, key: str) -> str:
    return item["assets"][key]["href"]


def raster_band(item: dict, key: str) -> dict:
    """The raster extension's first band of an asset: nodata, data_type, scale, offset (absent ones omitted)."""
    bands = item["assets"][key].get("raster:bands") or [{}]
    return bands[0]


def when(item: dict) -> dt.datetime:
    return dt.datetime.fromisoformat(item["properties"]["datetime"].replace("Z", "+00:00")).astimezone(UTC)


def _in_ring(lon: float, lat: float, ring: list) -> bool:
    """Whether a place lies inside a ring of [lon, lat] vertices, by counting crossings
    of a ray toward the east; in the plane of longitude and latitude, which suits a
    footprint of a few hundred km that does not cross the antimeridian."""
    inside = False
    for (x1, y1), (x2, y2) in zip(ring, ring[1:] + ring[:1]):
        if (y1 > lat) != (y2 > lat) and lon < x1 + (lat - y1) * (x2 - x1) / (y2 - y1):
            inside = not inside
    return inside


def contains(geometry: dict, lon: float, lat: float) -> bool:
    """Whether a GeoJSON Polygon or MultiPolygon contains a place (holes excepted)."""
    polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"] if geometry["type"] == "MultiPolygon" else []
    for rings in polygons:
        rings = [[tuple(p[:2]) for p in r] for r in rings]
        if _in_ring(lon, lat, rings[0]) and not any(_in_ring(lon, lat, h) for h in rings[1:]):
            return True
    return False


def covers(item: dict, west: float, south: float, east: float, north: float, n: int = 4) -> bool:
    """Whether an item's footprint contains a box: its corners and n points along each edge."""
    edge = [k / n for k in range(n + 1)]
    places = ([(west + t * (east - west), south) for t in edge] + [(west + t * (east - west), north) for t in edge]
              + [(west, south + t * (north - south)) for t in edge] + [(east, south + t * (north - south)) for t in edge])
    return all(contains(item["geometry"], lon, lat) for lon, lat in places)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("lon", type=float)
    parser.add_argument("lat", type=float)
    parser.add_argument("--days", type=float, default=30)
    parser.add_argument("--cloud", type=float, default=100, help="the largest cloud cover of the tile, percent")
    args = parser.parse_args(argv)
    end = dt.datetime.now(UTC)
    body = search_body([SENTINEL_2_L2A["collection"]], intersects=point(args.lon, args.lat),
                       start=end - dt.timedelta(days=args.days), end=end,
                       query={"eo:cloud_cover": {"lte": args.cloud}}, sortby=[("eo:cloud_cover", "asc")])
    items = search(EARTH_SEARCH, body)
    for item in items:
        p = item["properties"]
        print(f"{item['id']:28s} {p['datetime'][:19]}  cloud {p.get('eo:cloud_cover', float('nan')):6.2f} %  "
              f"baseline {p.get('s2:processing_baseline')}  {p.get('grid:code')}")
    print(f"{len(items)} scene(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
