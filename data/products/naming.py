"""The names of places in the data library's products, after the policy in naming.json (ADR 0069).

Places are named by Natural Earth's long names, except where a name is in dispute or out of date:
then by the short name of the UN's M49, and for territories M49 does not list, or lists under a name
in dispute, by the usage of the WTO, the OECD, the World Bank and the IMF, whose statistics the
infrastructure draws. A map draws no name for a territory, or a sea, whose name or status is in
dispute; a territory without a name in that usage is given none. Names current before are aliases,
so that data of one's own written with them still resolve. An area Natural Earth lists as in dispute
or breakaway is drawn in no claimant's colour and without a name; a place in it has no country, and
no place in it, or in a territory not labelled, is listed as a national capital.

    from naming import POLICY
    POLICY.code(properties)       the code a Natural Earth feature is known by
    POLICY.name(code, default)    the name to show (None: none)
    POLICY.labelled(code)         whether a map may draw the name
    POLICY.local(code, names)     a feature's names by language, after the policy
    POLICY.resolve(text)          the code a name or an alias stands for, or None
    POLICY.sea_labelled(name)     whether a map may draw a sea's name
    POLICY.place(name, kind, codes, lon, lat, areas)
                                  a populated place's kind, and whether it keeps its country
    AREAS.contains(lon, lat)      whether a point lies in an area in dispute (None if not fetched)

    python data/products/naming.py      (lists the policy)
"""

from __future__ import annotations

import json
from pathlib import Path

FILE = Path(__file__).resolve().with_name("naming.json")
DISPUTED = Path(__file__).resolve().parents[2] / "data" / "cache" / "natural-earth-50m" / "ne_50m_admin_0_breakaway_disputed_areas.geojson"
CAPITALS = ("Admin-0 capital", "Admin-0 capital alt", "Admin-0 region capital")


class Areas:
    """Natural Earth's areas in dispute or breakaway, as polygons in longitude and latitude, for
    asking whether a point lies in one (by the even-odd rule over each polygon's rings). Their names
    are not kept: they are not drawn."""

    def __init__(self, features: list[dict]):
        self.polygons = []
        for f in features:
            g = f["geometry"]
            for poly in g["coordinates"] if g["type"] == "MultiPolygon" else [g["coordinates"]]:
                xs = [x for ring in poly for x, _ in ring]
                ys = [y for ring in poly for _, y in ring]
                self.polygons.append(((min(xs), min(ys), max(xs), max(ys)), poly))

    @classmethod
    def load(cls, path: Path = DISPUTED) -> "Areas | None":
        return cls(json.loads(path.read_text(encoding="utf-8"))["features"]) if path.exists() else None

    def contains(self, lon: float, lat: float) -> bool:
        for (x0, y0, x1, y1), poly in self.polygons:
            if not (x0 <= lon <= x1 and y0 <= lat <= y1):
                continue
            inside = False
            for ring in poly:
                for (ax, ay), (bx, by) in zip(ring, ring[1:] + ring[:1]):
                    if (ay > lat) != (by > lat) and lon < (bx - ax) * (lat - ay) / (by - ay) + ax:
                        inside = not inside
            if inside:
                return True
        return False


class Naming:
    def __init__(self, data: dict):
        self.data = data
        self.names: dict[str, dict] = data.get("names", {})
        self.seas: dict[str, dict] = data.get("seas", {})
        self.places: dict[str, dict] = data.get("places", {})
        self.aliases: dict[str, str] = {}
        for code, entry in self.names.items():
            for alias in [entry.get("name"), *entry.get("aliases", [])]:
                if alias:
                    self.aliases.setdefault(alias.casefold(), code)

    @classmethod
    def load(cls, path: Path = FILE) -> "Naming":
        return cls(json.loads(path.read_text(encoding="utf-8")))

    @staticmethod
    def code(properties: dict) -> str | None:
        """ISO_A3, else ISO_A3_EH, else ADM0_A3 (Natural Earth's own code for a territory without one)."""
        for key in ("ISO_A3", "ISO_A3_EH", "ADM0_A3"):
            value = properties.get(key)
            if value and value != "-99" and len(str(value)) == 3:
                return str(value)
        return None

    def entry(self, *codes: str | None) -> dict | None:
        for c in codes:
            if c and c in self.names:
                return self.names[c]
        return None

    def name(self, code: str | None, default: str | None, *more: str | None) -> str | None:
        e = self.entry(code, *more)
        return default if e is None else e.get("name")

    def labelled(self, code: str | None, *more: str | None) -> bool:
        e = self.entry(code, *more)
        return True if e is None else e.get("label", True) is not False and e.get("name") is not None

    def local(self, code: str | None, names: dict, *more: str | None) -> dict:
        """A feature's names by language ({"default", "en", "zh", …}) after the policy: unchanged
        where the policy says nothing; the policy's name for default and English otherwise, and its
        own names by language where it gives them (the others then fall back to English); none at
        all where the policy gives no name."""
        e = self.entry(code, *more)
        if e is None:
            return names
        if e.get("name") is None:
            return {}
        out = dict(names) if "local" not in e else {}
        out.update({"default": e["name"], "en": e["name"]})
        if "long" in names:
            out["long"] = e["name"]
        for k in ("formal",):
            out.pop(k, None)
        out.update(e.get("local", {}))
        return out

    def resolve(self, text: str) -> str | None:
        return self.aliases.get((text or "").strip().casefold())

    def sea_labelled(self, name: str | None) -> bool:
        e = self.seas.get(name or "")
        return True if e is None else e.get("label", True) is not False

    def place(self, name: str, kind: str, codes: tuple, lon: float, lat: float, areas: "Areas | None") -> tuple[str, bool]:
        """A populated place's kind (Natural Earth's FEATURECLA) and whether it keeps its country:
        as the policy's places say where they name it; otherwise no country in an area in dispute,
        and a populated place, not a national capital, there or in a territory not labelled."""
        e = self.places.get(name)
        if e is not None:
            return e.get("kind", kind), not ("country" in e and e["country"] is None)
        inside = areas is not None and areas.contains(lon, lat)
        if kind in CAPITALS and (inside or not self.labelled(*codes)):
            kind = "Populated place"
        return kind, not inside


POLICY = Naming.load()
AREAS = Areas.load()


if __name__ == "__main__":
    for line in POLICY.data["policy"]:
        print("-", line)
    for code, e in POLICY.names.items():
        print(f"{code}: {e.get('name')!r}{'' if e.get('label', True) else ', not labelled'}")
    for name, e in POLICY.seas.items():
        print(f"sea {name!r}{'' if e.get('label', True) else ', not labelled'}")
    for name, e in POLICY.places.items():
        print(f"place {name!r}: {e.get('kind')}, {'no country' if e.get('country', 1) is None else 'its country'}")
    print(f"areas in dispute: {len(AREAS.polygons) if AREAS else 'not fetched'} polygons")
