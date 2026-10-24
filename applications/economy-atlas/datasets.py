"""Datasets of one's own for the Economy Atlas (ADR 0046).

A dataset is given as CSV or JSON and kept as a product of the same schema as the data library's,
so that the engine reads it unchanged and every view and plug-in works on it as on the library's
own:

Flows between places, kept as ``trade-flows/1`` (the map, its strip of years, and the plug-ins):

    places.csv   id,name,latitude,longitude[,group]      degrees on WGS 84; optional for ids that are
                                                         the ISO alpha-3 codes or the names of the
                                                         library's economies
    flows.csv    from,to,value[,year]                    ids; a flow without a year belongs to year 0
    or one JSON  {"title", "unit", "places": [{id, name, latitude, longitude, group}],
                  "flows": [{from, to, value, year}]}

Currencies held and exchanged, kept as ``currency-spacetime/1`` (the spacetime view):

    rates.csv    month,currency,rate,policy              month YYYY-MM; rate: units of the currency for
                                                         one unit of the base; policy: its holding rate,
                                                         per cent a year; the base's rows give its policy
                                                         rate (its rate is 1 or empty)
    or one JSON  {"title", "base", "rows": [{month, currency, rate, policy}]}

Kept under data/cache/economy-atlas/datasets/<id>/: ``product.json`` and ``about.json``. The
library's own products are listed beside them as ``library:trade`` and ``library:currencies``.

    python applications/economy-atlas/datasets.py add flows "My flows" places.csv flows.csv
    python applications/economy-atlas/datasets.py list
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import io
import json
import math
import re
import shutil
import sys
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
PRODUCTS = ROOT / "data" / "cache" / "products"
TRADE = PRODUCTS / "trade-flows.json"
CURRENCIES = PRODUCTS / "currency-spacetime.json"
STORE = ROOT / "data" / "cache" / "economy-atlas" / "datasets"
LIMIT_PLACES, LIMIT_FLOWS = 5000, 2_000_000
MONTH = re.compile(r"^\d{4}-(0[1-9]|1[0-2])$")
SAME_POINT_DEG = 1e-6


class DatasetError(ValueError):
    """A dataset that cannot be read, with what to change."""


# ---- reading what is given -------------------------------------------------------------------------------


def _rows(text: str) -> list[dict]:
    """CSV rows as dicts with lower-case, trimmed column names (comma, semicolon or tab separated)."""
    text = text.lstrip("﻿")
    try:
        dialect = csv.Sniffer().sniff(text[:4096], delimiters=",;\t")
    except csv.Error:
        dialect = csv.excel
    reader = csv.DictReader(io.StringIO(text), dialect=dialect)
    if not reader.fieldnames:
        return []
    names = [(n or "").strip().lower() for n in reader.fieldnames]
    out = []
    for row in reader:
        out.append({names[i]: (v or "").strip() for i, v in enumerate(row.get(k) for k in reader.fieldnames)})
    return out


def _number(value, what: str) -> float:
    try:
        x = float(str(value).replace(" ", "").replace("_", ""))
    except ValueError:
        raise DatasetError(f"{what}: {value!r} is not a number") from None
    if not math.isfinite(x):
        raise DatasetError(f"{what}: {value!r} is not a finite number")
    return x


def _files(files: list[dict]) -> list[tuple[str, object]]:
    """Each given file as (name, parsed): a JSON value, or CSV rows."""
    out = []
    for f in files:
        name, text = str(f.get("name") or "data"), str(f.get("text") or "")
        if not text.strip():
            continue
        if name.lower().endswith(".json") or text.lstrip()[:1] in "{[":
            try:
                out.append((name, json.loads(text)))
            except ValueError as error:
                raise DatasetError(f"{name}: not JSON ({error})") from None
        else:
            out.append((name, _rows(text)))
    if not out:
        raise DatasetError("no file with anything in it")
    return out


# ---- the library's economies, for ids given without places ------------------------------------------------


def library_places() -> tuple[dict, dict]:
    """The library's economies by ISO alpha-3 code and by name (case-folded)."""
    try:
        economies = json.loads(TRADE.read_text(encoding="utf-8"))["economies"]
    except (OSError, ValueError, KeyError):
        return {}, {}
    by_code = {e["iso3"]: e for e in economies}
    by_name = {e["name"].casefold(): e for e in economies}
    return by_code, by_name


# ---- flows ---------------------------------------------------------------------------------------------------


def make_flows(title: str, files: list[dict]) -> dict:
    """A trade-flows/1 product from places and flows given as CSV or JSON."""
    places: dict[str, dict] = {}
    flows: list[tuple[str, str, float, int]] = []
    unit = ""
    for name, data in _files(files):
        if isinstance(data, dict):
            title = title or str(data.get("title") or "")
            unit = unit or str(data.get("unit") or "")
            for p in data.get("places") or []:
                _place(places, p, name)
            for f in data.get("flows") or []:
                flows.append(_flow(f, name))
        elif isinstance(data, list):
            if not data:
                continue
            columns = set(data[0])
            if {"from", "to", "value"} <= columns:
                flows += [_flow(r, name) for r in data]
            elif {"id", "latitude", "longitude"} <= columns:
                for r in data:
                    _place(places, r, name)
            else:
                raise DatasetError(f"{name}: columns from,to,value (flows) or id,name,latitude,longitude (places) are needed; "
                                   f"found {', '.join(sorted(columns)) or 'none'}")
    if not flows:
        raise DatasetError("no flows: a file with columns from,to,value is needed")
    # Ids without a place of their own are the library's economies, by code or by name, each spelling
    # of one economy (CHN, chn, China) becoming its ISO alpha-3 code.
    by_code, by_name = library_places()
    canon: dict[str, str] = {}
    unknown = []
    for ident in sorted({f[0] for f in flows} | {f[1] for f in flows}):
        if ident in places:
            canon[ident] = ident
            continue
        e = by_code.get(ident.upper()) or by_name.get(ident.casefold())
        if e is None:
            unknown.append(ident)
            continue
        canon[ident] = e["iso3"]
        places.setdefault(e["iso3"], {"id": e["iso3"], "name": e["name"], "latitude": e["latitude"], "longitude": e["longitude"],
                                      "group": e.get("continent", "")})
    if unknown:
        shown = ", ".join(unknown[:6]) + (" …" if len(unknown) > 6 else "")
        raise DatasetError(f"no place for {len(unknown)} id(s): {shown}. Give them in a places file "
                           "(id,name,latitude,longitude), or use ISO alpha-3 codes or the names of the library's economies")
    flows = [(canon[a], canon[b], v, y) for a, b, v, y in flows]
    ids = sorted(set(canon.values()))
    if len(ids) < 2:
        raise DatasetError("flows need at least two places")
    if len(ids) > LIMIT_PLACES or len(flows) > LIMIT_FLOWS:
        raise DatasetError(f"at most {LIMIT_PLACES} places and {LIMIT_FLOWS} flows")
    # Two places at one point would make the map's triangulation degenerate.
    seen: dict[tuple, str] = {}
    for i in ids:
        p = places[i]
        key = (round(p["latitude"] / SAME_POINT_DEG), round(p["longitude"] / SAME_POINT_DEG))
        if key in seen:
            raise DatasetError(f"{seen[key]} and {i} are at the same point ({p['latitude']}, {p['longitude']}): move one, or merge them")
        seen[key] = i
    index = {i: k for k, i in enumerate(ids)}
    by_year: dict[int, dict[tuple[int, int], float]] = {}
    for a, b, v, y in flows:
        if a == b or v == 0:
            continue
        cell = by_year.setdefault(y, {})
        cell[(index[a], index[b])] = cell.get((index[a], index[b]), 0.0) + v
    if not by_year:
        raise DatasetError("every flow is zero or from a place to itself")
    years = sorted(by_year)
    return {
        "schema": "trade-flows/1",
        "title": title or "Flows of one's own",
        "domain": "economics",
        "unit": unit or "units of the data",
        "years": years,
        "economies": [{"iso3": i, "name": places[i]["name"], "continent": places[i].get("group", ""),
                       "latitude": places[i]["latitude"], "longitude": places[i]["longitude"]} for i in ids],
        "flows": {str(y): [[a, b, round(v, 6)] for (a, b), v in sorted(by_year[y].items())] for y in years},
        "source": {"citation": f"{title or 'Flows of one own'}: data given to the Economy Atlas", "credit": "the user"},
    }


def _place(places: dict, p: dict, name: str) -> None:
    ident = str(p.get("id") or "").strip()
    if not ident:
        raise DatasetError(f"{name}: a place without an id")
    lat = _number(p.get("latitude"), f"{name}, {ident}'s latitude")
    lon = _number(p.get("longitude"), f"{name}, {ident}'s longitude")
    if not (-90 <= lat <= 90 and -180 <= lon <= 360):
        raise DatasetError(f"{name}: {ident} is at ({lat}, {lon}), not on the sphere")
    if lon > 180:
        lon -= 360
    places[ident] = {"id": ident, "name": str(p.get("name") or ident), "latitude": lat, "longitude": lon,
                     "group": str(p.get("group") or "")}


def _flow(f: dict, name: str) -> tuple[str, str, float, int]:
    a, b = str(f.get("from") or "").strip(), str(f.get("to") or "").strip()
    if not a or not b:
        raise DatasetError(f"{name}: a flow without from or to")
    v = _number(f.get("value"), f"{name}, the flow {a} → {b}")
    if v < 0:
        raise DatasetError(f"{name}: the flow {a} → {b} is negative; give it the other way")
    y = f.get("year")
    year = int(_number(y, f"{name}, the year of {a} → {b}")) if str(y if y is not None else "").strip() else 0
    return a, b, v, year


# ---- currencies ----------------------------------------------------------------------------------------------


def make_currencies(title: str, files: list[dict], base: str = "") -> dict:
    """A currency-spacetime/1 product from monthly rates against a base and policy rates."""
    rows: list[dict] = []
    for name, data in _files(files):
        if isinstance(data, dict):
            title = title or str(data.get("title") or "")
            base = base or str(data.get("base") or "")
            rows += [{**r, "_file": name} for r in data.get("rows") or []]
        elif isinstance(data, list):
            if data and not {"month", "currency", "policy"} <= set(data[0]):
                raise DatasetError(f"{name}: columns month,currency,rate,policy are needed")
            rows += [{**r, "_file": name} for r in data]
    if not rows:
        raise DatasetError("no rows: columns month,currency,rate,policy are needed")
    rates: dict[str, dict[str, float]] = {}
    policy: dict[str, dict[str, float]] = {}
    ones: dict[str, int] = {}
    for r in rows:
        month, code = str(r.get("month") or "").strip(), str(r.get("currency") or "").strip().upper()
        where = f"{r['_file']}, {code} in {month}"
        if not MONTH.match(month):
            raise DatasetError(f"{r['_file']}: {month!r} is not a month (YYYY-MM)")
        if not code:
            raise DatasetError(f"{r['_file']}: a row without a currency")
        if str(r.get("policy") or "").strip():
            policy.setdefault(code, {})[month] = _number(r["policy"], f"{where}, the policy rate")
        rate = str(r.get("rate") or "").strip()
        if rate:
            x = _number(rate, f"{where}, the rate")
            if x <= 0:
                raise DatasetError(f"{where}: a rate is positive")
            rates.setdefault(code, {})[month] = x
            if x == 1:
                ones[code] = ones.get(code, 0) + 1
    codes = sorted(set(rates) | set(policy))
    if not base:
        # The base is the currency whose rate is always 1, or that has a policy rate and no rate.
        candidates = [c for c in codes if c not in rates or ones.get(c, 0) == len(rates[c])]
        if len(candidates) != 1:
            raise DatasetError("say which currency is the base (a field 'base'), or give its rows with rate 1 or empty")
        base = candidates[0]
    base = base.upper()
    if base not in policy:
        raise DatasetError(f"the base {base} has no policy rate")
    if len(codes) < 2:
        raise DatasetError("a record needs the base and at least one other currency")
    months = sorted({m for series in list(rates.values()) + list(policy.values()) for m in series})
    first, last = months[0], months[-1]
    months = _months_between(first, last)
    ordered = [base] + [c for c in codes if c != base]
    return {
        "schema": "currency-spacetime/1",
        "title": title or "Currencies of one's own",
        "domain": "economics",
        "base": base,
        "months": months,
        "time": "a month is its middle, in decimal years: year + (month − 0.5)/12",
        "currencies": [{"code": c, "name": c, "area": ""} for c in ordered],
        "rates": {c: [rates.get(c, {}).get(m) for m in months] for c in ordered if c != base},
        "rate_unit": f"units of the currency per {base}",
        "policy": {c: [policy.get(c, {}).get(m) for m in months] for c in ordered},
        "policy_unit": "per cent a year",
        "source": {"citation": f"{title or 'Currencies of one own'}: data given to the Economy Atlas", "credit": "the user"},
    }


def _months_between(first: str, last: str) -> list[str]:
    y, m = int(first[:4]), int(first[5:7])
    out = []
    while f"{y:04d}-{m:02d}" <= last:
        out.append(f"{y:04d}-{m:02d}")
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return out


# ---- the store -----------------------------------------------------------------------------------------------


class Datasets:
    """The datasets kept here, and the library's two products beside them."""

    def __init__(self, store: Path = STORE):
        self.store = Path(store)

    def list(self) -> list[dict]:
        out = []
        if TRADE.exists():
            out.append({"id": "library:trade", "kind": "flows", "title": "World trade in goods (OECD BIMTS)", "library": True,
                        **_summary(TRADE)})
        if CURRENCIES.exists():
            out.append({"id": "library:currencies", "kind": "currencies", "title": "Twenty currencies against the euro (ECB, BIS)",
                        "library": True, **_summary(CURRENCIES)})
        if self.store.exists():
            for about in sorted(self.store.glob("*/about.json")):
                try:
                    out.append(json.loads(about.read_text(encoding="utf-8")))
                except (OSError, ValueError):
                    continue
        return out

    def product(self, ident: str | None, kind: str) -> Path:
        """The product file of a dataset (the library's when none is named)."""
        if not ident or ident == f"library:{'trade' if kind == 'flows' else 'currencies'}":
            path = TRADE if kind == "flows" else CURRENCIES
            if not path.exists():
                raise DatasetError("the library's product is not built: run the recipe")
            return path
        if not re.fullmatch(r"[0-9a-f]{16}", ident):
            raise DatasetError(f"no dataset {ident!r}")
        about = self.store / ident / "about.json"
        if not about.exists():
            raise DatasetError(f"no dataset {ident!r}")
        if json.loads(about.read_text(encoding="utf-8")).get("kind") != kind:
            raise DatasetError(f"the dataset {ident} is not of {kind}")
        return self.store / ident / "product.json"

    def add(self, kind: str, title: str, files: list[dict], base: str = "") -> dict:
        if kind == "flows":
            product = make_flows(title, files)
        elif kind == "currencies":
            product = make_currencies(title, files, base)
        else:
            raise DatasetError("a dataset is of flows or of currencies")
        text = json.dumps(product, ensure_ascii=False, separators=(",", ":"))
        ident = hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]
        folder = self.store / ident
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "product.json").write_text(text + "\n", encoding="utf-8")
        about = {"id": ident, "kind": kind, "title": product["title"], "library": False,
                 "added": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), **_summary(folder / "product.json")}
        (folder / "about.json").write_text(json.dumps(about, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
        return about

    def remove(self, ident: str) -> bool:
        if not re.fullmatch(r"[0-9a-f]{16}", ident) or not (self.store / ident).exists():
            return False
        shutil.rmtree(self.store / ident)
        return True


def _summary(path: Path) -> dict:
    """What the page needs to offer a dataset: its years and places, or its months and currencies."""
    p = json.loads(path.read_text(encoding="utf-8"))
    if p.get("schema", "").startswith("trade-flows"):
        return {"years": p["years"], "places": len(p["economies"]), "unit": p.get("unit", ""),
                "flows": sum(len(v) for v in p["flows"].values())}
    return {"months": [p["months"][0], p["months"][-1]], "base": p.get("base", "EUR"),
            "currencies": [c["code"] for c in p["currencies"]]}


def main(argv: list[str]) -> int:
    store = Datasets()
    if argv[:1] == ["list"]:
        for d in store.list():
            print(json.dumps(d, ensure_ascii=False))
        return 0
    if len(argv) >= 4 and argv[0] == "add":
        kind, title, paths = argv[1], argv[2], argv[3:]
        files = [{"name": Path(p).name, "text": Path(p).read_text(encoding="utf-8")} for p in paths]
        try:
            print(json.dumps(store.add(kind, title, files), ensure_ascii=False))
        except DatasetError as error:
            print(f"not added: {error}", file=sys.stderr)
            return 1
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
