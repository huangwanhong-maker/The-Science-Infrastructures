"""Tables for the Economy Atlas's toolbox (ADR 0048): the data library's products as tables of
columns, and tables of one's own.

A table is named columns of equal length, each of numbers (null where missing) or of text:

    {"id", "title", "source", "rows", "columns": [{"name", "kind": "number"|"text", "label"?}],
     "data": {name: [values]}}

The library's tables, made from the products the atlas already fetches, and kept under
data/cache/economy-atlas/tables/ until their sources change:

    economies     each economy in each year, 2000–2024: its exports, imports and net inflow of goods
                  (OECD BIMTS, 2010–2024) beside the World Bank's indicators fetched (GDP per capita,
                  population, life expectancy, fertility, urban share, forest share, access to
                  electricity, CO2 and greenhouse gases per person, the Gini index, density), with the
                  World Bank's region and income group; a panel for the econometrics
    currencies    each currency in each month: its rate against the euro, its log, and its policy rate
    currency-rates  the same, wide: one column of rates and one of policy rates for each currency

Tables of one's own are CSV files (comma, semicolon or tab separated, a header row), kept as the
datasets of kind "table" (datasets.py); a column is numbers when every filled cell is a number.

    python applications/economy-atlas/tables.py list
    python applications/economy-atlas/tables.py make     (the recipe's step: fetch, then make)
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
CACHE = ROOT / "data" / "cache"
PRODUCTS = CACHE / "products"
KEPT = CACHE / "economy-atlas" / "tables"
WDI = CACHE / "worldbank-indicators"
sys.path.insert(0, str(APP))
import datasets  # noqa: E402
from naming import POLICY  # noqa: E402  (the names of places, ADR 0069; datasets put data/products on the path)

# The World Bank's profiles the economies table joins, as columns (data/registry/worldbank-indicators.json).
INDICATORS = [
    ("gdp-per-capita", "gdp_per_capita", "GDP per capita, current US$"),
    ("population", "population", "population"),
    ("population-density", "population_density", "people per km²"),
    ("life-expectancy", "life_expectancy", "life expectancy at birth, years"),
    ("fertility", "fertility", "births per woman"),
    ("urban-share", "urban_share", "urban population, % of total"),
    ("forest-share", "forest_share", "forest area, % of land"),
    ("electricity-access", "electricity_access", "access to electricity, % of population"),
    ("co2-per-capita", "co2_per_capita", "CO2 emissions, t per person"),
    ("ghg-per-capita", "ghg_per_capita", "greenhouse gases, t CO2e per person"),
    ("gini", "gini", "Gini index (World Bank estimate)"),
    ("gdp-growth", "gdp_growth", "GDP growth, annual %, constant prices"),
]


def _stamp(paths: list[Path]) -> str:
    return ";".join(f"{p.name}:{p.stat().st_mtime_ns}" for p in paths if p.exists())


def _kept(ident: str, sources: list[Path], make) -> dict | None:
    """A library table, made again only when its sources change."""
    if not any(p.exists() for p in sources):
        return None
    stamp = _stamp(sources)
    path = KEPT / f"{ident}.json"
    if path.exists():
        try:
            t = json.loads(path.read_text(encoding="utf-8"))
            if t.get("stamp") == stamp:
                return t
        except (OSError, ValueError):
            pass
    t = make()
    t["stamp"] = stamp
    KEPT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(t, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    return t


def _table(ident: str, title: str, source: str, columns: list[tuple[str, str, str]], rows: list[dict]) -> dict:
    return {"id": ident, "title": title, "source": source, "rows": len(rows),
            "columns": [{"name": n, "kind": k, **({"label": lab} if lab else {})} for n, k, lab in columns],
            "data": {n: [r.get(n) for r in rows] for n, _, _ in columns}}


def _economies() -> dict:
    trade = json.loads((PRODUCTS / "trade-flows.json").read_text(encoding="utf-8")) if (PRODUCTS / "trade-flows.json").exists() else None
    cells: dict[tuple[str, int], dict] = {}
    names: dict[str, str] = {}
    region: dict[str, str] = {}
    income: dict[str, str] = {}
    if trade:
        ids = [e["iso3"] for e in trade["economies"]]
        for e in trade["economies"]:
            names[e["iso3"]] = e["name"]
        for year, flows in trade["flows"].items():
            y = int(year)
            for i, j, v in flows:
                a = cells.setdefault((ids[int(i)], y), {})
                b = cells.setdefault((ids[int(j)], y), {})
                a["exports"] = a.get("exports", 0.0) + v
                b["imports"] = b.get("imports", 0.0) + v
    present = []
    for profile, column, _ in INDICATORS:
        folder = WDI / profile
        if not (folder / "values.json").exists():
            continue
        present.append(column)
        try:
            meta = json.loads((folder / "economies.json").read_text(encoding="utf-8"))[1]
            for e in meta:
                if e.get("region", {}).get("value") != "Aggregates":
                    region.setdefault(e["id"], e["region"]["value"].strip())
                    income.setdefault(e["id"], e["incomeLevel"]["value"].strip())
                    names.setdefault(e["id"], e["name"])
        except (OSError, ValueError, IndexError, KeyError):
            pass
        for r in json.loads((folder / "values.json").read_text(encoding="utf-8"))[1] or []:
            code, value = r.get("countryiso3code"), r.get("value")
            if not code or value is None or code not in region:
                continue
            cells.setdefault((code, int(r["date"])), {})[column] = float(value)
    rows = []
    for (code, year), c in sorted(cells.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        row = {"iso3": code, "name": POLICY.name(code, names.get(code, code)) or code, "region": region.get(code, ""), "income": income.get(code, ""), "year": year, **c}
        if "exports" in c or "imports" in c:
            row["exports"], row["imports"] = c.get("exports", 0.0), c.get("imports", 0.0)
            row["balance"] = row["imports"] - row["exports"]
            row["trade"] = row["imports"] + row["exports"]
        rows.append(row)
    columns = [("iso3", "text", "ISO alpha-3"), ("name", "text", ""), ("region", "text", "World Bank region"), ("income", "text", "World Bank income group"),
               ("year", "number", "")]
    if trade:
        columns += [("exports", "number", "goods exported, millions of US$ (OECD BIMTS)"), ("imports", "number", "goods imported, millions of US$"),
                    ("balance", "number", "net inflow: imports − exports, millions of US$"), ("trade", "number", "exports + imports, millions of US$")]
    columns += [(col, "number", label) for _, col, label in INDICATORS if col in present]
    return _table("economies", "Economies by year: trade and the World Bank's indicators",
                  "OECD BIMTS (CC BY 4.0), World Bank World Development Indicators (CC BY 4.0)", columns, rows)


def _currencies(wide: bool) -> dict:
    p = json.loads((PRODUCTS / "currency-spacetime.json").read_text(encoding="utf-8"))
    months, base = p["months"], p.get("base", "EUR")
    time = [int(m[:4]) + (int(m[5:7]) - 0.5) / 12 for m in months]
    codes = [c["code"] for c in p["currencies"]]
    if wide:
        columns = [("month", "text", ""), ("time", "number", "decimal year")]
        data_rows = [{"month": m, "time": t} for m, t in zip(months, time)]
        for c in codes:
            if c != base:
                columns.append((c, "number", f"units of {c} per {base}"))
                for r, v in zip(data_rows, p["rates"][c]):
                    r[c] = v
            columns.append((f"{c}_policy", "number", f"{c}'s policy rate, % a year"))
            for r, v in zip(data_rows, p["policy"][c]):
                r[f"{c}_policy"] = v
        return _table("currency-rates", f"Currencies by month, wide: rates against the {base} and policy rates",
                      "ECB (rates, the euro's deposit facility rate), BIS (policy rates)", columns, data_rows)
    rows = []
    for c in codes:
        for k, m in enumerate(months):
            rate = 1.0 if c == base else p["rates"][c][k]
            rows.append({"currency": c, "month": m, "time": time[k], "rate": rate, "log_rate": math.log(rate) if rate else None,
                         "policy": p["policy"][c][k]})
    return _table("currencies", f"Currencies by month: rates against the {base} and policy rates",
                  "ECB (rates, the euro's deposit facility rate), BIS (policy rates)",
                  [("currency", "text", ""), ("month", "text", ""), ("time", "number", "decimal year"), ("rate", "number", f"units per {base}"),
                   ("log_rate", "number", ""), ("policy", "number", "policy rate, % a year")], rows)


def library() -> list[dict]:
    out = []
    t = _kept("economies", [PRODUCTS / "trade-flows.json", *[WDI / p / "values.json" for p, _, _ in INDICATORS]], _economies)
    if t:
        out.append(t)
    for ident, wide in (("currencies", False), ("currency-rates", True)):
        t = _kept(ident, [PRODUCTS / "currency-spacetime.json"], lambda w=wide: _currencies(w))
        if t:
            out.append(t)
    return out


def own() -> list[dict]:
    """The tables of one's own, kept as datasets of kind "table"."""
    out = []
    store = datasets.Datasets()
    for d in store.list():
        if d.get("kind") == "table":
            try:
                t = json.loads(store.product(d["id"], "table").read_text(encoding="utf-8"))
                out.append({**t, "id": d["id"], "own": True})
            except (OSError, ValueError, datasets.DatasetError):
                continue
    return out


def listing() -> list[dict]:
    """Every table, without its data."""
    return [{k: v for k, v in t.items() if k not in ("data", "stamp")} for t in library() + own()]


def table(ident: str) -> dict:
    for t in library() + own():
        if t["id"] == ident:
            return {k: v for k, v in t.items() if k != "stamp"}
    raise datasets.DatasetError(f"no table {ident!r}")


def make() -> int:
    """The recipe's step: fetch the World Bank's profiles not yet in the cache, then make the tables."""
    import subprocess
    fetch = ROOT / "data" / "acquisition" / "fetch.py"
    failed = 0
    for profile, _, _ in INDICATORS:
        if (WDI / profile / "values.json").exists():
            continue
        print(f"fetching the World Bank's {profile}", flush=True)
        if subprocess.run([sys.executable, str(fetch), "worldbank-indicators", "--profile", profile], cwd=ROOT).returncode:
            failed += 1
    for t in library():
        print(t["id"], t["rows"], "rows", len(t["columns"]), "columns", flush=True)
    return 1 if failed == len(INDICATORS) else 0


if __name__ == "__main__":
    if sys.argv[1:2] == ["list"]:
        for t in listing():
            print(t["id"], t["rows"], [c["name"] for c in t["columns"]])
    elif sys.argv[1:2] == ["make"]:
        sys.exit(make())
    else:
        print(__doc__)
