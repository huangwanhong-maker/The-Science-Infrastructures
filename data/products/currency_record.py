"""Build the currency-spacetime product: currencies held through time and exchanged, month by month.

Each currency is an agent of the economics domain's spacetime complex: it is held from month to
month at its central bank's policy rate (BIS, end of month, per cent a year) and exchanged with
the euro at the month's average reference rate (ECB, units of the currency per euro). A currency
without a rate or a policy rate in some month keeps null there; the application chooses the
months and currencies with complete records.

The euro's rate is the ECB's deposit facility rate in force at the end of the month (ECB, from its
dates of change), not the BIS's series for the euro area. Since 2015 the Eurosystem has supplied
liquidity in excess, so that euro money has earned the deposit facility rate, the floor of the
market rates; the BIS's series follows the main refinancing rate, 25 to 50 basis points above it,
until August 2024, and the deposit facility rate from September 2024, when the ECB's new
operational framework took effect. Held at the refinancing rate, the euro would earn about half a
per cent a year that no holder of euros earned, which appears, for instance, as a spurious excess
return of the euro over the Danish krone, pegged to it.

The product keeps what the sources say, in their units; the conversion to transports (the log of
the rate, and log(1 + r/100) per year held) is the domain's and the application's.

Output, in data/cache/products/: currency-spacetime.json

Usage::

    python data/products/currency_record.py
"""

from __future__ import annotations

import calendar
import csv
import json
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from product_index import OUT, register  # noqa: E402

PROFILE = "major"
# ISO 4217 currency -> BIS reference area of its central bank (ISO 3166-1 alpha-2; XM the euro area).
AREA = {"EUR": "XM", "USD": "US", "JPY": "JP", "GBP": "GB", "CHF": "CH", "SEK": "SE", "NOK": "NO", "DKK": "DK",
        "PLN": "PL", "CZK": "CZ", "HUF": "HU", "AUD": "AU", "CAD": "CA", "NZD": "NZ", "CNY": "CN", "INR": "IN",
        "KRW": "KR", "BRL": "BR", "MXN": "MX", "ZAR": "ZA", "TRY": "TR"}
NAME = {"EUR": "euro", "USD": "US dollar", "JPY": "Japanese yen", "GBP": "pound sterling", "CHF": "Swiss franc",
        "SEK": "Swedish krona", "NOK": "Norwegian krone", "DKK": "Danish krone", "PLN": "Polish zloty",
        "CZK": "Czech koruna", "HUF": "Hungarian forint", "AUD": "Australian dollar", "CAD": "Canadian dollar",
        "NZD": "New Zealand dollar", "CNY": "Chinese yuan renminbi", "INR": "Indian rupee", "KRW": "South Korean won",
        "BRL": "Brazilian real", "MXN": "Mexican peso", "ZAR": "South African rand", "TRY": "Turkish lira"}


def read_rates(rows) -> dict[str, dict[str, float]]:
    """The ECB's monthly averages: currency -> month -> units per euro."""
    out: dict[str, dict[str, float]] = {}
    for r in rows:
        if r.get("CURRENCY_DENOM") != "EUR" or r.get("EXR_SUFFIX") != "A" or r.get("OBS_VALUE") in (None, ""):
            continue
        out.setdefault(r["CURRENCY"], {})[r["TIME_PERIOD"]] = float(r["OBS_VALUE"])
    return out


def read_policy(rows) -> dict[str, dict[str, float]]:
    """The BIS's policy rates: area -> month -> per cent a year."""
    out: dict[str, dict[str, float]] = {}
    for r in rows:
        if r.get("FREQ") != "M" or r.get("OBS_VALUE") in (None, "", "NaN"):
            continue
        out.setdefault(r["REF_AREA"], {})[r["TIME_PERIOD"]] = float(r["OBS_VALUE"])
    return out


def read_deposit_rate(rows) -> list[tuple[str, float]]:
    """The ECB's deposit facility rate: (date from which it applies, per cent a year), in order."""
    out = [(r["TIME_PERIOD"], float(r["OBS_VALUE"])) for r in rows if r.get("OBS_VALUE") not in (None, "", "NaN")]
    return sorted(out)


def in_force(changes: list[tuple[str, float]], month: str) -> float | None:
    """The rate in force at the end of a month (YYYY-MM), from its dates of change."""
    y, m = int(month[:4]), int(month[5:7])
    end = f"{month}-{calendar.monthrange(y, m)[1]:02d}"
    rate = None
    for date, value in changes:
        if date > end:
            break
        rate = value
    return rate


def months_between(first: str, last: str) -> list[str]:
    y, m = int(first[:4]), int(first[5:7])
    out = []
    while f"{y:04d}-{m:02d}" <= last:
        out.append(f"{y:04d}-{m:02d}")
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return out


def build() -> dict:
    ecb = DATA / "cache" / "ecb-exr" / PROFILE
    bis = DATA / "cache" / "bis-policy-rates" / PROFILE
    with (ecb / "exr.csv").open(encoding="utf-8", newline="") as fh:
        rates = read_rates(csv.DictReader(fh))
    with (bis / "policy.csv").open(encoding="utf-8", newline="") as fh:
        policy = read_policy(csv.DictReader(fh))
    with (ecb / "dfr.csv").open(encoding="utf-8", newline="") as fh:
        deposit = read_deposit_rate(csv.DictReader(fh))
    observed = sorted({m for series in rates.values() for m in series})
    months = months_between(observed[0], observed[-1])
    currencies = ["EUR"] + sorted(c for c in rates if c in AREA)
    product = {
        "schema": "currency-spacetime/1",
        "title": "Currencies held at their policy rates and exchanged with the euro, monthly",
        "domain": "economics",
        "base": "EUR",
        "months": months,
        "time": "a month is its middle, in decimal years: year + (month − 0.5)/12",
        "currencies": [{"code": c, "name": NAME.get(c, c), "area": AREA[c]} for c in currencies],
        "rates": {c: [rates[c].get(m) for m in months] for c in currencies if c != "EUR"},
        "rate_unit": "units of the currency per euro, monthly average of the ECB's daily reference rates",
        "policy": {c: [in_force(deposit, m) if c == "EUR" else policy.get(AREA[c], {}).get(m) for m in months] for c in currencies},
        "policy_unit": "per cent a year, the central bank's policy rate at the end of the month; the euro's, the ECB's deposit facility rate",
        "source": {
            "rates": {"dataset": "ecb-exr", "profile": PROFILE, "credit": "ECB",
                      "retrieved_at": json.loads((ecb / "retrieval.json").read_text(encoding="utf-8"))["retrieved_at"]},
            "euro": {"dataset": "ecb-exr", "series": "FM.B.U2.EUR.4F.KR.DFR.LEV", "credit": "ECB"},
            "policy": {"dataset": "bis-policy-rates", "profile": PROFILE, "credit": "BIS",
                       "retrieved_at": json.loads((bis / "retrieval.json").read_text(encoding="utf-8"))["retrieved_at"]},
            "licence": "ECB: reuse with the source acknowledged; BIS: use unrestricted with the BIS cited as the source",
        },
    }
    path = OUT / "currency-spacetime.json"
    OUT.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(product, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    complete = {c: sum(1 for i, m in enumerate(months) if product["policy"][c][i] is not None and (c == "EUR" or product["rates"][c][i] is not None))
                for c in currencies}
    register("currencies", {"name": "currency-spacetime", "file": path.name, "title": product["title"], "domain": "economics",
                            "months": [months[0], months[-1]], "currencies": currencies, "complete_months": complete,
                            "source": {"rates": "ECB", "policy": "BIS", "euro": "ECB deposit facility rate"}})
    return product


def main(argv: list[str]) -> int:
    for folder, file in (("ecb-exr", "dfr.csv"), ("bis-policy-rates", "policy.csv")):
        if not (DATA / "cache" / folder / PROFILE / file).exists():
            print(f"nothing fetched; run: python data/acquisition/fetch.py {folder}", file=sys.stderr)
            return 1
    p = build()
    gaps = [c["code"] for c in p["currencies"] if any(v is None for v in p["policy"][c["code"]])]
    print(f"currency-spacetime: {len(p['currencies'])} currencies, {p['months'][0]} to {p['months'][-1]}; "
          f"policy rates missing in some months for: {', '.join(gaps) or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
