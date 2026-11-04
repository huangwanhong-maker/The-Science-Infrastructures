"""The counters' weeks around their mean week, and their directions: what the field theory's Euclidean measure is
calibrated to, and what its current is tested against.

The shared problem (../2026-11-05-a-district-economy-the-problem) calibrates each counter's mean week, the week's
constituents of its tide (results/sensors.json). This script reads the same counts (the data library's street-counts
product, hourly, UTC) over the same reference year, the 52 weeks from Monday 27 October 2025 00:00 civil, and writes
for the study's field theory:

1. **The fluctuations.** At each counter c and hour t of the year, the count y_c(t) against its mean week
   m_c(h(t)) at the civil hour of the week h(t). The district's common relative deviation
   eps(t) = sum_c (y_c - m_c) / sum_c m_c over the counters of the district's box that report at t, its autocovariance
   at lags 0..336 h; and the local residual chi_c(t) = y_c(t) - (1 + eps(t)) m_c(h(t)), its variance by counter, its
   lag-0 correlation between counters, and its pooled autocorrelation. Public holidays are kept: a week that is not
   typical is what the fluctuations describe.
2. **The directions.** The counts by direction (the product's d1 and d2, with their compass names), averaged by hour
   of the civil week over the reference year: the net current past each counter of the district's box.

Usage (from the repository's root):

    .venv/Scripts/python.exe research/experiments/2026-11-05-a-district-economy-as-a-field-theory/counters.py

It writes results/counters.json, derived from the City of Melbourne's counts (the City's open data, CC BY 4.0).
"""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROBLEM = ROOT / "research" / "experiments" / "2026-11-05-a-district-economy-the-problem" / "results"
WEEK = 168
LAGS = 337


def utc(s: str) -> dt.datetime:
    return dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=dt.timezone.utc)


def civil_hours(t0: dt.datetime, n: int, offsets: list) -> np.ndarray:
    """The civil hour of the week (Monday 00:00 first) of each of n hours from t0 (UTC), by the record's offsets."""
    marks = [(utc(a), float(b)) for a, b in offsets]
    out = np.empty(n, dtype=np.int64)
    k = 0
    for i in range(n):
        t = t0 + dt.timedelta(hours=i)
        while k + 1 < len(marks) and marks[k + 1][0] <= t:
            k += 1
        c = t + dt.timedelta(hours=marks[k][1])
        out[i] = c.weekday() * 24 + c.hour
    return out


def main() -> None:
    inputs = json.loads((PROBLEM / "problem-inputs.json").read_text(encoding="utf-8"))
    sensors = json.loads((PROBLEM / "sensors.json").read_text(encoding="utf-8"))
    box = inputs["district"]["district_box"]
    start = utc(inputs["reference"]["start_utc"])
    hours = int(inputs["reference"]["weeks"]) * WEEK
    where = {s["id"]: s for s in inputs["sensors"]}
    rows = []
    for s in sensors["sensors"]:
        meta = where[s["id"]]
        inside = box["x0"] <= s["x"] <= box["x1"] and box["y0"] <= s["y"] <= box["y1"]
        series = json.loads((ROOT / meta["series"]).read_text(encoding="utf-8"))
        t0 = utc(series["t0_utc"])
        i0 = int((start - t0).total_seconds() // 3600)
        y = np.full(hours, np.nan)
        counts = series["counts"]
        for i in range(hours):
            j = i0 + i
            if 0 <= j < len(counts) and counts[j] is not None:
                y[i] = counts[j]
        h = civil_hours(start, hours, series["offsets"])
        mw = np.asarray(s["mean_week"], dtype=float)
        row = {"id": s["id"], "name": s["name"], "inside": inside, "y": y, "h": h, "m": mw[h]}
        d = series.get("directions")
        if inside and d and meta.get("directions"):
            dt0 = utc(d["t0_utc"])
            k0 = int((start - dt0).total_seconds() // 3600)
            dirs = {}
            for key in ("d1", "d2"):
                v = np.full(hours, np.nan)
                arr = d[key]
                for i in range(hours):
                    j = k0 + i
                    if 0 <= j < len(arr) and arr[j] is not None:
                        v[i] = arr[j]
                dirs[key] = v
            ok = np.isfinite(dirs["d1"]) & np.isfinite(dirs["d2"])
            weeks = {}
            for key in ("d1", "d2"):
                acc = np.zeros(WEEK)
                cnt = np.zeros(WEEK)
                np.add.at(acc, h[ok], dirs[key][ok])
                np.add.at(cnt, h[ok], 1)
                weeks[key] = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
            row["directions"] = {"names": meta["directions"], "hours": int(ok.sum()),
                                 "d1": [round(float(x), 2) for x in weeks["d1"]], "d2": [round(float(x), 2) for x in weeks["d2"]]}
        rows.append(row)

    inside = [r for r in rows if r["inside"]]
    Y = np.array([r["y"] for r in inside])
    M = np.array([r["m"] for r in inside])
    have = np.isfinite(Y)
    # The district's common relative deviation, over the counters that report at each hour.
    num = np.where(have, Y - M, 0.0).sum(0)
    den = np.where(have, M, 0.0).sum(0)
    eps = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)
    good = np.isfinite(eps)
    e = np.where(good, eps, 0.0)
    acov = []
    for lag in range(LAGS):
        a, b = e[: hours - lag], e[lag:]
        g = good[: hours - lag] & good[lag:]
        acov.append(float(np.mean(a[g] * b[g])) if g.any() else float("nan"))
    # The covariance of the common deviation between the hours of one week, over the year's 52 weeks (pairwise
    # complete): what a week's value is uncertain by when every counter of the district moves together.
    W = eps.reshape(-1, WEEK)
    weekly = np.full((WEEK, WEEK), np.nan)
    for a in range(WEEK):
        for b in range(a, WEEK):
            g = np.isfinite(W[:, a]) & np.isfinite(W[:, b])
            if g.sum() >= 10:
                weekly[a, b] = weekly[b, a] = float(np.mean(W[g, a] * W[g, b]))
    chi = Y - (1 + eps[None, :]) * M
    okc = have & good[None, :]
    var_chi = [float(np.nanmean(np.where(okc[c], chi[c] ** 2, np.nan))) for c in range(len(inside))]
    var_res = [float(np.nanmean(np.where(have[c], (Y[c] - M[c]) ** 2, np.nan))) for c in range(len(inside))]
    var_mean = [float(np.nanmean(np.where(have[c], (M[c] - np.nanmean(M[c])) ** 2, np.nan))) for c in range(len(inside))]
    n = len(inside)
    corr = np.full((n, n), np.nan)
    for a in range(n):
        for b in range(n):
            g = okc[a] & okc[b]
            if g.sum() > 100:
                x, z = chi[a][g], chi[b][g]
                corr[a, b] = float(np.mean(x * z) / np.sqrt(np.mean(x * x) * np.mean(z * z)))
    # The local residual's pooled autocorrelation, lags 0..167 h.
    auto = []
    for lag in range(WEEK):
        s, w = 0.0, 0.0
        for c in range(n):
            g = okc[c][: hours - lag] & okc[c][lag:]
            if g.sum() > 100:
                x, z = chi[c][: hours - lag][g], chi[c][lag:][g]
                s += np.mean(x * z) / np.mean(chi[c][okc[c]] ** 2)
                w += 1
        auto.append(s / w if w else float("nan"))
    out = {
        "reference": inputs["reference"],
        "credit": "City of Melbourne, Pedestrian Counting System (the City's open data, CC BY 4.0), the data library's street-counts product",
        "district_counters": [r["id"] for r in inside],
        "common": {"definition": "eps(t) = sum_c (y_c - m_c) / sum_c m_c over the district's counters reporting at t",
                   "hours": int(good.sum()), "mean": float(np.mean(e[good])), "sd": float(np.std(e[good])),
                   "autocovariance": [round(x, 7) for x in acov],
                   "series_by_hour_of_week_sd": [round(float(np.std(e[good & (inside[0]["h"] == k)])), 5) for k in range(WEEK)],
                   "weekly_covariance": [[round(float(x), 6) for x in row] for row in weekly],
                   "weeks": int(W.shape[0])},
        "local": {"definition": "chi_c(t) = y_c(t) - (1 + eps(t)) m_c(h(t))",
                  "variance": [round(v, 2) for v in var_chi], "residual_variance": [round(v, 2) for v in var_res],
                  "mean_week_variance": [round(v, 2) for v in var_mean],
                  "correlation": [[None if not np.isfinite(x) else round(float(x), 4) for x in row] for row in corr],
                  "autocorrelation": [round(float(x), 5) for x in auto]},
        "directions": {str(r["id"]): r["directions"] for r in inside if "directions" in r},
        "coverage": {str(r["id"]): round(float(np.isfinite(r["y"]).mean()), 4) for r in rows},
    }
    path = HERE / "results" / "counters.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"{n} counters of the district's box; common deviation sd {out['common']['sd']:.3f}; "
          f"local share of the residual variance {np.nanmedian(np.array(var_chi) / np.array(var_res)):.3f} (median); "
          f"{len(out['directions'])} with directions; wrote {path}")


if __name__ == "__main__":
    main()
