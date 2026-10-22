"""A minimal client for ERDDAP's griddap service, standard library only.

ERDDAP (Simons 2022, NOAA) serves gridded data through URLs of the form::

    <service>/griddap/<dataset>.nc?var1[axis constraints],var2[...]

where each axis is constrained as ``[(start):stride:(stop)]``, with values in
parentheses and the stride counted in indices. This module reads a dataset's
axes, chooses the indices that a retrieval asks for, and builds the request.
Downloads go through the caller, so that provenance is written in one place.

Reference: https://coastwatch.pfeg.noaa.gov/erddap/griddap/documentation.html
"""

from __future__ import annotations

import csv
import datetime as dt
import io
import json
import urllib.error
import urllib.parse
import urllib.request

AGENT = {"User-Agent": "laboratory-data-fetch/1 (+https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory)"}


class ErddapError(RuntimeError):
    pass


def _get(url: str, timeout: int = 120) -> bytes:
    request = urllib.request.Request(url, headers=AGENT)
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.read()
    except urllib.error.HTTPError as error:
        detail = error.read().decode("utf-8", "replace")[:400]
        raise ErddapError(f"HTTP {error.code} for {url}\n{detail}") from None


def info(service: str, dataset: str) -> dict:
    """Axes and variables of a dataset, with their attributes."""
    table = json.loads(_get(f"{service}/info/{dataset}/index.json"))["table"]
    out = {"axes": {}, "variables": {}, "global": {}}
    for kind, name, attribute, _, value in table["rows"]:
        if kind == "dimension":
            out["axes"][name] = {"order": len(out["axes"]), "summary": value}
        elif kind == "variable":
            out["variables"][name] = {"axes": value}
        elif kind == "attribute":
            target = (out["global"] if name == "NC_GLOBAL" else
                      out["axes"].get(name) or out["variables"].get(name))
            if target is not None:
                target[attribute] = value
    return out


def parse_time(text: str) -> dt.datetime:
    return dt.datetime.fromisoformat(text.replace("Z", "+00:00")).astimezone(dt.timezone.utc)


def iso(t: dt.datetime) -> str:
    return t.astimezone(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def time_values(service: str, dataset: str, start: dt.datetime, stop: dt.datetime) -> list[dt.datetime]:
    """The values of the time axis between start and stop, inclusive."""
    query = f"time[({iso(start)}):1:({iso(stop)})]"
    text = _get(f"{service}/griddap/{dataset}.csv?{urllib.parse.quote(query, safe='():,')}").decode()
    rows = list(csv.reader(io.StringIO(text)))[2:]          # header, units
    return [parse_time(r[0]) for r in rows if r]


def select_times(values: list[dt.datetime], step_hours: int) -> tuple[dt.datetime, dt.datetime, int, int]:
    """First and last chosen time, the index stride, and the count.

    Times are taken every step_hours, aligned to 00:00 UTC (the synoptic hours
    for steps of 3, 6 or 12 h). If no time of the axis falls on that alignment,
    as for daily means stamped at 12:00, the first time in the window sets it.
    The chosen times must be evenly spaced in the axis, so that a single strided
    request returns exactly them.
    """
    if not values:
        raise ErddapError("the dataset has no time in the requested window")
    step = dt.timedelta(hours=step_hours)

    def phase(t: dt.datetime) -> dt.timedelta:
        return (t - t.replace(hour=0, minute=0, second=0, microsecond=0)) % step

    target = dt.timedelta(0) if any(phase(t) == dt.timedelta(0) for t in values) else phase(values[0])
    chosen = [i for i, t in enumerate(values) if phase(t) == target]
    if not chosen:
        raise ErddapError(f"no time on a {step_hours} h step in the requested window")
    if len(chosen) == 1:
        return values[chosen[0]], values[chosen[0]], 1, 1
    stride = chosen[1] - chosen[0]
    if any(b - a != stride for a, b in zip(chosen, chosen[1:])):
        raise ErddapError(f"the times on a {step_hours} h step are not evenly spaced in the axis; "
                          "choose another step or a shorter window")
    return values[chosen[0]], values[chosen[-1]], stride, len(chosen)


def constraint(spec: dict) -> str:
    """'[(a):s:(b)]' or '[(a)]' from {'from', 'to', 'stride'} (values or ISO times)."""
    a = spec["from"]
    if spec.get("to") is None:
        return f"[({a})]"
    return f"[({a}):{spec.get('stride', 1)}:({spec['to']})]"


def request_url(service: str, dataset: str, variables: list[str], axes: list[dict], fmt: str = "nc") -> str:
    """The griddap URL for variables sharing the given ordered axis constraints."""
    box = "".join(constraint(a) for a in axes)
    query = ",".join(f"{v}{box}" for v in variables)
    return f"{service}/griddap/{dataset}.{fmt}?{urllib.parse.quote(query, safe='():,.-')}"
