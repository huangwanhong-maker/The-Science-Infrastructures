"""A minimal client for the NetCDF Subset Service (NCSS) of THREDDS data servers.

Unidata's THREDDS servers publish numerical weather prediction models (GFS,
NOAA's AI models AIGFS and AIGEFS, and others) as "best time series", which
join past analyses to the latest forecast. NCSS cuts a subset from such a
collection and returns it as a netCDF classic file::

    <service>/<dataset>?var=a&var=b&vertCoord=50000&time_start=...&time_end=...
        &timeStride=2&horizStride=4&accept=netcdf3

Variables requested together must share their coordinate system (time axis
and vertical level), so a manifest groups them, one request per group.

Reference: https://docs.unidata.ucar.edu/tds/current/userguide/netcdf_subset_service_ref.html
"""

from __future__ import annotations

import datetime as dt
import urllib.parse


def iso(t: dt.datetime) -> str:
    return t.astimezone(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def align(t: dt.datetime, step_hours: int, up: bool = False) -> dt.datetime:
    """t moved to a multiple of step_hours since 00:00 UTC (down, or up)."""
    midnight = t.replace(hour=0, minute=0, second=0, microsecond=0)
    step = dt.timedelta(hours=step_hours)
    k = (t - midnight) // step
    aligned = midnight + k * step
    if up and aligned < t:
        aligned += step
    return aligned


def request_url(service: str, dataset: str, variables: list[str], start: dt.datetime, end: dt.datetime,
                time_stride: int = 1, horiz_stride: int = 1, vert_coord: float | None = None) -> str:
    query = [("var", v) for v in variables]
    query += [("time_start", iso(start)), ("time_end", iso(end)), ("timeStride", str(time_stride)),
              ("horizStride", str(horiz_stride))]
    if vert_coord is not None:
        query.append(("vertCoord", f"{vert_coord:g}"))
    query.append(("accept", "netcdf3"))
    return f"{service.rstrip('/')}/{dataset}?{urllib.parse.urlencode(query)}"
