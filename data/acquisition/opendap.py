"""A minimal client for OPeNDAP servers (DAP2), and their responses as netCDF files.

An OPeNDAP server publishes each dataset at a URL, with three views of it:

    <dataset>.das                   attributes of each variable
    <dataset>.dds                   structure: variables, their types and shapes
    <dataset>.dods?v[a:s:b][...]    values of a hyperslab (start a, stride s,
                                    last index b, per dimension), in binary

The binary response repeats the structure of the subset as text, then the line
"Data:", then the values in XDR (big-endian): for each array, its length twice,
then its values, with 16-bit integers widened to 32 bits and bytes padded to a
multiple of four. A Grid is its array followed by its coordinate maps.

Many archives in the Earth sciences are served this way (THREDDS, Hyrax,
GrADS Data Server). Here a subset of one variable becomes a netCDF classic file
with the variable's attributes and coordinates, so that the product builders
read it like any other.

Reference: Gallagher, J., et al. (2007), "The Data Access Protocol, DAP 2.0",
NASA ESE-RFC-004.1.1.
"""

from __future__ import annotations

import re
import struct
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "formats"))
import netcdf3  # noqa: E402

AGENT = {"User-Agent": "laboratory-data-fetch/1 (+https://github.com/huangwanhong-maker/Integrative_Astronomy_Geography_Oceanography_Laboratory)"}
# DAP2 base type -> (XDR dtype on the wire, bytes per value, dtype written to netCDF)
TYPES = {"Byte": (">u1", 1, "int16"), "Int16": (">i4", 4, "int16"), "UInt16": (">u4", 4, "int32"),
         "Int32": (">i4", 4, "int32"), "UInt32": (">u4", 4, "float64"),
         "Float32": (">f4", 4, "float32"), "Float64": (">f8", 8, "float64")}
ATTRIBUTE_TYPES = {"Byte": np.uint8, "Int16": np.int16, "UInt16": np.uint16, "Int32": np.int32,
                   "UInt32": np.uint32, "Float32": np.float32, "Float64": np.float64}
DECLARATION = re.compile(r"\b(Byte|Int16|UInt16|Int32|UInt32|Float32|Float64|String|Url)\s+([\w.%-]+)\s*((?:\[[^\]]*\]\s*)*);")
DIMENSION = re.compile(r"\[\s*(?:([\w.%-]+)\s*=\s*)?(\d+)\s*\]")
HORIZONTAL = {"lat", "latitude", "lon", "longitude", "x", "y"}


class DapError(RuntimeError):
    pass


def _get(url: str, timeout: int = 120) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=AGENT), timeout=timeout) as r:
        return r.read()


def declarations(text: str) -> list[tuple[str, str, list[tuple[str, int]]]]:
    """Base-type declarations of a DDS, in order: (type, name, [(dimension, length)])."""
    return [(t, name, [(d or name, int(n)) for d, n in DIMENSION.findall(dims)])
            for t, name, dims in DECLARATION.findall(text)]


def dds(url: str) -> dict:
    """The variables of a dataset: name -> (type, [(dimension, length)]); a Grid gives its array."""
    out = {}
    for t, name, dims in declarations(_get(url + ".dds").decode("utf-8")):
        out.setdefault(name, (t, dims))
    return out


def _attribute_values(kind: str, raw: str):
    if kind in ("String", "Url"):
        values = [s.encode("latin-1", "backslashreplace").decode("unicode_escape")
                  for s in re.findall(r'"((?:[^"\\]|\\.)*)"', raw)]
        return values[0] if len(values) == 1 else "\n".join(values)
    numbers = [x.strip() for x in raw.split(",") if x.strip()]
    dtype = ATTRIBUTE_TYPES[kind]
    a = np.array([float(x) if kind.startswith("Float") else int(x) for x in numbers], dtype=dtype)
    return a[0] if a.size == 1 else a


def das(url: str) -> dict:
    """Attributes: variable name (or NC_GLOBAL) -> {attribute: value}, values typed as declared."""
    text = _get(url + ".das").decode("utf-8")
    out, stack = {}, []
    token = re.compile(r'\s*(?:(\})|([\w.%-]+)\s*\{|(\w+)\s+([\w.%-]+)\s+((?:"(?:[^"\\]|\\.)*"|[^;"])*);)', re.S)
    pos = text.find("{") + 1
    while pos < len(text):
        m = token.match(text, pos)
        if not m:
            break
        pos = m.end()
        if m.group(1):
            if stack:
                stack.pop()
        elif m.group(2):
            stack.append(m.group(2))
            out.setdefault(m.group(2), {})
        elif m.group(3) in ATTRIBUTE_TYPES or m.group(3) in ("String", "Url"):
            out.setdefault(stack[-1] if stack else "NC_GLOBAL", {})[m.group(4)] = _attribute_values(m.group(3), m.group(5))
    return out


def hyperslab(dims: list[tuple[str, int]], stride: int) -> str:
    """The constraint that takes every `stride`-th point along the horizontal dimensions, all of the others."""
    return "".join(f"[0:{stride if d.lower() in HORIZONTAL else 1}:{n - 1}]" for d, n in dims)


def request_url(url: str, variable: str, constraint: str) -> str:
    return f"{url}.dods?{urllib.parse.quote(variable + constraint, safe=':,')}"


def decode(data: bytes) -> list[tuple[str, str, list[tuple[str, int]], np.ndarray]]:
    """The arrays of a binary (.dods) response, in order: (type, name, dims, values)."""
    mark = data.find(b"\nData:\n")
    if mark < 0:
        raise DapError(data[:400].decode("utf-8", "replace"))
    pos = mark + len(b"\nData:\n")
    out = []
    for t, name, dims in declarations(data[:mark].decode("utf-8")):
        if t not in TYPES:
            raise DapError(f"{name}: values of type {t} are not supported")
        wire, size, _ = TYPES[t]
        if dims:
            n, n2 = struct.unpack_from(">II", data, pos)
            if n != n2:
                raise DapError(f"{name}: inconsistent lengths {n} and {n2}")
            pos += 8
        else:
            n = 1
        values = np.frombuffer(data, wire, n, pos)
        pos += n * size + ((4 - n % 4) % 4 if size == 1 else 0)
        shape = tuple(length for _, length in dims) or ()
        out.append((t, name, dims, values.reshape(shape)))
    return out


def coordinate(url: str, name: str, timeout: int = 120) -> np.ndarray:
    """The values of a one-dimensional coordinate variable."""
    arrays = decode(_get(request_url(url, name, ""), timeout))
    return arrays[0][3].astype(np.float64)


CF_TIME = re.compile(r"\s*(seconds|minutes|hours|days)\s+since\s+(\d{4}-\d{1,2}-\d{1,2})(?:[ T](\d{1,2}:\d{2}(?::\d{2}(?:\.\d*)?)?))?", re.I)
SECONDS = {"seconds": 1, "minutes": 60, "hours": 3600, "days": 86400}


def cf_epoch_seconds(values: np.ndarray, units: str) -> np.ndarray:
    """Times in seconds since 1970-01-01 UTC from CF values and units ('hours since 2026-10-12 12:00:00 UTC')."""
    import datetime as dt
    m = CF_TIME.match(units)
    if not m:
        raise DapError(f"time units {units!r} are not understood")
    unit, day, clock = m.groups()
    y, mo, d = (int(x) for x in day.split("-"))
    hh, mm, ss = 0, 0, 0.0
    if clock:
        parts = clock.split(":")
        hh, mm = int(parts[0]), int(parts[1])
        ss = float(parts[2]) if len(parts) > 2 else 0.0
    epoch = dt.datetime(y, mo, d, hh, mm, tzinfo=dt.timezone.utc).timestamp() + ss
    return epoch + np.asarray(values, np.float64) * SECONDS[unit.lower()]


def runs(indices: list[int], longest: int) -> list[tuple[int, int, int]]:
    """Indices as hyperslab runs (start, stride, last) of at most `longest` indices each."""
    out, i = [], 0
    while i < len(indices):
        stride = indices[i + 1] - indices[i] if i + 1 < len(indices) else 1
        j = i + 1
        while j < len(indices) and j - i < longest and indices[j] - indices[j - 1] == stride:
            j += 1
        out.append((indices[i], max(1, stride), indices[j - 1]))
        i = j
    return out


def fetch_slabs(tasks: list[tuple[str, str, str]], workers: int = 4, timeout: int = 600) -> list[np.ndarray]:
    """The values of each (dataset url, variable, constraint), requested in parallel.

    A server that opens one file per time step answers a long request slowly;
    several short ones at once share the work (about three times faster with
    four, on the HYCOM server). Results are in the order of the tasks.
    """
    from concurrent.futures import ThreadPoolExecutor

    def one(task: tuple[str, str, str]) -> np.ndarray:
        url, variable, constraint = task
        for attempt in range(3):
            try:
                arrays = decode(_get(request_url(url, variable, constraint), timeout))
                return next(a for _, name, _, a in arrays if name == variable)
            except (OSError, DapError):
                if attempt == 2:
                    raise
        raise DapError("unreachable")

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        return list(pool.map(one, tasks))


def fetch_to_netcdf(url: str, variable: str, stride: int, dest: Path, timeout: int = 900) -> dict:
    """A strided subset of one variable (and its coordinates) as a netCDF classic file.

    Returns a description of the subset for provenance.
    """
    structure = dds(url)
    if variable not in structure:
        raise DapError(f"no variable {variable}; the dataset has {', '.join(structure)}")
    attributes = das(url)
    t, dims = structure[variable]
    request = request_url(url, variable, hyperslab(dims, stride))
    arrays = decode(_get(request, timeout))
    dimensions, variables = {}, {}
    for kind, name, adims, values in arrays:
        for d, n in adims:
            dimensions[d] = n
        attrs = {k: v for k, v in attributes.get(name, {}).items() if not k.startswith("_Chunk") and k != "dimensions"}
        variables[name] = ([d for d, _ in adims], values.astype(TYPES[kind][2]), attrs)
    glob = {k: v for k, v in attributes.get("NC_GLOBAL", {}).items() if isinstance(v, str)}
    glob["source"] = f"{url} (OPeNDAP), subset {variable}{hyperslab(dims, stride)}"
    netcdf3.write(dest, dimensions, variables, glob)
    return {"url": request, "variable": variable, "shape": [n for n in dimensions.values()],
            "dimensions": dimensions, "stride": stride}
