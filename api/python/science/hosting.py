"""The hosting of plug-ins that run in a server (ADR 0029), for any application
and for the API: Python scripts, imported and run on numpy, and C++ programs,
compiled and run as processes behind one file contract. An application binds a
`Host` to its own plugins/ directory (applications/earth-atlas/py/host.py) and
runs jobs in threads of its server (`Jobs`); the API runs each job in a process
of its own, one at a time and with a time limit (`Runs`).

A field result is written as a product under data/cache/plugins/<id>/ and named
by its URL under /data/plugins/, so that a page reads it as it reads any product
and the API serves it as a field; sketches, tables and notes are returned as
they are.

The file contract of a process plug-in (every line type is documented in
api/cpp/plugin.hpp, whose classes a C++ plug-in derives from): the host writes a
directory with input.txt, one float32 file per input variable and one
tab-separated file per observation set, runs the program with --run <dir>, and
reads output.txt and its float32 files back: fields, sketches, tables and notes.
`--describe` prints the plug-in's description as JSON, its class among them. The
generic parts of the classes run here, in plugins.py's code for every language:
a ForecastModel's recorded fields are replayed through plugins.ForecastModel, so
that its verification is the Python one; a DataAssimilation's analysis through
plugins.DataAssimilation, for its increment and the sketch of its observations.
The host compiles with the compiler at hand, or, on Windows without one, with
g++ in WSL. Only files under the application's plugins/ directory are run.

    python -m science.hosting describe applications/earth-atlas plugins/vorticity.py
    python -m science.hosting run applications/earth-atlas plugins/vorticity.py [need=product|level ...] [name=value ...]
    python -m science.hosting job          one run for the API: a request as JSON on stdin, progress and results as JSON lines
"""

from __future__ import annotations

import collections
import datetime as dt
import importlib.util
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import uuid
from pathlib import Path

import numpy as np

from . import fields
from . import plugins as api
from .local import Busy, Local

ROOT = Path(__file__).resolve().parents[3]
# The data library's cache, which SCIENCE_CACHE may move, as for the checks; and
# where runs write their products, which SCIENCE_PLUGIN_OUTPUT moves (the API's
# checks do, so as not to replace a person's own results). An application's binding
# passes the paths its server serves instead.
CACHE = Path(os.environ.get("SCIENCE_CACHE") or ROOT / "data" / "cache")
PRODUCTS = CACHE / "products"
OUT = Path(os.environ.get("SCIENCE_PLUGIN_OUTPUT") or CACHE / "plugins")
UTC = dt.timezone.utc

# The Earth Atlas's order of models (src/layers/catalogue.js); an unknown model ranks last.
MODEL_ORDER = ["GFS", "AIGFS", "AIGEFS mean", "Pangu-Weather (run here)", "Pangu-Weather hindcast (run here)",
               "Shallow water (ours)", "Barotropic (ours)"]


def model_rank(model: str) -> int:
    return MODEL_ORDER.index(model) if model in MODEL_ORDER else len(MODEL_ORDER)


# The C++ classes of plug-ins, a header of the API's own (api/cpp/plugin.hpp), as
# the Python ones are plugins.py; a plug-in is compiled with -I beside it.
HEADER = Path(__file__).resolve().parents[2] / "cpp" / "plugin.hpp"
NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)     # no console flashes up when a server on Windows runs one


class Compiler:
    """A C++20 compiler: a command on this computer (c++, g++ or clang++), or, on
    Windows without one, g++ in WSL, which builds a Linux program that is run
    with `wsl -e`, its paths translated (/mnt/c/... for C:\\...)."""

    def __init__(self, command: list[str], wsl: bool = False):
        self.command, self.wsl = list(command), wsl

    def __repr__(self) -> str:
        return " ".join(self.command)

    @property
    def suffix(self) -> str:
        """The ending of the programs it builds: .wsl for WSL's, .exe for Windows's, none on Linux and macOS."""
        return ".wsl" if self.wsl else (".exe" if platform.system() == "Windows" else "")

    def path(self, path: Path) -> str:
        """A path as the compiler and the programs it builds see it."""
        return wsl_path(path) if self.wsl else str(path)

    def runner(self) -> list[str]:
        """What precedes a program it built on a command line."""
        return ["wsl", "-e"] if self.wsl else []


_wsl = {"checked": False, "compiler": None, "mounts": {}}
_wsl_lock = threading.Lock()


def compiler() -> Compiler | None:
    """The compiler at hand: a native one first; on Windows, else, g++ (or c++ or
    clang++) in WSL, looked for once (WSL may take seconds to start); or None."""
    native = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if native:
        return Compiler([native])
    if platform.system() != "Windows" or not shutil.which("wsl"):
        return None
    with _wsl_lock:
        if not _wsl["checked"]:
            _wsl["checked"] = True
            try:
                done = subprocess.run(["wsl", "-e", "sh", "-c", "command -v g++ || command -v c++ || command -v clang++"],
                                      capture_output=True, text=True, timeout=60, creationflags=NO_WINDOW)
                found = done.stdout.strip().splitlines()
                if done.returncode == 0 and found and found[0].startswith("/"):
                    _wsl["compiler"] = Compiler(["wsl", "-e", found[0]], wsl=True)
            except (OSError, subprocess.SubprocessError):
                pass
        return _wsl["compiler"]


def wsl_path(path: Path) -> str:
    """A Windows path as WSL sees it: the drive's mount point (asked of wslpath
    once per drive, /mnt/c by default) and the rest with forward slashes."""
    p = Path(path).resolve()
    drive = p.drive.upper()
    if len(drive) != 2 or drive[1] != ":":
        raise RuntimeError(f"{p} is not on a lettered drive, which WSL mounts")
    with _wsl_lock:
        if drive not in _wsl["mounts"]:
            done = subprocess.run(["wsl", "-e", "wslpath", "-a", f"{drive}/"], capture_output=True, text=True, timeout=60,
                                  creationflags=NO_WINDOW)
            mount = done.stdout.strip().rstrip("/")
            if done.returncode != 0 or not mount.startswith("/"):
                raise RuntimeError(f"WSL does not mount {drive} ({done.stderr.strip() or 'wslpath failed'})")
            _wsl["mounts"][drive] = mount
        mount = _wsl["mounts"][drive]
    return f"{mount}/{p.relative_to(p.anchor).as_posix()}"


class Library:
    """The field products on disk, by name: the library's, and those plug-ins made
    before; and the index's other products, for the context of a plug-in."""

    def __init__(self, products: Path = PRODUCTS, out: Path = OUT, root: Path = ROOT):
        self.products = Path(products)
        self.out = Path(out)
        self.meta: dict[str, dict] = {}
        self.cache: dict[str, fields.Field] = {}
        self.local = Local(root, cache=self.products.parent)
        index_path = self.products / "index.json"
        index = json.loads(index_path.read_text(encoding="utf-8")) if index_path.exists() else {}
        for entry in index.get("fields", []):
            path = self.products / entry["file"]
            if path.exists():
                self.meta[entry["name"]] = json.loads(path.read_text(encoding="utf-8"))
                self.meta[entry["name"]]["_path"] = str(path)
        for path in self.out.glob("*/fields-*.json"):
            meta = json.loads(path.read_text(encoding="utf-8"))
            meta["_path"] = str(path)
            self.meta[meta["name"]] = meta

    def fields(self) -> list[dict]:
        return list(self.meta.values())

    def field(self, name: str) -> fields.Field:
        if name not in self.cache:
            self.cache[name] = fields.read(Path(self.meta[name]["_path"]))
        return self.cache[name]

    def list(self, kind: str) -> list[dict]:
        return self.local.list(kind)

    def events(self, name: str) -> dict | None:
        try:
            return self.local.events(name)
        except LookupError:
            return None

    def observations(self, name: str) -> dict | None:
        try:
            return self.local.observations(name)
        except LookupError:
            return None


# ---- plug-ins: Python modules and C++ programs ------------------------------------------------


class ProcessPlugin:
    """A plug-in that is a program, compiled from C++ when its binary is older
    than its source or the header (api/cpp/plugin.hpp). It describes itself with
    --describe, its class among the rest, and runs on files (the contract is in
    the header); what it writes is read into the results a Python plug-in
    returns, and the generic parts of its class are added by plugins.py's own code."""

    def __init__(self, source: Path, out: Path = OUT, binaries: Path | None = None):
        self.source = Path(source)
        self.out = Path(out)
        self.cxx = compiler()
        # The programs are build products, kept apart from the results (a check may write its results elsewhere).
        self.binaries = Path(binaries) if binaries else self.out / "bin"
        self.binary = self.binaries / f"{self.source.stem}{self.cxx.suffix if self.cxx else ''}"
        self.spec = self.describe()
        self.id = self.spec["id"]
        self.title = self.spec.get("title", self.source.stem)
        self.description = self.spec.get("description", "")
        self.domain = self.spec.get("domain", "meteorology")
        self.parameters = self.spec.get("parameters", {})
        self.kind = self.spec.get("class", "Diagnostic")

    @property
    def product_name(self) -> str:
        return f"plugin-{self.id}"

    def needs(self) -> list[dict]:
        return self.spec.get("needs", [])

    def compile(self) -> None:
        cxx = self.cxx
        if not cxx:
            raise RuntimeError("no C++ compiler on this computer (c++, g++ or clang++, or on Windows g++ in WSL, is needed to run a C++ plug-in)")
        newest = max(self.source.stat().st_mtime, HEADER.stat().st_mtime if HEADER.exists() else 0)
        if self.binary.exists() and self.binary.stat().st_mtime >= newest:
            return
        self.binary.parent.mkdir(parents=True, exist_ok=True)
        building = self.binary.with_name(f"{self.binary.name}.{uuid.uuid4().hex[:8]}.part")     # two runs may build at once
        done = subprocess.run([*cxx.command, "-std=c++20", "-O2", "-I", cxx.path(HEADER.parent), "-o", cxx.path(building), cxx.path(self.source)],
                              cwd=self.source.parent, capture_output=True, text=True, encoding="utf-8", errors="replace",
                              creationflags=NO_WINDOW)
        if done.returncode != 0:
            building.unlink(missing_ok=True)
            raise RuntimeError(f"{self.source.name} did not compile ({cxx}):\n{done.stderr[-1500:]}")
        os.replace(building, self.binary)

    def command(self, *args: str) -> list[str]:
        return [*self.cxx.runner(), self.cxx.path(self.binary), *args]

    def describe(self) -> dict:
        self.compile()
        done = subprocess.run(self.command("--describe"), capture_output=True, text=True, encoding="utf-8", errors="replace",
                              creationflags=NO_WINDOW)
        if done.returncode != 0:
            raise RuntimeError(f"{self.source.name} --describe failed:\n{done.stderr[-800:]}")
        return json.loads(done.stdout)

    def run(self, inputs: dict, params: dict, context) -> list[dict]:
        needs = self.needs()
        first = inputs.get(needs[0]["name"]) if needs else None
        if not isinstance(first, api.Input):
            raise ValueError(f"{self.source.name}: the first need must be a field, and given")
        # The described defaults beneath the values given, so that the program and the replay count alike.
        params = {**{k: p.get("default") for k, p in self.parameters.items()}, **{k: v for k, v in (params or {}).items() if v is not None}}
        (self.out / "tmp").mkdir(parents=True, exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix="run-", dir=str(self.out / "tmp")))
        try:
            self._write_input(work, inputs, params, context)
            context.progress(0.05, "running")
            self._execute(work, context)
            out = self._read_output(work, first, context)
            return self._results(out, inputs, params, context)
        finally:
            shutil.rmtree(work, ignore_errors=True)

    # ---- the files of a run (api/cpp/plugin.hpp) ----

    def _write_input(self, work: Path, inputs: dict, params: dict, context) -> None:
        lines = [f"time\t{round(context.time.timestamp() * 1000)}"]       # the page's whole ms, kept whole
        if context.pinned:
            lines.append(f"pinned\t{float(context.pinned[0])}\t{float(context.pinned[1])}")
        for name, value in params.items():
            lines.append(f"param\t{_cell(name)}\t{_cell(value)}")
        blocks = 0
        for need_name, inp in inputs.items():
            if inp is None:                       # an optional need left without a source
                continue
            if isinstance(inp, api.Observed):
                file = f"in.{need_name}.tsv"
                _write_rows(work / file, inp.observations or {})
                lines.append("\t".join(["observations", _cell(need_name), _cell(inp.name), _cell(inp.model), _cell(inp.variable or ""), file]))
                continue
            g, n = inp.product.grid, len(inp.product.times)
            lines.append("\t".join(map(str, ["need", need_name, inp.product.name, _cell(inp.model), inp.level or "", max(1, n),
                                              g["lon0"], g["dlon"], g["nlon"], g["lat0"], g["dlat"], g["nlat"], int(bool(g.get("periodic")))])))
            lines.append("\t".join(["times", need_name, *[fields.iso(t) for t in inp.product.times]]))
            for var, meta in zip(inp.names, inp.variables):
                data = np.stack([inp.product.slice(var, k) for k in range(max(1, n))]).astype("<f4")
                file = f"in.{blocks}.f32"
                blocks += 1
                data.tofile(work / file)
                lines.append("\t".join(["data", need_name, var, file, _cell(meta.get("standard_name") or ""),
                                        _cell(meta.get("units") or ""), _cell(meta.get("long_name") or "")]))
        (work / "input.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    def _execute(self, work: Path, context) -> None:
        """Runs the program on the directory, passing its progress lines on."""
        process = subprocess.Popen(self.command("--run", self.cxx.path(work)), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, encoding="utf-8", errors="replace", creationflags=NO_WINDOW)
        errors = collections.deque(maxlen=200)
        drain = threading.Thread(target=lambda: errors.extend(process.stderr), daemon=True)
        drain.start()
        for line in process.stdout:
            parts = line.rstrip("\r\n").split("\t")
            if parts[0] == "progress" and len(parts) > 1:
                try:
                    fraction = min(1.0, max(0.0, float(parts[1])))
                except ValueError:
                    continue
                context.progress(0.05 + 0.85 * fraction, parts[2] if len(parts) > 2 else "")
        process.wait()
        drain.join(5)
        if process.returncode != 0:
            raise RuntimeError(f"{self.source.name} failed:\n{''.join(errors)[-1500:].strip()}")

    def _read_output(self, work: Path, first: api.Input, context) -> dict:
        """output.txt as data: the fields (with their slices), the members of the
        classes (start, verify, observations), and the results that need nothing
        more (sketches, tables, notes), in the order written."""
        out = {"class": self.kind, "start": None, "grid": None, "times": None, "fields": [], "vectors": [],
               "show": dict(self.spec.get("show") or {}), "verify": None, "observations": [], "results": []}
        sketch = table = None

        def current_sketch():
            nonlocal sketch
            if sketch is None:
                sketch = api.sketch("")
                out["results"].append(sketch)
            return sketch
        for line in (work / "output.txt").read_text(encoding="utf-8").splitlines():
            kind, *rest = line.split("\t")
            pad = rest + [""] * 8
            if kind == "class":
                out["class"] = pad[0]
            elif kind == "start":
                out["start"] = float(pad[0])
            elif kind == "grid":
                out["grid"] = {"lon0": float(pad[0]), "dlon": float(pad[1]), "nlon": int(pad[2]), "lat0": float(pad[3]),
                               "dlat": float(pad[4]), "nlat": int(pad[5]), "periodic": pad[6] == "1"}
            elif kind == "times":
                out["times"] = [fields.parse_time(t) for t in rest if t]
            elif kind == "field":
                name, long_name, standard_name, units, level, difference, file = pad[:7]
                entry = {"name": name, "long_name": long_name, "standard_name": standard_name or None, "units": units,
                         "level": level or first.level, "_file": file}
                if difference == "1":
                    entry["difference"] = True
                out["fields"].append(entry)
            elif kind == "vector":
                out["vectors"].append(dict(zip(["name", "u", "v", "units", "standard_name", "level"], pad[:6])))
            elif kind == "show" and len(rest) >= 2:
                out["show"][rest[0]] = rest[1]
            elif kind == "verify":
                out["verify"] = dict(zip(["output", "need", "label", "unit"], pad[:4]))
            elif kind == "observation":
                o = {"lon": float(pad[0]), "lat": float(pad[1]), "value": _number(pad[2])}
                if math.isfinite(_number(pad[3])):
                    o["error"] = _number(pad[3])
                if pad[4]:
                    o["id"] = pad[4]
                out["observations"].append(o)
            elif kind == "sketch":
                sketch = api.sketch(pad[0])
                out["results"].append(sketch)
            elif kind == "polyline":
                current_sketch()["polylines"].append({"coords": _numbers(pad[4]), "colour": pad[0], "width": _number(pad[1]),
                                                      "alpha": _number(pad[2]), "dashed": pad[3] == "1"})
            elif kind == "marker":
                current_sketch()["markers"].append({"lon": _number(pad[0]), "lat": _number(pad[1]), "colour": pad[2], "size": _number(pad[3])})
            elif kind == "label":
                current_sketch()["labels"].append({"lon": _number(pad[0]), "lat": _number(pad[1]), "text": pad[2]})
            elif kind == "moving":
                current_sketch()["moving"].append({"times": [int(t) for t in _numbers(pad[2])], "coords": _numbers(pad[3]),
                                                   "colour": pad[0], "size": _number(pad[1])})
            elif kind == "table":
                table = api.table(pad[0], rest[2:], [], pad[1])
                out["results"].append(table)
            elif kind == "row":
                if table is None:
                    table = api.table("", [], [])
                    out["results"].append(table)
                table["rows"].append(rest)
            elif kind == "note":
                out["results"].append(api.note(pad[0]))
            elif kind == "error":
                raise RuntimeError(f"{self.source.name}: {pad[0] or 'the program reported an error'}")
        # The fields lie on the grid and at the times written; else (the first contract) on the first need's.
        out["grid"] = out["grid"] or first.product.grid
        if out["times"] is None:
            out["times"] = first.product.times or [context.time]
        elif not out["times"]:
            out["times"] = [context.time]
        g = out["grid"]
        for entry in out["fields"]:
            raw = np.fromfile(work / entry.pop("_file"), dtype="<f4")
            if raw.size != len(out["times"]) * g["nlat"] * g["nlon"]:
                raise RuntimeError(f"{self.source.name}: {entry['name']} has {raw.size} values, not {len(out['times'])} × {g['nlat']} × {g['nlon']}")
            entry["slices"] = list(raw.reshape(len(out["times"]), g["nlat"], g["nlon"]))
        return out

    def _results(self, out: dict, inputs: dict, params: dict, context) -> list[dict]:
        first = inputs[self.needs()[0]["name"]]
        if out["class"] == "ForecastModel":
            return [*self._forecast(out, inputs, params, context), *out["results"]]
        if out["class"] == "DataAssimilation":
            return [*self._analysis(out, inputs, params, context), *out["results"]]
        made = []
        if out["fields"]:
            product = api.make_product(self.product_name, self.title, out["grid"], out["times"], out["fields"], model=f"{self.title} (plug-in)",
                                       domain=self.domain, vectors=out["vectors"],
                                       source={"dataset": first.meta["name"], "provider": f"computed by a C++ program in the server from {first.model}",
                                               "credit": f"{self.title} (plug-in), from {first.model}"})
            made.append(api.field(product, out["show"] or {"surface": out["fields"][0]["name"]}))
        return [*made, *out["results"]]

    def _forecast(self, out: dict, inputs: dict, params: dict, context) -> list[dict]:
        if not out["fields"]:
            raise RuntimeError(f"{self.source.name} recorded no fields")
        replay = _Replayed(self, out)
        results = replay.run(inputs, params, context)
        made = results[0]["product"]
        if replay.taken != len(out["times"]) or [fields.iso(t) for t in made.times] != [fields.iso(t) for t in out["times"]]:
            raise RuntimeError(f"{self.source.name} recorded {len(out['times'])} times, and the loop of plugins.ForecastModel takes "
                               f"{replay.taken}: are its steps those of the header's ForecastModel?")
        made.meta["source"]["provider"] = made.meta["source"]["provider"].replace("computed in the server", "computed by a C++ program in the server", 1)
        return results

    def _analysis(self, out: dict, inputs: dict, params: dict, context) -> list[dict]:
        if not out["fields"]:
            raise RuntimeError(f"{self.source.name} returned no analysis")
        results = _Analysed(self, out).run(inputs, params, context)
        made = results[0]["product"]
        made.meta["source"]["provider"] = made.meta["source"]["provider"].replace("analysed in the server", "analysed by a C++ program in the server", 1)
        return results


class _Replayed(api.ForecastModel):
    """A forecast that a program ran, replayed through plugins.ForecastModel's own
    loop: initial() and advance() do nothing, and record() hands back the fields
    the program recorded, in order. The loop counts steps and records as the
    header's does (rounding half to even, both), so the times agree; and the
    verification against the input's later times beside persistence, its table
    and the product are made by the one implementation Python plug-ins use."""

    def __init__(self, plugin: ProcessPlugin, out: dict):
        self.id, self.title, self.description, self.domain = plugin.id, plugin.title, plugin.description, plugin.domain
        self.parameters = plugin.parameters
        self.dt = float(plugin.spec.get("dt", 1800))
        self.every = float(plugin.spec.get("every", 6))
        self._needs, self._out, self._verify = plugin.needs(), out, out["verify"] or plugin.spec.get("verify")
        self.taken = 0

    def needs(self):
        return self._needs

    def grid(self):
        return self._out["grid"]

    def outputs(self):
        return [{k: v for k, v in f.items() if k != "slices"} for f in self._out["fields"]]

    def vectors(self):
        return self._out["vectors"]

    def show(self):
        return self._out["show"] or super().show()

    def verify(self):
        return self._verify if self._verify and self._verify.get("output") else None

    def initial(self, inputs, index, params, context):
        start = self._out["start"]
        if start is not None and abs(fields.ms_of(inputs[self._needs[0]["name"]].product.times[index]) - start) > 1000:
            raise RuntimeError(f"{self.id}: the program started at {fields.iso(dt.datetime.fromtimestamp(start / 1000, UTC))}, "
                               f"and the loop at {fields.iso(inputs[self._needs[0]['name']].product.times[index])}")
        return None

    def advance(self, state, seconds, params):
        pass

    def record(self, state, params):
        k = self.taken
        if k >= len(self._out["times"]):
            raise RuntimeError(f"{self.id}: the program recorded fewer times than the loop takes")
        self.taken += 1
        return {f["name"]: f["slices"][k] for f in self._out["fields"]}

    def breathe(self, context, fraction, text=""):
        context.progress(0.9 + 0.1 * fraction, text)


class _Analysed(api.DataAssimilation):
    """An analysis that a program made, through plugins.DataAssimilation's own run:
    observations() and analyse() hand back what the program used and found, and the
    class makes the background at the instant, the increment, the product and the
    sketch of the observations."""

    def __init__(self, plugin: ProcessPlugin, out: dict):
        self.id, self.title, self.description, self.domain = plugin.id, plugin.title, plugin.description, plugin.domain
        self.parameters = plugin.parameters
        self._needs, self._out = plugin.needs(), out

    def needs(self):
        return self._needs

    def observations(self, inputs, when, params, context):
        start = self._out["start"]
        if start is not None and abs(fields.ms_of(when) - start) > 1000:
            raise RuntimeError(f"{self.id}: the program analysed at {fields.iso(dt.datetime.fromtimestamp(start / 1000, UTC))}, "
                               f"and the class at {fields.iso(when)}")
        return self._out["observations"]

    def analyse(self, background, grid, observations, params, context):
        analysis = self._out["fields"][0]["slices"][0]
        if analysis.shape != np.shape(background):
            raise RuntimeError(f"{self.id}: the analysis is {analysis.shape}, the background {np.shape(background)}")
        return analysis


def _cell(value) -> str:
    """A value as one cell of the files: tabs and line breaks as spaces; None and NaN empty."""
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return ""
    if isinstance(value, float):
        return repr(value)
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _number(text: str) -> float:
    try:
        return float(text)
    except ValueError:
        return float("nan")


def _numbers(text: str) -> list[float]:
    return [float(x) for x in text.split()]


def _write_rows(path: Path, observations: dict) -> None:
    """An observation set's rows as a tab-separated file, the first line the column names."""
    names = list(observations.get("fields") or [])
    rows = observations.get("rows") or observations.get("readings") or []
    if rows and isinstance(rows[0], dict):          # readings as objects: their keys are the columns
        names = names or list(dict.fromkeys(k for r in rows for k in r))
        rows = [[r.get(k) for k in names] for r in rows]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(_cell(n) for n in names) + "\n")
        for row in rows:
            f.write("\t".join(_cell(c) for c in row) + "\n")


# ---- needs resolved to the library's products ------------------------------------------------


def resolve(plugin, library: Library) -> list[dict]:
    """The products that can serve each need: one candidate per product and level
    carrying every standard name asked, the plug-in's own output excepted. A need
    of kind 'observations' is served by the observation sets of the library's
    index that carry its variable; an optional need may be left without a source
    (the candidate 'none', keyed ''). As the page's host resolves them (src/plugins/host.js)."""
    out = []
    for need in plugin.needs():
        if need.get("kind") == "observations":
            candidates = [{"key": e["name"], "product": e["name"], "model": e.get("title") or e["name"], "kind": "observations"}
                          for e in library.list("observations")
                          if not need.get("variable") or need["variable"] in _variable_names(e.get("variables"))]
            if need.get("optional"):
                candidates.append({"key": "", "model": "none"})
            out.append({**need, "candidates": candidates})
            continue
        standards = need.get("standard_names") or [need.get("standard_name")]
        candidates = []
        for name, meta in library.meta.items():
            if name == plugin.product_name or not meta.get("variables"):
                continue
            if need.get("times") and len(meta.get("times", [])) < need["times"]:
                continue
            plain = [v for v in meta["variables"] if "realization" not in (v.get("cell_methods") or "")]
            levels = [need["level"]] if need.get("level") else sorted({fields.level_of(v) for v in plain})
            for level in levels:
                variables = [next((v for v in plain if v.get("standard_name") == std and fields.level_of(v) == (level or "")), None) for std in standards]
                if not all(variables):
                    continue
                model = meta.get("model") or meta.get("title")
                candidates.append({"key": f"{name}|{level or ''}", "product": name, "level": level or None,
                                   "model": f"{model}, {level}" if level and not need.get("level") else model, "rank": model_rank(model),
                                   "names": [v["name"] for v in variables], "variables": variables, "meta": meta})
        candidates.sort(key=lambda c: (c["rank"], c["model"]))
        if need.get("optional"):
            candidates.append({"key": "", "model": "none"})
        out.append({**need, "candidates": candidates})
    return out


def _variable_names(variables) -> list[str]:
    """The variables an index entry lists, by name (as the observation sets list them) or as objects."""
    return [v if isinstance(v, str) else (v.get("standard_name") or v.get("name")) for v in variables or []]


class Host:
    """The plug-ins of one application: the files under its plugins/ directory."""

    def __init__(self, app: Path, products: Path = PRODUCTS, out: Path = OUT, binaries: Path | None = None):
        self.app = Path(app).resolve()
        self.plugins_dir = self.app / "plugins"
        self.products = Path(products)
        self.out = Path(out)
        self.binaries = Path(binaries) if binaries else None      # where C++ plug-ins are built; else out/bin

    def library(self) -> Library:
        return Library(self.products, self.out)

    def plugin_path(self, file: str) -> Path:
        """The file under plugins/, checked to be there."""
        path = (self.app / file).resolve()
        if self.plugins_dir.resolve() not in path.parents or not path.is_file():
            raise ValueError(f"{file} is not a plug-in of {self.app.name} (only files under plugins/ are run)")
        return path

    def files(self) -> list[dict]:
        """The plug-ins this server can run: Python and C++ files under plugins/."""
        out = []
        for path in sorted(self.plugins_dir.rglob("*")) if self.plugins_dir.is_dir() else []:
            if path.suffix == ".py" and path.is_file():
                out.append({"file": path.relative_to(self.app).as_posix(), "runtime": "python"})
            elif path.suffix in (".cpp", ".cc", ".cxx") and path.is_file():
                out.append({"file": path.relative_to(self.app).as_posix(), "runtime": "cpp",
                            **({} if compiler() else {"unavailable": "no C++ compiler on this computer (nor g++ in WSL)"})})
        return out

    def load_plugin(self, file: str):
        """The plug-in in `file`, and its runtime: 'python' or 'cpp'."""
        path = self.plugin_path(file)
        if path.suffix == ".py":
            spec = importlib.util.spec_from_file_location(f"plugin_{path.stem.replace('-', '_')}", path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            if not hasattr(module, "create"):
                raise ValueError(f"{file} does not define create(api)")
            return module.create(api), "python"
        if path.suffix in (".cpp", ".cc", ".cxx"):
            return ProcessPlugin(path, self.out, self.binaries), "cpp"
        raise ValueError(f"{file}: the server runs .py and .cpp plug-ins; the page runs .js")

    def describe(self, file: str, library: Library | None = None) -> dict:
        """What a page shows for a plug-in before it runs."""
        library = library or self.library()
        plugin, runtime = self.load_plugin(file)
        needs = resolve(plugin, library)
        chosen, taken = {}, set()
        for need in needs:
            pick = next((c for c in need["candidates"] if not need.get("distinct") or c["key"] not in taken), None) or (need["candidates"][0] if need["candidates"] else None)
            chosen[need["name"]] = pick["key"] if pick else None
            if pick:
                taken.add(pick["key"])
        return {
            "id": plugin.id, "title": plugin.title, "description": plugin.description, "runtime": runtime,
            "parameters": plugin.parameters,
            "needs": [{"name": n["name"], "label": n.get("label", n["name"]), "chosen": chosen[n["name"]],
                       **{k: n[k] for k in ("kind", "variable", "optional") if n.get(k)},
                       "candidates": [{"key": c["key"], "model": c["model"]} for c in n["candidates"]]} for n in needs],
        }

    def run(self, file: str, chosen: dict, params: dict, context: api.Context, library: Library | None = None) -> tuple[object, list[dict]]:
        """Runs a plug-in on the chosen sources and returns (plugin, results); fields are written as products."""
        plugin, _ = self.load_plugin(file)
        return plugin, self.run_plugin(plugin, chosen, params, context, library)

    def run_plugin(self, plugin, chosen: dict, params: dict, context: api.Context, library: Library | None = None) -> list[dict]:
        """Runs a plug-in already loaded (a Python one, or a ProcessPlugin) on the chosen
        sources and returns its results; fields are written as products under out/<id>/."""
        library = library or self.library()
        inputs = {}
        for need in resolve(plugin, library):
            c = next((x for x in need["candidates"] if x["key"] == (chosen or {}).get(need["name"])), None)
            if not c and not need.get("optional"):
                raise ValueError(f"nothing in the library serves \"{need.get('label', need['name'])}\"")
            if not c or c["key"] == "":
                inputs[need["name"]] = None
            elif c.get("kind") == "observations":
                observed = library.observations(c["product"])
                if observed is None:
                    raise ValueError(f"the observation set {c['product']!r} is listed, but cannot be read")
                inputs[need["name"]] = api.Observed(observed, c["product"], c["model"], need.get("variable"))
            else:
                inputs[need["name"]] = api.Input(library.field(c["product"]), c["names"], c["variables"], c["meta"], c["model"], c["level"])
        context._fields, context._field = library.fields, library.field
        context._list, context._events, context._observations = library.list, library.events, library.observations
        results = plugin.run(inputs, params, context) or []
        out_dir = self.out / plugin.id
        shutil.rmtree(out_dir, ignore_errors=True)
        written = []
        for r in results:
            if r["kind"] == "field":
                path = fields.write(r["product"], out_dir)
                written.append({"kind": "field", "url": f"data/plugins/{plugin.id}/{path.name}", "field": r["product"].name,
                                "show": r.get("show", {})})
            else:
                written.append(r)
        return written


def defaults(description: dict, chosen: dict | None = None, params: dict | None = None) -> tuple[dict, dict]:
    """The sources and parameters of a run: those given, and the described defaults for the rest."""
    picked = {n["name"]: n["chosen"] for n in description["needs"]}
    picked.update({k: v for k, v in (chosen or {}).items() if v is not None})      # '' chooses none, for an optional need
    values = {k: v.get("default") for k, v in (description.get("parameters") or {}).items()}
    values.update(params or {})
    return picked, values


def context_time(context: dict) -> dt.datetime:
    """The instant of a run's context: ms since 1970 (as the page sends it), ISO 8601, or now."""
    t = (context or {}).get("time")
    if t is None or t == "":
        return dt.datetime.now(UTC)
    if isinstance(t, (int, float)):
        return dt.datetime.fromtimestamp(float(t) / 1000, UTC)
    from .local import instant
    return instant(t)


class Jobs:
    """Runs of plug-ins for an application's server, one thread each; the page polls their state."""

    def __init__(self, host: Host):
        self.host = host
        self.jobs: dict[str, dict] = {}
        self.lock = threading.Lock()

    def start(self, file: str, chosen: dict, params: dict, context: dict) -> dict:
        self.host.plugin_path(file)
        ident = uuid.uuid4().hex[:12]
        job = {"id": ident, "file": file, "state": "running", "progress": 0.0, "text": "", "message": "",
               "started": time.time(), "finished": None, "results": []}
        with self.lock:
            self.jobs[ident] = job
        threading.Thread(target=self._run, args=(job, chosen, params, context), daemon=True).start()
        return self.view(job)

    def status(self, ident: str) -> dict:
        with self.lock:
            job = self.jobs.get(ident)
        if not job:
            raise ValueError(f"no job {ident!r}")
        return self.view(job)

    @staticmethod
    def view(job: dict) -> dict:
        return {k: v for k, v in job.items() if k != "file"}

    def _run(self, job: dict, chosen: dict, params: dict, context: dict) -> None:
        def progress(fraction, text=""):
            job["progress"], job["text"] = float(fraction), text
        ctx = api.Context(time=context_time(context), pinned=(context or {}).get("pinned"), progress=progress)
        try:
            _, job["results"] = self.host.run(job["file"], chosen or {}, params or {}, ctx)
            job["state"] = "done"
        except Exception as error:         # reported to the page, which shows it
            job["state"], job["message"] = "failed", f"{type(error).__name__}: {error}"
            traceback.print_exc()
        finally:
            job["finished"] = time.time()


class Runs:
    """Runs of plug-ins for the API: each in a process of its own (python -m
    science.hosting job), one at a time, stopped at a time limit. A run computes
    and writes only its own products, under data/cache/plugins/<id>/."""

    def __init__(self, root: Path = ROOT, *, seconds: float = 600.0, keep: int = 100, python: str | None = None):
        self.root = Path(root).resolve()
        self.seconds = seconds
        self.python = python or sys.executable
        self.jobs: collections.OrderedDict[str, dict] = collections.OrderedDict()
        self.keep = keep
        self.lock = threading.Lock()

    def applications(self) -> dict[str, Host]:
        base = self.root / "applications"
        return {p.name: Host(p) for p in sorted(base.iterdir()) if (p / "plugins").is_dir()} if base.is_dir() else {}

    def host(self, app: str) -> Host:
        hosts = self.applications()
        if app not in hosts:
            raise LookupError(f"no application {app!r} with plug-ins; there are {', '.join(hosts) or 'none'}")
        return hosts[app]

    def start(self, app: str, file: str, chosen: dict | None, params: dict | None, context: dict | None) -> dict:
        host = self.host(app)
        host.plugin_path(file)
        with self.lock:
            running = next((j for j in self.jobs.values() if j["state"] in ("starting", "running")), None)
            if running:
                raise Busy(f"a plug-in is running ({running['app']}/{running['file']}, job {running['id']}); "
                           "one runs at a time: try again when it is done")
            ident = uuid.uuid4().hex[:12]
            job = {"id": ident, "app": app, "file": file, "state": "starting", "progress": 0.0, "text": "", "message": "",
                   "started": time.time(), "finished": None, "limit_seconds": self.seconds, "results": []}
            self.jobs[ident] = job
            while len(self.jobs) > self.keep:
                self.jobs.popitem(last=False)
        request = {"app": str(host.app), "file": file, "chosen": chosen or {}, "params": params or {}, "context": context or {}}
        threading.Thread(target=self._run, args=(job, request), daemon=True).start()
        return self.view(job)

    def status(self, ident: str) -> dict:
        with self.lock:
            job = self.jobs.get(ident)
        if not job:
            raise LookupError(f"no job {ident!r}")
        return self.view(job)

    @staticmethod
    def view(job: dict) -> dict:
        return dict(job)

    def _run(self, job: dict, request: dict) -> None:
        env = {**os.environ, "PYTHONPATH": os.pathsep.join([str(Path(__file__).resolve().parents[1]), os.environ.get("PYTHONPATH", "")]),
               "PYTHONIOENCODING": "utf-8"}
        try:
            process = subprocess.Popen([self.python, "-m", "science.hosting", "job"], cwd=str(self.root), env=env,
                                       stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       text=True, encoding="utf-8", creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
        except OSError as error:
            job.update(state="failed", message=f"the run could not start: {error}", finished=time.time())
            return
        errors = collections.deque(maxlen=60)
        threading.Thread(target=lambda: [errors.append(line) for line in process.stderr], daemon=True).start()
        stopped = threading.Event()

        def stop():
            stopped.set()
            process.kill()
        timer = threading.Timer(self.seconds, stop)
        timer.start()
        job["state"] = "running"
        try:
            process.stdin.write(json.dumps(request))
            process.stdin.close()
            for line in process.stdout:
                try:
                    message = json.loads(line)
                except ValueError:
                    continue
                if "progress" in message:
                    job["progress"], job["text"] = float(message["progress"]), message.get("text", "")
                if "results" in message:
                    job["results"] = message["results"]
                if "error" in message:
                    job["message"] = message["error"]
            process.wait()
        finally:
            timer.cancel()
        if stopped.is_set():
            job.update(state="failed", message=f"the run took longer than {self.seconds:g} s and was stopped")
        elif process.returncode == 0 and not job["message"]:
            job.update(state="done", progress=1.0)
        else:
            job.update(state="failed", message=job["message"] or "".join(errors)[-1500:].strip() or f"the run ended with code {process.returncode}")
        job["finished"] = time.time()


def job_main() -> int:
    """One run for the API: the request on stdin; progress, then results, as JSON lines on stdout."""
    protocol = sys.stdout
    sys.stdout = sys.stderr          # what a plug-in prints is kept out of the protocol

    def say(message: dict) -> None:
        protocol.write(json.dumps(message, ensure_ascii=False, allow_nan=False, default=str) + "\n")
        protocol.flush()
    try:
        request = json.loads(sys.stdin.read())
        host = Host(Path(request["app"]))
        library = host.library()
        description = host.describe(request["file"], library)
        chosen, params = defaults(description, request.get("chosen"), request.get("params"))
        context = request.get("context") or {}
        ctx = api.Context(time=context_time(context), pinned=context.get("pinned"),
                          progress=lambda fraction, text="": say({"progress": float(fraction), "text": text}))
        _, results = host.run(request["file"], chosen, params, ctx, library)
        say({"results": _finite(results)})
        return 0
    except Exception as error:      # reported in the job
        traceback.print_exc()
        say({"error": f"{type(error).__name__}: {error}"})
        return 1


def _finite(value):
    """JSON-safe: NaN and infinities as None, numpy scalars and arrays as numbers and lists."""
    if isinstance(value, dict):
        return {k: _finite(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite(v) for v in value]
    if isinstance(value, np.ndarray):
        return _finite(value.tolist())
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def main(argv: list[str]) -> int:
    if argv[:1] == ["job"]:
        return job_main()
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")       # a Windows console's code page cannot print every label (ε²)
    if len(argv) < 3 or argv[0] not in ("describe", "run"):
        print(__doc__)
        return 2
    host = Host(Path(argv[1]))
    if argv[0] == "describe":
        print(json.dumps(host.describe(argv[2]), indent=1, ensure_ascii=False))
        return 0
    description = host.describe(argv[2])
    chosen, params = defaults(description)
    for item in argv[3:]:
        key, _, value = item.partition("=")
        if key in chosen:
            chosen[key] = value
        else:
            params[key] = float(value) if value.replace(".", "", 1).lstrip("-").isdigit() else value
    _, results = host.run(argv[2], chosen, params, api.Context(time=dt.datetime.now(UTC), progress=lambda f, t="": print(f"  {f:4.0%} {t}")))
    for r in results:
        print(r.get("url") or r.get("title") or r.get("text") or r["kind"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
