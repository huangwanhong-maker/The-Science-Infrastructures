"""Plug-ins over networks of places (ADR 0046): scripts in Python and programs in
C++ that analyse the flows between places and answer what an application draws
on its map. They belong to no application: the Economy Atlas binds a `Host` to
its own plugins/ directory, as the Earth Atlas binds hosting.py to its own.

A plug-in is given a network:

    {"title", "year", "unit",
     "places": [{"id", "name", "latitude", "longitude", "group"}],    degrees on WGS 84
     "flows": [{"from", "to", "value", "kind"?}],                      place ids, the network's unit
     "decomposition": {...} | null}                                     the application's own reading, passed through

and answers results, data the application draws:

    values  a number for each place, coloured on its country or its circle (values());
    arcs    flows between places, drawn as arcs (arcs());
    sketch  polylines, markers and labels in longitude and latitude (plugins.sketch);
    table   rows in the panel (plugins.table);
    note    text in the panel (plugins.note).

A Python script exports

    def create(api):
        class Mine(api.NetworkAnalysis):
            id = "mine"; title = "..."; params = [api.param("k", "k", 3, kind="integer")]
            def run(self, network, params, context):
                return [api.values("balance", network.balance(), scale="diverging")]
        return Mine()

where `api` is this module; a C++ program derives from science::network::Analysis
(api/cpp/network_plugin.hpp) and ends with SCIENCE_NETWORK_PLUGIN(Mine).

Both languages run as processes behind one file contract, so that a script can
neither stop the server nor outlive its time: the host writes input.json,
{"network": ..., "params": ...}, into a directory of its own and runs the plug-in
with --run <dir> (a script through `python -m science.networks run <script> <dir>`);
the plug-in reports progress on standard output as lines "progress\\t<fraction>\\t<text>",
as hosting.py's C++ plug-ins do, and writes output.json, {"results": [...]}.
--describe (`python -m science.networks describe <script>`) prints its description.
Parameters reach the plug-in defaulted from its description and typed by kind.

    python -m science.networks describe <script>
    python -m science.networks run <script> <dir>
    python -m science.networks try <plugins dir> <file> <network.json> [-I<include dir> ...] [name=value ...]
"""

from __future__ import annotations

import collections
import importlib.util
import json
import math
import os
import re
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

from .plugins import note, sketch, table  # noqa: F401  (results shared with the field plug-ins)

EARTH_RADIUS_KM = 6371.0088          # the mean radius (IUGG), as the economics domain's places take it
PYTHON_PACKAGE = Path(__file__).resolve().parents[1]   # api/python
ROOT = PYTHON_PACKAGE.parents[1]
CPP = ROOT / "api" / "cpp"
HEADER = CPP / "network_plugin.hpp"
KINDS = ("values", "arcs", "sketch", "table", "note")
PARAM_KINDS = ("number", "integer", "choice", "boolean")
INCLUDE = re.compile(r'^\s*#\s*include\s*[<"]([^>"]+)[>"]', re.M)


# ---- the network --------------------------------------------------------------------------------


class Network:
    """The places and the flows between them, as a plug-in reads them: flows as
    index arrays `source`, `target` and their `value`, places by `index[id]`."""

    def __init__(self, data: dict):
        self.title = str(data.get("title") or "")
        year = data.get("year")
        self.year = int(year) if year not in (None, "") else None
        self.unit = str(data.get("unit") or "")
        self.places = []
        for p in data.get("places") or []:
            if p.get("id") in (None, ""):
                raise ValueError("a place without an id")
            self.places.append({"id": str(p["id"]), "name": str(p.get("name") or p["id"]),
                                "latitude": float(p["latitude"]), "longitude": float(p["longitude"]),
                                "group": str(p.get("group") or "")})
        self.ids = [p["id"] for p in self.places]
        self.index = {pid: i for i, pid in enumerate(self.ids)}
        if len(self.index) != len(self.ids):
            twice = [pid for pid, c in collections.Counter(self.ids).items() if c > 1]
            raise ValueError(f"two places share the id {twice[0]!r}")
        self.n = len(self.places)
        source, target, value, kinds = [], [], [], []
        for f in data.get("flows") or []:
            for end in ("from", "to"):
                if str(f.get(end)) not in self.index:
                    raise ValueError(f"a flow {end} {f.get(end)!r}: no place has that id")
            source.append(self.index[str(f["from"])])
            target.append(self.index[str(f["to"])])
            value.append(float(f["value"]))
            kinds.append(f.get("kind"))
        self.source = np.asarray(source, dtype=np.int64)
        self.target = np.asarray(target, dtype=np.int64)
        self.value = np.asarray(value, dtype=float)
        self.kinds = kinds if any(k is not None for k in kinds) else None
        self.decomposition = data.get("decomposition")

    def place(self, pid: str) -> dict:
        return self.places[self.index[pid]]

    def matrix(self) -> np.ndarray:
        """F[i, j], the flow from place i to place j, flows between the same two places added."""
        F = np.zeros((self.n, self.n))
        np.add.at(F, (self.source, self.target), self.value)
        return F

    def exports(self) -> np.ndarray:
        return np.bincount(self.source, weights=self.value, minlength=self.n)

    def imports(self) -> np.ndarray:
        return np.bincount(self.target, weights=self.value, minlength=self.n)

    def balance(self) -> np.ndarray:
        """The net inflow of each place, imports − exports (the domain's ∂F)."""
        return self.imports() - self.exports()

    def distances_km(self) -> np.ndarray:
        """Great-circle distances between every two places on the sphere of mean radius, km.
        The angle is atan2(|u × v|, u · v), exact for near and antipodal places alike."""
        lat = np.radians([p["latitude"] for p in self.places])
        lon = np.radians([p["longitude"] for p in self.places])
        u = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=1)
        dot = np.clip(u @ u.T, -1.0, 1.0)
        cross = np.linalg.norm(np.cross(u[:, None, :], u[None, :, :]), axis=2)
        return EARTH_RADIUS_KM * np.arctan2(cross, dot)


# ---- results and parameters ---------------------------------------------------------------------


def values(title: str, values, unit: str = "", scale: str = "diverging", note: str = "", ids: list | None = None) -> dict:
    """A number for each place, by id: {id: number}, or a list aligned with the places (with
    `ids`, or completed by the runner, which knows them). NaN and None are left out.
    `scale`: "diverging" about 0 (drawn and sent alike), or "sequential" (more is darker)."""
    if scale not in ("diverging", "sequential"):
        raise ValueError(f"a scale is 'diverging' or 'sequential', not {scale!r}")
    if not isinstance(values, dict):
        values = list(np.asarray(values, dtype=float).ravel()) if not isinstance(values, list) else values
        if ids is not None:
            values = dict(zip(ids, values))
    return {"kind": "values", "title": title, "unit": unit, "scale": scale, "note": note, "values": values}


def arcs(title: str, arcs: list, unit: str = "", note: str = "") -> dict:
    """Flows between places, drawn as arcs: [{from, to, value}], ids of places."""
    out = []
    for a in arcs:
        v = float(a["value"])
        if math.isfinite(v):
            out.append({"from": str(a["from"]), "to": str(a["to"]), "value": v})
    return {"kind": "arcs", "title": title, "unit": unit, "note": note, "arcs": out}


def param(name: str, label: str, default, kind: str = "number", minimum=None, maximum=None, step=None,
          options=None, help: str = "") -> dict:
    """A parameter of a plug-in, from which the application builds its form. `options`, for a
    choice: values, or {value, label}."""
    if kind not in PARAM_KINDS:
        raise ValueError(f"a parameter's kind is one of {', '.join(PARAM_KINDS)}, not {kind!r}")
    if kind == "choice" and not options:
        raise ValueError(f"the choice {name!r} has no options")
    out = {"name": name, "label": label, "kind": kind, "default": default}
    for key, v in (("min", minimum), ("max", maximum), ("step", step), ("options", options)):
        if v is not None:
            out[key] = v
    if help:
        out["help"] = help
    return out


def typed(described: list, given: dict | None) -> dict:
    """The parameters a run receives: each described one, given or else its default, as its kind:
    numbers as floats within their bounds, integers as ints, booleans, choices among their options."""
    given = given or {}
    out = {}
    for p in described or []:
        name, kind = p["name"], p.get("kind", "number")
        v = given.get(name)
        if v is None or v == "":
            v = p.get("default")
        try:
            if kind == "boolean":
                v = v if isinstance(v, bool) else str(v).strip().lower() in ("1", "true", "yes", "on")
            elif kind == "choice":
                allowed = [o["value"] if isinstance(o, dict) else o for o in p.get("options", [])]
                v = next((o for o in allowed if str(o) == str(v)), None)
                if v is None:
                    raise ValueError(f"{given.get(name)!r} is not one of {allowed}")
            else:
                x = float(v)
                if not math.isfinite(x):
                    raise ValueError(f"{v!r} is not a finite number")
                if kind == "integer":
                    if x != math.floor(x):
                        raise ValueError(f"{v!r} is not a whole number")
                    x = int(x)
                if p.get("min") is not None and x < p["min"]:
                    raise ValueError(f"{x} is less than its least, {p['min']}")
                if p.get("max") is not None and x > p["max"]:
                    raise ValueError(f"{x} is more than its most, {p['max']}")
                v = x
        except (TypeError, ValueError) as error:
            raise ValueError(f"the parameter {name}: {error}") from None
        out[name] = v
    return out


class Context:
    """What a run sees besides its network and parameters: its progress, reported to the host."""

    def __init__(self, progress=None):
        self._progress = progress or (lambda fraction, text="": None)

    def progress(self, fraction: float, text: str = "") -> None:
        self._progress(fraction, text)


class NetworkAnalysis:
    """An analysis of a network of places: identity, parameters, and run(), which answers results."""

    id = ""
    title = ""
    description = ""
    params: list = []

    def describe(self) -> dict:
        ident = self.id or type(self).__name__.lower()
        return {"kind": "network-analysis", "id": ident, "title": self.title or ident, "description": self.description,
                "params": list(self.params), "language": "python"}

    def run(self, network: Network, params: dict, context: Context) -> list[dict]:
        raise NotImplementedError


# ---- results as JSON ----------------------------------------------------------------------------


def plain(value):
    """A value as JSON holds it: numpy's scalars and arrays as Python's, NaN and the infinities as null."""
    if isinstance(value, np.ndarray):
        return [plain(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value


def finish(results, network: Network) -> list[dict]:
    """The results of a run, checked against its network: each of a known kind, values completed
    by the places' ids and keyed by them, arcs between places of the network."""
    if isinstance(results, dict):
        results = [results]
    out = []
    for r in results or []:
        if not isinstance(r, dict) or r.get("kind") not in KINDS:
            raise ValueError(f"a result is a dict whose kind is one of {', '.join(KINDS)}: {str(r)[:80]}")
        r = dict(r)
        if r["kind"] == "values":
            vals = r["values"]
            if not isinstance(vals, dict):
                vals = list(vals)
                if len(vals) != network.n:
                    raise ValueError(f"the values {r.get('title')!r} are {len(vals)}, for {network.n} places")
                vals = dict(zip(network.ids, vals))
            kept = {}
            for pid, v in vals.items():
                if str(pid) not in network.index:
                    raise ValueError(f"the values {r.get('title')!r} name {pid!r}, which no place is")
                if v is not None and math.isfinite(float(v)):
                    kept[str(pid)] = float(v)
            r["values"] = kept
        elif r["kind"] == "arcs":
            for a in r["arcs"]:
                for end in ("from", "to"):
                    if a[end] not in network.index:
                        raise ValueError(f"an arc of {r.get('title')!r} {end} {a[end]!r}, which no place is")
        out.append(plain(r))
    return out


# ---- the runner: a script in a process of its own ----------------------------------------------


def load(script: Path):
    """The plug-in a script makes: its create(api), with this module as `api`."""
    script = Path(script)
    spec = importlib.util.spec_from_file_location(f"network_plugin_{script.stem.replace('-', '_')}", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "create"):
        raise ValueError(f"{script.name} does not define create(api)")
    from science import networks as api
    plugin = module.create(api)
    if not hasattr(plugin, "run") or not hasattr(plugin, "describe"):
        raise ValueError(f"{script.name}: create(api) answers no plug-in (an instance of api.NetworkAnalysis)")
    return plugin


def say_progress(fraction: float, text: str = "") -> None:
    fraction = min(1.0, max(0.0, float(fraction)))
    print(f"progress\t{fraction:.4f}\t{' '.join(str(text).split())}", flush=True)


def run_script(script: Path, work: Path) -> None:
    """One run of a script: input.json in, output.json out (the process contract)."""
    request = json.loads((Path(work) / "input.json").read_text(encoding="utf-8"))
    plugin = load(script)
    network = Network(request["network"])
    params = typed(plugin.describe().get("params", []), request.get("params"))
    results = finish(plugin.run(network, params, Context(say_progress)), network)
    (Path(work) / "output.json").write_text(json.dumps({"results": results}, ensure_ascii=False, allow_nan=False), encoding="utf-8")


# ---- the host -----------------------------------------------------------------------------------


def include_closure(source: Path, dirs: list[Path]) -> list[Path]:
    """The source and the headers it includes, transitively, found beside the including file or
    in the include directories; the standard library's are not followed."""
    seen: list[Path] = []
    todo = [Path(source).resolve()]
    while todo:
        f = todo.pop()
        if f in seen:
            continue
        seen.append(f)
        for name in INCLUDE.findall(f.read_text(encoding="utf-8", errors="replace")):
            for d in [f.parent, *dirs]:
                candidate = (d / name)
                if candidate.is_file():
                    todo.append(candidate.resolve())
                    break
    return seen


class Host:
    """The plug-ins of one application's plugins/ directory: listed, described, compiled (C++)
    and run, each run a process in a directory of its own under `work`. C++ programs are built
    with the compiler hosting.py finds (natively, or g++ in WSL) against `include_dirs` and
    api/cpp, into work/bin/, again when the source or any header it includes is newer."""

    def __init__(self, plugins: Path, work: Path, include_dirs=(), python: str | None = None):
        self.plugins = Path(plugins).resolve()
        self.work = Path(work)
        self.include_dirs = [Path(d) for d in include_dirs] + [CPP]
        self.python = python or sys.executable
        self._described: dict[Path, tuple] = {}
        self._lock = threading.Lock()

    def files(self) -> list[dict]:
        if not self.plugins.is_dir():
            return []
        return [{"file": f"plugins/{p.name}", "language": "python" if p.suffix == ".py" else "cpp"}
                for p in sorted(self.plugins.iterdir())
                if p.is_file() and p.suffix in (".py", ".cpp") and not p.name.startswith("_")]

    def plugin_path(self, file: str) -> Path:
        """The plug-in named by `file` ("plugins/x.py" or "x.py"), refused when it is not a .py or
        .cpp file directly in plugins/."""
        name = str(file).replace("\\", "/")
        name = name[len("plugins/"):] if name.startswith("plugins/") else name
        path = (self.plugins / name).resolve()
        if path.parent != self.plugins or path.suffix not in (".py", ".cpp") or not path.is_file():
            raise ValueError(f"{file}: not a plug-in of {self.plugins.name}/ (a .py or .cpp file directly in it)")
        return path

    # ---- C++ ----

    def _compiler(self):
        from . import hosting
        cxx = hosting.compiler()
        if not cxx:
            raise RuntimeError("no C++ compiler on this computer (c++, g++ or clang++, or on Windows g++ in WSL, is needed to run a C++ plug-in)")
        return cxx

    def binary(self, path: Path) -> tuple:
        """The compiled program of a C++ plug-in, built when missing or older than its sources."""
        from .hosting import NO_WINDOW
        cxx = self._compiler()
        binary = self.work / "bin" / f"{path.stem}{cxx.suffix}"
        newest = max(f.stat().st_mtime for f in include_closure(path, self.include_dirs))
        with self._lock:
            if binary.exists() and binary.stat().st_mtime >= newest:
                return cxx, binary
            binary.parent.mkdir(parents=True, exist_ok=True)
            building = binary.with_name(f"{binary.name}.{uuid.uuid4().hex[:8]}.part")
            includes = [arg for d in self.include_dirs if d.is_dir() for arg in ("-I", cxx.path(d))]
            done = subprocess.run([*cxx.command, "-std=c++20", "-O2", *includes, "-o", cxx.path(building), cxx.path(path)],
                                  cwd=path.parent, capture_output=True, text=True, encoding="utf-8", errors="replace",
                                  creationflags=NO_WINDOW)
            if done.returncode != 0:
                building.unlink(missing_ok=True)
                raise RuntimeError(f"{path.name} did not compile ({cxx}):\n{done.stderr[-2000:]}")
            os.replace(building, binary)
        return cxx, binary

    # ---- describing and running ----

    def _environment(self) -> dict:
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join([str(PYTHON_PACKAGE), *([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])])
        env["PYTHONIOENCODING"] = "utf-8"
        return env

    def _command(self, path: Path, *args) -> tuple[list[str], object]:
        """The command of a description (args "--describe") or a run (args "--run", <dir>), and the
        compiler of a C++ plug-in (None for a script)."""
        if path.suffix == ".py":
            tail = ["describe", str(path)] if args[0] == "--describe" else ["run", str(path), str(args[1])]
            return [self.python, "-m", "science.networks", *tail], None
        cxx, binary = self.binary(path)
        tail = ["--describe"] if args[0] == "--describe" else ["--run", cxx.path(Path(args[1]))]
        return [*cxx.runner(), cxx.path(binary), *tail], cxx

    def describe(self, file: str) -> dict:
        """The plug-in's description, its file added; kept until its source (or, for C++, a header
        it includes) changes."""
        from .hosting import NO_WINDOW
        path = self.plugin_path(file)
        key = (path.stat().st_mtime_ns,) if path.suffix == ".py" else \
            tuple(f.stat().st_mtime_ns for f in include_closure(path, self.include_dirs))
        cached = self._described.get(path)
        if cached and cached[0] == key:
            return dict(cached[1])
        command, _ = self._command(path, "--describe")
        done = subprocess.run(command, cwd=ROOT, env=self._environment(), capture_output=True, text=True, encoding="utf-8",
                              errors="replace", timeout=120, creationflags=NO_WINDOW)
        if done.returncode != 0:
            raise RuntimeError(f"{path.name} --describe failed:\n{done.stderr[-1500:].strip()}")
        description = json.loads(done.stdout)
        description["file"] = f"plugins/{path.name}"
        description.setdefault("language", "python" if path.suffix == ".py" else "cpp")
        description.setdefault("kind", "network-analysis")
        description.setdefault("params", [])
        self._described[path] = (key, description)
        return dict(description)

    def run(self, file: str, network: dict, params: dict | None = None, progress=None, cancel: threading.Event | None = None,
            seconds: float = 600.0) -> list[dict]:
        """The results of one run of a plug-in on a network: `progress(fraction, text)` follows it,
        setting `cancel` stops it, and it is stopped after `seconds`."""
        from .hosting import NO_WINDOW
        path = self.plugin_path(file)
        description = self.describe(file)
        chosen = typed(description.get("params", []), params)
        (self.work / "tmp").mkdir(parents=True, exist_ok=True)
        work = Path(tempfile.mkdtemp(prefix="run-", dir=str(self.work / "tmp")))
        try:
            (work / "input.json").write_text(json.dumps({"network": plain(network), "params": chosen}, ensure_ascii=False,
                                                        allow_nan=False), encoding="utf-8")
            command, cxx = self._command(path, "--run", work)
            process = subprocess.Popen(command, cwd=ROOT, env=self._environment(), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       text=True, encoding="utf-8", errors="replace", creationflags=NO_WINDOW)
            errors = collections.deque(maxlen=400)
            drain = threading.Thread(target=lambda: errors.extend(process.stderr), daemon=True)
            drain.start()
            stopped: dict = {}

            def watch() -> None:
                deadline = time.monotonic() + seconds
                while process.poll() is None:
                    why = "cancelled" if cancel is not None and cancel.is_set() else \
                        f"stopped after {seconds:g} s" if time.monotonic() > deadline else None
                    if why:
                        stopped["why"] = why
                        process.kill()
                        if cxx is not None and cxx.wsl:   # killing wsl.exe need not end the program inside it
                            subprocess.run(["wsl", "-e", "pkill", "-f", cxx.path(work)], capture_output=True, creationflags=NO_WINDOW)
                        return
                    time.sleep(0.1)

            threading.Thread(target=watch, daemon=True).start()
            for line in process.stdout:
                parts = line.rstrip("\r\n").split("\t")
                if parts[0] == "progress" and len(parts) > 1 and progress:
                    try:
                        progress(min(1.0, max(0.0, float(parts[1]))), parts[2] if len(parts) > 2 else "")
                    except ValueError:
                        pass
            process.wait()
            drain.join(5)
            if stopped:
                raise RuntimeError(f"{path.name}: {stopped['why']}")
            if process.returncode != 0:
                raise RuntimeError(f"{path.name} failed:\n{''.join(errors)[-1500:].strip()}")
            out = json.loads((work / "output.json").read_text(encoding="utf-8"))
            return out.get("results", [])
        finally:
            shutil.rmtree(work, ignore_errors=True)


# ---- the command line ---------------------------------------------------------------------------


def summary(results: list[dict]) -> str:
    """A few lines on each result, for the command line."""
    lines = []
    for r in results:
        kind = r.get("kind")
        if kind == "values":
            top = sorted(r["values"].items(), key=lambda kv: -kv[1])
            lines.append(f"values  {r.get('title')} ({r.get('scale')}, {len(top)} places): highest "
                         + ", ".join(f"{k} {v:.4g}" for k, v in top[:5]) + "; lowest " + ", ".join(f"{k} {v:.4g}" for k, v in top[-3:]))
        elif kind == "arcs":
            lines.append(f"arcs    {r.get('title')}: {len(r['arcs'])}, e.g. "
                         + ", ".join(f"{a['from']}→{a['to']} {a['value']:.4g}" for a in r["arcs"][:4]))
        elif kind == "table":
            lines.append(f"table   {r.get('title')}: {', '.join(map(str, r['columns']))}")
            lines += ["          " + " | ".join(f"{c:.4g}" if isinstance(c, float) else str(c) for c in row) for row in r["rows"][:8]]
        elif kind == "sketch":
            lines.append(f"sketch  {r.get('title')}: {len(r.get('markers', []))} markers, {len(r.get('labels', []))} labels")
        else:
            lines.append(f"note    {r.get('text')}")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    for stream in (sys.stdout, sys.stderr):       # names and arrows print whatever the console's code page
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    if len(argv) >= 2 and argv[0] == "describe":
        print(json.dumps(plain(load(Path(argv[1])).describe()), ensure_ascii=False))
        return 0
    if len(argv) >= 3 and argv[0] == "run":
        try:
            run_script(Path(argv[1]), Path(argv[2]))
        except Exception as error:                  # the host reads the reason from the last lines
            traceback.print_exc()
            print(f"{type(error).__name__}: {error}", file=sys.stderr, flush=True)
            return 1
        return 0
    if len(argv) >= 4 and argv[0] == "try":
        plugins, file, network_file = Path(argv[1]), argv[2], Path(argv[3])
        params = dict(a.split("=", 1) for a in argv[4:] if "=" in a and not a.startswith("-I"))
        includes = [Path(a[2:]) for a in argv[4:] if a.startswith("-I")]
        cache = Path(os.environ.get("SCIENCE_CACHE") or ROOT / "data" / "cache")
        host = Host(plugins, cache / "network-plugins", includes)
        started = time.time()
        results = host.run(file, json.loads(network_file.read_text(encoding="utf-8")), params,
                           progress=lambda f, t: print(f"  {100 * f:5.1f}%  {t}", file=sys.stderr))
        print(summary(results))
        print(f"({time.time() - started:.1f} s)", file=sys.stderr)
        return 0
    print(__doc__.split("\n\n")[-1], file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
