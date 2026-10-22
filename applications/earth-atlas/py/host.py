"""The host of plug-ins that run in the server (ADR 0029): Python scripts,
imported into this process and run on numpy, and C++ programs, compiled and run
as processes behind one file contract. The page cannot tell them from its own
JavaScript plug-ins: the application's server exposes

    GET  /plugin-host/describe?file=plugins/x.py   the plug-in's identity, parameters and needs, with the sources that serve each
    POST /plugin-host/run  {file, chosen, params, context}   starts a job and returns {job}
    GET  /plugin-host/job/<id>                      the job's state, progress and results

and the host.js of the page loads, runs and shows the results as for any plug-in.
A field result is written as a product under data/cache/plugins/<id>/ and named
by its URL under /data/plugins/; sketches, tables and notes are returned as they are.

The file contract of a process plug-in (see plugins/cpp/plugin.hpp): the host
writes a directory with input.txt and one float32 file per input variable, runs
the program with --run <dir>, and reads output.txt and its float32 files back.
`--describe` prints the plug-in's description as JSON. Only files under the
application's plugins/ directory are run.

    python applications/earth-atlas/py/host.py describe plugins/vorticity.py       from the command line
    python applications/earth-atlas/py/host.py run plugins/vorticity.py [need=product|level ...] [name=value ...]
"""

from __future__ import annotations

import datetime as dt
import importlib.util
import json
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fields  # noqa: E402
import plugins as api  # noqa: E402

APP = HERE.parent
ROOT = APP.parents[1]
PLUGINS_DIR = APP / "plugins"
PRODUCTS = ROOT / "data" / "cache" / "products"
OUT = ROOT / "data" / "cache" / "plugins"
UTC = dt.timezone.utc

# The panel's order of models (src/layers/catalogue.js); an unknown model ranks last.
MODEL_ORDER = ["GFS", "AIGFS", "AIGEFS mean", "Pangu-Weather (run here)", "Pangu-Weather hindcast (run here)",
               "Shallow water (ours)", "Barotropic (ours)"]


def model_rank(model: str) -> int:
    return MODEL_ORDER.index(model) if model in MODEL_ORDER else len(MODEL_ORDER)


class Library:
    """The field products on disk, by name: the library's, and those plug-ins made before."""

    def __init__(self, products: Path = PRODUCTS):
        self.products = products
        self.meta: dict[str, dict] = {}
        self.cache: dict[str, fields.Field] = {}
        index_path = products / "index.json"
        index = json.loads(index_path.read_text(encoding="utf-8")) if index_path.exists() else {}
        for entry in index.get("fields", []):
            path = products / entry["file"]
            if path.exists():
                self.meta[entry["name"]] = json.loads(path.read_text(encoding="utf-8"))
                self.meta[entry["name"]]["_path"] = str(path)
        for path in OUT.glob("*/fields-*.json"):
            meta = json.loads(path.read_text(encoding="utf-8"))
            meta["_path"] = str(path)
            self.meta[meta["name"]] = meta

    def fields(self) -> list[dict]:
        return list(self.meta.values())

    def field(self, name: str) -> fields.Field:
        if name not in self.cache:
            self.cache[name] = fields.read(Path(self.meta[name]["_path"]))
        return self.cache[name]


# ---- plug-ins: Python modules and C++ programs ------------------------------------------------


def plugin_path(file: str) -> Path:
    """The file under plugins/, checked to be there."""
    path = (APP / file).resolve()
    if PLUGINS_DIR.resolve() not in path.parents or not path.is_file():
        raise ValueError(f"{file} is not a plug-in of this application (only files under plugins/ are run)")
    return path


def load_plugin(file: str):
    """The plug-in in `file`, and its runtime: 'python' or 'cpp'."""
    path = plugin_path(file)
    if path.suffix == ".py":
        spec = importlib.util.spec_from_file_location(f"plugin_{path.stem.replace('-', '_')}", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        if not hasattr(module, "create"):
            raise ValueError(f"{file} does not define create(api)")
        return module.create(api), "python"
    if path.suffix in (".cpp", ".cc", ".cxx"):
        return ProcessPlugin(path), "cpp"
    raise ValueError(f"{file}: the server runs .py and .cpp plug-ins; the page runs .js")


class ProcessPlugin:
    """A plug-in that is a program, compiled from C++ when its binary is older
    than its source or the header. It describes itself with --describe and runs
    as a Diagnostic on files: outputs on the first need's grid, at its times."""

    def __init__(self, source: Path):
        self.source = source
        exe = ".exe" if platform.system() == "Windows" else ""
        self.binary = OUT / "bin" / f"{source.stem}{exe}"
        self.spec = self.describe()
        self.id = self.spec["id"]
        self.title = self.spec.get("title", source.stem)
        self.description = self.spec.get("description", "")
        self.domain = self.spec.get("domain", "meteorology")
        self.parameters = self.spec.get("parameters", {})

    @property
    def product_name(self) -> str:
        return f"plugin-{self.id}"

    def needs(self) -> list[dict]:
        return self.spec.get("needs", [])

    def compile(self) -> None:
        compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
        if not compiler:
            raise RuntimeError("no C++ compiler on this computer (c++, g++ or clang++ is needed to run a C++ plug-in)")
        header = self.source.parent / "cpp" / "plugin.hpp"
        newest = max(self.source.stat().st_mtime, header.stat().st_mtime if header.exists() else 0)
        if self.binary.exists() and self.binary.stat().st_mtime >= newest:
            return
        self.binary.parent.mkdir(parents=True, exist_ok=True)
        done = subprocess.run([compiler, "-std=c++20", "-O2", "-o", str(self.binary), str(self.source)],
                              cwd=self.source.parent, capture_output=True, text=True)
        if done.returncode != 0:
            raise RuntimeError(f"{self.source.name} did not compile:\n{done.stderr[-1500:]}")

    def describe(self) -> dict:
        self.compile()
        done = subprocess.run([str(self.binary), "--describe"], capture_output=True, text=True)
        if done.returncode != 0:
            raise RuntimeError(f"{self.source.name} --describe failed:\n{done.stderr[-800:]}")
        return json.loads(done.stdout)

    def run(self, inputs: dict, params: dict, context) -> list[dict]:
        first_name = self.needs()[0]["name"]
        first = inputs[first_name]
        grid, times = first.product.grid, first.product.times or [context.time]
        work = Path(tempfile.mkdtemp(prefix="run-", dir=str(OUT / "tmp" if (OUT / "tmp").mkdir(parents=True, exist_ok=True) is None else OUT / "tmp")))
        try:
            lines = [f"time\t{int(context.time.timestamp() * 1000)}"]
            if context.pinned:
                lines.append(f"pinned\t{context.pinned[0]}\t{context.pinned[1]}")
            for name, value in params.items():
                lines.append(f"param\t{name}\t{value}")
            for need_name, inp in inputs.items():
                g, n = inp.product.grid, len(inp.product.times)
                lines.append("\t".join(map(str, ["need", need_name, inp.product.name, inp.model, inp.level or "", max(1, n),
                                                  g["lon0"], g["dlon"], g["nlon"], g["lat0"], g["dlat"], g["nlat"], int(bool(g.get("periodic")))])))
                lines.append("\t".join(["times", need_name, *[fields.iso(t) for t in inp.product.times]]))
                for var in inp.names:
                    data = np.stack([inp.product.slice(var, k) for k in range(max(1, n))]).astype("<f4")
                    file = f"{need_name}.{len([l for l in lines if l.startswith('data')])}.f32"
                    data.tofile(work / file)
                    lines.append(f"data\t{need_name}\t{var}\t{file}")
            (work / "input.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
            context.progress(0.1, "running")
            done = subprocess.run([str(self.binary), "--run", str(work)], capture_output=True, text=True)
            if done.returncode != 0:
                raise RuntimeError(f"{self.source.name} failed:\n{done.stderr[-1500:]}")
            outputs, show, notes = [], dict(self.spec.get("show", {})), []
            for line in (work / "output.txt").read_text(encoding="utf-8").splitlines():
                parts = line.split("\t")
                if parts[0] == "field":
                    _, name, long_name, standard_name, units, level, difference, file = parts
                    raw = np.fromfile(work / file, dtype="<f4").reshape(len(times), grid["nlat"], grid["nlon"])
                    entry = {"name": name, "long_name": long_name, "standard_name": standard_name or None, "units": units,
                             "level": level or first.level, "slices": list(raw)}
                    if difference == "1":
                        entry["difference"] = True
                    outputs.append(entry)
                elif parts[0] == "show" and len(parts) == 3:
                    show[parts[1]] = parts[2]
                elif parts[0] == "note":
                    notes.append(api.note(parts[1]))
                elif parts[0] == "error":
                    raise RuntimeError(parts[1])
            made = api.make_product(self.product_name, self.title, grid, times, outputs, model=f"{self.title} (plug-in)", domain=self.domain,
                                    source={"dataset": first.meta["name"], "provider": f"computed by a C++ program in the server from {first.model}",
                                            "credit": f"{self.title} (plug-in), from {first.model}"})
            return [api.field(made, show), *notes]
        finally:
            shutil.rmtree(work, ignore_errors=True)


# ---- needs resolved to the library's products ------------------------------------------------


def resolve(plugin, library: Library) -> list[dict]:
    """The products that can serve each need: one candidate per product and level
    carrying every standard name asked, the plug-in's own output excepted."""
    out = []
    for need in plugin.needs():
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
        out.append({**need, "candidates": candidates})
    return out


def describe(file: str, library: Library | None = None) -> dict:
    """What the page shows for a plug-in before it runs."""
    library = library or Library()
    plugin, runtime = load_plugin(file)
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
                   "candidates": [{"key": c["key"], "model": c["model"]} for c in n["candidates"]]} for n in needs],
    }


def run(file: str, chosen: dict, params: dict, context: api.Context, library: Library | None = None) -> tuple[object, list[dict]]:
    """Runs a plug-in on the chosen sources and returns (plugin, results); fields are written as products."""
    library = library or Library()
    plugin, _ = load_plugin(file)
    inputs = {}
    for need in resolve(plugin, library):
        c = next((x for x in need["candidates"] if x["key"] == chosen.get(need["name"])), None)
        if not c:
            raise ValueError(f"nothing in the library serves \"{need.get('label', need['name'])}\"")
        inputs[need["name"]] = api.Input(library.field(c["product"]), c["names"], c["variables"], c["meta"], c["model"], c["level"])
    context._fields = library.fields
    context._field = library.field
    results = plugin.run(inputs, params, context) or []
    out_dir = OUT / plugin.id
    shutil.rmtree(out_dir, ignore_errors=True)
    written = []
    for r in results:
        if r["kind"] == "field":
            path = fields.write(r["product"], out_dir)
            written.append({"kind": "field", "url": f"data/plugins/{plugin.id}/{path.name}", "show": r.get("show", {})})
        else:
            written.append(r)
    return plugin, written


class Jobs:
    """Runs of plug-ins for the server, one thread each; the page polls their state."""

    def __init__(self):
        self.jobs: dict[str, dict] = {}
        self.lock = threading.Lock()

    def start(self, file: str, chosen: dict, params: dict, context: dict) -> dict:
        plugin_path(file)
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
        when = dt.datetime.fromtimestamp(float(context.get("time", time.time() * 1000)) / 1000, UTC)
        ctx = api.Context(time=when, pinned=context.get("pinned"), progress=progress)
        try:
            _, job["results"] = run(job["file"], chosen or {}, params or {}, ctx)
            job["state"] = "done"
        except Exception as error:         # reported to the page, which shows it
            job["state"], job["message"] = "failed", f"{type(error).__name__}: {error}"
            traceback.print_exc()
        finally:
            job["finished"] = time.time()


def main(argv: list[str]) -> int:
    if len(argv) < 2 or argv[0] not in ("describe", "run"):
        print(__doc__)
        return 2
    if argv[0] == "describe":
        print(json.dumps(describe(argv[1]), indent=1, ensure_ascii=False))
        return 0
    d = describe(argv[1])
    chosen = {n["name"]: n["chosen"] for n in d["needs"]}
    params = {k: v["default"] for k, v in d["parameters"].items()}
    for item in argv[2:]:
        key, _, value = item.partition("=")
        if key in chosen:
            chosen[key] = value
        else:
            params[key] = float(value) if value.replace(".", "", 1).lstrip("-").isdigit() else value
    plugin, results = run(argv[1], chosen, params, api.Context(time=dt.datetime.now(UTC), progress=lambda f, t="": print(f"  {f:4.0%} {t}")))
    for r in results:
        print(r.get("url") or r.get("title") or r.get("text") or r["kind"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
