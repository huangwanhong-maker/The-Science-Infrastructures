"""Serve the Economy Atlas on this computer.

The page is served from this directory, and the code it shares from the repository's
``domains/``, ``foundations/`` and ``visualization/`` under the same paths, with the data cache
at ``/data/``. What it shows is computed by the engine (``engine.py``: the economics domain's
C++, on physis's discrete exterior calculus and morphe's complexes) from the data library's
products. The protocol is in the README, §3:

    GET  /api/engine                     whether the engine is built and current, and the limits of a request
    GET  /api/data                       which products are built, and the datasets of one's own
    POST /api/map                        {"year", "dataset"?}: flows in a year, their exchange, drive and circulation
    POST /api/timeline                   {"dataset"?}: the same shares for every year
    POST /api/spacetime                  {"first", "last", "currencies"?, "dataset"?}: currencies held and exchanged
    POST /api/production                 {"table", "value_added_coefficients"?, "final_demand"?}: a Leontief model
    POST /api/circuits                   {"network"}: a transfer network's deficits, prices and circulation
    POST /api/cancel                     stop the computations in progress

Datasets of one's own (``datasets.py``: flows between places, or currencies held and exchanged, as
CSV or JSON) and plug-ins (Python scripts and C++ programs in ``plugins/``, run by the API's
``science.networks`` on the network of a dataset in a year, ADR 0046):

    GET    /api/datasets                 the library's products and the datasets kept here
    POST   /api/datasets                 {"kind": "flows"|"currencies", "title", "files": [{name, text}], "base"?}
    DELETE /api/datasets/<id>            forget one
    GET    /api/plugins                  the plug-ins in plugins/, with their descriptions once known
    GET    /api/plugins/describe?file=   one plug-in's description (a C++ one is compiled first)
    POST   /api/plugins/run              {"file", "params", "dataset"?, "year"?}: starts a job; returns {"job"}
    GET    /api/plugins/job/<id>         its state, progress and results
    POST   /api/plugins/job/<id>/cancel  stop it

The engine and the data are made by the application's recipe (``recipe.json``), which the page
may run:

    GET  /setup                          its steps, and whether each is done
    POST /setup/run[?again=1]            run them (again)

Answers for a year or a span are kept, keyed by the engine's build and the product, so that each
is computed once. Only 127.0.0.1 is bound unless ``--host`` says otherwise; on a public address
the recipe, the datasets and the plug-ins need the server's token (``--token``), and requests are
held within smaller limits.

Usage::

    python applications/economy-atlas/serve.py [--host 127.0.0.1] [--port 8910] [--open]

then open http://127.0.0.1:8910/ in a browser.
"""

from __future__ import annotations

import argparse
import hashlib
import http.server
import json
import sys
import threading
import time
import urllib.parse
import uuid
import webbrowser
from pathlib import Path

APP = Path(__file__).resolve().parent
ROOT = APP.parents[1]
MOUNTS = {
    "/domains/": ROOT / "domains",
    "/foundations/": ROOT / "foundations",
    "/visualization/": ROOT / "visualization",
    "/data/": ROOT / "data" / "cache",
}
TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "text/javascript; charset=utf-8",
    ".mjs": "text/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json",
    ".geojson": "application/geo+json",
    ".md": "text/markdown; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
}
PRODUCTS = ROOT / "data" / "cache" / "products"
TRADE = PRODUCTS / "trade-flows.json"
CURRENCIES = PRODUCTS / "currency-spacetime.json"
KEPT = ROOT / "data" / "cache" / "economy-atlas" / "answers"

sys.path.insert(0, str(APP))
sys.path.insert(0, str(ROOT / "api" / "python"))
import datasets  # noqa: E402
import engine  # noqa: E402

RECIPE = None
WRITES = {"allowed": True, "token": None}
LIMITS = {"bytes": 1024 * 1024, "upload": 48 * 1024 * 1024,
          "seconds": {"map": 120.0, "timeline": 600.0, "spacetime": 120.0, "production": 30.0, "circuits": 30.0, "plugin": 600.0}}
RUNNING: set = set()
SLOTS = threading.BoundedSemaphore(3)
OPS = {"/api/map": "map", "/api/timeline": "timeline", "/api/spacetime": "spacetime", "/api/production": "production",
       "/api/circuits": "circuits"}
DATASETS = datasets.Datasets()
PLUGINS = APP / "plugins"
PLUGIN_WORK = ROOT / "data" / "cache" / "economy-atlas" / "plugins"
_host = {"host": None}


def plugin_host():
    """The API's host of plug-ins over networks of places, bound to this atlas's plugins/."""
    if _host["host"] is None:
        from science import networks
        _host["host"] = networks.Host(PLUGINS, PLUGIN_WORK, include_dirs=engine.include_dirs())
    return _host["host"]


def recipe_job():
    global RECIPE
    if RECIPE is None:
        sys.path.insert(0, str(ROOT / "data" / "acquisition"))
        import recipe
        RECIPE = recipe.Job(APP / "recipe.json", ROOT)
    return RECIPE


def write_refusal(handler) -> str | None:
    if not WRITES["allowed"]:
        return "this server does not run its recipe: start it with --token to allow it with a token"
    token = WRITES["token"]
    if token and handler.headers.get("Authorization", "") != f"Bearer {token}":
        return "running the recipe on this server needs its token"
    return None


def data_state() -> dict:
    """The products built, what the page may ask of them, and the datasets of one's own."""
    out: dict = {"trade": None, "currencies": None, "datasets": DATASETS.list()}
    try:
        index = json.loads((PRODUCTS / "index.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        index = {}
    if TRADE.exists():
        entry = next((e for e in index.get("trade", []) if e["name"] == "trade-flows"), None)
        out["trade"] = entry or {"name": "trade-flows"}
    if CURRENCIES.exists():
        entry = next((e for e in index.get("currencies", []) if e["name"] == "currency-spacetime"), None)
        out["currencies"] = entry or {"name": "currency-spacetime"}
    return out


def answer_key(op: str, request: dict, product: Path | None) -> str:
    h = hashlib.sha256(json.dumps({"engine": engine.read_build() or {}, "op": op, "request": request}, sort_keys=True).encode("utf-8"))
    if product is not None:
        stat = product.stat()
        h.update(f"{product.name}:{stat.st_size}:{stat.st_mtime_ns}".encode())
    return h.hexdigest()[:32]


def compute(op: str, body: dict) -> tuple[str | dict, int]:
    """The engine's answer to a request of the page (JSON text), or a message and its status."""
    toolchain = (engine.read_build() or {}).get("toolchain") or "native"
    request: dict = {"op": op}
    product = None
    dataset = body.get("dataset") if isinstance(body.get("dataset"), str) else None
    try:
        if op in ("map", "timeline"):
            product = DATASETS.product(dataset, "flows")
            request["flows"] = engine.path_for(toolchain, product)
            if op == "map" and isinstance(body.get("year"), int):
                request["year"] = body["year"]
        elif op == "spacetime":
            product = DATASETS.product(dataset, "currencies")
            request["record"] = engine.path_for(toolchain, product)
            for key in ("first", "last"):
                if isinstance(body.get(key), str):
                    request[key] = body[key]
            if isinstance(body.get("currencies"), list):
                request["currencies"] = [c for c in body["currencies"] if isinstance(c, str)]
        elif op == "production":
            if not isinstance(body.get("table"), dict):
                return {"message": "a production request needs a table"}, 400
            request["table"] = body["table"]
            for key in ("value_added_coefficients", "final_demand"):
                if isinstance(body.get(key), list):
                    request[key] = body[key]
        elif op == "circuits":
            if not isinstance(body.get("network"), dict):
                return {"message": "a circuits request needs a network"}, 400
            request["network"] = body["network"]
    except datasets.DatasetError as error:
        return {"message": str(error)}, 404
    state = engine.status()
    if not state["current"]:
        return {"message": "the engine is not built from the present sources: build it first", "engine": state}, 503
    keyed = {k: v for k, v in request.items() if k not in ("flows", "record")}
    kept = KEPT / f"{answer_key(op, keyed, product)}.json" if product is not None else None
    if kept is not None and kept.exists():
        return kept.read_text(encoding="utf-8"), 200
    if not SLOTS.acquire(timeout=120):
        return {"message": "the engine is busy with other requests; try again in a moment"}, 409
    run = None
    try:
        run = engine.Run(json.dumps(request))
        RUNNING.add(run)
        code, out, err = run.execute(timeout=LIMITS["seconds"][op])
    finally:
        RUNNING.discard(run)
        SLOTS.release()
    try:
        answer = json.loads(out) if out.strip() else None
    except ValueError:
        answer = None
    if isinstance(answer, dict) and "error" in answer:
        return {"message": str(answer["error"])}, 400
    if code != 0 or answer is None:
        message = err.replace("\0", "").strip().splitlines()[-1] if err.strip() else "the engine failed"
        return {"message": message}, 499 if code < 0 else 500
    if kept is not None:
        try:
            KEPT.mkdir(parents=True, exist_ok=True)
            kept.write_text(out, encoding="utf-8")
        except OSError:
            pass
    return out, 200


def network_of(dataset: str | None, year: int | None) -> dict:
    """The network a plug-in is given: the dataset's places, its gross flows in the year, and the
    engine's decomposition of them (the map's answer), every place named by its id."""
    product = json.loads(DATASETS.product(dataset, "flows").read_text(encoding="utf-8"))
    out, code = compute("map", {"dataset": dataset, **({"year": year} if year is not None else {})})
    if isinstance(out, dict):
        raise RuntimeError(out.get("message", "the map could not be computed"))
    A = json.loads(out)
    ids = [e["iso3"] for e in product["economies"]]
    return {
        "title": product.get("title", ""), "year": A["year"], "unit": product.get("unit", ""),
        "places": [{"id": e["iso3"], "name": e["name"], "latitude": e["latitude"], "longitude": e["longitude"],
                    "group": e.get("continent", "")} for e in product["economies"]],
        "flows": [{"from": ids[int(i)], "to": ids[int(j)], "value": v} for i, j, v in product["flows"][str(A["year"])]],
        "decomposition": {
            "attraction": {e["iso3"]: e["attraction"] for e in A["economies"]},
            "legs": [{"from": ids[l[0]], "to": ids[l[1]], "net": l[2], "driven": l[3], "local": l[4], "global": l[5]} for l in A["legs"]],
            "triangles": [{"places": [ids[i] for i in t["places"]], "circulation": t["curl"]} for t in A["triangles"]],
            "totals": A["totals"],
        },
    }


class PluginJobs:
    """Plug-in runs in threads of this server, each followed by the page until it ends."""

    KEEP = 50

    def __init__(self):
        self.jobs: dict[str, dict] = {}
        self.descriptions: dict[str, dict] = {}
        self.lock = threading.Lock()

    def listing(self, host) -> list[dict]:
        """The plug-ins, a Python one described at once (a fraction of a second), a C++ one when it is
        first opened (its first compiling takes half a minute)."""
        out = []
        for f in host.files():
            if f["language"] == "python" and f["file"] not in self.descriptions:
                try:
                    self.describe(host, f["file"])
                except (ValueError, RuntimeError, OSError) as error:
                    f = {**f, "error": str(error)}
            out.append({**f, **({"description": self.descriptions[f["file"]]} if f["file"] in self.descriptions else {})})
        return out

    def describe(self, host, file: str) -> dict:
        d = host.describe(file)
        self.descriptions[file] = d
        return d

    def start(self, host, file: str, params: dict, dataset, year) -> dict:
        ident = uuid.uuid4().hex[:12]
        job = {"id": ident, "file": file, "state": "running", "progress": 0.0, "text": "preparing the network",
               "started": time.time(), "cancel": threading.Event()}
        with self.lock:
            self.jobs[ident] = job
            for old in sorted(self.jobs.values(), key=lambda j: j["started"])[:-self.KEEP]:
                self.jobs.pop(old["id"], None)
        threading.Thread(target=self._run, args=(host, job, params, dataset, year), daemon=True).start()
        return {"job": ident}

    def _run(self, host, job: dict, params: dict, dataset, year) -> None:
        def progress(fraction: float, text: str = "") -> None:
            job["progress"], job["text"] = float(fraction), text or job["text"]
        try:
            network = network_of(dataset, year)
            job["text"] = "running"
            job["results"] = host.run(job["file"], network, params, progress=progress, cancel=job["cancel"],
                                      seconds=LIMITS["seconds"]["plugin"])
            job["year"], job["state"], job["progress"] = network["year"], "done", 1.0
        except Exception as error:  # a plug-in's failure is reported to the page, never raised in the server
            job["state"] = "cancelled" if job["cancel"].is_set() else "failed"
            job["message"] = str(error)
        job["seconds"] = round(time.time() - job["started"], 2)

    def status(self, ident: str) -> tuple[dict, int]:
        job = self.jobs.get(ident)
        if job is None:
            return {"message": f"no job {ident}"}, 404
        return {k: v for k, v in job.items() if k != "cancel"}, 200

    def cancel(self, ident: str) -> tuple[dict, int]:
        job = self.jobs.get(ident)
        if job is None:
            return {"message": f"no job {ident}"}, 404
        job["cancel"].set()
        return {"cancelled": ident}, 200


PLUGIN_JOBS = PluginJobs()


class Handler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {**http.server.SimpleHTTPRequestHandler.extensions_map, **TYPES}

    def translate_path(self, path: str) -> str:
        path = urllib.parse.unquote(path.split("?", 1)[0].split("#", 1)[0])
        prefix = next((p for p in MOUNTS if path.startswith(p)), None)
        base, rel = (MOUNTS[prefix], path[len(prefix):]) if prefix else (APP, path.lstrip("/"))
        target = (base / rel).resolve()
        if target != base and base not in target.parents:
            return str(base / "__outside__")
        return str(target)

    def do_GET(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path == "/api/engine":
            return self.reply_json({**engine.status(), "limits": LIMITS, "running": len(RUNNING)}, 200)
        if url.path == "/api/data":
            return self.reply_json(data_state(), 200)
        if url.path == "/setup":
            return self.reply_json(*self.recipe_request(url))
        if url.path == "/api/datasets":
            return self.datasets_request("GET", url)
        if url.path.startswith("/api/plugins"):
            return self.plugins_request("GET", url)
        return super().do_GET()

    def do_DELETE(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path.startswith("/api/datasets/"):
            return self.datasets_request("DELETE", url)
        self.send_error(404)

    def do_POST(self) -> None:
        url = urllib.parse.urlparse(self.path)
        if url.path == "/setup/run":
            return self.reply_json(*self.recipe_request(url))
        if url.path == "/api/cancel":
            runs = list(RUNNING)
            for run in runs:
                run.cancel()
            return self.reply_json({"cancelled": len(runs)}, 200)
        if url.path in OPS:
            return self.engine_request(OPS[url.path])
        if url.path == "/api/datasets":
            return self.datasets_request("POST", url)
        if url.path.startswith("/api/plugins/"):
            return self.plugins_request("POST", url)
        self.send_error(404)

    def read_body(self, limit: int | None = None) -> dict | None:
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > (limit or LIMITS["bytes"]):
            self.reply_json({"message": "the request is too large"}, 413)
            return None
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
            if not isinstance(body, dict):
                raise ValueError("a request is a JSON object")
            return body
        except ValueError as error:
            self.reply_json({"message": f"not a request: {error}"}, 400)
            return None

    def engine_request(self, op: str) -> None:
        body = self.read_body()
        if body is None:
            return
        out, code = compute(op, body)
        if isinstance(out, dict):
            return self.reply_json(out, code)
        self.reply_text(out, code)

    def datasets_request(self, method: str, url) -> None:
        if method == "GET":
            return self.reply_json({"datasets": DATASETS.list()}, 200)
        refusal = write_refusal(self)
        if refusal:
            return self.reply_json({"message": refusal}, 403)
        if method == "DELETE":
            ident = url.path[len("/api/datasets/"):]
            return self.reply_json({"removed": DATASETS.remove(ident)}, 200)
        body = self.read_body(LIMITS["upload"])
        if body is None:
            return
        files = body.get("files")
        if not isinstance(files, list) or not all(isinstance(f, dict) for f in files):
            return self.reply_json({"message": "a dataset is given as files: [{name, text}]"}, 400)
        try:
            about = DATASETS.add(str(body.get("kind") or ""), str(body.get("title") or ""), files, str(body.get("base") or ""))
        except datasets.DatasetError as error:
            return self.reply_json({"message": str(error)}, 400)
        self.reply_json(about, 200)

    def plugins_request(self, method: str, url) -> None:
        try:
            host = plugin_host()
        except ImportError as error:
            return self.reply_json({"message": f"the API's plug-in host is missing: {error}"}, 503)
        path = url.path
        if method == "GET" and path == "/api/plugins":
            return self.reply_json({"plugins": PLUGIN_JOBS.listing(host)}, 200)
        if method == "GET" and path == "/api/plugins/describe":
            file = dict(urllib.parse.parse_qsl(url.query)).get("file", "")
            try:
                return self.reply_json(PLUGIN_JOBS.describe(host, file), 200)
            except (ValueError, RuntimeError, OSError) as error:
                return self.reply_json({"message": str(error), "file": file}, 400)
        if method == "GET" and path.startswith("/api/plugins/job/"):
            return self.reply_json(*PLUGIN_JOBS.status(path[len("/api/plugins/job/"):]))
        refusal = write_refusal(self)
        if refusal:
            return self.reply_json({"message": refusal.replace("runs its recipe", "runs plug-ins")}, 403)
        if method == "POST" and path == "/api/plugins/run":
            body = self.read_body()
            if body is None:
                return
            if not isinstance(body.get("file"), str):
                return self.reply_json({"message": "name the plug-in's file"}, 400)
            year = body.get("year") if isinstance(body.get("year"), int) else None
            params = body.get("params") if isinstance(body.get("params"), dict) else {}
            return self.reply_json(PLUGIN_JOBS.start(host, body["file"], params, body.get("dataset"), year), 200)
        if method == "POST" and path.startswith("/api/plugins/job/") and path.endswith("/cancel"):
            return self.reply_json(*PLUGIN_JOBS.cancel(path[len("/api/plugins/job/"):-len("/cancel")]))
        self.send_error(404)

    def recipe_request(self, url) -> tuple[dict, int]:
        refusal = write_refusal(self)
        policy = {"writable": refusal is None, "refusal": refusal}
        try:
            if url.path == "/setup":
                return {**recipe_job().view(), **policy}, 200
            if refusal:
                return {"message": refusal, **policy}, 403
            query = dict(urllib.parse.parse_qsl(url.query))
            return {**recipe_job().start(None, again=query.get("again") == "1"), **policy}, 200
        except RuntimeError as error:
            return {"message": str(error)}, 409

    def reply_json(self, body: dict, code: int) -> None:
        self.reply_text(json.dumps(body, ensure_ascii=False), code)

    def reply_text(self, text: str, code: int) -> None:
        data = text.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def end_headers(self) -> None:
        self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def log_message(self, format: str, *args) -> None:
        if args and str(args[1]).startswith(("4", "5")):
            super().log_message(format, *args)


class Server(http.server.ThreadingHTTPServer):
    request_queue_size = 64
    daemon_threads = True


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="Serve the Economy Atlas locally.")
    parser.add_argument("--host", default="127.0.0.1", help="the address to bind; 0.0.0.0 for every interface")
    parser.add_argument("--port", type=int, default=8910)
    parser.add_argument("--open", action="store_true", help="open the page in the default browser")
    parser.add_argument("--token", help="on a public address, the token that allows running the recipe")
    args = parser.parse_args(argv)
    public = args.host not in ("127.0.0.1", "localhost", "::1")
    if public:
        WRITES.update(allowed=bool(args.token), token=args.token)
        LIMITS.update(bytes=256 * 1024, upload=4 * 1024 * 1024,
                      seconds={"map": 60.0, "timeline": 300.0, "spacetime": 60.0, "production": 10.0, "circuits": 10.0, "plugin": 120.0})
    if not engine.status()["current"]:
        print("The engine is not built from the present sources. Build it from the page, or run:\n"
              "  python data/acquisition/recipe.py applications/economy-atlas/recipe.json", file=sys.stderr)
    server = Server((args.host, args.port), Handler)
    url = f"http://127.0.0.1:{args.port}/"
    where = url if not public else f"{url} and on {args.host}"
    print(f"Economy Atlas at {where}  (Ctrl+C to stop)", flush=True)
    if args.open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
